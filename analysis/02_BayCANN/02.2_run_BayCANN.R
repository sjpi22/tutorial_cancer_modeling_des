###########################  BayCANN   #########################################
#
#  Objective: Script to perform an emulator-based Bayesian calibration
########################### <<<<<>>>>> #########################################

# Sources: Jalal H, Trikalinos TA, Alarid-Escudero F. BayCANN: Streamlining 
# Bayesian Calibration With Artificial Neural Network Metamodeling. Front 
# Physiol. 2021 May 25;12:662314. doi: 10.3389/fphys.2021.662314. PMID: 
#   34113262; PMCID: PMC8185956.

# Pineda-Antunez C, Seguin C, van Duuren LA, Knudsen AB, Davidi B, de Lima PN, 
# Rutter C, Kuntz KM, Lansdorp-Vogelaar I, Collier N, Ozik J, Alarid-Escudero F. 
# Emulator-based Bayesian calibration of the CISNET colorectal cancer models. 
# medRxiv [Preprint]. 2024 Feb 5:2023.02.27.23286525. doi: 
#   10.1101/2023.02.27.23286525. PMID: 36909607; PMCID: PMC10002763.

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(keras3)   # Install tensorflow beforehand
library(tfruns)
library(rstan)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(doParallel)
library(patchwork)
library(MASS)
library(bestNormalize)
library(data.table)
library(assertthat)
rstan_options(auto_write = TRUE)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 file paths

######## 2.1.1 directories
outpath <- "output/calibration/BayCANN"
outpath_plots <- file.path(outpath, "plots")
outpath_runs <- file.path(outpath, "runs")

######## 2.1.2 inputs
file_params <- "data/calibration_params.rds"
file_model <- file.path("analysis/02_BayCANN/model.R")
file_sample <- file.path(outpath, "sample_BayCANN.rds")
file_stan <- "stan/post_multi_perceptron_flex.stan"

######## 2.1.3 outputs
file_keras_model <- file.path(outpath, "model_ANN.keras") # File path for the compiled model
file_test_outputs <- file.path(outpath, "test_ANN_outputs.rds") # Scaled test outputs for validating STAN model
file_fig_history <- file.path(outpath_plots, "fig_ANN_history.png")
file_fig_valid <- file.path(outpath_plots, "fig_ANN_validation.png")
file_fig_posterior <- file.path(outpath_plots, "fig_posterior_chains_BayCANN.png")
file_fig_corr <- file.path(outpath_plots, "fig_posterior_corr_BayCANN.png")
file_fig_prior_post <- file.path(outpath_plots, "fig_prior_posterior_BayCANN.png")
file_stan_model <- file.path(outpath, "model_stan_BayCANN.rds")
file_posterior <- file.path(outpath, "calibrated_posteriors_BayCANN.csv")
file_baycann_params <- file.path(outpath, "parameters_BayCANN.rds")


###### 2.2 parameters for data preparation
v_groups_categorical <- c("stage_distr") # Target groups that are categorical (percentages sum to 1)
scale_type <- 3  # 1: for scale from -1 to 1; 2: for standardization ; 3: for scale from 0 to 1
train_split <- 0.8 # Train vs. test
valid_split <- 0.8 # Train vs. validation within train set for hyperparameter tuning
Normalize_outputs   <- FALSE  # TRUE if we want to normalize outputs 
Scale_inputs        <- TRUE   # TRUE if we want to scale inputs
Scale_outputs       <- TRUE   # TRUE if we want to scale outputs 
Remove_outliers     <- FALSE  # TRUE if we want to remove outliers
Standardize_targets <- FALSE  # TRUE if we want to standardize targets
Saved_data          <- FALSE  # TRUE if we want to saved the data
Selected_targets    <- FALSE  # TRUE if we want to use an specific list of calibration targets

###### 2.3 parameters for ANN
verbose          <- 0
n_batch_size     <- 128
n_epochs         <- 1000
patience         <- 100
init_W = NULL

# Other options to initialize weights
#init_W=initializer_random_uniform(minval = -0.7, maxval = 0.7, seed = 2312)   ###initialization of weights with uniform distribution
#init_W=initializer_random_normal(mean = 0, stddev = 0.1, seed = 2312)  ###initialization of weights with normal distribution

###### 2.4 parameters for ANN hyperparameter tuning
v_hidden_layers <- seq(6)
v_hidden_nodes <- c(32, 64, 96, 128, 256)
v_dropout <- c(TRUE, FALSE)
v_dropout_rate <- c(0.25)
v_activation_fun <- c("relu", "tanh", "sigmoid")
p_hp_sample      <- 0.1 # Proportion of hyperparameter possibilities to sample

###### 2.5 parameters for Bayesian calibration
n_iter   <- 300000
n_thin   <- 100
n_chains <- 4

###### 2.6 parameters for validating Stan representation of ANN
validate_stan_ANN <- TRUE
eps_stan <- 1e-4


#### 3. Pre-processing actions  ===========================================

# Create directory for plots if it does not exist
outpath_split <- unlist(strsplit(outpath_plots, split = "/"))
for (i in 1:length(outpath_split)) {
  dir.create(do.call(file.path, as.list(outpath_split[1:i])))
}

# Load model and calibration parameters
l_params_calib <- readRDS(file_params)

# Set columns to scale (currently only include targets that are not categorical)
# Select all columns if desired
scale_cols <- which(!l_params_calib$df_true_targets$target_group %in% v_groups_categorical)

# Create data frame of output layer parameters (activation functions, weights, loss function) based on target and function groups
df_fn_grps <- l_params_calib$df_true_targets %>%
  mutate(target_groups = factor(target_groups, levels = unique(l_params_calib$df_true_targets$target_groups))) %>% # To preserve order of groups
  group_by(target_groups) %>%
  summarize(n_targets = n(), .groups = "drop") %>%
  mutate(
    fn_grp = 1:n(),
    activation = ifelse( # Assign activation functions for output layer
      target_groups %in% v_groups_categorical, "softmax", # Categorical targets that sum to 1
      "sigmoid")) %>%
  mutate(
    loss_fn = case_when( # Assign loss functions
      activation %in% "softmax" ~ "categorical_crossentropy",
      activation %in% "sigmoid" ~ "mean_squared_error"), 
    loss_weight = 1, # Assign loss weights
    metric = case_when( # Assign evaluation metrics
      target_groups %in% v_groups_categorical ~ "accuracy", 
      .default="mae"
    ))

# Set seed
seed <- l_params_calib$l_params_all$seed
set.seed(seed)

# Set cores
options(mc.cores = parallel::detectCores() - l_params_calib$n_cores_reserved_local)


#### 4. Load the training and test data for the simulations =======================

BayCANN_sample <- readRDS(file_sample)

###### 4.1 Input Parameters ####
data_sim_param <- as.matrix(BayCANN_sample$m_param_samp)

###### 4.2 Model Outputs ####
data_sim_target <- as.matrix(BayCANN_sample$m_calib_outputs)


#### 5. Remove outliers from model outputs ####

if (Remove_outliers) {
  vec_out <- outlier_vector(data_sim_target)
  data_sim_param  <- data_sim_param[!vec_out,]
  data_sim_target <- data_sim_target[!vec_out,]
}


#### 6. Train/test partition ======================================================

# Create train-test split
train.rows <- sample.int(nrow(data_sim_param), 
                         size = train_split * nrow(data_sim_param),
                         replace = FALSE)

data_sim_param_train  <- data_sim_param[train.rows,]
data_sim_param_test   <- data_sim_param[-train.rows,]

data_sim_target_train <- data_sim_target[train.rows,]
data_sim_target_test  <- data_sim_target[-train.rows,]

prepared_data <- prepare_data(xtrain = data_sim_param_train,
                              ytrain = data_sim_target_train,
                              xtest  = data_sim_param_test,
                              ytest  = data_sim_target_test,
                              scale  = scale_type,
                              scale_cols = scale_cols)

list2env(prepared_data, envir = .GlobalEnv)

# Reshape outputs for keras
ytrain_scaled_reshape <- list()
ytest_scaled_reshape <- list()
for (grp in unique(df_fn_grps$target_groups)) {
  ytrain_scaled_reshape <- c(ytrain_scaled_reshape, 
                             list(ytrain_scaled[, l_params_calib$df_true_targets$target_groups %in% grp]))
  ytest_scaled_reshape <- c(ytest_scaled_reshape, 
                            list(ytest_scaled[, l_params_calib$df_true_targets$target_groups %in% grp]))
}


####  7. Load and scale targets and their SE ============================================

# Load targets and SE
true_targets_mean  <- l_params_calib$df_true_targets$targets
true_targets_se   <- l_params_calib$df_true_targets$se

# Scale targets and SE
if (scale_type==1) {
  y_targets <- 2 * (true_targets_mean - ymins) / (ymaxs - ymins) - 1   ## range from -1 to 1
  y_targets_se <- 2 * (true_targets_se) / (ymaxs - ymins)
}

if (scale_type==2) {
  y_targets <- (true_targets_mean - ymeans)/ysds   ## Standardization
  y_targets_se <-(true_targets_se)/ysds
}

if (scale_type==3) {
  y_targets <- (true_targets_mean - ymins) / (ymaxs - ymins)   ## range from 0 to 1
  y_targets_se <- (true_targets_se) / (ymaxs - ymins)
}

# Rescale only relevant cols
y_targets_final <- true_targets_mean
y_targets_final[scale_cols] <- y_targets[scale_cols]
y_targets <- y_targets_final

y_targets_se_final <- true_targets_se
y_targets_se_final[scale_cols] <- y_targets_se[scale_cols]
y_targets_se <- y_targets_se_final

y_targets <- t(as.matrix(y_targets))
y_targets_se <- t(as.matrix(y_targets_se))   


#### 8. Keras Section BayCANN ==============================================

###### 8.1 Hyperparameter tuning  ####

# Set hyperparameter flags
l_hyperparams <- list(
  n_hidden_layers = v_hidden_layers,
  n_hidden_nodes = v_hidden_nodes,
  dropout = v_dropout,
  dropout_rate = v_dropout_rate,
  activation_fun = v_activation_fun)

# Set validation data to NULL to use random training split
l_validation_data <- NULL

# Run hyperparameter tuning
runs <- tuning_run(
  file_model, 
  sample = p_hp_sample, 
  runs_dir = outpath_runs, 
  flags = l_hyperparams,
  confirm = TRUE, # Set false to override confirmation prompt
)

###### 8.2 Train best model  ####

# See runs in directory
ls_runs_df <- ls_runs(runs_dir = outpath_runs)
runs <- ls_runs_df # If interrupted

# Select best run
best_run <- runs[order(runs$metric_val_loss), ][1,]
l_hyperparams_best <- as.list(best_run)[paste0("flag_", names(l_hyperparams))]
names(l_hyperparams_best) <- names(l_hyperparams)

# Set validation data as test data
l_validation_data <- list(xtest_scaled, ytest_scaled_reshape)

# Train model from best run
keras.time <- proc.time()
run <- training_run(file_model, 
                    flags = l_hyperparams_best,
                    run_dir = file.path(outpath_runs, "best"))
t_training <- proc.time() - keras.time #keras ann fitting time
t_training <- t_training[3]/60
print(t_training)

# Save model
save_model(model, file_keras_model, overwrite = TRUE)  # Save the model and reload
model <- load_model(file_keras_model)

# Plot loss function and accuracy function
png(filename = file_fig_history)
plot(history)   
dev.off()

###### 8.3 Evaluate model  ####

# Model performance evaluation
acc_err <- model %>% evaluate(xtest_scaled, ytest_scaled_reshape) 

# Save predictions on test data
pred <- model %>% predict(xtest_scaled)
ytest_scaled_pred <- data.frame(pred)
colnames(ytest_scaled_pred) <- y_names
saveRDS(ytest_scaled_pred, file = file_test_outputs) 

# Plot ANN validation
ann_valid <- rbind(data.frame(sim = 1:n_test, ytest_scaled, type = "model"),
                   data.frame(sim = 1:n_test, ytest_scaled_pred, type = "pred"))
ann_valid_transpose <- ann_valid %>%
  pivot_longer(cols = -c(sim, type)) %>%
  pivot_wider(id_cols = c(sim, name), names_from = type, values_from = value)

n_partition <- round(sqrt(n_outputs))
n_part_bach <- floor(n_outputs/n_partition)

ann_valid_transpose <- arrange(ann_valid_transpose, desc(name))

ann_plot_full <- ggplot(data = ann_valid_transpose, aes(x = model, y = pred)) +
  geom_point(alpha = 0.5, color = "tomato") +
  facet_wrap(~name, scales="free", nrow=n_partition) +
  xlab("Model outputs (scaled)") +
  ylab("ANN predictions (scaled)") +
  theme_bw()

ggsave(filename = file_fig_valid,
       ann_plot_full,
       width = 2*n_part_bach, height = 2*n_partition)


#### 9. Stan ==============================================================

# Convert best hyperparameters to global environment
list2env(l_hyperparams_best, envir = .GlobalEnv)

# Get output layer activation type
df_fn_grps <- df_fn_grps %>%
  mutate(stan_activation_num = case_when(
    activation %in% "linear" ~ 0,
    activation %in% "sigmoid" ~ 1,
    activation %in% "exponential" ~ 2,
    activation %in% "softmax" ~ 3,
    .default = NA)
  )

stan_activation_num <- df_fn_grps$stan_activation_num
if (sum(is.na(stan_activation_num) > 0)) {
  stop("Activation function not recognized")
}

# Set activation function for Stan
if (activation_fun == "relu") {
  hidden_activation_num = 1
} else if (activation_fun == "tanh") {
  hidden_activation_num = 2
} else {
  hidden_activation_num = 3
}

# Get ANN weights
weights <- get_weights(model) 

# Pass the weights and biases to Stan for Bayesian calibration
weight_first <- weights[[1]]
beta_first <- 1 %*% weights[[2]]

weight_middle <- array(0, c(n_hidden_layers-1, n_hidden_nodes, n_hidden_nodes))
beta_middle <- array(0, c(n_hidden_layers-1, 1, n_hidden_nodes))
for (l in 1:(n_hidden_layers-1)){ # Hidden layers
  weight_middle[l,,] <- weights[[l*2+1]]
  beta_middle[l,,] <- weights[[l*2+2]]
}

weight_last <- array(0, c(n_hidden_nodes, length(y_targets)))
beta_last <- array(0, c(1, length(y_targets)))
v_group_sizes <- c()
idx_start <- 1
for (l in n_hidden_layers:(length(weights)/2-1)){ # Output layers
  group_size <- length(weights[[l*2+2]])
  v_group_sizes <- c(v_group_sizes, group_size) # Count number of units for group size
  idx_end <- idx_start+group_size-1
  weight_last[,idx_start:idx_end] <- weights[[l*2+1]]
  beta_last[,idx_start:idx_end] <- weights[[l*2+2]]
  idx_start <- idx_end+1
}

stan.dat=list(
  num_hidden_nodes = n_hidden_nodes,
  num_hidden_layers = n_hidden_layers,
  num_inputs = n_inputs,
  num_outputs = n_outputs,
  num_targets = 1,
  num_groups = length(v_group_sizes),
  group_sizes = v_group_sizes,
  hidden_activation = hidden_activation_num,
  activation_type = stan_activation_num,
  y_targets = y_targets,
  y_targets_se = y_targets_se,
  beta_first = beta_first,
  beta_middle = beta_middle,
  beta_last = beta_last,
  weight_first = weight_first,
  weight_middle = weight_middle,
  weight_last = weight_last)

# Verify that the Stan file produces the correct data
if (validate_stan_ANN) {
  # Load ANN outputs for validating STAN
  ytest_scaled_pred <- readRDS(file_test_outputs)
  
  # Create R representation of Stan functions
  expose_stan_functions(stanmodel = file_stan)
  
  # Calculate ANN predictions with Stan functions
  stan.dat_test <- stan.dat[c("beta_first", "beta_middle", "beta_last", "weight_first",
                              "weight_middle", "weight_last", "num_hidden_layers", "num_hidden_nodes",
                              "group_sizes", "hidden_activation", "activation_type")]
  stan.dat_test <- c(X = list(t(xtest_scaled[1,])), stan.dat_test)
  stan.dat_test$beta_middle <- lapply(1:nrow(stan.dat_test$beta_middle), FUN = function(i) t(stan.dat_test$beta_middle[i,,]))
  stan.dat_test$weight_middle <- lapply(1:nrow(stan.dat_test$weight_middle), FUN = function(i) stan.dat_test$weight_middle[i,,])
  stan_test_output <- matrix(0, nrow = nrow(xtest_scaled), ncol = n_outputs)
  for (i in 1:nrow(xtest_scaled)) {
    stan.dat_test$X <- t(xtest_scaled[i,])
    stan_test_output[i,] <- as.vector(do.call(calculate_alpha, stan.dat_test))
  }
  
  # Count errors beyond epsilon margin
  assert_that(sum(abs(stan_test_output - ytest_scaled_pred) > eps_stan) == 0)
}

# Run stan file
stan.time <- proc.time()
m <- stan(file = file_stan,
          data = stan.dat,
          iter = n_iter,
          chains = n_chains,
          thin = n_thin,
          pars = c("Xq"),
          warmup = floor(n_iter/2),   ## (cp)
          seed = seed) # For reproducibility. R's set.seed() will not work for stan
t_calibration <- proc.time() - stan.time # stan sampling time
t_calibration <- t_calibration[3] / 60
summary(m)

param_names    <- colnames(data_sim_param)

names(m)[1:n_inputs] <- param_names

saveRDS(m, file_stan_model)

m <- readRDS(file_stan_model) 
param_names    <- colnames(data_sim_param)
names(m)[1:n_inputs] <- param_names

###### 9.1 Stan Diagnose  ----

stan_trace(m, pars=param_names, inc_warmup = FALSE)

stan_plot(m, pars=param_names, point_est = "mean", show_density = TRUE, fill_color = "maroon", ncol=2)

stan_hist(m, pars=param_names, inc_warmup = FALSE)

standensity <- stan_dens(m, pars=param_names, inc_warmup = FALSE, separate_chains=TRUE)
ggsave(filename = file_fig_posterior,
       standensity,
       width = 24, height = 16)

stan_dens(m, pars=param_names, inc_warmup = FALSE, separate_chains=FALSE)

stan_ac(m, pars=param_names, inc_warmup = FALSE, separate_chains=TRUE)

stan_rhat(m,pars=param_names)          # Rhat statistic 
stan_par(m,par=param_names[1])         # Mean metrop. acceptances, sample step size
stan_ess(m,pars=param_names)           # Effective sample size / Sample size
stan_mcse(m,pars=param_names)          # Monte Carlo SE / Posterior SD
stan_diag(m,)


#### STEP 5: Rescaling posterior distribution ####

###### 9.2 Stan extraction  ----

params <- rstan::extract(m, permuted=TRUE, inc_warmup = FALSE)
lp <- params$lp__
Xq <- params$Xq
Xq_df = as.data.frame(Xq)

# Scale the posteriors
if (Scale_inputs) {
  Xq_unscaled <- unscale_data(Xq_df, vec.mins = xmins, vec.maxs = xmaxs, vec.means = xmeans, vec.sds = xsds, type = scale_type)
} else {
  Xq_unscaled <- Xq_df
}


Xq_lp <- cbind(Xq_unscaled, lp)

# Save the unscaled posterior samples
write.csv(Xq_lp,
          file = file_posterior,
          row.names = FALSE)

cal_mdl_1 <- file_posterior

### Load ANN posterior
Xq_lp <- read.csv(file = cal_mdl_1)
n_col <- ncol(Xq_lp)
lp <- Xq_lp[, n_col]
Xq_unscaled <- Xq_lp[, -n_col]
map_baycann <- Xq_unscaled[which.max(lp), ]     ### MAP for first BayCANN model
df_post_ann <- read.csv(file = cal_mdl_1)[, -n_col]
colnames(df_post_ann) <- x_names


## Correlation graph
library(GGally)

df_post <- data.frame(Xq_unscaled)
colnames(df_post) <- x_names
df_post_long <- reshape2::melt(df_post,
                               variable.name = "Parameter")

df_post_long$Parameter <- factor(df_post_long$Parameter,
                                 levels = levels(df_post_long$Parameter),
                                 ordered = TRUE)

gg_calib_post_pair_corr <- GGally::ggpairs(df_post,
                                           upper = list(continuous = wrap("cor",
                                                                          color = "black",
                                                                          size = 5)),
                                           diag = list(continuous = wrap("barDiag",
                                                                         alpha = 0.8)),
                                           lower = list(continuous = wrap("points",
                                                                          alpha = 0.3,
                                                                          size = 0.5)),
                                           labeller = "label_parsed") +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size=6),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "white",
                                        color = "white"),
        strip.text = element_text(hjust = 0))
gg_calib_post_pair_corr

ggsave(filename = file_fig_corr,
       gg_calib_post_pair_corr,
       width = 36, height = 24)


#### Prior and prior graph
n_samp <- min(1000, nrow(data_sim_param_train))
df_samp_prior <- melt(cbind(Distribution = "Prior",
                            as.data.frame(data_sim_param_train[1:n_samp, ])),
                      variable.name = "Parameter")

df_samp_post_ann   <- melt(cbind(Distribution = "Posterior BayCANN",
                                 as.data.frame(df_post_ann[1:n_samp, ])),
                           variable.name = "Parameter")


df_samp_prior_post <- rbind(df_samp_prior,
                            df_samp_post_ann)
df_samp_prior_post$Distribution <- ordered(df_samp_prior_post$Distribution,
                                           levels = c("Prior",
                                                      "Posterior BayCANN"))



df_samp_prior_post$Parameter <- factor(df_samp_prior_post$Parameter,
                                       levels = levels(df_samp_prior_post$Parameter),
                                       ordered = TRUE)

df_maps_n_true_params <- data.frame(Type = ordered(rep(c( "BayCANN MAP"), each = n_inputs),
                                                   levels = c("BayCANN MAP")),
                                    value = c(t(map_baycann)))
df_maps_n_true_params


### Plot priors and ANN posteriors
df_maps_n_true_params$Parameter<-as.factor(x_names)


library(dampack)

gg_prior_post <- ggplot(df_samp_prior_post,
                         aes(x = value, y = ..density.., fill = Distribution)) +
  facet_wrap(~Parameter, scales = "free",
             ncol = 5,
             labeller = label_parsed) +
  scale_x_continuous(breaks = number_ticks(5)) +
  scale_color_manual("", values = c("black", "navy blue", "tomato","green")) +
  geom_density(alpha=0.5) +
  theme_bw(base_size = 16) +
  guides(fill = guide_legend(title = "", order = 1),
         linetype = guide_legend(title = "", order = 2),
         color = guide_legend(title = "", order = 2)) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin=margin(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill = "white",
                                        color = "white"),
        strip.text = element_text(hjust = 0))
gg_prior_post
ggsave(gg_prior_post,
       filename = file_fig_prior_post,
       width = 36, height = 24)


# 10. Save BayCANN parameters  -------------------------------------------------

## Save BayCANN parameters
param_BayCANN <- list(train_split = train_split,
                      n_batch_size = n_batch_size,
                      n_epochs = n_epochs,
                      patience = patience,
                      l_hyperparams_best = l_hyperparams_best,
                      init_W = init_W,
                      scale_type = scale_type,
                      df_fn_grps = df_fn_grps,
                      scale_cols = scale_cols,
                      verbose = verbose,
                      Normalize_outputs = Normalize_outputs,
                      Scale_inputs = Scale_inputs,
                      Scale_outputs = Scale_outputs,
                      Remove_outliers = Remove_outliers,
                      Standardize_targets = Standardize_targets,
                      Saved_data = Saved_data,
                      Selected_targets = Selected_targets,
                      n_chains = n_chains,
                      n_iter = n_iter,
                      n_thin = n_thin,
                      file_sample = file_sample,
                      file_keras_model = file_keras_model,
                      file_test_outputs = file_test_outputs,
                      file_stan = file_stan,
                      t_training = t_training,
                      t_calibration = t_calibration,
                      acc_err = acc_err
)

saveRDS(param_BayCANN, file = file_baycann_params)

