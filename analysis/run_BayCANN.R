###########################  BayCANN   #########################################
#
#  Objective: Script to perform an emulator-based Bayesian calibration on SimCRC
########################### <<<<<>>>>> ##############################################

# Sources: Jalal H, Trikalinos TA, Alarid-Escudero F. BayCANN: Streamlining 
# Bayesian Calibration With Artificial Neural Network Metamodeling. Front 
# Physiol. 2021 May 25;12:662314. doi: 10.3389/fphys.2021.662314. PMID: 
#   34113262; PMCID: PMC8185956.

# Pineda-Antunez C, Seguin C, van Duuren LA, Knudsen AB, Davidi B, de Lima PN, 
# Rutter C, Kuntz KM, Lansdorp-Vogelaar I, Collier N, Ozik J, Alarid-Escudero F. 
# Emulator-based Bayesian calibration of the CISNET colorectal cancer models. 
# medRxiv [Preprint]. 2024 Feb 5:2023.02.27.23286525. doi: 
#   10.1101/2023.02.27.23286525. PMID: 36909607; PMCID: PMC10002763.


#### 1.Libraries and functions  ==================================================

library(keras)   #Install previously tensorflow
library(rstan)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(doParallel)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


###### 1.1 Load baycann functions =================================================
#* Clean environment
rm(list = ls())
source("R/03_calibration_general_functions.R")
source("R/03a_baycann_functions.R")

#### 2. General parameters ========================================================

###### 2.1 parameters for Data preparation 
scale_type <- 1  ## 1: for scale from -1 to 1; 2: for standardization ; 3: for scale from 0 to 1
seed <- 123
train_split <- 0.8
sample_file <- "data/calibration_sample.RData"
targets_files <- list(prevalence = list(
  a = "data/prevalence_lesion_a.csv",
  b = "data/prevalence_lesion_b.csv"),
  incidence = "data/incidence_cancer.csv",
  stage_distr = "data/stage_distr.csv")

###### 2.1 parameters for ANN 
verbose          <- 0
n_batch_size     <- 2000
n_epochs         <- 15000
patience         <- 10000
n_hidden_nodes   <- 360
n_hidden_layers  <- 4
activation_fun   <- 'tanh'
init_W = NULL

# Other options to initialize weights
#init_W=initializer_random_uniform(minval = -0.7, maxval = 0.7,seed = 2312)   ###initialization of weights with uniform distribution
#init_W=initializer_random_normal(mean = 0, stddev = 0.1, seed = 2312)  ###initialization of weights with normal distribution


###### 2.2 parameters for Bayesian calibration
n_iter   <- 300000
n_thin   <- 100
n_chains <- 4

#### 3. Pre-processing actions  ===========================================

set.seed(seed)

Normalize_inputs    <- FALSE  # TRUE if we want to normalize inputs
Normalize_outputs   <- FALSE  # TRUE if we want to normalize outputs 
Scale_inputs        <- TRUE   # TRUE if we want to scale inputs
Scale_outputs       <- TRUE   # TRUE if we want to scale outputs 
Remove_outliers     <- FALSE  # TRUE if we want to remove outliers
Standardize_targets <- FALSE  # TRUE if we want to standardize targets
Saved_data          <- FALSE  # TRUE if we want to saved the data
Selected_targets    <- FALSE  # TRUE if we want to use an specific list of calibration targets


#### 4. Load the training and test data for the simulations =======================

load(sample_file)

###### 4.1 Input Parameters ####

data_sim_param <- as.matrix(m_param_samp)

###### 4.2 Model Outputs ####

data_sim_target <- as.matrix(out_calib_targets)

#### 5. Removing outliers from model outputs ####

if (Remove_outliers) {
  vec_out <- outlier_vector(data_sim_target)
  data_sim_param  <- data_sim_param[!vec_out,]
  data_sim_target <- data_sim_target[!vec_out,]
}

#### STEP 2 in paper: Splitting and rescaling data ####

#### 7. Train/test partition ======================================================

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
                              scale  = scale_type)

list2env(prepared_data, envir = .GlobalEnv)


####  9. Load the targets and their se ============================================

l_true_targets <- recursive_read_csv(targets_files)
data_true_targets <- reshape_calib_targets(l_true_targets, output_se = TRUE)

true_targets_mean  <- data_true_targets$v_targets
true_targets_se   <- data_true_targets$v_se


#### 10. Scale the targets and their SE  ####

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
  y_targets_se <-(true_targets_se) / (ymaxs - ymins)
}

y_targets <- t(as.matrix(y_targets))
y_targets_se <- t(as.matrix(y_targets_se))   

#### STEP 3 in paper: Artificial Neural Network ####

#### 11. Keras Section BayCANN SimCRC ==============================================

# File name of keras model
path_keras_model <- paste0("output/model_keras_BayCANN.h5")    ##File path for the compiled model

model <- keras_model_sequential()
mdl_string <- paste("model %>% layer_dense(units = n_hidden_nodes, kernel_initializer=init_W, activation = activation_fun, input_shape = n_inputs) %>%",
                    paste(rep(x = "layer_dense(units = n_hidden_nodes, activation = activation_fun) %>%",
                              (n_hidden_layers)), collapse = " "),
                    "layer_dense(units = n_outputs)")

eval(parse(text = mdl_string))

summary(model)

model %>% compile(
  loss = 'mean_squared_error',
  optimizer = 'adam'  ,
  metrics = list('mae',"accuracy")
)

keras.time <- proc.time()

history <- model %>% fit(
  xtrain_scaled, ytrain_scaled,
  epochs = n_epochs,
  batch_size = n_batch_size,
  validation_data = list(xtest_scaled, ytest_scaled),
  verbose = verbose,
  callback_early_stopping(
    monitor = "val_loss",
    patience = patience,
    verbose = 0,
    restore_best_weights = TRUE
  )
)

t_training <- proc.time() - keras.time #keras ann fitting time

t_training <- t_training[3]/60

acc_err<-model%>%evaluate(xtest_scaled,ytest_scaled) # Model performance evaluation
metric_loss      <- acc_err[1]
metric_mae       <- acc_err[2]
metric_accuracy  <- acc_err[3]

# save_model_hdf5(model,path_keras_model)  #Save the model
save_model_tf(model, path_keras_model)
model <- load_model_hdf5(path_keras_model)


###### 11.4 History Graph ####

plot(history)   #Plot loss function and accuracy function

###### 11.5 Prediction Graph  ####

pred <- model %>% predict(xtest_scaled)
ytest_scaled_pred <- data.frame(pred)
colnames(ytest_scaled_pred) <- y_names
head(ytest_scaled_pred)    #

ann_valid <- rbind(data.frame(sim = 1:n_test, ytest_scaled, type = "model"),
                   data.frame(sim = 1:n_test, ytest_scaled_pred, type = "pred"))
ann_valid_transpose <- ann_valid %>%
  pivot_longer(cols = -c(sim, type)) %>%
  pivot_wider(id_cols = c(sim, name), names_from = type, values_from = value)

##Partition of validation data for Graph (4 parts)
n_partition <-4
n_part_bach <-floor(n_outputs/n_partition)

ann_valid_transpose <- arrange(ann_valid_transpose,desc(name))

ann_valid_transpose1 <- ann_valid_transpose[(1):(n_part_bach*n_test),]
ann_valid_transpose2 <- ann_valid_transpose[(n_part_bach*n_test+1):(2*n_part_bach*n_test),]
ann_valid_transpose3 <- ann_valid_transpose[(2*n_part_bach*n_test+1):(3*n_part_bach*n_test),]
ann_valid_transpose4 <- ann_valid_transpose[(3*n_part_bach*n_test+1):dim(ann_valid_transpose)[1],]

#part 1
ggplot(data = ann_valid_transpose1, aes(x = model, y = pred)) +
  geom_point(alpha = 0.5, color = "tomato") +
  facet_wrap(~name, scales="free",  ncol = 7) +
  xlab("Model outputs (scaled)") +
  ylab("ANN predictions (scaled)") +
  #coord_equal() +
  theme_bw()

#part 2
ggplot(data = ann_valid_transpose2, aes(x = model, y = pred)) +
  geom_point(alpha = 0.5, color = "tomato") +
  facet_wrap(~name, scales="free", ncol = 7) +
  xlab("Model outputs (scaled)") +
  ylab("ANN predictions (scaled)") +
  #coord_equal() +
  theme_bw()

#part 3
ggplot(data = ann_valid_transpose3, aes(x = model, y = pred)) +
  geom_point(alpha = 0.5, color = "tomato") +
  facet_wrap(~name, scales="free", ncol = 7) +
  xlab("Model outputs (scaled)") +
  ylab("ANN predictions (scaled)") +
  #coord_equal() +
  theme_bw()

#part 4
ggplot(data = ann_valid_transpose4, aes(x = model, y = pred)) +
  geom_point(alpha = 0.5, color = "tomato") +
  facet_wrap(~name, scales="free", ncol = 7) +
  xlab("Model outputs (scaled)") +
  ylab("ANN predictions (scaled)") +
  #coord_equal() +
  theme_bw()

#### STEP 4 in paper: Bayesian calibration ####

#### 12. Stan section ==============================================================

path_posterior <- "output/calibrated_posteriors_BayCANN.csv"

weights <- get_weights(model) #get ANN weights

n_hidden_layers <- length(weights)/2-2    #Removing bias layers and input and output layers
n_hidden_nodes  <- ncol(weights[[1]])     #Get number of hidden nodes from the firs layers

# pass the weights and biases to Stan for Bayesian calibration
n_layers <- length(weights)
weight_first <- weights[[1]]
beta_first <- 1 %*% weights[[2]]
weight_last <- weights[[n_layers-1]]
beta_last <- 1 %*% weights[[n_layers]]
weight_middle <- array(0, c(n_hidden_layers, n_hidden_nodes, n_hidden_nodes))
beta_middle <- array(0, c(n_hidden_layers, 1, n_hidden_nodes))
for (l in 1:n_hidden_layers){
  weight_middle[l,,] <- weights[[l*2+1]]
  beta_middle[l,,] <- weights[[l*2+2]]
}

stan.dat=list(
  num_hidden_nodes = n_hidden_nodes,
  num_hidden_layers= n_hidden_layers,
  num_inputs=n_inputs,
  num_outputs=n_outputs,
  num_targets=1,
  y_targets = y_targets,
  y_targets_se = y_targets_se,
  beta_first = beta_first,
  beta_middle = beta_middle,
  beta_last = beta_last,
  weight_first = weight_first,
  weight_middle = weight_middle,
  weight_last = weight_last)

# Select the stan file based on data transformation

if (Normalize_inputs) {
  file_perceptron <- "stan/post_multi_perceptron_normal.stan"
} else {
  file_perceptron <- "stan/post_multi_perceptron.stan"  
}

# Run stan file
stan.time <- proc.time()
m <- stan(file = file_perceptron,
          data = stan.dat,
          iter = n_iter,
          chains = n_chains,
          thin = n_thin,
          pars = c("Xq"),
          warmup = floor(n_iter/2),   ## (cp)
          seed = seed) #for reproducibility. R's set.seed() will not work for stan
t_calibration <- proc.time() - stan.time # stan sampling time
t_calibration <- t_calibration[3] / 60
summary(m)

path_stan_model <- paste0("output/model_stan_BayCANN.rds")

param_names    <- colnames(data_sim_param)

names(m)[1:n_inputs] <- param_names

saveRDS(m, path_stan_model)

m <- readRDS(path_stan_model) 
param_names    <- colnames(data_sim_param)
names(m)[1:n_inputs] <- param_names

###### 12.1 Stan Diagnose  ----

stan_trace(m,pars=param_names,inc_warmup = FALSE)

stan_plot(m,pars=param_names, point_est = "mean", show_density = TRUE, fill_color = "maroon", ncol=2)

stan_hist(m,pars=param_names, inc_warmup = FALSE)

standensity <- stan_dens(m,pars=param_names, inc_warmup = FALSE, separate_chains=TRUE)
ggsave(filename = paste0("output/fig_posterior_distribution_chains_BayCANN.png"),
       standensity,
       width = 24, height = 16)

stan_dens(m,pars=param_names, inc_warmup = FALSE, separate_chains=FALSE)

stan_ac(m,pars=param_names[1:15], inc_warmup = FALSE, separate_chains=TRUE)
stan_ac(m,pars=param_names[16:30], inc_warmup = FALSE, separate_chains=TRUE)

stan_rhat(m,pars=param_names)          # Rhat statistic 
stan_par(m,par=param_names[1])         # Mean metrop. acceptances, sample step size
stan_ess(m,pars=param_names)           # Effective sample size / Sample size
stan_mcse(m,pars=param_names)          # Monte Carlo SE / Posterior SD
stan_diag(m,)


#### STEP 5 in paper: Rescaling posterior distribution ####

###### 12.2 Stan extraction  ----

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
          file = path_posterior,
          row.names = FALSE)

cal_mdl_1 <- path_posterior
### Load ANN posterior
Xq_lp <- read.csv(file = cal_mdl_1)
n_col <- ncol(Xq_lp)
lp <- Xq_lp[, n_col]
Xq_unscaled <- Xq_lp[, -n_col]
map_baycann <- Xq_unscaled[which.max(lp), ]     ### MAP for first BayCANN model
df_post_ann <- read.csv(file = cal_mdl_1)[, -n_col]
colnames(df_post_ann) <- x_names


##correlation graph
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

ggsave(filename = paste0("output/fig_posterior_distribution_pairwise_corr_BayCANN_version.png"),
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
       filename = paste0("output/fig_BayCANN-prior-posterior.pdf"),
       width = 36, height = 24)
ggsave(gg_prior_post,
       filename = paste0("output/fig_BayCANN-prior-posterior.png"),
       width = 36, height = 24)


# 14. Save BayCANN parameters  -------------------------------------------------

##Save BayCANN parameters
path_baycann_params <- paste0("output/parameters_BayCANN.RData")
param_BayCANN <- list(scale_type,
                      scale_type, 
                      verbose,
                      n_batch_size,
                      n_chains,
                      n_epochs,
                      patience,
                      n_hidden_nodes,
                      n_hidden_layers,
                      activation_fun,
                      init_W,
                      n_iter,
                      n_thin,
                      Normalize_inputs,
                      Normalize_outputs,
                      Scale_inputs,
                      Scale_outputs,
                      Remove_outliers, 
                      Standardize_targets,
                      Saved_data,
                      Selected_targets,
                      sample_file,
                      targets_files,
                      path_keras_model,
                      t_training,
                      metric_loss,
                      metric_mae,
                      metric_accuracy, 
                      path_posterior, 
                      file_perceptron, 
                      t_calibration,
                      path_stan_model,
                      path_baycann_params
)

save(param_BayCANN, file = path_baycann_params)

