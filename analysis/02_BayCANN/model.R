# library(keras3)

# Hyperparameter flags ---------------------------------------------------

FLAGS <- flags(
  flag_integer("n_hidden_layers", 2),
  flag_integer("n_hidden_nodes", 64),
  flag_boolean("dropout", TRUE),
  flag_numeric("dropout_rate", 0.3),
  flag_string("activation_fun", "relu")
)

# Data Preparation ---------------------------------------------------

# --- Assumes data already loaded in global environment


# Define Model --------------------------------------------------------------

# Input layer
input <- layer_input(shape = ncol(data_sim_param_train))

# Add hidden layers and dropout if specified
for (i in 1:FLAGS$n_hidden_layers) {
  if (i==1) {
    x <- layer_dense(input, units = FLAGS$n_hidden_nodes, activation = FLAGS$activation_fun)
  } else {
    x <- layer_dense(x, units = FLAGS$n_hidden_nodes, activation = FLAGS$activation_fun)
  }
  if (FLAGS$dropout) {
    x <- layer_dropout(x, rate = FLAGS$dropout_rate)
  }
}

# Create output branches and loss
output <- list()
loss_list <- list()
loss_weights <- list()
metrics_list <- list()
for (i in 1:nrow(df_fn_grps)) {
  output[[i]] <- layer_dense(x, 
                             units = max(1, ncol(ytrain_scaled_reshape[[i]])), 
                             activation = df_fn_grps$activation[i])
  loss_list[[i]] <- df_fn_grps$loss_fn[i]
  loss_weights[[i]] <- df_fn_grps$loss_weight[i]
  metrics_list[[i]] <- df_fn_grps$metric[i]
}

# Define the model
model <- keras_model(inputs = input, outputs = output)

# View the model summary
summary(model)

# Compile the model with both MSE and Categorical Cross-Entropy losses
# loss_list <- list(loss_mean_squared_error(), loss_mean_squared_error())
# for (i in unique(data_true_targets$fn_grp[data_true_targets$fn_grp != 0])) {
#   loss_list <- c(loss_list, loss_categorical_crossentropy())
# }

model %>% compile(
  loss = loss_list,
  loss_weights = loss_weights, 
  optimizer = optimizer_adam(),
  metrics = metrics_list
)


# Training & Evaluation ----------------------------------------------------

history <- model %>% fit(
  x = xtrain_scaled, 
  y = ytrain_scaled_reshape,
  epochs = n_epochs,
  batch_size = n_batch_size,
  validation_split = 0.2,
  validation_data = l_validation_data,
  verbose = verbose,
  callbacks = callback_early_stopping(
    monitor = "val_loss",
    patience = patience,
    verbose = 0,
    restore_best_weights = TRUE
  )
)
