functions {
 matrix calculate_alpha(matrix X, matrix beta_first, matrix[] beta_middle, matrix beta_last, matrix weight_first, matrix[] weight_middle, matrix weight_last, int num_hidden_layers, int num_hidden_nodes, int[] group_sizes, int hidden_activation, int[] activation_type){
    matrix[rows(X), cols(beta_first)] layer_values_first;
    matrix[rows(X), cols(beta_first)] layer_values_middle[num_hidden_layers];
    matrix[rows(X), cols(beta_last)] alpha;
    // Forward pass through hidden layers
    if (hidden_activation == 1) {
      layer_values_first = fmax(beta_first + X * weight_first, 0);  
      layer_values_middle[1] = fmax(beta_middle[1] + layer_values_first * weight_middle[1], 0);
      for(i in 2:(num_hidden_layers-1)){
        layer_values_middle[i] = fmax(beta_middle[i] + layer_values_middle[i-1] * weight_middle[i], 0);
      }
    } else if (hidden_activation == 2) {
      layer_values_first = tanh(beta_first + X * weight_first);  
      layer_values_middle[1] = tanh(beta_middle[1] + layer_values_first * weight_middle[1]);
      for(i in 2:(num_hidden_layers-1)){
        layer_values_middle[i] = tanh(beta_middle[i] + layer_values_middle[i-1] * weight_middle[i]);
      }
    } else {
      layer_values_first = inv_logit(beta_first + X * weight_first);  
      layer_values_middle[1] = inv_logit(beta_middle[1] + layer_values_first * weight_middle[1]);
      for(i in 2:(num_hidden_layers-1)){
        layer_values_middle[i] = inv_logit(beta_middle[i] + layer_values_middle[i-1] * weight_middle[i]);
      }
    }
    
    // Calculate initial alpha value
    alpha = beta_last + layer_values_middle[num_hidden_layers-1] * weight_last;
    
    // Apply activation functions
    int idx_start = 1;
    for (k in 1:size(group_sizes)) {
      // Get group size and index end
      int group_size = group_sizes[k];
      int idx_end = idx_start + group_size - 1;

      // Apply activation functions based on specified type
      if (activation_type[k] == 1) { // Sigmoid
        alpha[:, idx_start:idx_end] = inv_logit(alpha[:, idx_start:idx_end]);
      } else if (activation_type[k] == 2) {  // Exponential
        alpha[:, idx_start:idx_end] = exp(alpha[:, idx_start:idx_end]);
      } else if (activation_type[k] == 3) {  // Softmax
        // Apply softmax to the softmax groups
        // Loop over each row and apply softmax
        for (n in 1:rows(X)) {
          alpha[n, idx_start:idx_end] = to_row_vector(softmax(to_vector(alpha[n, idx_start:idx_end])));  // Convert row slice to vector and apply softmax
        }
      }
      idx_start = idx_end + 1;
    }
    
    return alpha;
  }
}

data {
  int<lower=0> num_targets;
  int<lower=0> num_inputs;
  int<lower=0> num_outputs;
  int<lower=0> num_hidden_nodes;
  int<lower=1> num_hidden_layers;
  int<lower=1> num_groups;
  int<lower=0> group_sizes[num_groups];
  int hidden_activation; // 1 for ReLU, 2 for tanh
  int activation_type[num_groups]; // 0 for none/linear, 1 for sigmoid, 2 for exponential, 3 for softmax
  matrix[num_targets,num_outputs] y_targets;
  matrix[num_targets,num_outputs] y_targets_se;
  matrix[num_inputs, num_hidden_nodes] weight_first;
  matrix[num_hidden_nodes, num_hidden_nodes] weight_middle[num_hidden_layers-1];
  matrix[num_hidden_nodes, num_outputs] weight_last;
  matrix[1, num_hidden_nodes] beta_first;
  matrix[1, num_hidden_nodes] beta_middle[num_hidden_layers-1];
  matrix[1, num_outputs] beta_last;
}

parameters {
  matrix<lower=0, upper=1>[num_targets,num_inputs] Xq;
}

model{
  matrix[1, num_outputs] alpha_post;
    alpha_post = calculate_alpha(Xq, beta_first, beta_middle, beta_last, weight_first, weight_middle, weight_last, num_hidden_layers, num_hidden_nodes, group_sizes, hidden_activation, activation_type);
    to_vector(y_targets) ~ normal(to_vector(alpha_post),to_vector(y_targets_se)); //get SE directly from data
    to_vector(Xq) ~ uniform(0,1); // assume uniform prior  
}


