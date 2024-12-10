# =================================================
# functions
prepare_data <- function(xtrain, ytrain, xtest, ytest, scale_type, scale_cols = NULL){
  y_names <- colnames(ytrain)
  x_names <- colnames(xtrain)
  n_train <- nrow(xtrain)
  n_test <- nrow(xtest)
  x <- rbind(xtrain, xtest)
  y <- rbind(ytrain, ytest)
  n <- nrow(x)
  n_inputs <- length(x_names)
  n_outputs <- length(y_names)
  
  # scale the PSA inputs and outputs
  if (length(scale_type)==1) {
    xresults <- scale_data(x, scale_type) 
    yresults <- scale_data(y, scale_type, scale_cols)  ##I think this is not adequate, because we are including the test data
  } else {
    xresults <- scale_data(x, scale_type[1]) 
    yresults <- scale_data(y, scale_type[2], scale_cols)
  }
  
  xscaled <- xresults$scaled_data 
  yscaled <- yresults$scaled_data 
  xmins <- xresults$vec.mins
  xmaxs <- xresults$vec.maxs
  xmeans<- xresults$vec.means  # for standardizing
  xsds  <- xresults$vec.sds    # for standardizing
  
  ymins <- yresults$vec.mins
  ymaxs <- yresults$vec.maxs
  ymeans<- yresults$vec.means
  ysds  <- yresults$vec.sds
  
  xtrain_scaled <- xscaled[1:n_train, ]
  ytrain_scaled <- yscaled[1:n_train, ]
  xtest_scaled  <- xscaled[(n_train+1):n, ]
  ytest_scaled  <- yscaled[(n_train+1):n, ]
  
  return(list(n_inputs = n_inputs,
              n_outputs = n_outputs,
              n_train = n_train,
              n_test = n_test,
              x_names = x_names, 
              y_names = y_names,
              xscaled = xscaled,
              yscaled = yscaled,
              xtrain_scaled = xtrain_scaled,
              ytrain_scaled = ytrain_scaled,
              xtest_scaled  = xtest_scaled ,
              ytest_scaled  = ytest_scaled,
              xmins = xmins,
              xmaxs = xmaxs,
              xmeans= xmeans,
              xsds  = xsds,
              ymins = ymins,
              ymaxs = ymaxs,
              ymeans= ymeans,
              ysds  = ysds
  ))
}

scale_data  <- function(unscaled_data, type, scale_cols=NULL){
  vec.maxs  <- apply(unscaled_data, 2, max) 
  vec.mins  <- apply(unscaled_data, 2, min)
  vec.means <- apply(unscaled_data, 2, mean)
  vec.sds   <- apply(unscaled_data, 2, sd)
  
  vec.quant <- apply(unscaled_data, 2, quantile)
  vec.q50   <- vec.quant[3,]
  vec.q25   <- vec.quant[2,]
  vec.q75   <- vec.quant[4,]
  vec.dist  <- apply(unscaled_data, 2, norm_vector)
  
  vec.ones  <- matrix(1, nrow = nrow(unscaled_data), 1)
  mat.maxs  <- vec.ones %*% vec.maxs
  mat.mins  <- vec.ones %*% vec.mins
  mat.means <- vec.ones %*% vec.means
  mat.sds   <- vec.ones %*% vec.sds
  mat.q50   <- vec.ones %*% vec.q50
  mat.q25   <- vec.ones %*% vec.q25
  mat.q75   <- vec.ones %*% vec.q75
  mat.dist  <- vec.ones %*% vec.dist
  
  
  if (type==1) {
    scaled_data <- 2 * (unscaled_data - mat.mins) / (mat.maxs - mat.mins) - 1   ###  range from -1 to 1
  }
  
  if (type==2) {
    scaled_data <- (unscaled_data - mat.means) / (mat.sds)                      ###  standardized
  }
  
  if (type==3) {
    scaled_data <- (unscaled_data - mat.mins) / (mat.maxs - mat.mins)           ### range from 0 to 1
  }
  
  if (!is.null(scale_cols)) { # Applies scaling only to the columns specified in scale_cols
    scaled_data_final <- unscaled_data
    scaled_data_final[, scale_cols] <- scaled_data[, scale_cols]
    scaled_data <- scaled_data_final
  }
  
  results <- list(scaled_data = scaled_data, vec.mins = vec.mins, vec.maxs = vec.maxs,vec.means=vec.means,vec.sds=vec.sds)
  return(results)
}

unscale_data <- function(scaled_data, vec.mins, vec.maxs,vec.means,vec.sds, type, scale_cols=NULL){
  vec.ones <- matrix(1, nrow = nrow(scaled_data), 1)
  mat.mins <- vec.ones %*% vec.mins
  mat.maxs <- vec.ones %*% vec.maxs
  mat.means<- vec.ones %*% vec.means
  mat.sds  <- vec.ones %*% vec.sds
  
  
  if (type==1) {
    unscaled_data <- (scaled_data + 1) * (mat.maxs - mat.mins) / 2 + mat.mins  ###  range from -1 to 1
  }
  
  if (type==2) {
    unscaled_data <- (scaled_data * mat.sds) + mat.means                     ###  standardized
  } 
  
  if (type==3) {
    unscaled_data <- scaled_data * (mat.maxs - mat.mins) + mat.mins          ### range from 0 to 1
  }
  
  if (!is.null(scale_cols)) {
    unscaled_data_final <- scaled_data
    unscaled_data_final[, scale_cols] <- unscaled_data[, scale_cols]
    unscaled_data <- unscaled_data_final
  }
  
  return(unscale_data= unscaled_data)
}


number_ticks <- function(n) {
  function(limits) {
    pretty(limits, n + 1)
  }
}

norm_vector <- function(x) {
  sqrt(sum(x^2))   ##Euclidean distance of vector
}



best_normal <- function(vec) {
  BN_obj  <- bestNormalize(vec)
  vec     <- predict(BN_obj)
  results <- list(vec = vec, BN_obj = BN_obj)
  return  (results)
}

best_normal_dataset <- function(data) { 
  
  inverse_dist <- vector(mode = "list", length = dim(data)[2])       #Variable that will contain the inverse transformation
  data_normal  <- data
  for (i in 1:dim(data)[2]) {
    norm_info <-best_normal(data[,i])
    data_normal [,i] <- norm_info[["vec"]]
    inverse_dist[[i]] <- norm_info[["BN_obj"]]
    print(i)
    vec_comp <- predict(inverse_dist[[i]], newdata = data_normal[,i], inverse = TRUE)
    print(all.equal(vec_comp, data[,i]))
  }
  results <- list(data_normal=data_normal, inverse_dist=inverse_dist)
  return(results)
}


#function to identify rows with outliers to be removed
outlier_vector <-function(df) {
  
  require(outliers)
  num_targets <- dim(df)[2]  
  num_obs     <- dim(df)[1]  
  Out_mat     <- matrix(nrow = num_obs, ncol=num_targets)
  
  for (i in 1:num_targets) {
    
    vec<-outlier(df[,i], logical = TRUE) 
    Out_mat[,i] <- vec
  }
  o_vector<-apply(Out_mat,1,max)
  return(o_vector)
  
}

quantile_vector <- function (df, prob)  {      #function to get vector for values over te quantile probs
  
  num_targets      <- dim(df)[2]  
  num_obs          <- dim(df)[1]  
  quantile_mat     <- matrix(nrow = num_obs, ncol=num_targets)
  
  for (i in 1:num_targets) {
    
    q<-quantile(df[,i], probs=prob) 
    quantile_mat[,i] <- (df[,i]>q)
  }
  q_vector<-apply(quantile_mat,1,max)
  return(q_vector)
}

normalize_predict <- function(df,list_transf,inverse= FALSE)  {
  
  if (is.matrix(df) | is.data.frame(df)) {
    lim=dim(df)[2]
    for (i in seq(1, lim)){
      df[,i] <- predict(list_transf[[i]], newdata = df[,i], inverse = inverse )
    }
  }
  if (is.vector(df)) {
    lim=length(df)  
    for (i in seq(1, lim)){
      df[i] <- predict(list_transf[[i]], newdata = df[i], inverse = inverse )
    }
  }
  
  return(df)
}