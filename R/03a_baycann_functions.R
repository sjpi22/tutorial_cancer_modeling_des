# =================================================
# Libraries
library(MASS)
library(bestNormalize)

# functions
prepare_data <- function(xtrain, ytrain, xtest, ytest, scale_type){
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
  xresults <- scale_data(x,scale_type) 
  yresults <- scale_data(y,scale_type)  ##I think this is not adequate, because we are including the test data
  
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

scale_data  <- function(unscaled_data,type){
  vec.maxs  <- apply(unscaled_data, 2, max) 
  vec.mins  <- apply(unscaled_data, 2, min)
  vec.means <- apply(unscaled_data, 2,mean)
  vec.sds   <- apply(unscaled_data, 2,sd)

  vec.quant <- apply(unscaled_data, 2,quantile)
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
  

  
  results <- list(scaled_data = scaled_data, vec.mins = vec.mins, vec.maxs = vec.maxs,vec.means=vec.means,vec.sds=vec.sds)
  return(results)
}

unscale_data <- function(scaled_data, vec.mins, vec.maxs,vec.means,vec.sds, type){
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


# Uncomment if an specific list of SimCRC parameters are required 
#Selected targets 
# Selected_SimCRC_Targets <- c("Prev0AdOrPreClin_age27",
# "Prev0AdOrPreClin_age52",
# "Prev0AdOrPreClin_age72",
# "Prev0AdOrPreClin_age97",
# "Prev1AdOrPreClin_age27",
# "Prev1AdOrPreClin_age52",
# "Prev1AdOrPreClin_age72",
# "Prev1AdOrPreClin_age97",
# "Prev2AdOrPreClin_age27",
# "Prev2AdOrPreClin_age52",
# "Prev2AdOrPreClin_age72",
# "Prev2AdOrPreClin_age97",
# "Prev3AdOrPreClin_age27",
# "Prev3AdOrPreClin_age52",
# "Prev3AdOrPreClin_age72",
# "Prev3AdOrPreClin_age97",
# "SizeLRGivenAdInP_wtd",
# "SizeHRGivenAdInP_wtd",
# "SizeLRGivenAdInD_wtd",
# "SizeHRGivenAdInD_wtd",
# "SizeLRGivenAdInR_wtd",
# "SizeHRGivenAdInR_wtd",
# "PrevPreclinical_ages40_49",
# "PrevPreclinical_ages50_59",
# "PrevPreclinical_ages60_69",
# "PrevPreclinical_ages70_79",
# "PrevPreclinical_ages80_89",
# "CRCincPer100K_P_ages20_39",
# "CRCincPer100K_P_ages40_49",
# "CRCincPer100K_P_ages50_59",
# "CRCincPer100K_P_ages60_69",
# "CRCincPer100K_P_ages70_79",
# "CRCincPer100K_P_ages80_99",
# "CRCincPer100K_D_ages20_39",
# "CRCincPer100K_D_ages40_49",
# "CRCincPer100K_D_ages50_59",
# "CRCincPer100K_D_ages60_69",
# "CRCincPer100K_D_ages70_79",
# "CRCincPer100K_D_ages80_99",
# "CRCincPer100K_R_ages20_39",
# "CRCincPer100K_R_ages40_49",
# "CRCincPer100K_R_ages50_59",
# "CRCincPer100K_R_ages60_69",
# "CRCincPer100K_R_ages70_79",
# "CRCincPer100K_R_ages80_99",
# "CRCincPer100K_ages20_39",
# "CRCincPer100K_ages40_49",
# "CRCincPer100K_ages50_59",
# "CRCincPer100K_ages60_69",
# "CRCincPer100K_ages70_79",
# "CRCincPer100K_ages80_99",
# "StageDistribution_P_S1",
# "StageDistribution_P_S3",
# "StageDistribution_P_S4",
# "StageDistribution_D_S1",
# "StageDistribution_D_S3",
# "StageDistribution_D_S4",
# "StageDistribution_R_S1",
# "StageDistribution_R_S3",
# "StageDistribution_R_S4")
# 

# Uncomment if an specific list of MISCAN parameters are required 
# Selected_MISCAN_R_targets <- c(
#   "inc_ad_40",
#   "inc_ad_45",
#   "inc_ad_50",
#   "inc_ad_55",
#   "inc_ad_60",
#   "inc_ad_65",
#   "inc_ad_70",
#   "inc_ad_75",
#   "inc_ad_80",
#   "inc_ad_85",
#   "ad_prev_32",
#   "ad_prev_37",
#   "ad_prev_42",
#   "ad_prev_47",
#   "ad_prev_52",
#   "ad_prev_57",
#   "ad_prev_62",
#   "ad_prev_67",
#   "ad_prev_72",
#   "va_prev_37",
#   "va_prev_42",
#   "va_prev_47",
#   "va_prev_52",
#   "va_prev_57",
#   "va_prev_62",
#   "va_prev_67",
#   "va_prev_72",
#   "ad_mult_1_37",
#   "ad_mult_1_42",
#   "ad_mult_1_47",
#   "ad_mult_1_52",
#   "ad_mult_1_57",
#   "ad_mult_1_62",
#   "ad_mult_1_67",
#   "ad_mult_1_72",
#   "ad_mult_2_37",
#   "ad_mult_2_42",
#   "ad_mult_2_47",
#   "ad_mult_2_52",
#   "ad_mult_2_57",
#   "ad_mult_2_62",
#   "ad_mult_2_67",
#   "ad_mult_2_72",
#   "ad_mult_3_37",
#   "ad_mult_3_42",
#   "ad_mult_3_47",
#   "ad_mult_3_52",
#   "ad_mult_3_57",
#   "ad_mult_3_62",
#   "ad_mult_3_67",
#   "ad_mult_3_72",
#   "va_mult_1_37",
#   "va_mult_1_42",
#   "va_mult_1_47",
#   "va_mult_1_52",
#   "va_mult_1_57",
#   "va_mult_1_62",
#   "va_mult_1_67",
#   "va_mult_1_72",
#   "va_mult_2_52",
#   "va_mult_2_57",
#   "va_mult_2_62",
#   "va_mult_2_67",
#   "va_mult_3_37",
#   "va_mult_3_42",
#   "va_mult_3_47",
#   "va_mult_3_52",
#   "va_mult_3_57",
#   "va_mult_3_62",
#   "ta_size_5_37",
#   "ta_size_5_42",
#   "ta_size_5_47",
#   "ta_size_5_52",
#   "ta_size_5_57",
#   "ta_size_5_62",
#   "ta_size_5_67",
#   "ta_size_5_72",
#   "ta_size_6_32",
#   "ta_size_6_37",
#   "ta_size_6_42",
#   "ta_size_6_47",
#   "ta_size_6_52",
#   "ta_size_6_57",
#   "ta_size_6_62",
#   "ta_size_6_67",
#   "ta_size_6_72",
#   "ta_size_10_32",
#   "ta_size_10_37",
#   "ta_size_10_42",
#   "ta_size_10_47",
#   "ta_size_10_52",
#   "ta_size_10_57",
#   "ta_size_10_62",
#   "ta_size_10_67",
#   "ta_size_10_72",
#   "va_size_5_52",
#   "va_size_5_57",
#   "va_size_5_62",
#   "va_size_5_67",
#   "va_size_6_52",
#   "va_size_6_57",
#   "va_size_6_62",
#   "va_size_6_67",
#   "va_size_6_72",
#   "va_size_10_37",
#   "va_size_10_47",
#   "va_size_10_52",
#   "va_size_10_57",
#   "va_size_10_62",
#   "va_size_10_67",
#   "va_size_10_72",
#   "ta_va_42",
#   "ta_va_47",
#   "ta_va_52",
#   "ta_va_57",
#   "ta_va_62",
#   "ta_va_67",
#   "ta_va_72",
#   "va10_hgd_52",
#   "va10_hgd_57",
#   "va10_hgd_67",
#   "va10_hgd_72"
# )



