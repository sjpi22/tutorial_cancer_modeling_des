#' Load default parameters
#'
#' \code{load_default_params} loads default parameters for the decision model 
#' and creates a list
#' 
#' @return 
#' A list of all parameters used for the decision model
#' 
#' @export
load_default_params <- function(file.mort = "data/background_mortality.xlsx",
                                file.surv = "data/relative_survival_cancer.csv"){
  
  #### General setup ####
  # Initialize list to store all parameters
  l_params_all <- list(
    seed        = 123,                # Random seed
    n_cohort    = 100000,             # Cohort size
    conf_level  = 0.95,               # Confidence level
    v_states    = c(0, 1, 2, 3, "D"), # Health states
    v_lesions   = c("a", "b"),        # Precancerous lesion types
    v_cancer    = c("i", "ii", "iii", "iv"), # Cancer stages
    v_D         = c("o", "c"),        # Death states
    v_strats    = c("None", "Screen1", "Screen2")  # CEA strategies
  )
  
  #### Load input data ####
  l_lifetables <- load_from_excel(file.mort) # Background mortality
  
  # Disease-specific relative survival
  if (!is.null(file.surv)) {
    surv_data <- read.csv(file = file.surv) %>%
      mutate(pct_died = 1 - relative_surv)
  }
  
  # Strategy parameters (sensitivity, treatment effect) @@@
  # Screening @@@
  
  #### Create useful variables ####
  # Maximum age
  max_age <- max(l_lifetables$female$age, l_lifetables$male$age) + 1
  
  # Precancerous lesion starter variables
  l_distr_lesions <- list(
    time_0_1 = list(distr = "weibull", params = list(shape = 2, scale = 75)),
    n_add = list(distr = "pois", params = list(lambda = 0.05)),
    time_1_1j = list(distr = "unif", params = "list(max = (time_0_Do - time_0_1))"),
    time_1j_2i = list(distr = "gamma", params = list(shape = 8, scale = 8)),
    size_1j_2 = list(distr = "unif", params = list(min = 10, max = 40))
  )
  
  # Cancer starter variable
  time_2i_2ii <- list(distr = "gamma", params = list(shape = 2, scale = 2))
  
  # Survival after cancer diagnosis
  if (!is.null(file.surv)) {
    l_distr_surv <- list()
    for (stg in sort(unique(surv_data$stage))) {
      temp_surv_data <- surv_data[surv_data$stage == stg, ]
      l_distr_surv[[stg]] <- list(distr = "empirical", params = list(xs = temp_surv_data$years_from_dx[-length(temp_surv_data$years_from_dx)], probs = pmax(diff(temp_surv_data$pct_died), 0), max_x = max_age))
    }
  }
  
  #### Update parameter list with distributions ####
  l_params_all <- within(l_params_all, {
    # Max age
    max_age <- max_age
    
    # Sex
    b_male <- list(distr = "binom", params = list(size = 1, prob = 0.5))
    
    # Time to death from other causes
    time_0_Do_male <- list(distr = "empirical", params = list(xs = l_lifetables$male$age, probs = l_lifetables$male$p_death , max_x = max_age))
    time_0_Do_female <- list(distr = "empirical", params = list(xs = l_lifetables$female$age, probs = l_lifetables$female$p_death , max_x = max_age))
    
    # Precancerous lesions
    for(lesiontype in v_lesions) {
      for(distr in names(l_distr_lesions)) {
        assign(paste0(distr, lesiontype), l_distr_lesions[[distr]])
      }
    }
    
    # Cancer stage progression
    for(i in seq(length(v_cancer)-1)) {
      assign(paste0("time_2", v_cancer[i], "_2", v_cancer[i+1]), time_2i_2ii)
    }
    
    # Save cancer stage progression variables
    vars_cancer <- paste0("time_2", v_cancer, "_2", lead(v_cancer))[-length(v_cancer)]
    
    # Stage at diagnosis
    stage_dx <- list(distr = "empirical", params = list(xs = seq(length(v_cancer)), probs = rep(1/length(v_cancer), length(v_cancer)), continuity_correction = NULL))
    
    # Time of diagnosis within stage
    p_time_2x_3 <- list(distr = "unif", params = list(min = 0, max = 1))
    
    # Survival after cancer diagnosis
    if (!is.null(file.surv)) {
      for (stg in sort(unique(surv_data$stage))) {
        assign(paste0("time_3", v_cancer[stg], "_Dc"), l_distr_surv[[stg]])
      }
    }
    
  })
  
  
  # Screening
  
  return(l_params_all)
}


#' Update parameters
#'
#' \code{update_param_list} is used to update list of all parameters with new 
#' values for specific parameters
#'
#' @param l_params_all List with all parameters of decision model
#' @param params_updated Parameters for which values need to be updated
#' 
#' @return A list with all parameters updated
#' @export
update_param_list <- function(l_params_all, params_updated){
  
  if (typeof(params_updated) != "list"){
    # Convert the named vector to a list
    params_updated <- split(unname(params_updated), names(params_updated)) 
  }
  # Update the values
  l_params_all <- modifyList(l_params_all, params_updated) 
  return(l_params_all)
}


#' Load list of Excel tables from list of filepaths
#'
#' \code{load_from_excel} loads and processes data from Excel sheets 
#' tables
#'
#' @param filepath File path for tables
#' 
#' @return A labeled list of processed tables
#' 
#' @import tidyverse
#' @import readxl
#' @import dplyr
#' 
#' @export
load_from_excel <- function(filepath) {
  if(!is.null(filepath)) {
    # Read sheet names
    sheet_names <- excel_sheets(filepath)
    
    # Read data
    l_lifetables <- list()
    for(nm in sheet_names) {
      df_lifetable <- read_excel(filepath, sheet = nm)
      l_lifetables <- c(l_lifetables, list(df_lifetable))
    }
    names(l_lifetables) <- sheet_names
    
    # Return list of data tables if there were >1 sheets, otherwise return singular data table
    if(length(sheet_names) > 1) return(l_lifetables)
    else return(l_lifetables[[1]])
  } else {
    # Return nothing if no filepath given
    return()
  }
}
