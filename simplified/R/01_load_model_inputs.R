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
    v_states    = c(0, 1, 2, "D"), # Health states
    v_cancer    = c("i", "ii"),       # Cancer stages in order
    v_D         = c("o", "c"),        # Death states
    v_strats    = c("None", "Screen") # CEA strategies
  )
  
  #### Load input data ####
  l_lifetables <- load_from_excel(file.mort) # Background mortality
  
  # Disease-specific relative survival
  if (!is.null(file.surv)) {
    surv_data <- read.csv(file = file.surv) %>%
      mutate(pct_died = 1 - surv)
    assertthat::validate_that(length(l_params_all$v_cancer) == length(unique(surv_data$stage)), 
                            msg = 'Number of cancer states in survival file not consistent with number of states expected in v_cancer')
  }
  
  #### Create useful variables ####
  # Maximum age
  max_age <- max(l_lifetables$female$age, l_lifetables$male$age) + 1
  
  # Survival after cancer diagnosis
  if (!is.null(file.surv)) {
    l_distr_surv <- list()
    for (i in sort(unique(surv_data$stage))) {
      # Filter survival data to stage at diagnosis
      temp_surv_data <- surv_data[surv_data$stage == i, ]
      
      # Calculate probability mass function from CDF
      # probs <- pmax(diff(temp_surv_data$pct_died), 0)
      probs <- diff(temp_surv_data$pct_died)
      probs <- c(probs, 1 - sum(probs))
      
      # Create distribution data
      l_distr_surv[[i]] <- list(distr = "empirical", 
                                params = list(xs = temp_surv_data$years_from_dx, 
                                              probs = probs, 
                                              max_x = max_age), 
                                src = "known",
                                description = paste("Survival from diagnosis at stage", i))
    }
  } else {
    # If no survival file uploaded, use true survival distribution of exponential from diagnosis
    l_distr_surv <- list()
    for (i in 1:length(l_params_all$v_cancer)) {    
      l_distr_surv[[i]] <- list(distr = NULL,
                                params = NULL,
                                src = "known",
                                description = paste("Survival from diagnosis at stage", i))
    }
  }
  
  #### Update parameter list with distributions ####
  l_params_all <- within(l_params_all, {
    # Max age
    max_age <- max_age
    
    # Sex
    b_male <- list(distr = "binom", params = list(size = 1, prob = 0.5), src = "known")
    
    # Time to death from other causes
    time_0_Do_male <- list(distr = "empirical", 
                           params = list(xs = l_lifetables$male$age, 
                                         probs = l_lifetables$male$p_death, 
                                         max_x = max_age), 
                           src = "known",
                           description = "Background mortality for males")
    
    time_0_Do_female <- list(distr = "empirical", 
                             params = list(xs = l_lifetables$female$age, 
                                           probs = l_lifetables$female$p_death,
                                           max_x = max_age),
                             src = "known",
                             description = "Background mortality for females")
    
    # Time to cancer onset
    time_0_1 = list(distr = "weibull", 
                    params = list(shape = 1, scale = 1), 
                    src = "unknown",
                    description = "Time from birth to cancer onset")
      
    # Cancer preclinical stage progression
    time_1i_1ii <- list(distr = "exp", 
                        params = list(rate = 1), 
                        src = "unknown",
                        description = "Time from early to late stage cancer")
    
    # Cancer symptomatic detection by stage
    time_1i_2 <- list(distr = "exp", 
                      params = list(rate = 1), 
                      src = "unknown",
                      description = "Time from early stage cancer onset to symptomatic detection")
    time_1ii_2 <- list(distr = "exp", 
                       params = list(rate = 1), 
                       src = "unknown",
                       description = "Time from late stage cancer onset to symptomatic detection")
    
    # Survival after cancer diagnosis
    time_2i_Dc <- l_distr_surv[[1]]
    time_2ii_Dc <- l_distr_surv[[2]]
    
  })
  
  # Screening
  
  return(l_params_all)
}


#' Update parameters from named list
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

#' Update parameters given in the specified parameter vector and mapping
#'
#' \code{update_param_from_map} is used to update list of all parameters with 
#' new values for specific parameters using a mapping dataframe
#'
#' @param l_params_all List with all parameters of decision model
#' @param v_params_update Vector of new parameter values
#' @param param_map Dataframe as created by make_param_map below with columns 
#' var_name (parameter label in l_params_all), var_distr (distribution of 
#' parameter), param_name (name of parameter in distribution function), 
#' param_index (index of parameter if a vector, otherwise 1), param_val 
#' (old saved value of parameter), and var_id (unique variable ID consisting of
#' var_name, param_name, and param_index concatenated)
#' 
#' @return A list with all parameters updated
#' @export
update_param_from_map <- function(l_params_all, v_params_update, param_map) {
  # Make copy of parameter list
  l_params_update <- copy(l_params_all)
  
  # Update parameters based on mapping
  assertthat::are_equal(length(v_params_update), nrow(param_map))
  for (i in seq(length(v_params_update))) {
    # Get value to update variable to
    val <- v_params_update[i]
    
    # Get parameter to update
    var_info <- unlist(param_map[i,])
    l_params_update[[var_info['var_name']]]$params[[var_info['param_name']]] <- val
  }
  
  return(l_params_update)
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


# Function to make dataframe of parameters
make_param_map <- function(l_params_all, src = 'unknown') {
  param_map <- lapply(names(l_params_all), function(x) {
    # Find distribution variables, which have 'params' item
    var <- l_params_all[[x]]
    if ('params' %in% names(var)) {
      # Only keep distributions that are unknown (not from data or assumption)
      if (var$src == src) {
        # Get params
        par <- var$params
        
        # Consolidate params and names
        df <- data.frame(var_name = x,
                         var_distr = var$distr,
                         param_name = names(par),
                         param_val = unname(unlist(par)),
                         var_description = var$description)
        
        return(df)
        
      }
    }
  })
  
  # Combine list into dataframe
  param_map <- rbindlist(param_map) %>%
    # Set variable identifier
    mutate(var_id = paste(var_name, param_name, sep = '.'))
  
  return(param_map)
}
