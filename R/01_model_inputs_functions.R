#' Load default parameters
#'
#' \code{load_default_params} loads default parameters for the decision model 
#' and creates a list
#' 
#' @return 
#' A list of all parameters used for the decision model
#' 
#' @export
load_model_params <- function(
    seed          = 123,                   # Random seed
    n_cohort      = 100000,                # Cohort size
    conf_level    = 0.95,                  # Confidence level
    sex           = c("male", "female"),   # Determines which mortality table(s) to use - male only, female only, or male and female
    p_male        = 0.5,                   # Proportion of males born in population - only used if sex is c("male", "female")
    lesion_state  = FALSE,                 # Indicator to include precancerous lesion state
    n_lesions_max = 20,                    # Maximum number of precancerous lesions
    v_cancer      = seq(4),                # Cancer stages to model in order of disease progression
    p_cancer      = NULL,                  # Probability of being diagnosed at each cancer stage; should be same length at v_cancer
    hr_cancer     = c(1, 1, 1, 1),         # Hazard ratios for progression from one preclinical cancer stage to the next, or detection for stage IV; set to NULL to update cancer progression variables directly
    v_cancer_surv = NULL,                  # Cancer stages from relative survival data corresponding to each modeled cancer stage; if NULL, assumed to be the same as v_cancer
    v_death       = c("o", "c"),           # Death causes
    file.distr    = NULL,                  # Path to distribution file if loading from file
    file.mort     = "data/background_mortality.xlsx",     # Path to background mortality data
    file.surv     = "data/relative_survival_cancer.csv"){ # Path to relative survival data
  
  #### General setup ####
  # Assign health states and reset max number of lesions depending on whether lesions are included
  if (lesion_state) {
    v_states <- c("H", "L", "P", "C", "D")
  } else {
    v_states <- c("H", "P", "C", "D")
    n_lesions_max = NULL
  }
  
  # Assign probability of diagnosis at each cancer stage
  if (is.null(p_cancer)) {
    p_cancer <- rep(1 / length(v_cancer), length(v_cancer))
  }
  
  # Assign relative survival cancer stages if NULL
  if (is.null(v_cancer_surv)) {
    v_cancer_surv <- v_cancer
  }
  
  # Initialize list to store all parameters
  l_params_model <- list(
    seed          = seed,
    n_cohort      = n_cohort,
    conf_level    = conf_level,
    sex           = sex,
    v_states      = v_states,
    v_cancer      = v_cancer,
    p_cancer      = p_cancer,
    diff_p_cancer = rep(0, length(v_cancer)),
    hr_cancer     = hr_cancer,
    v_death       = v_death,
    n_lesions_max = n_lesions_max
  )
  
  #### Load input data ####
  l_lifetables <- load_lifetables(file.mort) # Background mortality
  
  #### Create background mortality distributions and get maximum age ####
  d_time_H_Do <- list()
  max_age <- 0
  for(label in names(l_lifetables)) {
    d_time_H_Do[[label]] <- with(l_params_model, {
      l_distr <- set_mort_distr(
        l_lifetables, 
        label
      )
    })
    max_age <- max(max_age, l_lifetables[[label]]$age)
  }
  max_age <- max_age + 1
  
  # If survival data filepath is given, load disease-specific relative survival 
  # and create distributions
  l_d_time_C_Dc <- list()
  if (!is.null(file.surv)) {
    df_surv <- load_surv_data(file.surv)
    for (i in unique(v_cancer_surv)) {
      # Create distribution data
      l_d_time_C_Dc[[i]] <- set_surv_distr(df_surv, i, max_age)
    }
  } else {
    # If no survival file, create placeholder for true survival distribution 
    # from diagnosis. Manually input parameters after running function in the 
    # form distr = <string> and params = <list of named parameters>
    for (i in unique(v_cancer_surv)) {
      l_d_time_C_Dc[[i]] <- list(distr = NULL,
                                 params = NULL,
                                 src = "known")
    }
  }
  
  #### Update parameter list with distributions ####
  l_params_update <- list()
  l_params_update <- within(l_params_update, {
    # Max age
    max_age <- max_age
    
    # Sex
    d_male <- list(distr = "binom", 
                   params = list(size = 1, prob = p_male), 
                   src = "known")
    
    # Time to death from other causes
    d_time_H_Do <- d_time_H_Do
    
    # Add next distributions depending on whether model includes precancerous lesion state
    if("L" %in% v_states) {
      # Time to precancerous lesion onset
      d_time_H_L <- list(distr = "weibull", 
                         params = list(shape = 1, scale = 1), 
                         src = "unknown")
      
      # Time from precancerous lesion onset to preclinical cancer onset
      d_time_L_P <- list(distr = "exp", 
                         params = list(rate = 1), 
                         src = "unknown")
      
      # Rate of lesion development after onset
      d_n_L <- list(distr = "pois", 
                    params = list(lambda = 1), 
                    src = "unknown")
      
      # Time to additional lesion development after onset
      d_time_L_Lj <- list(distr = "unif",
                          params = list(min = 0, max = 1),
                          src = "assumed")
    } else {
      # Time to preclinical cancer onset
      d_time_H_P = list(distr = "weibull", 
                        params = list(shape = 1, scale = 1), 
                        src = "unknown")
    }
    
    # Time to next stage of preclinical cancer
    for (i in 1:(length(v_cancer) - 1)) {
      assign(paste0("d_time_P", i, "_P", i + 1), 
             list(distr = "exp", 
                  params = list(rate = 1), 
                  src = ifelse(i == 1 | is.null(hr_cancer), "unknown", "calculated")))
    }
    
    # Cancer symptomatic detection by stage
    for (i in 1:length(v_cancer)) {
      assign(paste0("d_time_P", i, "_C", i), 
             list(distr = "exp", 
                  params = list(rate = 1), 
                  src = ifelse(is.null(hr_cancer), "unknown", "calculated")))
    }
    
    # Assign distributions for time to death from cancer by stage at detection
    for (i in 1:length(v_cancer)) {
      assign(paste0("d_time_C", i, "_Dc"), l_d_time_C_Dc[[v_cancer_surv[i]]])
    }
  })
  
  # Remove variable for looping
  l_params_update$i <- NULL
  
  # Add updated variables to l_params_model
  l_params_model <- c(l_params_model, rev(l_params_update))
  
  # Load distribution file if given
  if (!is.null(file.distr)) {
    # Determine file extension
    file_distr_ext <- tools::file_ext(file.distr)
    
    # Read file based on extension
    if (file_distr_ext == "csv") {
      df_distr <- read.csv(file.distr)
    } else if (file_distr_ext %in% c("xls", "xlsx")) {
      df_distr <- read_excel(file.distr, sheet = 1)  # Read first sheet by default
    } else {
      stop("Unsupported file type: ", file_ext)
    }
    
    # Update parameters with distributions from file
    l_params_model <- update_param_from_map(l_params_model, param_map = df_distr, update_distr = T)
  }
  
  return(l_params_model)
}


#' Update parameters from named list
#'
#' \code{update_param_list} is used to update list of all parameters with new 
#' values for specific parameters
#'
#' @param l_params_model List with all parameters of decision model
#' @param params_updated Parameters for which values need to be updated
#' 
#' @return A list with all parameters updated
#' @export
update_param_list <- function(l_params_model, params_updated){
  
  if (typeof(params_updated) != "list"){
    # Convert the named vector to a list
    params_updated <- split(unname(params_updated), names(params_updated)) 
  }
  # Update the values
  l_params_model <- modifyList(l_params_model, params_updated) 
  return(l_params_model)
}

#' Update parameters given in the specified parameter vector and mapping
#'
#' \code{update_param_from_map} is used to update list of all parameters with 
#' new values for specific parameters using a mapping dataframe
#'
#' @param l_params_model List with all parameters of decision model
#' @param v_params_update Vector of new parameter values (not used if 
#' \code{update_distr} = T)
#' @param param_map Dataframe as created by make_param_map below with columns 
#' var_name (parameter label in l_params_model), var_distr (distribution of 
#' parameter), param_name (name of parameter in distribution function), 
#' param_index (index of parameter if a vector, otherwise 1), param_val 
#' (old saved value of parameter), and var_id (unique variable ID consisting of
#' var_name, param_name, and param_index concatenated). If \code{update_distr} 
#' = T, should also include a column param_val with the updated parameter values
#' @param update_distr Logical indicating whether to update distribution type as well
#' 
#' @return A list with all parameters updated
#' @export
update_param_from_map <- function(l_params_model, v_params_update, param_map, 
                                  update_distr = F) {
  # Update parameter values but not distributions
  if (update_distr == F) {
    # Update parameters based on mapping
    assertthat::are_equal(length(v_params_update), nrow(param_map))
    if ("idx" %in% names(param_map)) { # Check if distribution variable has nested distributions
      for (i in seq(length(v_params_update))) {
        # Get value to update variable to
        val <- unlist(v_params_update)[i]
        
        # Update parameter
        if (is.na(param_map$idx[i])) {
          l_params_model[[param_map$var_name[i]]]$params[[param_map$param_name[i]]] <- val
        } else {
          if (is.na(param_map$var_distr[i])) {
            l_params_model[[param_map$var_name[i]]][param_map$idx[i]] <- val
          } else {
            l_params_model[[param_map$var_name[i]]][[param_map$idx[i]]]$params[[param_map$param_name[i]]] <- val
          }
        }
      }
    } else { # If no nested distributions, run same updating code without extra checks
      for (i in seq(length(v_params_update))) {
        # Get value to update variable to
        val <- unlist(v_params_update)[i]
        
        # Update parameter
        l_params_model[[param_map$var_name[i]]]$params[[param_map$param_name[i]]] <- val
      }
    }
  } else { # Update parameter values and distributions
    if ("idx" %in% names(param_map)) { # Check if distribution variable has nested distributions
      # Group parameters by variable name and index
      df_distr_grouped <- param_map %>%
        group_by(var_name, idx) %>%
        summarise(
          distr = first(var_distr),  # Get the distribution type
          params = list(setNames(param_val, param_name)),  # Create a named list for params
          .groups = "drop"
        ) 
      
      # Loop over grouped parameters and update l_params_model for variables in the list
      for (i in seq(nrow(df_distr_grouped))) {
        if (df_distr_grouped$var_name[i] %in% names(l_params_model)) {
          if (is.na(df_distr_grouped$distr[i])) {
            l_params_model[[df_distr_grouped$var_name[i]]][df_distr_grouped$idx[i]] <- unname(df_distr_grouped$params[[i]][[1]])
          } else {
            # Get list of of values to update
            l_distr_updates <- list(distr = df_distr_grouped$distr[i],
                                    params = as.list(df_distr_grouped$params[[i]]))
            
            # Check if distribution has an index within a nested list and update the distribution
            if (is.na(df_distr_grouped$idx[i])) {
              l_params_model[[df_distr_grouped$var_name[i]]]$distr <- l_distr_updates$distr
              l_params_model[[df_distr_grouped$var_name[i]]]$params <- l_distr_updates$params
            } else {
              l_params_model[[df_distr_grouped$var_name[i]]][[df_distr_grouped$idx[[i]]]]$distr <- l_distr_updates$distr
              l_params_model[[df_distr_grouped$var_name[i]]][[df_distr_grouped$idx[[i]]]]$params <- l_distr_updates$params
            }
          }
        }
      }
    } else {
      # Group parameters by variable name
      df_distr_grouped <- param_map %>%
        group_by(var_name) %>%
        summarise(
          distr = first(var_distr),  # Get the distribution type
          params = list(setNames(param_val, param_name)),  # Create a named list for params
          .groups = "drop"
        ) 
      
      # Loop over grouped parameters and update l_params_model for variables in the list
      for (i in seq(nrow(df_distr_grouped))) { # If no nested distributions, run same updating code without extra checks
        if (df_distr_grouped$var_name[i] %in% names(l_params_model)) {
          # Get list of of values to update
          l_distr_updates <- list(distr = df_distr_grouped$distr[i],
                                  params = as.list(df_distr_grouped$params[[i]]))
          
          # Check if distribution has an index within a nested list and update the distribution
          l_params_model[[df_distr_grouped$var_name[i]]]$distr <- l_distr_updates$distr
          l_params_model[[df_distr_grouped$var_name[i]]]$params <- l_distr_updates$params
        }
      }
    }
  }
  
  # If using hazard ratios for cancer progression variables, update variable parameterse
  if (!is.null(l_params_model$hr_cancer)) {
    l_params_model <- modifyList(l_params_model, update_cancer_variables(l_params_model))
  }
  
  return(l_params_model)
}


# Function to make dataframe of parameters
make_param_map <- function(l_params_model, src = 'unknown') {
  param_map <- lapply(names(l_params_model), function(x) {
    # Find distribution variables, which have 'params' item
    var <- l_params_model[[x]]
    if ('params' %in% names(var)) {
      # Only keep distributions that are unknown (not from data or assumption)
      if (var$src %in% src) {
        # Get params
        par <- var$params
        
        # Consolidate params and names
        df <- data.frame(var_name = x,
                         var_distr = var$distr,
                         param_name = names(par),
                         param_val = unname(unlist(par)))
        
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


# Custom function to update parameters for cancer progression and detection 
# based on hazard ratios and proportion diagnosed at each stage
update_cancer_variables <- function(l_params_model) {
  l_vars_update <- with(l_params_model, {
    # Initialize list of updated variables
    l_vars_update <- list()
    
    # Assign distributions for preclinical cancer progression by stage
    if (length(v_cancer) > 2) {
      for (i in 2:(length(v_cancer) - 1)) {
        # Get prior stage progression distribution
        if (i == 2) {
          d_prior <- get(paste0("d_time_P", i - 1, "_P", i))
        } else {
          d_prior <- l_vars_update[[paste0("d_time_P", i - 1, "_P", i)]]
        }
        
        # Calculate current rate
        rate_progress <- hr_cancer[i - 1] * d_prior$params$rate
        
        # Update distribution
        l_vars_update[[paste0("d_time_P", i, "_P", i + 1)]] <- list(
          params = list(rate = rate_progress))
      }
    }
    
    # Calculate shifted cancer stage distribution
    p_cancer_shift <- p_cancer + diff_p_cancer
    
    # Add default distributions for cancer symptomatic detection by stage
    for (i in 1:(length(v_cancer) - 1)) {
      # Get probability of being diagnosed at current stage conditional on not being diagnosed earlier
      p_conditional <- p_cancer_shift[i] / sum(p_cancer_shift[i:length(v_cancer)])
      
      # Get stage distribution
      if (i == 1) {
        d_stage <- get(paste0("d_time_P", i, "_P", i + 1))
      } else {
        d_stage <- l_vars_update[[paste0("d_time_P", i, "_P", i + 1)]]
      }
      
      # Calculate rate consistent with conditional probability of diagnosis at stage
      r_consistent <- p_conditional * d_stage$params$rate / (1 - p_conditional)
      
      # Assign rate
      l_vars_update[[paste0("d_time_P", i, "_C", i)]] <- list(
        params = list(rate =  r_consistent))
    }
    
    # Calculate rate of detection for final stage as hazard ratio applied to sum of prior stage rates
    rate_detect <- hr_cancer[length(v_cancer)] * (d_prior$params$rate + d_stage$params$rate)
    l_vars_update[[paste0("d_time_P", length(v_cancer), "_C", length(v_cancer))]] <- list(
      params = list(rate = rate_detect))
    
    return(l_vars_update)
  })
  
  return(l_vars_update)
}