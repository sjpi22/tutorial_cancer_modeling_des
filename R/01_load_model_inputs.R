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
    v_cancer      = seq(4),                # Cancer stages in order
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
  
  # Initialize list to store all parameters
  l_params_all <- list(
    seed          = seed,
    n_cohort      = n_cohort,
    conf_level    = conf_level,
    sex           = sex,
    v_states      = v_states,
    v_cancer      = v_cancer,
    v_death       = v_death,
    n_lesions_max = n_lesions_max
  )
  
  #### Load input data ####
  l_lifetables <- load_lifetables(file.mort) # Background mortality
  
  #### Create background mortality distributions and get maximum age ####
  d_time_H_Do <- list()
  max_age <- 0
  for(label in names(l_lifetables)) {
    d_time_H_Do[[label]] <- with(l_params_all, {
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
  if (!is.null(file.surv)) {
    surv_data <- load_surv_data(file.surv)
    d_time_C_Dc <- list()
    for (i in v_cancer) {
      # Filter survival data to stage at diagnosis
      temp_surv_data <- surv_data[surv_data$stage == i, ]
      
      # Calculate probability mass function from CDF
      probs <- diff(temp_surv_data$pct_died)
      probs <- c(probs, 1 - sum(probs))
      
      # Create distribution data
      d_time_C_Dc[[i]] <- list(distr = "empirical", 
                               params = list(xs = temp_surv_data$years_from_dx, 
                                             probs = probs, 
                                             max_x = max_age), 
                               src = "known")
    }
  } else {
    # If no survival file, create placeholder for true survival distribution 
    # from diagnosis. Manually input parameters after running function in the 
    # form distr = <string> and params = <list of named parameters>
    d_time_C_Dc <- list()
    for (i in v_cancer) {    
      d_time_C_Dc[[i]] <- list(distr = NULL,
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
    
    # Add default distributions for cancer progression by stage
    for (i in 1:length(v_cancer)) {
      if (i < length(v_cancer)) {
        # Time to next stage of preclinical cancer
        assign(paste0("d_time_P", v_cancer[i], "_P", v_cancer[i+1]), 
               list(distr = "exp", 
                    params = list(rate = 1), 
                    src = "unknown"))
      }
      
      # Cancer symptomatic detection by stage
      assign(paste0("d_time_P", v_cancer[i], "_C"), 
             list(distr = "exp", 
                  params = list(rate = 1), 
                  src = "unknown"))
    }
    
    # Survival after cancer diagnosis by stage
    d_time_C_Dc <- d_time_C_Dc
  })
  
  # Remove variable for looping
  l_params_update$i <- NULL
  
  # Add updated variables to l_params_all
  l_params_all <- c(l_params_all, rev(l_params_update))
  
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
    
    # Group parameters by variable name and index
    df_distr_grouped <- df_distr %>%
      group_by(var_name, idx) %>%
      summarise(
        distr = first(var_distr),  # Get the distribution type
        params = list(setNames(param_val, param_name)),  # Create a named list for params
        .groups = "drop"
      ) 
    
    # Loop over grouped parameters and update l_params_all for variables in the list
    for (i in seq(nrow(df_distr_grouped))) {
      if (df_distr_grouped$var_name[i] %in% names(l_params_all)) {
        # Get list of of values to update
        l_distr_updates <- list(distr = df_distr_grouped$distr[i],
                                params = as.list(df_distr_grouped$params[[i]]))
        
        # Check if distribution has an index within a nested list and update the distribution
        if (is.na(df_distr_grouped$idx[i])) {
          l_params_all[[df_distr_grouped$var_name[i]]]$distr <- l_distr_updates$distr
          l_params_all[[df_distr_grouped$var_name[i]]]$params <- l_distr_updates$params
        } else {
          l_params_all[[df_distr_grouped$var_name[i]]][[df_distr_grouped$idx[[i]]]]$distr <- l_distr_updates$distr
          l_params_all[[df_distr_grouped$var_name[i]]][[df_distr_grouped$idx[[i]]]]$params <- l_distr_updates$params
        }
      }
    }
  }
  
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
    val <- unlist(v_params_update)[i]
    
    # Update parameter
    if (is.na(param_map$idx[i])) {
      l_params_update[[param_map$var_name[i]]]$params[[param_map$param_name[i]]] <- val
    } else {
      l_params_update[[param_map$var_name[i]]][[param_map$idx[i]]]$params[[param_map$param_name[i]]] <- val
    }
  }
  
  return(l_params_update)
}


# Function to make dataframe of parameters
make_param_map <- function(l_params_all, src = 'unknown') {
  param_map <- lapply(names(l_params_all), function(x) {
    # Find distribution variables, which have 'params' item
    var <- l_params_all[[x]]
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
