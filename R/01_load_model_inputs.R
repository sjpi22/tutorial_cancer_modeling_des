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
                                file.incidence = NULL,
                                file.prevalence = NULL,
                                file.surv = NULL){
  
  #### General setup ####
  # Initialize list to store all parameters
  l_params_all <- list() 
  l_params_all <- within(l_params_all, {
    seed        <- 123                # Random seed
    n_cohort    <- 100000             # Cohort size
    conf_level  <- 0.95               # Confidence level
    v_states    <- c(0, 1, 2, 3, "D") # Health states
    v_lesions   <- c("a", "b")        # Precancerous lesion types
    v_cancer    <- c("i", "ii", "iii", "iv")  # Cancer stages
    v_D         <- c("o", "c")        # Death states
    v_strats    <- c("None", "Screen1", "Screen2")  # CEA strategies
    p_male      <- 0.5                # Probability of being born male
  })
  
  #### Load input data ####
  l_lifetables <- load_from_excel(file.mort) # Background mortality
  # Disease-specific relative survival @@@
  # Strategy parameters (sensitivity, treatment effect) @@@
  # Screening @@@
  
  #### Update parameter list with input data ####
  l_params_all <- within(l_params_all, {
    probs_death_male1 <- l_lifetables$male$p_death 
    probs_death_male0 <- l_lifetables$female$p_death
  })
  
  ### Create table of health states and descriptions
  
  # df_states <- data.frame(
  #   state = v_states,
  #   description = c("Healthy", "Precancerous lesion", "Cancer (preclinical)", "Cancer (clinical)", "Dead"),
  #   preceded_by = lag(v_states),
  #   followed_by = lead(v_states)
  # )
  # df_states$substates = c(list(NULL), list(v_1), list(v_2), list(v_3), list(v_D))
  
  ### Create table with variables and their probability distributions
  df_vars <- setDT(list(
    vartype = character(0), # Type of variable ("bin" for binary, "time" for time to event, "p" for probability)
    from_state = character(0), # Starting state
    to_state = character(0), # Ending state
    subtype = character(0), # Subtype within states (e.g., lesion type; NA if none)
    vargroup = character(0), 
    distr = character(0), 
    params = c())
  )
  
  # Baseline characteristics (demographics, risk factors, time to death from other causes)
  df_vars <- add_row(df_vars,
      varname = "cat_male", 
      description = "Binary variable indicating whether patient is male (1) or female (0)", 
      vargroup = "baseline", 
      distr = "binom", 
      params = list(list(size = 1, prob = "p_male"))) %>%
    # add_row(
    #   varname = "time_0_Do_male1", 
    #   description = "Time from birth to death from other causes (male)", 
    #   vargroup = "baseline", 
    #   distr = "empirical", 
    #   params = list(list(xs = l_lifetables$male$age, probs = l_lifetables$male$p_death))) %>%
    # add_row(
    #   varname = "time_0_Do_male0", 
    #   description = "Time from birth to death from other causes (female)", 
    #   vargroup = "baseline", 
    #   distr = "empirical", 
    #   params = list(list(xs = l_lifetables$female$age, probs = l_lifetables$female$p_death))) %>%
    add_row(
      varname = "time_0_Do", 
      description = "Time from birth to death from other causes (female)", 
      vargroup = "baseline", 
      distr = "empirical", 
      params = list(list(xs = l_lifetables$female$age, probs = l_lifetables$female$p_death, max_x = l_params_all$max_age)))
    
  # Precancerous lesions: Every variable name should end in the lesion type
  for (lesiontype in l_params_all$v_lesions) {
    # Initialize temporary data table for lesion variables
    temp_df_vars <- setDT(list(
      varname = character(0), 
      subtype = character(0),
      varID = character(0),
      description = character(0), 
      vargroup = character(0), 
      distr = character(0), 
      params = c())
    )
    
    # Add rows for lesion variables
    temp_df_vars <- add_row(
      temp_df_vars,
      varname = "time_0_1", 
      description = paste("Time from birth to first precancerous lesion of type", lesiontype),
      vargroup = "lesion pt-level", 
      distr = "weibull", 
      params = list(list(shape = 2, scale = 75))
    ) %>%
      add_row(varname = "n_add", 
              description = paste("Number of additional precancerous lesions of type", lesiontype, "in next max_age years"),
              vargroup = "lesion pt-level", 
              distr = "pois", 
              params = list(list(lambda = 0.05 * l_params_all$max_age))) %>%
      add_row(varname = "time_1_1j", 
              description = paste("Time from first precancerous lesion of type", lesiontype, "to additional lesion of same type"),
              vargroup = "lesion lesion-level", 
              distr = "unif", 
              params = list(list(max = l_params_all$max_age))) %>%
      add_row(varname = "time_1j_2i", 
              description = paste("Time from onset of precancerous lesion of type", lesiontype, "to conversion to cancer"),
              vargroup = "lesion lesion-level", 
              distr = "gamma", 
              params = list(list(shape = 8, scale = 8))) %>%
      add_row(varname = "size_1j_2", 
              description = paste("Size of precancerous lesion of type", lesiontype, "at conversion to cancer"),
              vargroup = "lesion lesion-level", 
              distr = "unif", 
              params = list(list(min = 10, max = 40))) %>%
      add_row(varname = "size_1j_t", 
              description = paste("Size of precancerous lesion of type", lesiontype, "at time t from lesion onset"),
              vargroup = "screening lesion", 
              distr = "direct", 
              params = list(list(paste0("exp(", paste("size_1j_2", lesiontype, sep="_"), "/", paste("time_1j_2", lesiontype, sep="_"), ")"))))
    
    # Set lesion type as subtype
    temp_df_vars$subtype = lesiontype
    
    # Add lesion type to vargroup
    temp_df_vars$vargroup = paste(temp_df_vars$vargroup, lesiontype)
    
    # Join with full data table of variables
    df_vars <- rbind(df_vars, temp_df_vars)
  }
  
  # Cancer
  for (i in seq(length(l_params_all$v_cancer)-1)) {
    df_vars <- add_row(
      df_vars,
      varname = paste0("time_2", l_params_all$v_cancer[i], "_2", l_params_all$v_cancer[i+1]), 
      description = paste("Time from preclinical cancer stage", toupper(l_params_all$v_cancer[i]), "to", toupper(l_params_all$v_cancer[i+1])),
      vargroup = "cancer", 
      distr = "gamma", 
      params = list(list(shape = 2, scale = 2))
    )
  }
  
  # Screening
  
  # Set ID var
  df_vars$varID <- paste0(df_vars$varname, ifelse(is.na(df_vars$subtype), "", paste0("_", df_vars$subtype)))
  
  l_params_all$df_vars <- df_vars
  
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
