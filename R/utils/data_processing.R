################################################################################
# Functions to load and process data
################################################################################

################### Background mortality ###################

#' Read lifetable text file from Human Mortality Database
#'
#' \code{read_lifetable} is used to load life table data from a Human Mortality 
#' Database .txt file into a data.table when called by 
#' \code{load_lifetables}
#'
#' @param filepath String with the location and name of the file with life table 
#' data
#' @param skip Number of rows to skip before reading data
#' @param header Whether data has header with variable names
#' 
#' @return A data.table object with the life table's raw data
#' 
#' @import data.table
#' 
#' @export
read_lifetable <- function(filepath, skip=2, header=TRUE) {
  lifetable <- data.table(read.table(filepath, skip=skip, header=header))
  return(lifetable)
}


#' Process lifetable data frame (Human Mortality Database format)
#'
#' \code{process_lifetable} processes the Human Mortality Database life tables 
#' to include additional variables used in downstream analysis. These include 
#' \code{Age} in numeric form
#' \code{yob} as year of birth
#' \code{S} as the probability of survival to the given age
#' \code{neg_log_S} as negative log survival (for exponential fit)
#' \code{log_hr} as the log of the hazard (for Gompertz fit)
#' \code{Fx} as the probability of death before the given age (complement of S)
#'
#' @param lifetable_df Dataframe or data table with life table data
#' 
#' @return A life table in the format of the inputted life table with additional
#' variables (listed in function description)
#' 
#' @import dplyr
#' 
#' @export
process_lifetable <- function(lifetable_df, radix=100000) {
  processed_lifetable_df <- lifetable_df  %>%
    mutate(Age = as.integer(gsub("\\+", "", as.character(Age))),
           S = lx / radix,
           neg_log_S = -log(lx / radix),
           Fx = 1 - S)
  
  # Process depending on whether Year is in data
  if ("Year" %in% colnames(processed_lifetable_df)) {
    processed_lifetable_df <- processed_lifetable_df %>%
      mutate(yob = Year - Age) %>%
      group_by(Year) %>%
      mutate(px = c(diff(Fx), 1 - Fx[n()])) %>% # Probability mass
      ungroup()
  } else {
    processed_lifetable_df <- processed_lifetable_df %>%
      mutate(px = c(diff(Fx), 1 - Fx[n()])) # Probability mass
  }
  
  return(processed_lifetable_df)
}


#' Load list of life tables from list of filepaths
#'
#' \code{load_lifetables} loads and processes Human Mortality Database life 
#' tables
#'
#' @param l_filepaths A list of file paths for life tables; items can 
#' optionally have names that will be carried down to the list of results
#' 
#' @return A labeled list of processed life tables
#' 
#' @import dplyr
#' 
#' @export
load_lifetables <- function(l_filepaths, skip=2, header=TRUE, radix=100000) {
  # Read lifetables and add to list
  l_lifetables <- c()
  for (filepath in l_filepaths) {
    lifetable <- read_lifetable(filepath=filepath, skip=skip, header=header) %>%
      process_lifetable(radix = radix)
    l_lifetables <- c(l_lifetables, list(lifetable))
  }
  
  # Set names of lifetables as names of input list
  if (!is.null(names(l_filepaths))) names(l_lifetables) <- names(l_filepaths)
  return(l_lifetables)
}


#' Set background mortality distribution for a specific year from a life table 
#'
#' \code{set_mort_distr} sets mortality distributions from Human Mortality Database lifetables
#' tables
#'
#' @param l_lifetables A list of life tables
#' @param label Label for desired life table (default first table in list)
#' @param year Year for life table
#' 
#' @return Distribution object with distribution type and parameters
#' 
#' @import dplyr
#' 
#' @export
# 
set_mort_distr <- function(l_lifetables, label = 1, year = NULL) {
  lifetable_dat <- l_lifetables[[label]]
  
  # Filter to year if necessary
  if ("Year" %in% colnames(lifetable_dat)) {
    if (is.null(year)) { # If null, use most recent year
      lifetable_dat <- lifetable_dat %>%
        filter(Year == max(year))
    } else {
      lifetable_dat <- lifetable_dat %>%
        filter(Year == year)
    }
  }
  
  # Create distribution
  distr <- list(distr = "empirical",
                params = list(xs = lifetable_dat$Age, 
                              probs = lifetable_dat$px,
                              max_x = max(lifetable_dat$Age) + 1),
                src = "known")
  return(distr)
}


################### Survival from diagnosis ###################

#' Read survival from diagnosis data downloaded from SEER
#'
#' \code{read_CRC_surv} reads survival from diagnosis data downloaded from SEER
#' 
#' @param filepath String with the location and name of the file with CRC
#' incidence data
#' 
#' @return A data frame with survival from diagnosis by stage
#' 
#' @import tidyverse
#' @import readxl
#' @import dplyr
#' 
#' @export
load_surv_data <- function(
    l_filepaths,
    version = c("SEER", "simulated"),
    max_age = 110
){
  if (version == "simulated") {
    # Read data and calculate CDF of death
    surv_data <- read.csv(l_filepaths) %>%
      mutate(pct_died = 1 - surv)
    
    # Set distribution for each stage
    l_distr_surv <- list()
    for (i in unique(surv_data$stage)) {
      distr_surv <- set_surv_distr(surv_data, i, max_age)
      l_distr_surv[[i]] <- distr_surv
    }
    
    # Return as data.table
    return(list(l_data_surv = surv_data, l_distr_surv = l_distr_surv))
  } else {
    # Read survival data tables and add to list
    l_data_surv <- c()
    l_distr_surv <- list()
    data_names <- c("relative_survival", "lb_95ci", "ub_95ci")
    
    for (filepath in l_filepaths) {
      # Read survival data and clean NAs
      surv_data <- read_csv(filepath, skip = 3) %>%
        drop_na() %>%
        slice(-1)
      surv_data[surv_data == "^"] <- NA
      
      # Change column names
      orig_cols <- names(surv_data)
      cols_add <- c("", paste("__", rep(c(1, 2, 3), (length(orig_cols) - 1)/3), sep = ""))
      
      new_cols <- paste(sub("\\...+", '', orig_cols), cols_add, sep = "")
      new_cols[1] <- "time_char"
      names(surv_data) <- new_cols
      
      # Change from wide to long
      surv_data_long <- surv_data %>%
        pivot_longer(cols = -"time_char", 
                     names_to = c("age_group", ".value"), 
                     names_sep = "__") %>%
        mutate(time_char = sub("Diagnosis", "0", time_char)) %>%
        mutate(time = as.integer(sub("(\\w+).*", "\\1", time_char)),
               age_min = as.integer(substr(age_group, 5, 7)),
               age_max = as.integer(substring(age_group, nchar(age_group)-1)),
               age_group = gsub("&lt; ", "<", age_group)) %>%
        mutate(age_group = gsub("Ages ", "", age_group)) %>%
        mutate(age_group = gsub(" Ages", "", age_group))
      
      # Correct NAs
      surv_data_long$age_min <- surv_data_long$age_min %>%
        replace_na(0)
      surv_data_long$age_max <- surv_data_long$age_max %>%
        replace_na(max_age)
      
      # Change column names for data and set to numeric
      names(surv_data_long)[3:5] <- data_names
      surv_data_long <- surv_data_long %>%
        filter(!is.na(relative_survival)) %>%
        mutate_at(data_names, as.numeric) %>%
        mutate(cdf = 1 - relative_survival/100) %>%
        group_by(age_group) %>%
        mutate(pmf = c(diff(cdf), 1 - max(cdf)))
      
      # Filter initially to all ages and set distribution@@@
      surv_data_all <- surv_data_long %>%
        filter(age_group == "All")  
      
      distr_surv <- list(distr = "empirical",
                         params = list(xs = surv_data_all$time,
                                       probs = surv_data_all$pmf,
                                       max_x = max_age),
                         src = "known")
      l_distr_surv <- c(l_distr_surv, list(distr_surv))
      
      l_data_surv <- c(l_data_surv, list(surv_data_long))
    }
    names(l_distr_surv) <- names(l_filepaths)
    names(l_data_surv) <- names(l_filepaths)
    
    # Return as data.table
    return(list(l_data_surv = l_data_surv, l_distr_surv = l_distr_surv))
  }
}


#' Set relative survival distribution for a stage of cancer
#'
#' \code{set_surv_distr} sets relative survival distributions for a stage of cancer.
#'
#' @param df_surv Data frame with columns \code{years_from_dx} for years from
#'   diagnosis, \code{pct_died} for cumulative percentage dead, and \code{stage}
#'   for stage of cancer
#' @param stage Stage of cancer
#' @param max_age Maximum survival age
#' 
#' @return Distribution object with distribution type and parameters
#' 
#' @export
# 
set_surv_distr <- function(df_surv, stage, max_age) {
  
  # Filter survival data to stage at diagnosis
  temp_df_surv <- df_surv[df_surv$stage == stage, ]
  
  # Calculate probability mass function from CDF
  probs <- diff(temp_df_surv$pct_died)
  probs <- c(probs, 1 - sum(probs))
  
  distr <- list(distr = "empirical", 
                params = list(xs = temp_df_surv$years_from_dx, 
                              probs = probs, 
                              max_x = max_age), 
                src = "known")
  return(distr)
}


################### Target data ###################

#' Load calibration targets
#'
#' @param filepath String with the location and name of the file with data
#' 
#' @return A data.table object with calibration targets
#' 
#' @import tidyverse
#' @import readxl
#' @import tools
#' 
#' @export
load_calibration_targets <- function(file_targets){
  # Read file
  if (file_ext(file_targets) %in% c("xls", "xlsx")) {
    df_targets <- read_excel(file_targets)  
  } else {
    df_targets <- read.csv(file_targets)
  }
  
  return(df_targets)
}