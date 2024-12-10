################################################################################
# Functions to load and process individual data input files
################################################################################

################### Background mortality ###################
#' Load list of life tables from list of filepaths
#'
#' \code{load_lifetables} loads and processes Excel workbook with sex-specific 
#' life tables
#'
#' @param l_filepaths A list of file paths for life tables; items can 
#' optionally have names that will be carried down to the list of results
#' 
#' @return A labeled list of processed life tables
#' 
#' @import tidyverse
#' @import readxl
#' 
#' @export
load_lifetables <- function(filepath) {
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


#' Set background mortality distribution for a specific year from a life table 
#'
#' \code{set_mort_distr} sets mortality distributions from a life table
#'
#' @param l_lifetables A list of life tables
#' @param label Label for desired life table (default first table in list)
#' 
#' @return Distribution object with distribution type and parameters
#' 
#' @import dplyr
#' 
#' @export
# 
set_mort_distr <- function(l_lifetables, label = 1) {
  
  lifetable_dat <- l_lifetables[[label]]
  
  distr <- list(distr = "empirical",
                params = list(xs = lifetable_dat$age, 
                              probs = lifetable_dat$p_death,
                              max_x = max(lifetable_dat$age) + 1),
                src = "known")
  return(distr)
}


################### Survival from diagnosis ###################

#' Load survival from diagnosis data
#'
#' \code{load_surv_data} reads survival from diagnosis data
#' 
#' @param filepath String with the location and name of the data file
#' @param max_age Maximum age to use for survival data
#' 
#' @return A data frame with survival from diagnosis by stage
#' 
#' @import tidyverse
#' @import readxl
#' @import dplyr
#' 
#' @export
load_surv_data <- function(
    filepath,
    max_age = 110
){
  surv_data <- read.csv(filepath) %>%
    mutate(pct_died = 1 - surv)
  return(surv_data)
}


################### Target data ###################

# Load preclinical cancer prevalence data
load_prevalence <- function(filepath, target_groups = NULL) {
  if (is.null(target_groups)) {
    target_groups <- "prevalence"
  }
  filedata <- read.csv(filepath) %>%
    mutate(target_index = (age_start + age_end)/2,
           target_groups = target_groups,
           target_names = paste(target_groups, age_start, age_end, sep="_"))
  return(filedata)
}

# Load cancer incidence data
load_incidence <- function(filepath, target_groups = NULL) {
  if (is.null(target_groups)) {
    target_groups <- "incidence"
  }
  filedata <- read.csv(filepath) %>%
    mutate(target_index = (age_start + age_end)/2,
           target_groups = target_groups,
           target_names = paste(target_groups, age_start, age_end, sep="_"))
  return(filedata)
}

# Load cancer stage distribution data
load_stage_distr <- function(filepath, target_groups = NULL) {
  if (is.null(target_groups)) {
    target_groups <- "stage_distr"
  }
  filedata <- read.csv(filepath) %>% 
    rename(target_index = stage_dx) %>%
    mutate(target_groups = target_groups,
           target_names = paste(target_groups, target_index, sep="_"))
  return(filedata)
}

#' Load calibration targets
#'
#' @param filepath String with the location and name of the file with data
#' 
#' @return A data.table object with calibration targets
#' 
#' @import tidyverse
#' @import readxl
#' 
#' @export
load_calibration_targets <- function(l_filepaths){
  # Read files into list of data files
  l_true_targets <- list()
  for (label in names(l_filepaths)) {
    filedata <- do.call(paste0("load_", label), 
                        list(l_filepaths[[label]], label)) %>%
      dplyr::select(any_of(c("target_names", "target_groups", 
                             "target_index",  
                             "age_start", "age_end",
                             "targets", "se", 
                             "ci_lb", "ci_ub", 
                             "total_atrisk", "n_events",
                             "sex", "lesion_type")))
    l_true_targets[[label]] <- filedata 
  }
  return(l_true_targets)
}