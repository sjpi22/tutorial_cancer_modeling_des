###########################  Create background mortality distributions  ##########################
#
#  Objective: Save artificial lifetable data
########################### <<<<<>>>>> #########################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(data.table)
library(VGAM)

###### 1.2 Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

# Output path
path_out <- "data/simulated"

# Distribution parameters
shape_male <- 0.000078
scale_male <- 0.08
multiplier_female <- 1/1.02
max_age <- 110
n_total <- 100000

# Plot to visualize
curve(pgompertz(x, shape = shape_male, scale = scale_male), 0, 110, ylim = c(0, 1)) # male
curve(pgompertz(x, shape = shape_male*multiplier_female, scale = scale_male*multiplier_female), 0, 110, ylim = c(0, 1)) # female


#### 3. Generate data ========================================================

# Create necessary variables
v_ages <- seq(0, max_age) # vector of ages
shape <- c(male = shape_male, female = shape_male*multiplier_female)
scale <- c(male = scale_male, female = scale_male*multiplier_female)

for (sex in c("male", "female")) {
  # Create mortality data table
  df_mort <- data.table(
    Age = v_ages,
    pct_dead = pgompertz(v_ages, shape = shape[sex], scale = scale[sex])
  )
  
  # Calculate probability of dying (qx) and number of survivors (lx)
  df_mort[, `:=` (qx = diff(c(pct_dead, 1)),
                  lx = n_total*(1-pct_dead))]
  
  # Remove pct_dead variable
  df_mort[, pct_dead := NULL]
  
  # Write to csv file
  write.table(df_mort, 
            file = file.path(path_out, paste0("background_mortality_", sex, ".txt")),
            row.names=FALSE)
}
