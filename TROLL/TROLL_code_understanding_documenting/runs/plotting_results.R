#code to plot results from tests regarding WTD
library(ggplot2)
library(gridExtra)
library(dplyr)


### PLOTTING RUNS SEPARATELY
# Define the plotting function
plot_time_series <- function(data, time_col, value_col, line_color = "blue") {
  ggplot(data, aes_string(x = time_col, y = value_col)) +
    geom_line(color = line_color, size = 1) +
    theme_minimal() +
    labs(title = paste(value_col, "Over Time"),
         x = "Iteration",
         y = value_col) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Function to read data and generate plots
plot_variables <- function(file_path, variables, ncol, nrow, time_col = "iter", line_color = "darkblue") {
  data <- read.table(file_path, header = TRUE)
  plots <- lapply(variables, function(col) plot_time_series(data, time_col, col, line_color))
  grid.arrange(grobs = plots, ncol = ncol, nrow = nrow)
}

#########
biogem_file_nowtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/no_WTD/(null)_0_sumstats.txt"
water_file_nowtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/no_WTD/(null)_0_water_balance.txt"

biogem_file_shallowwtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/shallow_WTD/(null)_0_sumstats.txt"
water_file_shallowwtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/shallow_WTD/(null)_0_water_balance.txt"

biogem_file_deepwtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/deep_WTD/(null)_0_"
water_file_deepwtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/deep_WTD/(null)_0_water_balance.txt"

biogem_shallow_30 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_30/shallow_WTD/(null)_0_sumstats.txt"

# Biogeochemical variables
plot_variables(biogem_file_shallowwtd, 
               c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                 'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall'),
               ncol = 3, nrow = 5)

# Water balance variables
plot_variables(water_file_shallowwtd, 
               c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4'),
               ncol = 2, nrow = 3)

plot_variables(biogem_file_deepwtd, 
               c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                 'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall'),
               ncol = 3, nrow = 5)

# Water balance variables
plot_variables(water_file_deepwtd, 
               c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4'),
               ncol = 2, nrow = 3)

plot_variables(biogem_file_nowtd, 
               c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                 'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall'),
               ncol = 3, nrow = 5)

# Water balance variables
plot_variables(water_file_nowtd, 
               c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4'),
               ncol = 2, nrow = 3)

plot_variables(biogem_shallow_30, 
              c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall'),
              ncol = 3, nrow = 5)

# Water balance variables
plot_variables(water_file_nowtd, 
               c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4'),
               ncol = 2, nrow = 3)


####################################
####################################
##### PLOTTING RUNS TOGETHER #######

# Function to read data
read_data <- function(base_path, condition, file_type) {
  file_path <- file.path(base_path, condition, paste0("(null)_0_", file_type, ".txt"))
  read.table(file_path, header = TRUE)
}


# Define paths 
##### Regular climate, redprec50 e redprec30
climate_paths = list(
  regular = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/",
  redprec50 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_50/",
  redprec30 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_30/"
  )

# Define the variables of interest
biogeochem_vars <- c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                     'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall')

water_soil_content <- c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4')

water_soil_potential <- c('SWP_0', 'SWP_1', 'SWP_2', 'SWP_3', 'SWP_4')

transpiration_layer <- c('transpitation_0', 'transpitation_1', 'transpitation_2', 'transpitation_3', 'transpitation_4')

waterflux_vars <- c('precipitation', 'interception', 'throughfall',	'runoff',	'leak',	'evaporation')

#### function with optional plotting nowtd
plot_time_series <- function(shallow, deep, nowtd = TRUE, time_col, value_col) {
  p <- ggplot()
  
  if (!is.null(nowtd)) {
    p <- p + geom_line(data = nowtd, aes_string(x = time_col, y = value_col), 
                       color = "darkblue", linetype = "solid", size = 0.7, alpha = 0.6)
  }
  
  p <- p +
    geom_line(data = shallow, aes_string(x = time_col, y = value_col), 
              color = "darkgreen", linetype = "solid", size = 0.7, alpha = 0.6) +
    geom_line(data = deep, aes_string(x = time_col, y = value_col), 
              color = "darkorange", linetype = "solid", size = 0.7, alpha = 0.6) +
    xlab("") +
    ylab(value_col) +
    theme_minimal()
  
  return(p)
}

get_variable_info <- function(variable_type) {
  switch(variable_type,
         "biogeochemical" = list(vars = biogeochem_vars, file_type = "sumstats"),
         "water soil content" = list(vars = water_soil_content, file_type = "water_balance"),
         "water soil potential" = list(vars = water_soil_potential, file_type = "water_balance"),
         "transpiration layer" = list(vars = transpiration_layer, file_type = "water_balance"),
         "water fluxes" = list(vars = waterflux_vars, file_type = "water_balance"),
         stop("Invalid variable type.")
  )
}

plot_selected_variables <- function(variable_type,base_path, include_nowtd = TRUE) {
  info <- get_variable_info(variable_type)
  vars <- info$vars
  file_type <- info$file_type
  
  shallow <- read_data(base_path, "shallow_WTD", file_type)
  deep <- read_data(base_path, "deep_WTD", file_type)
  
  if (include_nowtd) {
    nowtd <- read_data(base_path, "no_WTD", file_type)
  } else {
    nowtd <- NULL
  }
  
  plots <- lapply(vars, function(var) plot_time_series(shallow, deep, nowtd, "iter", var))
  grid.arrange(grobs = plots, ncol = 3, nrow = ceiling(length(vars) / 3))
}

## Here you can select which climate and it will give you the correspondent
# path and you also choose if you want to include nowtd or not (by default, nowtd = TRUE)
plot_selected_variables("biogeochemical", climate_paths$regular)
plot_selected_variables("water soil content", climate_paths$regular)         
plot_selected_variables("water soil potential", climate_paths$regular)        
plot_selected_variables("transpiration layer", climate_paths$regular)
plot_selected_variables("water fluxes", climate_paths$regular)

plot_selected_variables("biogeochemical", climate_paths$redprec50)
plot_selected_variables("water soil content", climate_paths$redprec50)         
plot_selected_variables("water soil potential", climate_paths$redprec50)        
plot_selected_variables("transpiration layer", climate_paths$redprec50)
plot_selected_variables("water fluxes", climate_paths$redprec50)


plot_selected_variables("biogeochemical", climate_paths$redprec30)
plot_selected_variables("water soil content", climate_paths$redprec30)         
plot_selected_variables("water soil potential", climate_paths$redprec30)        
plot_selected_variables("transpiration layer", climate_paths$redprec30)
plot_selected_variables("water fluxes", climate_paths$redprec30)



##############################
## PLOTTING FOR EACH WTD CONDITION
## DIFFERENT LINES FOR EACH CLIMATE

####################################
####################################
##### PLOTTING RUNS TOGETHER #######

# Function to read data
read_data <- function(base_path, condition, file_type) {
  file_path <- file.path(base_path, condition, paste0("(null)_0_", file_type, ".txt"))
  read.table(file_path, header = TRUE)
}

# Define paths 
##### Regular climate, redprec50 e redprec30
climate_paths = list(
  regular = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/",
  redprec50 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_50/",
  redprec30 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_30/"
)

# Define the variables of interest
biogeochem_vars <- c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                     'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall')

water_soil_content <- c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4')

water_soil_potential <- c('SWP_0', 'SWP_1', 'SWP_2', 'SWP_3', 'SWP_4')

transpiration_layer <- c('transpitation_0', 'transpitation_1', 'transpitation_2', 'transpitation_3', 'transpitation_4')

waterflux_vars <- c('precipitation', 'interception', 'throughfall', 'runoff', 'leak', 'evaporation')

#### function with optional plotting nowtd
plot_time_series_combined <- function(data_list, condition_names, time_col, value_col) {
  p <- ggplot()
  
  for (i in seq_along(data_list)) {
    p <- p + geom_line(data = data_list[[i]], aes_string(x = time_col, y = value_col, color = shQuote(condition_names[i])), 
                       linetype = "solid", size = 0.7, alpha = 0.9)
  }
  
  p <- p +
    xlab("") +
    ylab(value_col) +
    theme_minimal() +
    scale_color_discrete(name = "Climate")
  
  return(p)
}

get_variable_info <- function(variable_type) {
  switch(variable_type,
         "biogeochemical" = list(vars = biogeochem_vars, file_type = "sumstats"),
         "water soil content" = list(vars = water_soil_content, file_type = "water_balance"),
         "water soil potential" = list(vars = water_soil_potential, file_type = "water_balance"),
         "transpiration layer" = list(vars = transpiration_layer, file_type = "water_balance"),
         "water fluxes" = list(vars = waterflux_vars, file_type = "water_balance"),
         stop("Invalid variable type.")
  )
}

plot_combined_climates <- function(variable_type) {
  info <- get_variable_info(variable_type)
  vars <- info$vars
  file_type <- info$file_type
  
  conditions <- c("shallow_WTD", "deep_WTD", "no_WTD")
  condition_names <- names(climate_paths)
  
  for (condition in conditions) {
    data_list <- lapply(climate_paths, function(path) read_data(path, condition, file_type))
    
    plots <- lapply(vars, function(var) plot_time_series_combined(data_list, condition_names, "iter", var))
    
    print(paste("Condition:", condition))
    grid.arrange(grobs = plots, ncol = 3, nrow = ceiling(length(vars) / 3))
  }
}

# Plotting for each condition (shallow, deep, no_WTD) with all climates combined
plot_combined_climates("biogeochemical")
plot_combined_climates("water soil content")
plot_combined_climates("water soil potential")
plot_combined_climates("transpiration layer")
plot_combined_climates("water fluxes")
