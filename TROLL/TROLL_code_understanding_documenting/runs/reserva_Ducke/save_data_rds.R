# ============================================================
# Load and aggregate TROLL output files
# ============================================================

# Load required packages
library(tidyverse)
library(fs)
library(dplyr)

# ------------------------------------------------------------
# 1. Define the main project directory
# ------------------------------------------------------------

project_dir <- path_expand(
  "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke"
)

# ------------------------------------------------------------
# 2. Define the simulation scenarios and their output folders
# ------------------------------------------------------------

runs <- tibble(
  scenario = c(
    "ducke_deep_clayey_regclim_ks14",
    "ducke_deep_clayey_redprec30_ks14",
    "ducke_deep_sandy_redprec30_ks14",
    "ducke_deep_sandy_regclim_ks14",
    "ducke_shallow_sandy_regclim_ks14",
    "ducke_shallow_sandy_redprec30_ks14",
    "ducke_shallow_clayey_redprec30_ks14",
    "ducke_shallow_clayey_regclim_ks14"
  ),
  output_dir = path(
    project_dir,
    c(
      "deepWT_clayeysoil/ducke_deep_clayey_regclim_ks14",
      "deepWT_clayeysoil/ducke_deep_clayey_redprec30_ks14",
      "deepWT_sandysoil/ducke_deep_sandy_redprec30_ks14",
      "deepWT_sandysoil/ducke_deep_sandy_regclim_ks14",
      "shallowWT_sandysoil/ducke_shallow_sandy_regclim_ks14",
      "shallowWT_sandysoil/ducke_shallow_sandy_redprec30_ks14",
      "shallowWT_clayeysoil/ducke_shallow_clayey_redprec30_ks14",
      "shallowWT_clayeysoil/ducke_shallow_clayey_regclim_ks14"
    )
  )  
) %>%
  mutate(
    final_pattern_file = path(output_dir, "output")
  )

# ------------------------------------------------------------
# 3. Function to find one output file
# ------------------------------------------------------------

find_output_file <- function(output_dir, pattern, file_label) {
  
  # Search for files that match the expected pattern
  files <- list.files(
    output_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  # Stop if the file does not exist
  if (length(files) == 0 || !file_exists(files[1])) {
    stop("Missing ", file_label, " file in: ", output_dir)
  }
  
  # Return the first matching file
  files[1]
}

# ------------------------------------------------------------
# 4. Function to remove empty columns created by read_tsv()
# ------------------------------------------------------------

clean_output_columns <- function(df) {
  
  # Remove unnamed columns such as ...44
  df %>%
    select(-matches("^\\.\\.\\.[0-9]+$"))
}

# ------------------------------------------------------------
# 5. Function to read one simulation folder
# ------------------------------------------------------------

read_simulation_data <- function(scenario, output_dir) {
  
  # Find water balance file
  file_wb <- find_output_file(
    output_dir = output_dir,
    pattern = "water_balance\\.txt$",
    file_label = "water balance"
  )
  
  # Find vertical water flux file
  file_vf <- find_output_file(
    output_dir = output_dir,
    pattern = "vertical_water_flux\\.txt$",
    file_label = "vertical water flux"
  )
  
  # Find biogeochemical summary file
  file_biog <- find_output_file(
    output_dir = output_dir,
    pattern = "sumstats\\.txt$",
    file_label = "biogeochemical"
  )
  
  # Read water balance output
  wb <- read_tsv(file_wb, show_col_types = FALSE) %>%
    clean_output_columns() %>%
    mutate(
      iter = as.numeric(iter),
      scenario = scenario,
      year = floor(iter / 365) + 1
    )
  
  # Read vertical flux output
  vf <- read_tsv(file_vf, show_col_types = FALSE) %>%
    clean_output_columns() %>%
    mutate(
      iter = as.numeric(iter),
      scenario = scenario,
      year = floor(iter / 365) + 1
    )
  
  # Read biogeochemical output
  biog <- read_tsv(file_biog, show_col_types = FALSE) %>%
    clean_output_columns() %>%
    mutate(
      iter = as.numeric(iter),
      scenario = scenario,
      year = floor(iter / 365) + 1
    )
  
  # Return the three outputs as a list
  list(
    wb = wb,
    vf = vf,
    biog = biog
  )
}

# ------------------------------------------------------------
# 6. Read all simulations
# ------------------------------------------------------------

# Apply the reading function to all scenarios
all_data <- map2(
  runs$scenario,
  runs$final_pattern_file,
  read_simulation_data
)

# Bind all scenarios into three daily dataframes
df_water_balance <- map_dfr(all_data, "wb")
df_vertical_flux <- map_dfr(all_data, "vf")
df_biogem        <- map_dfr(all_data, "biog")

# Remove temporary object
rm(all_data)

# ------------------------------------------------------------
# 7. Create cache directory and save daily dataframes
# ------------------------------------------------------------

cache_dir <- path(project_dir, "_rds")
dir_create(cache_dir)

saveRDS(df_water_balance, path(cache_dir, "df_water_balance_daily.rds"))
saveRDS(df_vertical_flux, path(cache_dir, "df_vertical_flux_daily.rds"))
saveRDS(df_biogem,        path(cache_dir, "df_biogem_daily.rds"))

# ============================================================
# Aggregate annual data
# ============================================================

# ------------------------------------------------------------
# 8. Water balance annual data
# ------------------------------------------------------------

# Identify water balance variables that are specific to soil layers.
# Examples:
# SWC_0, SWC_1, SWP_0, SWP_1, transpiration_0, transpiration_1
layer_vars_wb <- names(df_water_balance) %>%
  str_subset("^(SWC|SWP|transpiration)_[0-9]+$")

# Convert water balance layer variables from wide to long format.
# This creates:
# - variable: SWC, SWP, or transpiration
# - layer: soil layer number
# - value: the value of the variable
df_water_balance_annual <- df_water_balance %>%
  select(
    scenario,
    iter,
    year,
    all_of(layer_vars_wb)
  ) %>%
  pivot_longer(
    cols = all_of(layer_vars_wb),
    names_to = c("variable", "layer"),
    names_pattern = "([A-Za-z]+)_([0-9]+)",
    values_to = "value"
  ) %>%
  mutate(
    layer = as.integer(layer)
  ) %>%
  group_by(
    scenario,
    year,
    variable,
    layer
  ) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sum_value  = sum(value, na.rm = TRUE),
    min_value  = min(value, na.rm = TRUE),
    max_value  = max(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    last_value = value[which.max(iter)],
    .groups = "drop"
  )

# Save annual water balance data
saveRDS(
  df_water_balance_annual,
  path(cache_dir, "df_water_balance_annual.rds")
)

# ------------------------------------------------------------
# 9. Vertical flux annual data
# ------------------------------------------------------------

# Identify vertical flux variables between soil layers.
# Examples:
# mean_flux_layers0_1
# mean_abs_flux_layers0_1
# gross_volumetric_change_layers0_1
# net_volumetric_change_layers0_1
interface_vars_vf <- names(df_vertical_flux) %>%
  str_subset("_layers[0-9]+_[0-9]+$")

# Stop if no vertical flux interface variables were found
if (length(interface_vars_vf) == 0) {
  stop("No vertical flux interface variables were found.")
}

# Convert vertical flux variables from wide to long format.
# This is the important corrected step.
#
# Example:
# gross_volumetric_change_layers0_1
#
# becomes:
# variable    = gross_volumetric_change
# layer_upper = 0
# layer_lower = 1
# interface   = 0_1
df_vertical_flux_daily_long <- df_vertical_flux %>%
  select(
    scenario,
    iter,
    year,
    all_of(interface_vars_vf)
  ) %>%
  pivot_longer(
    cols = all_of(interface_vars_vf),
    names_to = c("variable", "layer_upper", "layer_lower"),
    names_pattern = "(.+)_layers([0-9]+)_([0-9]+)",
    values_to = "value"
  ) %>%
  mutate(
    layer_upper = as.integer(layer_upper),
    layer_lower = as.integer(layer_lower),
    interface = paste0(layer_upper, "_", layer_lower)
  )

# Aggregate vertical flux data by year.
# The column 'variable' keeps the name of the flux variable.
df_vertical_flux_annual <- df_vertical_flux_daily_long %>%
  group_by(
    scenario,
    year,
    variable,
    layer_upper,
    layer_lower,
    interface
  ) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sum_value  = sum(value, na.rm = TRUE),
    min_value  = min(value, na.rm = TRUE),
    max_value  = max(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    last_value = value[which.max(iter)],
    .groups = "drop"
  )

# Save vertical flux data
saveRDS(
  df_vertical_flux_daily_long,
  path(cache_dir, "df_vertical_flux_daily_long.rds")
)

saveRDS(
  df_vertical_flux_annual,
  path(cache_dir, "df_vertical_flux_annual.rds")
)

# ------------------------------------------------------------
# 10. Biogeochemical annual data
# ------------------------------------------------------------

# Identify biogeochemical variables.
# These are all numeric variables except identifiers.
biogeochemical_vars <- names(df_biogem) %>%
  setdiff(c("scenario", "iter", "year"))

biogeochemical_vars <- biogeochemical_vars[
  map_lgl(df_biogem[biogeochemical_vars], is.numeric)
]

# Convert biogeochemical variables from wide to long format.
# This creates:
# - variable: name of the biogeochemical variable
# - value: value of the variable
df_biogem_annual <- df_biogem %>%
  select(
    scenario,
    iter,
    year,
    all_of(biogeochemical_vars)
  ) %>%
  pivot_longer(
    cols = all_of(biogeochemical_vars),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(
    scenario,
    year,
    variable
  ) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sum_value  = sum(value, na.rm = TRUE),
    min_value  = min(value, na.rm = TRUE),
    max_value  = max(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    last_value = value[which.max(iter)],
    .groups = "drop"
  )

# Save annual biogeochemical data
saveRDS(
  df_biogem_annual,
  path(cache_dir, "df_biogeochemical_annual.rds")
)

# ------------------------------------------------------------
# 11. Simple checks
# ------------------------------------------------------------

# Check the variables saved in the vertical flux annual dataframe
print("Vertical flux variables:")
print(unique(df_vertical_flux_annual$variable))

# Check the interfaces saved in the vertical flux annual dataframe
print("Vertical flux interfaces:")
print(unique(df_vertical_flux_annual$interface))

# Check the first rows of the annual vertical flux dataframe
print("Preview of df_vertical_flux_annual:")
print(head(df_vertical_flux_annual))

# Check where the files were saved
print("RDS files saved in:")
print(cache_dir)


# # Loading required packages
# library(tidyverse)
# library(data.table)
# library(fs)
# 
# project_dir <- path_expand(
#   "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke"
# )
# 
# runs <- tibble(
#   scenario = c(
#     "ducke_deep_clayey_regclim_ks12",
#     "ducke_deep_clayey_redprec",
#     "ducke_shallow_sandy_regclim_ks12",
#     "ducke_shallow_sandy_redprec"
#   ),
#   output_dir = path(
#     project_dir,
#     c(
#       "deepWT_clayeysoil/ducke_deep_clayey_regclim_ks12/output",
#       "deepWT_clayeysoil/ducke_deep_clayey_redprec/output",
#       "shallowWT_sandysoil/ducke_shallow_sandy_regclim_ks12/output",
#       "shallowWT_sandysoil/ducke_shallow_sandy_redprec/output"
#     )
#   )
# )
# read_simulation_data <- function(scenario, output_dir) {
#   
#   file_wb   <- list.files(output_dir, pattern = "water_balance\\.txt$", full.names = TRUE)[1]
#   file_vf   <- list.files(output_dir, pattern = "vertical_water_flux\\.txt$", full.names = TRUE)[1]
#   file_biog <- list.files(output_dir, pattern = "sumstats\\.txt$", full.names = TRUE)[1]
#   
#   if (any(is.na(c(file_wb, file_vf, file_biog)))) {
#     stop("Missing output files in: ", output_dir)
#   }
#   
#   list(
#     wb = read_tsv(file_wb, show_col_types = FALSE) %>%
#       mutate(scenario = scenario, year = floor(iter / 365) + 1),
#     
#     vf = read_tsv(file_vf, show_col_types = FALSE) %>%
#       mutate(scenario = scenario, year = floor(iter / 365) + 1),
#     
#     biog = read_tsv(file_biog, show_col_types = FALSE) %>%
#       mutate(scenario = scenario, year = floor(iter / 365) + 1)
#   )
# }
# 
# # Apply the function to all folders and bind them together
# # map extracts the list, map_dfr binds the rows together into a single dataframe
# all_data <- map2(runs$scenario, runs$output_dir, read_simulation_data)
# 
# df_water_balance <- map_dfr(all_data, "wb")
# df_vertical_flux <- map_dfr(all_data, "vf")
# df_biogem <- map_dfr(all_data, "biog")
# 
# rm(all_data)
# 
# # Identify groups of variables available in the outputs
# 
# wb_vars <- names(df_water_balance)
# vf_vars <- names(df_vertical_flux)
# biog_vars <- names(df_biogem)
# 
# # Layer variables in water balance:
# # SWC_i, SWP_i, transpiration_i
# layer_vars_wb <- wb_vars %>%
#   str_subset("^(SWC|SWP|transpiration)_[0-9]+$")
# 
# # Whole-system water balance variables
# water_balance_vars <- wb_vars %>%
#   setdiff(c("scenario", "year")) %>%
#   setdiff(layer_vars_wb)
# 
# # Interface variables in vertical flux
# interface_vars_vf <- vf_vars %>%
#   str_subset("layers[0-9]+_[0-9]+")
# 
# # Select all biogeochemical variables, excluding identifiers
# biogeochemical_vars <- biog_vars %>%
#   setdiff(c(
#     "iter",
#     "scenario",
#     "year"
#   ))
# 
# # Discomment if you have new data
# #Saving data as .rds
# cache_dir <- path(project_dir, "_rds")
# # dir_create(cache_dir)
# saveRDS(df_water_balance, path(cache_dir, "df_water_balance_daily.rds"))
# saveRDS(df_vertical_flux, path(cache_dir, "df_vertical_flux_daily.rds"))
# saveRDS(df_biogem, path(cache_dir, "df_biogem_daily.rds"))
# 
# 
# 
# ##################
# # Aggregate by year
# 
# # Identify layer variables in the water balance output
# layer_vars_wb <- names(df_water_balance) %>%
#   str_subset("^(SWC|SWP|transpiration)_[0-9]+$")
# layer_vars_wb
# 
# # Convert layer variables to long format and aggregate by year
# df_layers_annual <- df_water_balance %>%
#   select(
#     scenario,
#     iter,
#     year,
#     all_of(layer_vars_wb)
#   ) %>%
#   pivot_longer(
#     cols = all_of(layer_vars_wb),
#     names_to = c("variable", "layer"),
#     names_pattern = "([A-Za-z]+)_([0-9]+)",
#     values_to = "value"
#   ) %>%
#   mutate(
#     layer = as.integer(layer)
#   ) %>%
#   group_by(scenario, year, variable, layer) %>%
#   summarise(
#     mean_value = mean(value, na.rm = TRUE),
#     min_value  = min(value, na.rm = TRUE),
#     max_value  = max(value, na.rm = TRUE),
#     sd_value   = sd(value, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# saveRDS(df_layers_annual, path(cache_dir, "df_water_balance_annual.rds"))
# 
# # Identify interface variables in the vertical flux output
# interface_vars_vf <- names(df_vertical_flux) %>%
#   str_subset("layers[0-9]+_[0-9]+")
# 
# interface_vars_vf
# 
# # Convert interface variables to long format
# df_interfaces_daily <- df_vertical_flux %>%
#   select(
#     scenario,
#     iter,
#     year,
#     all_of(interface_vars_vf)
#   ) %>%
#   pivot_longer(
#     cols = all_of(interface_vars_vf),
#     names_to = c("metric", "layer_upper", "layer_lower"),
#     names_pattern = "(.+)_layers([0-9]+)_([0-9]+)",
#     values_to = "value"
#   ) %>%
#   mutate(
#     layer_upper = as.integer(layer_upper),
#     layer_lower = as.integer(layer_lower),
#     interface = paste0(layer_upper, "_", layer_lower)
#   )
# 
# # Aggregate by year
# df_interfaces_annual <- df_interfaces_daily %>%
#   group_by(
#     scenario,
#     year,
#     metric,
#     layer_upper,
#     layer_lower,
#     interface
#   ) %>%
#   summarise(
#     mean_value = mean(value, na.rm = TRUE),
#     sum_value  = sum(value, na.rm = TRUE),
#     min_value  = min(value, na.rm = TRUE),
#     max_value  = max(value, na.rm = TRUE),
#     sd_value   = sd(value, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# saveRDS(df_interfaces_annual, path(cache_dir, "df_vertical_flux_annual.rds"))
# 
# df_biogem_annual <- df_biogem %>%
#   select(
#     scenario,
#     iter,
#     year,
#     all_of(biogeochemical_vars)
#   ) %>%
#   pivot_longer(
#     cols = all_of(biogeochemical_vars),
#     names_to = "variable",
#     values_to = "value"
#   ) %>%
#   group_by(
#     scenario,
#     year,
#     variable
#   ) %>%
#   summarise(
#     mean_value = mean(value, na.rm = TRUE),
#     sum_value  = sum(value, na.rm = TRUE),
#     min_value  = min(value, na.rm = TRUE),
#     max_value  = max(value, na.rm = TRUE),
#     sd_value   = sd(value, na.rm = TRUE),
#     last_value = value[which.max(iter)],
#     .groups = "drop"
#   )
# 
# saveRDS(df_biogem_annual, path(cache_dir, "df_biogeochemical_annual.rds"))
# 
