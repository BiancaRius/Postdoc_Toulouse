# ============================================================
# Load and aggregate TROLL output files
# ============================================================

# Load required packages
library(tidyverse)
library(fs)

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
    "ducke_shallow_sandy_regclim_ks14",
    "ducke_shallow_sandy_redprec30_ks14"
  ),
  output_dir = path(
    project_dir,
    c(
      "deepWT_clayeysoil/ducke_deep_clayey_regclim_ks14/output",
      "deepWT_clayeysoil/ducke_deep_clayey_redprec30_ks14/output",
      "deepWT_sandysoil/ducke_deep_sandy_redprec30_ks14/output",
      "shallowWT_sandysoil/ducke_shallow_sandy_regclim_ks14/output",
      "shallowWT_sandysoil/ducke_shallow_sandy_redprec30_ks14/output"
    )
  )
)

# ------------------------------------------------------------
# 3. Helper function to safely find one output file
# ------------------------------------------------------------

find_output_file <- function(output_dir, pattern, file_label) {
  
  # Look for files matching the expected pattern
  files <- list.files(
    output_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  # Stop the script if the file is missing
  if (length(files) == 0 || !file_exists(files[1])) {
    stop("Missing ", file_label, " file in: ", output_dir)
  }
  
  # Return the first matching file
  files[1]
}

# ------------------------------------------------------------
# 4. Function to read one simulation folder
# ------------------------------------------------------------

read_simulation_data <- function(scenario, output_dir) {
  
  # Find the three main output files
  file_wb <- find_output_file(
    output_dir,
    pattern = "water_balance\\.txt$",
    file_label = "water balance"
  )
  
  file_vf <- find_output_file(
    output_dir,
    pattern = "vertical_water_flux\\.txt$",
    file_label = "vertical water flux"
  )
  
  file_biog <- find_output_file(
    output_dir,
    pattern = "sumstats\\.txt$",
    file_label = "biogeochemical"
  )
  
  # Read water balance data
  wb <- read_tsv(file_wb, show_col_types = FALSE) %>%
    mutate(
      scenario = scenario,
      iter = as.numeric(iter),
      year = floor(iter / 365) + 1
    )
  
  # Read vertical flux data
  vf <- read_tsv(file_vf, show_col_types = FALSE) %>%
    mutate(
      scenario = scenario,
      iter = as.numeric(iter),
      year = floor(iter / 365) + 1
    )
  
  # Read biogeochemical data
  biog <- read_tsv(file_biog, show_col_types = FALSE) %>%
    mutate(
      scenario = scenario,
      iter = as.numeric(iter),
      year = floor(iter / 365) + 1
    )
  
  # Return the three dataframes as a list
  list(
    wb = wb,
    vf = vf,
    biog = biog
  )
}

# ------------------------------------------------------------
# 5. Read all simulations
# ------------------------------------------------------------

# Apply the reading function to all scenarios
all_data <- map2(
  runs$scenario,
  runs$output_dir,
  read_simulation_data
)

# Bind all scenarios into three daily dataframes
df_water_balance <- map_dfr(all_data, "wb")
df_vertical_flux <- map_dfr(all_data, "vf")
df_biogem        <- map_dfr(all_data, "biog")

# Remove temporary list to save memory
rm(all_data)

# ------------------------------------------------------------
# 6. Create cache directory and save daily data
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
# 7. Water balance: layer variables
# ------------------------------------------------------------

# Identify variables with layer information:
# Examples: SWC_0, SWP_0, transpiration_0
layer_vars_wb <- names(df_water_balance) %>%
  str_subset("^(SWC|SWP|transpiration)_[0-9]+$")

# Convert layer variables from wide to long format
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

saveRDS(
  df_water_balance_annual,
  path(cache_dir, "df_water_balance_annual.rds")
)

# ------------------------------------------------------------
# 8. Vertical flux: interface variables
# ------------------------------------------------------------

# Identify vertical flux variables that contain soil interface information.
# Examples:
# mean_flux_layers0_1
# mean_abs_flux_layers0_1
# mean_delta_swp_layers0_1
# mean_ks_harmonic_layers0_1
#
# The important part is that these names have:
# <variable name>_layers<upper layer>_<lower layer>
interface_vars_vf <- names(df_vertical_flux) %>%
  str_subset("_layers[0-9]+_[0-9]+$")

# Stop if no interface variables were found
if (length(interface_vars_vf) == 0) {
  stop(
    "No vertical flux interface variables were found. ",
    "Check the column names in df_vertical_flux."
  )
}

# Convert interface variables from wide to long format.
# This is the corrected part:
# - variable stores the name of the flux variable
# - layer_upper stores the upper soil layer
# - layer_lower stores the lower soil layer
df_vertical_flux_interfaces_daily <- df_vertical_flux %>%
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
    interface = paste0(layer_upper, "_", layer_lower),
    
    # Keep this extra column because some plotting functions
    # may still expect a column called metric.
    metric = variable,
    
    # This identifies that these variables are layer-interface fluxes.
    flux_group = "interface"
  )

# ------------------------------------------------------------
# 9. Vertical flux: non-interface variables, if they exist
# ------------------------------------------------------------

# Some vertical flux outputs may also contain variables without layers.
# Examples could be gross_volumetric_change, net_volumetric_change,
# directionality, or similar variables.
other_vf_vars <- names(df_vertical_flux) %>%
  setdiff(c("scenario", "iter", "year", interface_vars_vf))

# Keep only numeric non-interface variables
other_vf_vars <- other_vf_vars[
  map_lgl(df_vertical_flux[other_vf_vars], is.numeric)
]

# Convert non-interface vertical flux variables to long format.
# If there are no such variables, create an empty dataframe.
if (length(other_vf_vars) > 0) {
  
  df_vertical_flux_other_daily <- df_vertical_flux %>%
    select(
      scenario,
      iter,
      year,
      all_of(other_vf_vars)
    ) %>%
    pivot_longer(
      cols = all_of(other_vf_vars),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      layer_upper = NA_integer_,
      layer_lower = NA_integer_,
      interface = NA_character_,
      metric = variable,
      flux_group = "whole_profile"
    )
  
} else {
  
  df_vertical_flux_other_daily <- tibble(
    scenario = character(),
    iter = numeric(),
    year = numeric(),
    variable = character(),
    value = numeric(),
    layer_upper = integer(),
    layer_lower = integer(),
    interface = character(),
    metric = character(),
    flux_group = character()
  )
}

# ------------------------------------------------------------
# 10. Combine all vertical flux variables and aggregate by year
# ------------------------------------------------------------

# Combine interface and non-interface vertical flux variables
df_vertical_flux_daily_long <- bind_rows(
  df_vertical_flux_interfaces_daily,
  df_vertical_flux_other_daily
)

# Aggregate vertical flux data by year
df_vertical_flux_annual <- df_vertical_flux_daily_long %>%
  group_by(
    scenario,
    year,
    variable,
    metric,
    flux_group,
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

saveRDS(
  df_vertical_flux_daily_long,
  path(cache_dir, "df_vertical_flux_daily_long.rds")
)

saveRDS(
  df_vertical_flux_annual,
  path(cache_dir, "df_vertical_flux_annual.rds")
)

# ------------------------------------------------------------
# 11. Biogeochemical variables
# ------------------------------------------------------------

# Identify all biogeochemical variables, excluding identifiers
biogeochemical_vars <- names(df_biogem) %>%
  setdiff(c("scenario", "iter", "year"))

# Keep only numeric variables
biogeochemical_vars <- biogeochemical_vars[
  map_lgl(df_biogem[biogeochemical_vars], is.numeric)
]

# Convert biogeochemical variables from wide to long format
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

saveRDS(
  df_biogem_annual,
  path(cache_dir, "df_biogeochemical_annual.rds")
)

# ------------------------------------------------------------
# 12. Quick checks
# ------------------------------------------------------------

# Check which vertical flux variables were saved
print("Vertical flux variables saved:")
print(unique(df_vertical_flux_annual$variable))

# Check the first rows of the annual vertical flux dataframe
print("Preview of annual vertical flux dataframe:")
print(head(df_vertical_flux_annual))

# Check output files
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
