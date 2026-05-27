# Loading required packages
library(tidyverse)
library(data.table) 
library(fs) 
base_dir <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/ksfloor/veg/deep_WT/sandy/longterm"

# List simulation folder names
folders <- c(
  "longterm_deep_sandy_ks1e_14_redprec"
)

read_simulation_data <- function(folder) {
  
  ks_string <- str_extract(folder, "1e_[0-9]+")
  ks_numeric <- as.numeric(str_replace(ks_string, "_", "-"))
  
  folder_path <- path(base_dir, folder, "output")
  
  file_wb <- list.files(folder_path, pattern = "water_balance\\.txt$", full.names = TRUE)
  file_vf <- list.files(folder_path, pattern = "vertical_water_flux\\.txt$", full.names = TRUE)
  file_biog <- list.files(folder_path, pattern = "sumstats\\.txt$", full.names = TRUE)
  
  if (length(file_wb) == 0 || length(file_vf) == 0 || length(file_biog) == 0) {
    stop(paste("Could not find the output files inside:", folder_path))
  }
  
  df_wb <- read_tsv(file_wb[1], show_col_types = FALSE) %>%
    mutate(
      scenario = folder,
      ksfloor_val = ks_numeric,
      ksfloor_factor = ks_string,
      year = floor(iter / 365) + 1
    )
  
  df_vf <- read_tsv(file_vf[1], show_col_types = FALSE) %>%
    mutate(
      scenario = folder,
      ksfloor_val = ks_numeric,
      ksfloor_factor = ks_string,
      year = floor(iter / 365) + 1
    )
  
  df_biog <- read_tsv(file_biog[1], show_col_types = FALSE) %>%
    mutate(
      scenario = folder,
      ksfloor_val = ks_numeric,
      ksfloor_factor = ks_string,
      year = floor(iter / 365) + 1
    )
  
  return(list(wb = df_wb, vf = df_vf, biog = df_biog))
}

# Apply the function to all folders and bind them together
# map extracts the list, map_dfr binds the rows together into a single dataframe
all_data <- map(folders, read_simulation_data)

# Creating the two main dataframes we will use in the analysis
df_water_balance <- map_dfr(all_data, "wb")
df_vertical_flux <- map_dfr(all_data, "vf")
df_biogem <- map_dfr(all_data, "biog")

# Clearing memory
rm(all_data)


# Identify groups of variables available in the outputs

wb_vars <- names(df_water_balance)
vf_vars <- names(df_vertical_flux)
biog_vars <- names(df_biogem)

# Layer variables in water balance:
# SWC_i, SWP_i, transpiration_i
layer_vars_wb <- wb_vars %>%
  str_subset("^(SWC|SWP|transpiration)_[0-9]+$")

# Whole-system water balance variables
water_balance_vars <- wb_vars %>%
  setdiff(c("scenario", "ksfloor_val", "ksfloor_factor", "year")) %>%
  setdiff(layer_vars_wb)

# Interface variables in vertical flux
interface_vars_vf <- vf_vars %>%
  str_subset("layers[0-9]+_[0-9]+")

# Select all biogeochemical variables, excluding identifiers
biogeochemical_vars <- biog_vars %>%
  setdiff(c(
    "iter",
    "scenario",
    "ksfloor_val",
    "ksfloor_factor",
    "year"
  ))

## Discomment if you have new data
#Saving data as .rds
# cache_dir <- path(base_dir, "_rds")
# dir_create(cache_dir)
# saveRDS(df_water_balance, path(cache_dir, "df_water_balance_daily.rds"))
# saveRDS(df_vertical_flux, path(cache_dir, "df_vertical_flux_daily.rds"))
# saveRDS(df_biogem, path(cache_dir, "df_biogem_daily.rds"))



##################
# Aggregate by year

# Identify layer variables in the water balance output
layer_vars_wb <- names(df_water_balance) %>%
  str_subset("^(SWC|SWP|transpiration)_[0-9]+$")
layer_vars_wb

# Convert layer variables to long format and aggregate by year
df_layers_annual <- df_water_balance %>%
  select(
    scenario,
    ksfloor_val,
    ksfloor_factor,
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
  group_by(scenario, ksfloor_val, ksfloor_factor, year, variable, layer) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    min_value  = min(value, na.rm = TRUE),
    max_value  = max(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(df_layers_annual, path(cache_dir, "df_water_balance_annual.rds"))

# Identify interface variables in the vertical flux output
interface_vars_vf <- names(df_vertical_flux) %>%
  str_subset("layers[0-9]+_[0-9]+")

interface_vars_vf

# Convert interface variables to long format
df_interfaces_daily <- df_vertical_flux %>%
  select(
    scenario,
    ksfloor_val,
    ksfloor_factor,
    iter,
    year,
    all_of(interface_vars_vf)
  ) %>%
  pivot_longer(
    cols = all_of(interface_vars_vf),
    names_to = c("metric", "layer_upper", "layer_lower"),
    names_pattern = "(.+)_layers([0-9]+)_([0-9]+)",
    values_to = "value"
  ) %>%
  mutate(
    layer_upper = as.integer(layer_upper),
    layer_lower = as.integer(layer_lower),
    interface = paste0(layer_upper, "_", layer_lower)
  )

# Aggregate by year
df_interfaces_annual <- df_interfaces_daily %>%
  group_by(
    scenario,
    ksfloor_val,
    ksfloor_factor,
    year,
    metric,
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
    .groups = "drop"
  )

saveRDS(df_interfaces_annual, path(cache_dir, "df_vertical_flux_annual.rds"))

df_biogem_annual <- df_biogem %>%
  select(
    scenario,
    ksfloor_val,
    ksfloor_factor,
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
    ksfloor_val,
    ksfloor_factor,
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

saveRDS(df_biogem_annual, path(cache_dir, "df_biogeochemical_annual.rds"))

