# ============================================================
# Load and aggregate TROLL output files
# Mimicking Reserva Ducke Topography (shallow WT sandy soil & 
# deep WT clayey soil)
# ============================================================

# Load required packages
library(tidyverse)
library(fs)
library(dplyr)

# ------------------------------------------------------------
# 1. Define the main project directory
# ------------------------------------------------------------

project_dir <- path_expand(
  "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_ducke_corrected_units/"
)

# ------------------------------------------------------------
# 2. Define the simulation scenarios and their output folders
# ------------------------------------------------------------

runs <- tibble(
  scenario = c(
    "deepWT_clayey_regclim",
    "deepWT_clayey_redprec",
    "shallowWT_sandy_regclim",
    "shallowWT_sandy_redprec"
  ),
  output_dir = path(
    project_dir,
    c(
      "deepWT_clayey_regclim",
      "deepWT_clayey_redprec",
      "shallowWT_sandy_regclim",
      "shallowWT_sandy_redprec"
    ),
    "output"
  )
)

# ------------------------------------------------------------
# 3. Function to find one required output file
# ------------------------------------------------------------

find_output_file <- function(output_dir, pattern, file_label) {
  
  files <- list.files(
    output_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0 || !file_exists(files[1])) {
    stop("Missing ", file_label, " file in: ", output_dir)
  }
  
  files[1]
}

# ------------------------------------------------------------
# 4. Function to find one optional output file
# ------------------------------------------------------------

find_optional_output_file <- function(output_dir, pattern) {
  
  files <- list.files(
    output_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0 || !file_exists(files[1])) {
    return(NA_character_)
  }
  
  files[1]
}

# ------------------------------------------------------------
# 5. Function to remove empty columns created by read_tsv()
# ------------------------------------------------------------

clean_output_columns <- function(df) {
  
  df %>%
    select(-matches("^\\.\\.\\.[0-9]+$"))
}

# ------------------------------------------------------------
# 6. Function to read one simulation folder
# ------------------------------------------------------------

read_simulation_data <- function(scenario, output_dir) {
  
  message("Reading scenario: ", scenario)
  
  # Required file: water balance
  file_wb <- find_output_file(
    output_dir = output_dir,
    pattern = "water_balance\\.txt$",
    file_label = "water balance"
  )
  
  # Optional file: vertical water flux
  # This file does not exist for the original TROLL code because
  # vertical water flux / capillary flux diagnostics were not represented.
  file_vf <- find_optional_output_file(
    output_dir = output_dir,
    pattern = "vertical_water_flux\\.txt$"
  )
  
  # Required file: biogeochemical summary
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
  
  # Read vertical flux output only if available
  if (!is.na(file_vf)) {
    
    vf <- read_tsv(file_vf, show_col_types = FALSE) %>%
      clean_output_columns() %>%
      mutate(
        iter = as.numeric(iter),
        scenario = scenario,
        year = floor(iter / 365) + 1
      )
    
  } else {
    
    message("No vertical water flux file for scenario: ", scenario)
    vf <- NULL
    
  }
  
  # Read biogeochemical output
  biog <- read_tsv(file_biog, show_col_types = FALSE) %>%
    clean_output_columns() %>%
    mutate(
      iter = as.numeric(iter),
      scenario = scenario,
      year = floor(iter / 365) + 1
    )
  
  list(
    wb = wb,
    vf = vf,
    biog = biog
  )
}

# ------------------------------------------------------------
# 7. Read all simulations
# ------------------------------------------------------------

all_data <- map2(
  runs$scenario,
  runs$output_dir,
  read_simulation_data
)

# Bind all scenarios into daily dataframes
df_water_balance <- map_dfr(all_data, "wb")
df_biogem        <- map_dfr(all_data, "biog")

# Bind vertical flux only for scenarios where it exists
df_vertical_flux <- map(all_data, "vf") %>%
  compact() %>%
  bind_rows()

# ------------------------------------------------------------
# Limit shallow WT scenarios to the first 150 years
# ------------------------------------------------------------

target_scenarios <- c(
  "shallowWT_sandy_regclim",
  "shallowWT_sandy_redprec"
)

max_year_target <- 150

keep_first_years_target <- function(df, target_scenarios, max_year_target) {
  
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  
  df %>%
    filter(
      !scenario %in% target_scenarios | year <= max_year_target
    )
}

df_water_balance <- keep_first_years_target(
  df_water_balance,
  target_scenarios,
  max_year_target
)

df_biogem <- keep_first_years_target(
  df_biogem,
  target_scenarios,
  max_year_target
)

df_vertical_flux <- keep_first_years_target(
  df_vertical_flux,
  target_scenarios,
  max_year_target
)
# Remove temporary object
rm(all_data)

# ------------------------------------------------------------
# 8. Create cache directory and save daily dataframes
# ------------------------------------------------------------

cache_dir <- path(project_dir, "_rds")
dir_create(cache_dir)

saveRDS(
  df_water_balance,
  path(cache_dir, "df_water_balance_daily.rds")
)

saveRDS(
  df_biogem,
  path(cache_dir, "df_biogem_daily.rds")
)

# Save vertical flux daily dataframe only if it exists
if (nrow(df_vertical_flux) > 0) {
  
  saveRDS(
    df_vertical_flux,
    path(cache_dir, "df_vertical_flux_daily.rds")
  )
  
} else {
  
  message("No vertical flux data available in any scenario.")
  
}

# ============================================================
# Aggregate annual data
# ============================================================

# ------------------------------------------------------------
# 9. Water balance annual data: layer-specific variables
# ------------------------------------------------------------

# Identify water balance variables that are specific to soil layers.
# Examples:
# SWC_0, SWC_1, SWP_0, SWP_1, transpiration_0, transpiration_1
layer_vars_wb <- names(df_water_balance) %>%
  str_subset("^(SWC|SWP|transpiration)_[0-9]+$")

if (length(layer_vars_wb) == 0) {
  stop("No layer-specific water balance variables were found.")
}

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
# 10. Water balance annual data: non-layer variables
# ------------------------------------------------------------

# This keeps variables such as precipitation, runoff, leak,
# evaporation, total transpiration, etc., if they exist.
water_balance_other_vars <- names(df_water_balance) %>%
  setdiff(c("scenario", "iter", "year", layer_vars_wb))

water_balance_other_vars <- water_balance_other_vars[
  map_lgl(df_water_balance[water_balance_other_vars], is.numeric)
]

if (length(water_balance_other_vars) > 0) {
  
  df_water_balance_other_annual <- df_water_balance %>%
    select(
      scenario,
      iter,
      year,
      all_of(water_balance_other_vars)
    ) %>%
    pivot_longer(
      cols = all_of(water_balance_other_vars),
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
    df_water_balance_other_annual,
    path(cache_dir, "df_water_balance_other_annual.rds")
  )
  
} else {
  
  message("No non-layer water balance variables were found.")
  
}

# ------------------------------------------------------------
# 11. Vertical flux annual data
# ------------------------------------------------------------

# Initialize objects as NULL in case vertical flux does not exist
df_vertical_flux_daily_long <- NULL
df_vertical_flux_annual <- NULL

if (nrow(df_vertical_flux) > 0) {
  
  # Identify vertical flux variables between soil layers.
  # Examples:
  # mean_flux_layers0_1
  # mean_abs_flux_layers0_1
  # gross_volumetric_change_layers0_1
  # net_volumetric_change_layers0_1
  interface_vars_vf <- names(df_vertical_flux) %>%
    str_subset("_layers[0-9]+_[0-9]+$")
  
  if (length(interface_vars_vf) > 0) {
    
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
    
    saveRDS(
      df_vertical_flux_daily_long,
      path(cache_dir, "df_vertical_flux_daily_long.rds")
    )
    
    saveRDS(
      df_vertical_flux_annual,
      path(cache_dir, "df_vertical_flux_annual.rds")
    )
    
  } else {
    
    message("Vertical flux file exists, but no interface variables were found.")
    
  }
  
} else {
  
  message("Skipping vertical flux annual aggregation because no vertical flux file was found.")
  
}

# ------------------------------------------------------------
# 12. Biogeochemical annual data
# ------------------------------------------------------------

# Identify biogeochemical variables.
# These are all numeric variables except identifiers.
biogeochemical_vars <- names(df_biogem) %>%
  setdiff(c("scenario", "iter", "year"))

biogeochemical_vars <- biogeochemical_vars[
  map_lgl(df_biogem[biogeochemical_vars], is.numeric)
]

if (length(biogeochemical_vars) == 0) {
  stop("No numeric biogeochemical variables were found.")
}

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
# 13. Simple checks
# ------------------------------------------------------------

print("Scenarios in water balance data:")
print(unique(df_water_balance$scenario))

print("Scenarios in biogeochemical data:")
print(unique(df_biogem$scenario))

print("Water balance layer variables:")
print(unique(df_water_balance_annual$variable))

if (exists("df_water_balance_other_annual")) {
  print("Water balance non-layer variables:")
  print(unique(df_water_balance_other_annual$variable))
}

print("Biogeochemical variables:")
print(unique(df_biogem_annual$variable))

if (!is.null(df_vertical_flux_annual)) {
  
  print("Scenarios in vertical flux data:")
  print(unique(df_vertical_flux_annual$scenario))
  
  print("Vertical flux variables:")
  print(unique(df_vertical_flux_annual$variable))
  
  print("Vertical flux interfaces:")
  print(unique(df_vertical_flux_annual$interface))
  
  print("Preview of df_vertical_flux_annual:")
  print(head(df_vertical_flux_annual))
  
} else {
  
  print("No vertical flux annual dataframe was created.")
  
}

print("RDS files saved in:")
print(cache_dir)
