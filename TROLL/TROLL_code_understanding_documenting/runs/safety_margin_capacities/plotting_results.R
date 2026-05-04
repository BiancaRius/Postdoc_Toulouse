# =========================================================
# TROLL biogeochemical, transpiration and water-flux comparison script
# =========================================================
# Purpose:
# This script compares TROLL simulations across different run branches
# using an experiment_table structure.
#
# It reads:
# - (null)_0_sumstats.txt        for biogeochemical variables
# - (null)_0_water_balance.txt   for transpiration and water-flux variables
#
# It plots:
# - biogeochemical variables: npp, gpp, agb, sum1, sum10, sum30, ba, ba10, litterfall
# - transpiration layers: transpitation_0 to transpitation_5
# - total transpiration across soil layers
# - water flux variables: precipitation, interception, throughfall, runoff, leak, evaporation
# - cumulative water fluxes
# - annual water-flux totals
#
# Important:
# The water_balance file is loaded only once and reused for both
# transpiration and water-flux analyses.
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
})

# =========================================================
# 1) USER SETTINGS
# =========================================================

project_dir <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting"

out_prefix <- "(null)"
start_date <- as.Date("2004-01-01")

# Variables to analyze from sumstats
biogeochemical_vars <- c(
  "npp", "gpp", "agb",
  "sum1", "sum10", "sum30",
  "ba", "ba10",
  "litterfall"
)

# Important: your current water_balance output seems to use "transpitation", not "transpiration".
transpiration_layers <- c(
  "transpitation_0",
  "transpitation_1",
  "transpitation_2",
  "transpitation_3",
  "transpitation_4",
  "transpitation_5"
)

# Water-balance flux variables to analyze
water_flux_vars <- c(
  "precipitation",
  "interception",
  "throughfall",
  "runoff",
  "leak",
  "evaporation"
)

# =========================================================
# 2) EXPERIMENT TABLE
# =========================================================

experiment_table <- tibble(
  family = c(
    "safety_margin_capacities", "safety_margin_capacities",
    "safety_margin_capacities", "safety_margin_capacities",
    "safety_margin_capacities", "safety_margin_capacities",
    "safety_margin_capacities", "safety_margin_capacities"
  ),
  vegetation = c(
    "veg", "veg",
    "veg", "veg",
    "veg", "veg",
    "veg", "veg"
  ),
  water_table = c(
    "shallowWT", "shallowWT",
    "deepWT", "deepWT",
    "shallowWT", "shallowWT",
    "deepWT", "deepWT"
  ),
  texture = c(
    "sandy", "sandy",
    "sandy", "sandy",
    "clayey", "clayey",
    "clayey", "clayey"
  ),
  experiment_name = c(
    "safetymargin_shallowWT_sandy_regclim",
    "safetymargin_shallowWT_sandy_redprec",
    "safetymargin_deepWT_sandy_regclim",
    "safetymargin_deepWT_sandy_redprec",
    "safetymargin_shallowWT_clayey_regclim",
    "safetymargin_shallowWT_clayey_redprec",
    "safetymargin_deepWT_clayey_regclim",
    "safetymargin_deepWT_clayey_redprec"
  )
) %>%
  mutate(
    scenario = case_when(
      str_detect(experiment_name, "^safetymargin_") ~
        str_replace(experiment_name, "^safetymargin_", "margin_"),
      TRUE ~ experiment_name
    ),
    run_dir = file.path(
      project_dir,
      "runs",
      family,
      vegetation,
      water_table,
      texture,
      experiment_name
    ),
    output_dir = file.path(run_dir, "output")
  )

# Optional: control plotting order
experiment_table <- experiment_table %>%
  mutate(
    scenario = factor(
      scenario,
      levels = scenario
    )
  )

# =========================================================
# 3) HELPER FUNCTIONS
# =========================================================

check_experiment_files <- function(experiment_table, file_type, out_prefix = "(null)") {
  
  files_df <- experiment_table %>%
    mutate(
      file_path = file.path(
        output_dir,
        paste0(out_prefix, "_0_", file_type, ".txt")
      )
    )
  
  missing_files <- files_df %>%
    filter(!file.exists(file_path))
  
  if (nrow(missing_files) > 0) {
    warning(
      "Missing files for file type '", file_type, "':\n",
      paste(missing_files$file_path, collapse = "\n")
    )
  }
  
  invisible(files_df)
}

read_troll_file <- function(output_dir, file_type, scenario, out_prefix = "(null)") {
  
  file_path <- file.path(
    output_dir,
    paste0(out_prefix, "_0_", file_type, ".txt")
  )
  
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }
  
  df <- tryCatch(
    {
      x <- read.table(
        file_path,
        header = TRUE,
        sep = "",
        check.names = FALSE,
        fill = TRUE
      )
      
      # Remove empty unnamed columns if present
      empty_idx <- which(names(x) == "")
      if (length(empty_idx) > 0) {
        names(x)[empty_idx] <- paste0("empty_col_", seq_along(empty_idx))
      }
      
      x <- x %>%
        select(-starts_with("empty_col_"))
      
      x$scenario <- as.character(scenario)
      x$file_type <- file_type
      x
    },
    error = function(e) {
      warning("Failed to read: ", file_path, " -> ", conditionMessage(e))
      NULL
    }
  )
  
  df
}

load_all_runs <- function(experiment_table,
                          file_type,
                          out_prefix = "(null)",
                          start_date = as.Date("2004-01-01")) {
  
  check_experiment_files(
    experiment_table = experiment_table,
    file_type = file_type,
    out_prefix = out_prefix
  )
  
  dfs <- purrr::map2(
    experiment_table$output_dir,
    experiment_table$scenario,
    ~ read_troll_file(
      output_dir = .x,
      file_type = file_type,
      scenario = .y,
      out_prefix = out_prefix
    )
  )
  
  all_data <- bind_rows(dfs[!sapply(dfs, is.null)])
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No valid data found for file type: ", file_type)
  }
  
  if (!"iter" %in% names(all_data)) {
    stop("Column 'iter' not found in file type: ", file_type)
  }
  
  all_data %>%
    mutate(
      iter = as.numeric(iter),
      date = start_date + iter,
      sim_year = as.numeric(date - start_date) / 365,
      scenario = factor(
        scenario,
        levels = as.character(experiment_table$scenario)
      )
    )
}

check_requested_variables <- function(data, variables, variable_group, file_type) {
  
  vars_present <- intersect(variables, names(data))
  vars_missing <- setdiff(variables, names(data))
  
  if (length(vars_missing) > 0) {
    warning(
      "These ", variable_group, " variables were not found in ", file_type, ":\n",
      paste(vars_missing, collapse = ", ")
    )
  }
  
  if (length(vars_present) == 0) {
    stop("None of the requested ", variable_group, " variables were found in ", file_type, ".")
  }
  
  vars_present
}

# =========================================================
# 4) DATA PREPARATION
# =========================================================

prepare_biogeochemical_data <- function(experiment_table,
                                        variables = biogeochemical_vars,
                                        out_prefix = "(null)",
                                        start_date = as.Date("2004-01-01")) {
  
  sumstats <- load_all_runs(
    experiment_table = experiment_table,
    file_type = "sumstats",
    out_prefix = out_prefix,
    start_date = start_date
  )
  
  vars_present <- check_requested_variables(
    data = sumstats,
    variables = variables,
    variable_group = "biogeochemical",
    file_type = "sumstats"
  )
  
  sumstats %>%
    mutate(across(all_of(vars_present), as.numeric)) %>%
    select(iter, date, sim_year, scenario, all_of(vars_present)) %>%
    pivot_longer(
      cols = all_of(vars_present),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      variable = factor(variable, levels = vars_present)
    )
}

load_water_balance_data <- function(experiment_table,
                                    out_prefix = "(null)",
                                    start_date = as.Date("2004-01-01")) {
  
  load_all_runs(
    experiment_table = experiment_table,
    file_type = "water_balance",
    out_prefix = out_prefix,
    start_date = start_date
  )
}

prepare_transpiration_data <- function(water_balance,
                                       variables = transpiration_layers) {
  
  vars_present <- check_requested_variables(
    data = water_balance,
    variables = variables,
    variable_group = "transpiration",
    file_type = "water_balance"
  )
  
  water_balance %>%
    mutate(across(all_of(vars_present), as.numeric)) %>%
    select(iter, date, sim_year, scenario, all_of(vars_present)) %>%
    pivot_longer(
      cols = all_of(vars_present),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      layer = str_extract(variable, "\\d+$"),
      layer = factor(layer, levels = as.character(0:9)),
      variable = factor(variable, levels = vars_present)
    )
}

prepare_total_transpiration_data <- function(water_balance,
                                             variables = transpiration_layers) {
  
  vars_present <- check_requested_variables(
    data = water_balance,
    variables = variables,
    variable_group = "transpiration",
    file_type = "water_balance"
  )
  
  water_balance %>%
    mutate(across(all_of(vars_present), as.numeric)) %>%
    mutate(
      total_transpiration = rowSums(
        across(all_of(vars_present)),
        na.rm = TRUE
      )
    ) %>%
    select(iter, date, sim_year, scenario, total_transpiration)
}

prepare_water_flux_data <- function(water_balance,
                                    variables = water_flux_vars) {
  
  vars_present <- check_requested_variables(
    data = water_balance,
    variables = variables,
    variable_group = "water-flux",
    file_type = "water_balance"
  )
  
  water_balance %>%
    mutate(across(all_of(vars_present), as.numeric)) %>%
    select(iter, date, sim_year, scenario, all_of(vars_present)) %>%
    pivot_longer(
      cols = all_of(vars_present),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      variable = factor(variable, levels = vars_present)
    )
}

prepare_cumulative_water_flux_data <- function(water_flux_long) {
  
  water_flux_long %>%
    arrange(scenario, variable, iter) %>%
    group_by(scenario, variable) %>%
    mutate(
      cumulative_value = cumsum(replace_na(value, 0))
    ) %>%
    ungroup()
}

prepare_annual_water_flux_data <- function(water_flux_long) {
  
  water_flux_long %>%
    mutate(
      sim_year_int = floor(sim_year) + 1
    ) %>%
    group_by(scenario, variable, sim_year_int) %>%
    summarise(
      annual_total = sum(value, na.rm = TRUE),
      annual_mean = mean(value, na.rm = TRUE),
      .groups = "drop"
    )
}

# =========================================================
# 5) PLOT FUNCTIONS
# =========================================================

plot_biogeochemical_grid <- function(biogeo_long) {
  
  ggplot(
    biogeo_long,
    aes(x = sim_year, y = value, color = scenario)
  ) +
    geom_line(linewidth = 0.7, alpha = 0.9) +
    facet_wrap(~ variable, scales = "free_y", ncol = 3) +
    labs(
      title = "Biogeochemical variables across experiments",
      x = "Simulation year",
      y = NULL,
      color = "Experiment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

plot_biogeochemical_faceted_by_experiment <- function(biogeo_long) {
  
  ggplot(
    biogeo_long,
    aes(x = sim_year, y = value)
  ) +
    geom_line(linewidth = 0.6) +
    facet_grid(variable ~ scenario, scales = "free_y") +
    labs(
      title = "Biogeochemical variables by experiment",
      x = "Simulation year",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold")
    )
}

plot_transpiration_layers <- function(transp_long) {
  
  ggplot(
    transp_long,
    aes(x = sim_year, y = value, color = scenario)
  ) +
    geom_line(linewidth = 0.7, alpha = 0.9) +
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    labs(
      title = "Transpiration by soil layer",
      x = "Simulation year",
      y = "Transpiration",
      color = "Experiment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

plot_transpiration_layers_by_experiment <- function(transp_long) {
  
  ggplot(
    transp_long,
    aes(x = sim_year, y = value)
  ) +
    geom_line(linewidth = 0.6) +
    facet_grid(layer ~ scenario, scales = "free_y") +
    labs(
      title = "Transpiration by layer and experiment",
      x = "Simulation year",
      y = "Transpiration"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold")
    )
}

plot_total_transpiration <- function(total_transp) {
  
  ggplot(
    total_transp,
    aes(x = sim_year, y = total_transpiration, color = scenario)
  ) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    labs(
      title = "Total transpiration across experiments",
      x = "Simulation year",
      y = "Total transpiration",
      color = "Experiment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
}

plot_total_transpiration_by_experiment <- function(total_transp) {
  
  ggplot(
    total_transp,
    aes(x = sim_year, y = total_transpiration)
  ) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ scenario, scales = "free_y", ncol = 2) +
    labs(
      title = "Total transpiration by experiment",
      x = "Simulation year",
      y = "Total transpiration"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

plot_water_fluxes <- function(water_flux_long) {
  
  ggplot(
    water_flux_long,
    aes(x = sim_year, y = value, color = scenario)
  ) +
    geom_line(linewidth = 0.6, alpha = 0.75) +
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    labs(
      title = "Daily water-balance fluxes across experiments",
      x = "Simulation year",
      y = "Flux",
      color = "Experiment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

plot_water_fluxes_by_experiment <- function(water_flux_long) {
  
  ggplot(
    water_flux_long,
    aes(x = sim_year, y = value)
  ) +
    geom_line(linewidth = 0.5, alpha = 0.75) +
    facet_grid(variable ~ scenario, scales = "free_y") +
    labs(
      title = "Daily water-balance fluxes by experiment",
      x = "Simulation year",
      y = "Flux"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold")
    )
}

plot_cumulative_water_fluxes <- function(cumulative_water_flux) {
  
  ggplot(
    cumulative_water_flux,
    aes(x = sim_year, y = cumulative_value, color = scenario)
  ) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    labs(
      title = "Cumulative water-balance fluxes across experiments",
      x = "Simulation year",
      y = "Cumulative flux",
      color = "Experiment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

plot_annual_water_fluxes <- function(annual_water_flux) {
  
  ggplot(
    annual_water_flux,
    aes(x = sim_year_int, y = annual_total, color = scenario)
  ) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    geom_point(size = 1.3, alpha = 0.9) +
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    labs(
      title = "Annual total water-balance fluxes across experiments",
      x = "Simulation year",
      y = "Annual total flux",
      color = "Experiment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

# =========================================================
# 6) SUMMARY TABLES
# =========================================================

summarise_biogeochemical <- function(biogeo_long) {
  
  biogeo_long %>%
    group_by(scenario, variable) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      min_value = min(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE),
      final_value = value[which.max(iter)],
      .groups = "drop"
    ) %>%
    arrange(variable, scenario)
}

summarise_transpiration_layers <- function(transp_long) {
  
  transp_long %>%
    group_by(scenario, layer, variable) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      min_value = min(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE),
      cumulative_value = sum(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(layer, scenario)
}

summarise_total_transpiration <- function(total_transp) {
  
  total_transp %>%
    group_by(scenario) %>%
    summarise(
      mean_total_transpiration = mean(total_transpiration, na.rm = TRUE),
      median_total_transpiration = median(total_transpiration, na.rm = TRUE),
      min_total_transpiration = min(total_transpiration, na.rm = TRUE),
      max_total_transpiration = max(total_transpiration, na.rm = TRUE),
      cumulative_total_transpiration = sum(total_transpiration, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scenario)
}

summarise_water_fluxes <- function(water_flux_long) {
  
  water_flux_long %>%
    group_by(scenario, variable) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      min_value = min(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE),
      cumulative_value = sum(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(variable, scenario)
}

# Optional diagnostic: compare cumulative precipitation in redprec vs regclim
summarise_precipitation_ratio <- function(water_flux_long) {
  
  water_flux_long %>%
    filter(variable == "precipitation") %>%
    mutate(
      climate = case_when(
        str_detect(as.character(scenario), "redprec") ~ "redprec",
        str_detect(as.character(scenario), "regclim") ~ "regclim",
        TRUE ~ "other"
      ),
      scenario_base = as.character(scenario) %>%
        str_replace("_redprec$", "") %>%
        str_replace("_regclim$", "")
    ) %>%
    group_by(scenario_base, climate) %>%
    summarise(
      cumulative_precipitation = sum(value, na.rm = TRUE),
      mean_precipitation = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = climate,
      values_from = c(cumulative_precipitation, mean_precipitation)
    ) %>%
    mutate(
      cumulative_redprec_over_regclim = cumulative_precipitation_redprec / cumulative_precipitation_regclim,
      mean_redprec_over_regclim = mean_precipitation_redprec / mean_precipitation_regclim
    )
}

# =========================================================
# 7) MAIN ANALYSIS FUNCTION
# =========================================================

analyze_biogeo_water_balance <- function(experiment_table,
                                         biogeochemical_variables = biogeochemical_vars,
                                         transpiration_variables = transpiration_layers,
                                         water_flux_variables = water_flux_vars,
                                         out_prefix = "(null)",
                                         start_date = as.Date("2004-01-01"),
                                         save_plots = FALSE,
                                         plot_dir = NULL) {
  
  # sumstats is loaded once here
  biogeo_long <- prepare_biogeochemical_data(
    experiment_table = experiment_table,
    variables = biogeochemical_variables,
    out_prefix = out_prefix,
    start_date = start_date
  )
  
  # water_balance is loaded once here and reused below
  water_balance <- load_water_balance_data(
    experiment_table = experiment_table,
    out_prefix = out_prefix,
    start_date = start_date
  )
  
  transp_long <- prepare_transpiration_data(
    water_balance = water_balance,
    variables = transpiration_variables
  )
  
  total_transp <- prepare_total_transpiration_data(
    water_balance = water_balance,
    variables = transpiration_variables
  )
  
  water_flux_long <- prepare_water_flux_data(
    water_balance = water_balance,
    variables = water_flux_variables
  )
  
  cumulative_water_flux <- prepare_cumulative_water_flux_data(water_flux_long)
  annual_water_flux <- prepare_annual_water_flux_data(water_flux_long)
  
  plots <- list(
    p_biogeo_grid = plot_biogeochemical_grid(biogeo_long),
    p_biogeo_by_experiment = plot_biogeochemical_faceted_by_experiment(biogeo_long),
    p_transp_layers = plot_transpiration_layers(transp_long),
    p_transp_layers_by_experiment = plot_transpiration_layers_by_experiment(transp_long),
    p_total_transpiration = plot_total_transpiration(total_transp),
    p_total_transpiration_by_experiment = plot_total_transpiration_by_experiment(total_transp),
    p_water_fluxes = plot_water_fluxes(water_flux_long),
    p_water_fluxes_by_experiment = plot_water_fluxes_by_experiment(water_flux_long),
    p_cumulative_water_fluxes = plot_cumulative_water_fluxes(cumulative_water_flux),
    p_annual_water_fluxes = plot_annual_water_fluxes(annual_water_flux)
  )
  
  summaries <- list(
    biogeochemical = summarise_biogeochemical(biogeo_long),
    transpiration_layers = summarise_transpiration_layers(transp_long),
    total_transpiration = summarise_total_transpiration(total_transp),
    water_fluxes = summarise_water_fluxes(water_flux_long),
    annual_water_fluxes = annual_water_flux,
    precipitation_ratio = summarise_precipitation_ratio(water_flux_long)
  )
  
  print(plots$p_biogeo_grid)
  print(plots$p_transp_layers)
  print(plots$p_total_transpiration)
  print(plots$p_water_fluxes)
  print(plots$p_cumulative_water_fluxes)
  
  if (save_plots) {
    
    if (is.null(plot_dir)) {
      plot_dir <- file.path(project_dir, "runs", "_plots_biogeo_water_balance_compare")
    }
    
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    ggsave(
      file.path(plot_dir, "biogeochemical_grid.png"),
      plots$p_biogeo_grid,
      width = 14,
      height = 10,
      dpi = 300
    )
    
    ggsave(
      file.path(plot_dir, "biogeochemical_by_experiment.png"),
      plots$p_biogeo_by_experiment,
      width = 18,
      height = 14,
      dpi = 300
    )
    
    ggsave(
      file.path(plot_dir, "transpiration_layers.png"),
      plots$p_transp_layers,
      width = 14,
      height = 9,
      dpi = 300
    )
    
    ggsave(
      file.path(plot_dir, "transpiration_layers_by_experiment.png"),
      plots$p_transp_layers_by_experiment,
      width = 18,
      height = 10,
      dpi = 300
    )
    
    ggsave(
      file.path(plot_dir, "total_transpiration.png"),
      plots$p_total_transpiration,
      width = 12,
      height = 7,
      dpi = 300
    )
    
    ggsave(
      file.path(plot_dir, "total_transpiration_by_experiment.png"),
      plots$p_total_transpiration_by_experiment,
      width = 14,
      height = 10,
      dpi = 300
    )
    
    ggsave(
      file.path(plot_dir, "water_fluxes.png"),
      plots$p_water_fluxes,
      width = 14,
      height = 10,
      dpi = 300
    )
    
    ggsave(
      file.path(plot_dir, "water_fluxes_by_experiment.png"),
      plots$p_water_fluxes_by_experiment,
      width = 18,
      height = 14,
      dpi = 300
    )
    
    ggsave(
      file.path(plot_dir, "cumulative_water_fluxes.png"),
      plots$p_cumulative_water_fluxes,
      width = 14,
      height = 10,
      dpi = 300
    )
    
    ggsave(
      file.path(plot_dir, "annual_water_fluxes.png"),
      plots$p_annual_water_fluxes,
      width = 14,
      height = 10,
      dpi = 300
    )
    
    message("Plots saved in: ", plot_dir)
  }
  
  invisible(list(
    biogeo_long = biogeo_long,
    water_balance = water_balance,
    transp_long = transp_long,
    total_transp = total_transp,
    water_flux_long = water_flux_long,
    cumulative_water_flux = cumulative_water_flux,
    annual_water_flux = annual_water_flux,
    plots = plots,
    summaries = summaries
  ))
}

# =========================================================
# 8) RUN ANALYSIS
# =========================================================

res_biogeo <- analyze_biogeo_water_balance(
  experiment_table = experiment_table,
  biogeochemical_variables = biogeochemical_vars,
  transpiration_variables = transpiration_layers,
  water_flux_variables = water_flux_vars,
  out_prefix = out_prefix,
  start_date = start_date,
  save_plots = FALSE
)

# =========================================================
# 9) ACCESS RESULTS
# =========================================================

# Main plots
# res_biogeo$plots$p_biogeo_grid
# res_biogeo$plots$p_biogeo_by_experiment
# res_biogeo$plots$p_transp_layers
# res_biogeo$plots$p_transp_layers_by_experiment
# res_biogeo$plots$p_total_transpiration
# res_biogeo$plots$p_total_transpiration_by_experiment
# res_biogeo$plots$p_water_fluxes
# res_biogeo$plots$p_water_fluxes_by_experiment
# res_biogeo$plots$p_cumulative_water_fluxes
# res_biogeo$plots$p_annual_water_fluxes

# Example: show water-flux plots
res_biogeo$plots$p_water_fluxes
res_biogeo$plots$p_cumulative_water_fluxes
res_biogeo$plots$p_annual_water_fluxes

# Summary tables
# res_biogeo$summaries$biogeochemical
# res_biogeo$summaries$transpiration_layers
# res_biogeo$summaries$total_transpiration
# res_biogeo$summaries$water_fluxes
# res_biogeo$summaries$annual_water_fluxes
# res_biogeo$summaries$precipitation_ratio

# Useful check: is redprec really about 50% of regclim?
res_biogeo$summaries$precipitation_ratio
