suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# =========================================================
# 1) PATH HELPER FUNCTIONS
# =========================================================

# Build the base branch path for a given experiment setup
make_branch_path <- function(project_dir, family, vegetation, water_table, texture) {
  file.path(
    project_dir,
    "runs",
    family,
    vegetation,
    water_table,
    texture
  )
}

# Build a configuration list for a single run
make_run_config <- function(
    project_dir,
    family,
    vegetation,
    water_table,
    texture,
    experiment_name,
    out_prefix = "(null)"
) {
  branch_path <- make_branch_path(
    project_dir = project_dir,
    family = family,
    vegetation = vegetation,
    water_table = water_table,
    texture = texture
  )
  
  run_dir <- file.path(branch_path, experiment_name)
  output_dir <- file.path(run_dir, "output")
  common_inputs_dir <- file.path(branch_path, "common_inputs")
  
  list(
    branch_path = branch_path,
    run_dir = run_dir,
    output_dir = output_dir,
    common_inputs_dir = common_inputs_dir,
    pedology_path = file.path(common_inputs_dir, "Paracou_input_pedology.txt"),
    global_path = file.path(common_inputs_dir, "Paracou_input_global.txt"),
    waterbal_path = file.path(output_dir, paste0(out_prefix, "_0_water_balance.txt")),
    figure_dir = file.path(run_dir, "_figures"),
    model_name = experiment_name
  )
}

# Build a dataframe for multiple runs under the same branch
make_files_df <- function(
    project_dir,
    family,
    vegetation,
    water_table,
    texture,
    experiment_names,
    out_prefix = "(null)"
) {
  branch_path <- make_branch_path(
    project_dir = project_dir,
    family = family,
    vegetation = vegetation,
    water_table = water_table,
    texture = texture
  )
  
  common_inputs_dir <- file.path(branch_path, "common_inputs")
  pedology_path <- file.path(common_inputs_dir, "Paracou_input_pedology.txt")
  
  tibble(
    model_name = experiment_names,
    run_dir = file.path(branch_path, experiment_names),
    experiment_dir = file.path(run_dir, "output"),
    pedology_path = pedology_path,
    waterbal_path = file.path(experiment_dir, paste0(out_prefix, "_0_water_balance.txt"))
  )
}

# Build a dataframe for runs that may belong to different branches
# Required columns in experiment_table:
# family, vegetation, water_table, texture, experiment_name
make_files_df_from_table <- function(project_dir, experiment_table, out_prefix = "(null)") {
  required_cols <- c("family", "vegetation", "water_table", "texture", "experiment_name")
  
  missing_cols <- setdiff(required_cols, names(experiment_table))
  if (length(missing_cols) > 0) {
    stop("experiment_table is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  
  experiment_table %>%
    mutate(
      branch_path = purrr::pmap_chr(
        list(family, vegetation, water_table, texture),
        ~ file.path(project_dir, "runs", ..1, ..2, ..3, ..4)
      ),
      run_dir = file.path(branch_path, experiment_name),
      experiment_dir = file.path(run_dir, "output"),
      pedology_path = file.path(branch_path, "common_inputs", "Paracou_input_pedology.txt"),
      waterbal_path = file.path(experiment_dir, paste0(out_prefix, "_0_water_balance.txt")),
      model_name = experiment_name
    )
}

# =========================================================
# 2) FILE READING FUNCTIONS
# =========================================================

# Read the water balance file and clean malformed trailing whitespace
read_waterbal_clean <- function(path) {
  tmp <- tempfile(fileext = ".tsv")
  
  lines <- readLines(path, warn = FALSE)
  lines <- sub("\r$", "", lines)
  lines <- sub("[\t ]+$", "", lines)
  
  writeLines(lines, tmp, useBytes = TRUE)
  
  readr::read_tsv(tmp, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::select(-dplyr::matches("^\\.\\.\\.[0-9]+$"))
}

# Rebuild depth intervals from the pedology file
build_depths_from_pedology <- function(pedology_path) {
  ped <- readr::read_tsv(pedology_path, show_col_types = FALSE, progress = FALSE)
  
  if (!("layer_thickness" %in% names(ped))) {
    stop("Column 'layer_thickness' not found in: ", pedology_path)
  }
  
  layer_thickness <- as.numeric(ped$layer_thickness)
  
  if (any(!is.finite(layer_thickness)) || any(layer_thickness <= 0, na.rm = TRUE)) {
    stop("Non-finite or non-positive layer_thickness in: ", pedology_path)
  }
  
  depth_max <- cumsum(layer_thickness)
  depth_min <- c(0, head(depth_max, -1))
  
  tibble(
    layer     = 0:(length(layer_thickness) - 1),
    depth_min = depth_min,
    depth_max = depth_max
  )
}

# =========================================================
# 3) DATA PROCESSING FUNCTIONS
# =========================================================

# Load all selected runs into a single long-format table
load_all_runs <- function(files_df, start_date) {
  purrr::map_dfr(seq_len(nrow(files_df)), function(i) {
    
    depths <- build_depths_from_pedology(files_df$pedology_path[i])
    wb <- read_waterbal_clean(files_df$waterbal_path[i])
    
    if (!"iter" %in% names(wb)) {
      stop("Column 'iter' not found in water balance file: ", files_df$waterbal_path[i])
    }
    
    wb %>%
      mutate(
        model_name = files_df$model_name[i],
        date = start_date + iter
      ) %>%
      select(model_name, date, matches("^(SWC|SWP)_[0-9]+$")) %>%
      pivot_longer(
        cols = -c(model_name, date),
        names_to = "var_layer",
        values_to = "value"
      ) %>%
      tidyr::extract(
        var_layer,
        into = c("variable", "layer"),
        regex = "^(SWC|SWP)_([0-9]+)$",
        convert = TRUE
      ) %>%
      left_join(depths, by = "layer")
  })
}

# Expand each layer to a regular depth grid for raster plotting
rasterize_profiles <- function(data_long_all, resolution, start_date, experiment_names) {
  data_long_all %>%
    mutate(
      depth_min = round(depth_min, 1),
      depth_max = round(depth_max, 1)
    ) %>%
    rowwise() %>%
    mutate(depth = list(seq(depth_min, depth_max, by = resolution))) %>%
    unnest(depth) %>%
    ungroup() %>%
    mutate(
      days_run = as.numeric(date - start_date),
      sim_year = days_run / 365,
      model_name = factor(model_name, levels = experiment_names)
    )
}

# Check whether required files exist
check_required_files <- function(files_df) {
  missing_pedology <- files_df$pedology_path[!file.exists(files_df$pedology_path)]
  missing_waterbal <- files_df$waterbal_path[!file.exists(files_df$waterbal_path)]
  
  if (length(missing_pedology) > 0) {
    stop("Missing pedology file(s):\n", paste(unique(missing_pedology), collapse = "\n"))
  }
  
  if (length(missing_waterbal) > 0) {
    stop("Missing water balance file(s):\n", paste(unique(missing_waterbal), collapse = "\n"))
  }
}

# =========================================================
# 4) PLOTTING FUNCTIONS
# =========================================================

# Plot SWC profiles for all selected runs
plot_swc_profiles <- function(raster_all, ncol_facets = 2) {
  raster_all %>%
    filter(variable == "SWC") %>%
    ggplot(aes(sim_year, depth, fill = value)) +
    geom_tile(width = 10 / 365) +
    scale_y_reverse() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(expression("SWC [" ~ m^3 ~ m^-3 ~ "]"), direction = -1) +
    facet_wrap(~ model_name, ncol = ncol_facets) +
    labs(x = "Time [years]", y = "Depth [m]", title = "SWC") +
    theme_bw() +
    theme(legend.position = "bottom")
}

# Plot SWP profiles for all selected runs with a shared log color scale
plot_swp_profiles <- function(raster_all, ncol_facets = 2) {
  swp_vals <- raster_all %>%
    filter(variable == "SWP") %>%
    transmute(val = -value) %>%
    filter(is.finite(val), val > 0)
  
  swp_min <- min(swp_vals$val, na.rm = TRUE)
  swp_max <- max(swp_vals$val, na.rm = TRUE)
  
  raster_all %>%
    filter(variable == "SWP") %>%
    mutate(
      val = -value,
      val = ifelse(val <= 0, NA, val)
    ) %>%
    ggplot(aes(sim_year, depth, fill = val)) +
    geom_tile(width = 10 / 365) +
    scale_y_reverse() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(
      expression("|SWP| [" ~ MPa ~ "]"),
      direction = 1,
      trans = "log",
      limits = c(swp_min, swp_max),
      na.value = "grey50"
    ) +
    facet_wrap(~ model_name, ncol = ncol_facets) +
    labs(x = "Time [years]", y = "Depth [m]", title = "SWP") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(angle = 45, hjust = 1)
    )
}

# Save plots to disk
save_profile_plots <- function(p_swc, p_swp, figure_dir, figure_suffix) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename = file.path(figure_dir, paste0("SWC_", figure_suffix, ".png")),
    plot = p_swc,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(figure_dir, paste0("SWP_", figure_suffix, ".png")),
    plot = p_swp,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  cat("Saved figures in:\n", figure_dir, "\n")
}

# =========================================================
# 5) MAIN FUNCTION: SAME BRANCH
# =========================================================

# Analyze one or multiple runs that belong to the same branch
# Example: all runs inside
# runs/theta_w_tests/veg/deepWT/sandy/
analyze_troll_runs <- function(
    project_dir,
    family,
    vegetation = c("veg", "noveg"),
    water_table = c("deepWT", "deepWT"),
    texture = c("sandy", "clayey"),
    experiment_names,
    mode = 0,                 # 0 = single run, 1 = multiple runs in same branch
    save_figures = FALSE,
    resolution = 0.1,
    start_date = as.Date("2004-01-01"),
    out_prefix = "(null)",
    ncol_facets = 2
) {
  vegetation <- match.arg(vegetation)
  water_table <- match.arg(water_table)
  texture <- match.arg(texture)
  
  if (!mode %in% c(0, 1)) {
    stop("'mode' must be 0 or 1.")
  }
  
  if (mode == 0 && length(experiment_names) != 1) {
    stop("For mode = 0, provide exactly one experiment name.")
  }
  
  if (mode == 1 && length(experiment_names) < 2) {
    stop("For mode = 1, provide at least two experiment names.")
  }
  
  if (mode == 0) {
    cfg <- make_run_config(
      project_dir = project_dir,
      family = family,
      vegetation = vegetation,
      water_table = water_table,
      texture = texture,
      experiment_name = experiment_names,
      out_prefix = out_prefix
    )
    
    files_df <- tibble(
      model_name = cfg$model_name,
      run_dir = cfg$run_dir,
      experiment_dir = cfg$output_dir,
      pedology_path = cfg$pedology_path,
      waterbal_path = cfg$waterbal_path
    )
    
    figure_dir <- cfg$figure_dir
    figure_suffix <- experiment_names
    
  } else {
    files_df <- make_files_df(
      project_dir = project_dir,
      family = family,
      vegetation = vegetation,
      water_table = water_table,
      texture = texture,
      experiment_names = experiment_names,
      out_prefix = out_prefix
    )
    
    figure_dir <- file.path(
      make_branch_path(project_dir, family, vegetation, water_table, texture),
      "_figures_compare"
    )
    
    figure_suffix <- paste0(family, "_", vegetation, "_", water_table, "_", texture)
  }
  
  check_required_files(files_df)
  
  data_long_all <- load_all_runs(files_df, start_date = start_date)
  
  raster_all <- rasterize_profiles(
    data_long_all = data_long_all,
    resolution = resolution,
    start_date = start_date,
    experiment_names = files_df$model_name
  )
  
  p_swc_all <- plot_swc_profiles(raster_all, ncol_facets = ncol_facets)
  p_swp_all <- plot_swp_profiles(raster_all, ncol_facets = ncol_facets)
  
  print(p_swc_all)
  print(p_swp_all)
  
  if (save_figures) {
    save_profile_plots(
      p_swc = p_swc_all,
      p_swp = p_swp_all,
      figure_dir = figure_dir,
      figure_suffix = figure_suffix
    )
  }
  
  invisible(list(
    files_df = files_df,
    data_long_all = data_long_all,
    raster_all = raster_all,
    p_swc_all = p_swc_all,
    p_swp_all = p_swp_all,
    figure_dir = figure_dir
  ))
}

# =========================================================
# 6) MAIN FUNCTION: MIXED BRANCHES
# =========================================================

# Analyze runs that may belong to different branches
# Example: compare deepWT vs shallowWT, or sandy vs clayey
#
# experiment_table must contain:
# family, vegetation, water_table, texture, experiment_name
analyze_troll_runs_table <- function(
    project_dir,
    experiment_table,
    save_figures = FALSE,
    resolution = 0.1,
    start_date = as.Date("2004-01-01"),
    out_prefix = "(null)",
    ncol_facets = 2,
    figure_dir = NULL,
    figure_suffix = "compare"
) {
  files_df <- make_files_df_from_table(
    project_dir = project_dir,
    experiment_table = experiment_table,
    out_prefix = out_prefix
  )
  
  check_required_files(files_df)
  
  data_long_all <- load_all_runs(files_df, start_date = start_date)
  
  raster_all <- rasterize_profiles(
    data_long_all = data_long_all,
    resolution = resolution,
    start_date = start_date,
    experiment_names = files_df$model_name
  )
  
  p_swc_all <- plot_swc_profiles(raster_all, ncol_facets = ncol_facets)
  p_swp_all <- plot_swp_profiles(raster_all, ncol_facets = ncol_facets)
  
  print(p_swc_all)
  print(p_swp_all)
  
  if (save_figures) {
    if (is.null(figure_dir)) {
      figure_dir <- file.path(project_dir, "runs", "_figures_compare")
    }
    
    save_profile_plots(
      p_swc = p_swc_all,
      p_swp = p_swp_all,
      figure_dir = figure_dir,
      figure_suffix = figure_suffix
    )
  }
  
  invisible(list(
    files_df = files_df,
    data_long_all = data_long_all,
    raster_all = raster_all,
    p_swc_all = p_swc_all,
    p_swp_all = p_swp_all
  ))
}

project_dir <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting"

#single run 
# res <- analyze_troll_runs(
#   project_dir = project_dir,
#   family = "theta_w_tests",
#   vegetation = "veg",
#   water_table = "deepWT",
#   texture = "sandy",
#   experiment_names = "theta1e-3_ksfloor1e-10_veg_redprec_deepwt_sandy",
#   mode = 0,
#   save_figures = FALSE
# )
# 
# #multiple run, same branch (EXMPLE)
# res2 <- analyze_troll_runs(
#   project_dir = project_dir,
#   family = "theta_w_tests",
#   vegetation = "veg",
#   water_table = "deepWT",
#   texture = "sandy",
#   experiment_names = c(
#     "theta1e-3_ksfloor1e-10_veg_redprec_deepwt_sandy",
#     "theta1e-4_ksfloor1e-10_veg_redprec_deepwt_sandy",
#     "theta1e-6_ksfloor1e-10_veg_redprec_deepwt_sandy"
#   ),
#   mode = 1,
#   save_figures = FALSE
# )

#mixed branches
# experiment_table <- tibble(
#   family = c("infiltration", "infiltration", "infiltration", "infiltration", "infiltration", "infiltration", "infiltration", "infiltration"),
#   vegetation = c("veg", "veg", "veg", "veg", "veg", "veg", "veg", "veg" ),
#   water_table = c("shallowWT", "shallowWT", "deepWT", "deepWT", "shallowWT", "shallowWT", "deepWT", "deepWT"),
#   texture = c("sandy", "clayey", "sandy",  "clayey", "sandy", "clayey", "sandy",  "clayey"),
#   experiment_name = c(
#     "inf_shallowWT_sandy_regclim",
#     "inf_shallowWT_clayey_regclim",
#     "inf_deepWT_sandy_regclim",
#     "inf_deepWT_clayey_regclim",
#     "inf_shallowWT_sandy_redprec",
#     "inf_shallowWT_clayey_redprec",
#     "inf_deepWT_sandy_redprec",
#     "inf_deepWT_clayey_redprec"
#   )
# )

# experiment_table <- tibble(
#   family = c("infiltration", "infiltration", "infiltration", "infiltration"),
#   vegetation = c("veg", "veg", "veg", "veg"),
#   water_table = c("deepWT", "deepWT", "deepWT", "deepWT"),
#   texture = c("sandy", "sandy", "sandy", "sandy"),
#   experiment_name = c(
#     "redprec_lowertheta",
#     "redprec_highertheta",
#     "regclim_highertheta_lowerks",
#     "redprec_highertheta_lowerks"
#   )
# )
experiment_table <- tibble(
  family = c("infiltration"),
  vegetation = c("veg"),
  water_table = c("deepWT"),
  texture = c("sandy"),
  experiment_name = c(
    "regclim_highertheta_lowerks"
  )
)


res3 <- analyze_troll_runs_table(
  project_dir = project_dir,
  experiment_table = experiment_table,
  save_figures = FALSE
)

# ============================================================
# # Compare SWC and SWP in the first soil layer: regclim vs redprec
# # ============================================================
# 
# library(tidyverse)
# library(lubridate)
# 
# # ============================================================
# # Create wide water_balance object from the runs used in res3
# # ============================================================
# 
# load_water_balance_wide <- function(files_df, start_date = as.Date("2004-01-01")) {
#   purrr::map_dfr(seq_len(nrow(files_df)), function(i) {
#     
#     wb <- read_waterbal_clean(files_df$waterbal_path[i])
#     
#     if (!"iter" %in% names(wb)) {
#       stop("Column 'iter' not found in water balance file: ", files_df$waterbal_path[i])
#     }
#     
#     wb %>%
#       mutate(
#         scenario = files_df$model_name[i],
#         date = start_date + iter,
#         sim_day = as.numeric(iter),
#         sim_year = sim_day / 365
#       )
#   })
# }
# 
# water_balance <- load_water_balance_wide(
#   files_df = res3$files_df,
#   start_date = as.Date("2004-01-01")
# )
# 
# names(water_balance)
# 
# # ------------------------------------------------------------
# # 1. Check if water_balance exists
# # ------------------------------------------------------------
# 
# if (!exists("water_balance")) {
#   stop("Object 'water_balance' not found. Load water_balance before running this block.")
# }
# 
# # ------------------------------------------------------------
# # 2. Define first-layer variables
# # ------------------------------------------------------------
# 
# first_layer_vars <- c("SWC_0", "SWP_0")
# 
# missing_first_layer_vars <- setdiff(first_layer_vars, names(water_balance))
# 
# if (length(missing_first_layer_vars) > 0) {
#   stop(
#     "These first-layer variables were not found in water_balance: ",
#     paste(missing_first_layer_vars, collapse = ", "),
#     "\nAvailable variables are:\n",
#     paste(names(water_balance), collapse = ", ")
#   )
# }
# 
# # ------------------------------------------------------------
# # 3. Prepare first-layer data
# # ------------------------------------------------------------
# 
# first_layer_water <- water_balance %>%
#   mutate(
#     scenario_chr = as.character(scenario),
#     climate = case_when(
#       str_detect(scenario_chr, "regclim") ~ "regclim",
#       str_detect(scenario_chr, "redprec") ~ "redprec",
#       TRUE ~ NA_character_
#     ),
#     scenario_base = scenario_chr %>%
#       str_remove("_regclim$") %>%
#       str_remove("_redprec$")
#   ) %>%
#   filter(!is.na(climate)) %>%
#   select(
#     scenario,
#     scenario_base,
#     climate,
#     date,
#     sim_day,
#     sim_year,
#     SWC_0,
#     SWP_0
#   ) %>%
#   pivot_longer(
#     cols = c(SWC_0, SWP_0),
#     names_to = "variable",
#     values_to = "value"
#   ) %>%
#   mutate(
#     variable = factor(variable, levels = c("SWC_0", "SWP_0")),
#     climate = factor(climate, levels = c("regclim", "redprec"))
#   )
# 
# # ------------------------------------------------------------
# # 4. Monthly means — lighter and clearer than daily plots
# # ------------------------------------------------------------
# 
# first_layer_monthly <- first_layer_water %>%
#   mutate(
#     year = year(date),
#     month = month(date),
#     month_date = as.Date(sprintf("%04d-%02d-01", year, month))
#   ) %>%
#   group_by(scenario_base, climate, month_date, variable) %>%
#   summarise(
#     value = mean(value, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # ------------------------------------------------------------
# # 5. Plot monthly SWC_0 and SWP_0: regclim vs redprec
# # ------------------------------------------------------------
# 
# p_first_layer_monthly <- ggplot(
#   first_layer_monthly,
#   aes(x = month_date, y = value, color = climate)
# ) +
#   geom_line(linewidth = 0.7) +
#   facet_grid(variable ~ scenario_base, scales = "free_y") +
#   theme_bw(base_size = 13) +
#   labs(
#     x = "Date",
#     y = "Value",
#     title = "First soil layer: SWC and SWP under regclim vs redprec",
#     color = "Climate"
#   )
# 
# print(p_first_layer_monthly)
# 
# # ------------------------------------------------------------
# # 6. Difference redprec - regclim
# # ------------------------------------------------------------
# 
# first_layer_diff <- first_layer_monthly %>%
#   pivot_wider(
#     names_from = climate,
#     values_from = value
#   ) %>%
#   mutate(
#     diff_redprec_minus_regclim = redprec - regclim,
#     relative_diff = redprec / regclim
#   )
# 
# # ------------------------------------------------------------
# # 7. Plot absolute difference
# # ------------------------------------------------------------
# 
# p_first_layer_diff <- ggplot(
#   first_layer_diff,
#   aes(x = month_date, y = diff_redprec_minus_regclim)
# ) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_line(linewidth = 0.7) +
#   facet_grid(variable ~ scenario_base, scales = "free_y") +
#   theme_bw(base_size = 13) +
#   labs(
#     x = "Date",
#     y = "redprec - regclim",
#     title = "Difference in first-layer SWC and SWP: redprec minus regclim"
#   )
# 
# print(p_first_layer_diff)
# 
# # ------------------------------------------------------------
# # 8. Summary table
# # ------------------------------------------------------------
# 
# first_layer_summary <- first_layer_diff %>%
#   group_by(scenario_base, variable) %>%
#   summarise(
#     mean_regclim = mean(regclim, na.rm = TRUE),
#     mean_redprec = mean(redprec, na.rm = TRUE),
#     mean_difference = mean(diff_redprec_minus_regclim, na.rm = TRUE),
#     min_difference = min(diff_redprec_minus_regclim, na.rm = TRUE),
#     max_difference = max(diff_redprec_minus_regclim, na.rm = TRUE),
#     mean_relative_diff = mean(relative_diff, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# print(first_layer_summary)
# # write.csv(first_layer_summary, "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/infiltration/veg/first_layer_summary.csv")
# 
# 
# # ============================================================
# # Diagnostic: does SWP_0 behave consistently with SWC_0?
# # ============================================================
# 
# first_layer_daily <- water_balance %>%
#   mutate(
#     scenario_chr = as.character(scenario),
#     climate = case_when(
#       str_detect(scenario_chr, "regclim") ~ "regclim",
#       str_detect(scenario_chr, "redprec") ~ "redprec",
#       TRUE ~ NA_character_
#     ),
#     scenario_base = scenario_chr %>%
#       str_remove("_regclim$") %>%
#       str_remove("_redprec$")
#   ) %>%
#   filter(!is.na(climate)) %>%
#   select(
#     scenario,
#     scenario_base,
#     climate,
#     date,
#     sim_day,
#     sim_year,
#     SWC_0,
#     SWP_0
#   )
# 
# # 1. Print full summary without truncation
# print(first_layer_summary, n = Inf, width = Inf)
# 
# 
# # 2. Scatterplot SWC_0 vs SWP_0
# p_swc_swp_relation <- ggplot(
#   first_layer_daily,
#   aes(x = SWC_0, y = SWP_0, color = climate)
# ) +
#   geom_point(alpha = 0.25, size = 0.7) +
#   facet_wrap(~ scenario_base, scales = "free") +
#   theme_bw(base_size = 13) +
#   labs(
#     x = "SWC_0",
#     y = "SWP_0",
#     title = "Relationship between first-layer SWC and SWP",
#     color = "Climate"
#   )
# 
# print(p_swc_swp_relation)
# 
# # 3. Correlation by scenario
# swc_swp_correlation <- first_layer_daily %>%
#   group_by(scenario_base, climate) %>%
#   summarise(
#     cor_swc_swp = cor(SWC_0, SWP_0, use = "complete.obs"),
#     mean_SWC_0 = mean(SWC_0, na.rm = TRUE),
#     mean_SWP_0 = mean(SWP_0, na.rm = TRUE),
#     min_SWP_0 = min(SWP_0, na.rm = TRUE),
#     max_SWP_0 = max(SWP_0, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# print(swc_swp_correlation, n = Inf, width = Inf)
# # write.csv(swc_swp_correlation,"/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/infiltration/veg/swc_swp_correlation.csv")
# 
# # ============================================================
# # Distribution of SWC_0 and SWP_0 by climate
# # ============================================================
# 
# first_layer_distribution <- first_layer_daily %>%
#   pivot_longer(
#     cols = c(SWC_0, SWP_0),
#     names_to = "variable",
#     values_to = "value"
#   ) %>%
#   group_by(scenario_base, climate, variable) %>%
#   summarise(
#     mean = mean(value, na.rm = TRUE),
#     median = median(value, na.rm = TRUE),
#     q01 = quantile(value, 0.01, na.rm = TRUE),
#     q05 = quantile(value, 0.05, na.rm = TRUE),
#     q10 = quantile(value, 0.10, na.rm = TRUE),
#     q25 = quantile(value, 0.25, na.rm = TRUE),
#     q75 = quantile(value, 0.75, na.rm = TRUE),
#     q90 = quantile(value, 0.90, na.rm = TRUE),
#     q95 = quantile(value, 0.95, na.rm = TRUE),
#     q99 = quantile(value, 0.99, na.rm = TRUE),
#     min = min(value, na.rm = TRUE),
#     max = max(value, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# print(first_layer_distribution, n = Inf, width = Inf)
# # write.csv(first_layer_distribution,"/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/infiltration/veg/firts_layer_distribution.csv")
# 
# # ============================================================
# # Frequency of dry SWP_0 conditions
# # ============================================================
# 
# swp_dry_frequency <- first_layer_daily %>%
#   group_by(scenario_base, climate) %>%
#   summarise(
#     n_days = n(),
#     
#     prop_SWP_below_500  = mean(SWP_0 < -500, na.rm = TRUE),
#     prop_SWP_below_1000 = mean(SWP_0 < -1000, na.rm = TRUE),
#     prop_SWP_below_2000 = mean(SWP_0 < -2000, na.rm = TRUE),
#     prop_SWP_below_4000 = mean(SWP_0 < -4000, na.rm = TRUE),
#     prop_SWP_below_6000 = mean(SWP_0 < -6000, na.rm = TRUE),
#     
#     .groups = "drop"
#   )
# 
# print(swp_dry_frequency, n = Inf, width = Inf)
# # write.csv(swp_dry_frequency,"/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/infiltration/veg/swp_dry_frequency.csv")
# 
# # ============================================================
# # Density plots: SWC_0 and SWP_0
# # ============================================================
# 
# first_layer_density <- first_layer_daily %>%
#   pivot_longer(
#     cols = c(SWC_0, SWP_0),
#     names_to = "variable",
#     values_to = "value"
#   )
# 
# p_first_layer_density <- ggplot(
#   first_layer_density,
#   aes(x = value, fill = climate)
# ) +
#   geom_density(alpha = 0.35) +
#   facet_grid(variable ~ scenario_base, scales = "free") +
#   theme_bw(base_size = 13) +
#   labs(
#     x = "Value",
#     y = "Density",
#     title = "Distribution of first-layer SWC and SWP",
#     fill = "Climate"
#   )
# 
# print(p_first_layer_density)
# 
# library(tidyverse)
# library(lubridate)
# 
# wb_all = water_balance
# wb_inf <- wb_all %>%
#   mutate(
#     date = as.Date("2004-01-01") + iter - 1,
#     year = year(date),
#     
#     # same unit as throughfall/runoff, probably m water per timestep
#     infiltration_real = pmax(throughfall - runoff, 0),
#     
#     # fractions only when there is water reaching the soil
#     infil_frac = if_else(throughfall > 0, infiltration_real / throughfall, NA_real_),
#     runoff_frac = if_else(throughfall > 0, runoff / throughfall, NA_real_)
#   )
# 
# inf_summary <- wb_inf %>%
#   group_by(scenario_base, climate) %>%
#   summarise(
#     total_throughfall = sum(throughfall, na.rm = TRUE),
#     total_infiltration = sum(infiltration_real, na.rm = TRUE),
#     total_runoff = sum(runoff, na.rm = TRUE),
#     
#     infil_frac_total = total_infiltration / total_throughfall,
#     runoff_frac_total = total_runoff / total_throughfall,
#     
#     mean_daily_throughfall = mean(throughfall, na.rm = TRUE),
#     mean_daily_infiltration = mean(infiltration_real, na.rm = TRUE),
#     mean_daily_runoff = mean(runoff, na.rm = TRUE),
#     
#     wet_days = sum(throughfall > 0, na.rm = TRUE),
#     runoff_days = sum(runoff > 0 & throughfall > 0, na.rm = TRUE),
#     frac_wet_days_with_runoff = runoff_days / wet_days,
#     
#     .groups = "drop"
#   )
# 
# print(inf_summary, n = Inf, width = Inf)
# 
# # ============================================================
# # 2 & 3. Prepare data for ALL layers
# # ============================================================
# 
# # We use regex to select all columns like SWC_0, SWC_1, SWP_0, SWP_1, etc.
# all_layers_water <- water_balance %>%
#   mutate(
#     scenario_chr = as.character(scenario),
#     climate = case_when(
#       str_detect(scenario_chr, "regclim") ~ "regclim",
#       str_detect(scenario_chr, "redprec") ~ "redprec",
#       TRUE ~ NA_character_
#     ),
#     scenario_base = scenario_chr %>%
#       str_remove("_regclim$") %>%
#       str_remove("_redprec$")
#   ) %>%
#   filter(!is.na(climate)) %>%
#   # Select metadata and ALL columns matching SWC_x or SWP_x
#   select(
#     scenario, scenario_base, climate, date, sim_day, sim_year,
#     matches("^(SWC|SWP)_[0-9]+$")
#   ) %>%
#   # This is the MAGIC STEP: 
#   # It takes "SWC_3", puts "SWC" in 'variable' and "3" in 'layer'
#   pivot_longer(
#     cols = matches("^(SWC|SWP)_[0-9]+$"),
#     names_to = c("variable", "layer"),
#     names_pattern = "^(SWC|SWP)_([0-9]+)$",
#     values_to = "value"
#   ) %>%
#   mutate(
#     variable = factor(variable, levels = c("SWC", "SWP")),
#     climate = factor(climate, levels = c("regclim", "redprec")),
#     layer = as.numeric(layer) # Make layer numeric so it sorts correctly 0, 1, 2...
#   )
# 
# # ============================================================
# # 4. Monthly means for ALL layers
# # ============================================================
# 
# all_layers_monthly <- all_layers_water %>%
#   mutate(
#     year = year(date),
#     month = month(date),
#     month_date = as.Date(sprintf("%04d-%02d-01", year, month))
#   ) %>%
#   # Notice we added 'layer' to the grouping
#   group_by(scenario_base, climate, layer, month_date, variable) %>%
#   summarise(
#     value = mean(value, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # ============================================================
# # 5. Plotting (Tip: Filter layers to avoid overcrowded plots)
# # ============================================================
# 
# # Plotting all 6 layers * 4 scenarios * 2 variables is too much for one figure.
# # Let's filter to plot specific layers of interest (e.g., surface, middle, bottom)
# layers_to_plot <- c(0, 2, 5) 
# 
# p_layers_monthly <- all_layers_monthly %>%
#   filter(layer %in% layers_to_plot) %>%
#   ggplot(aes(x = month_date, y = value, color = climate)) +
#   geom_line(linewidth = 0.7) +
#   # Facet by Variable (Row) and Scenario + Layer (Column)
#   facet_grid(variable ~ scenario_base + paste("Layer", layer), scales = "free_y") +
#   theme_bw(base_size = 11) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(
#     x = "Date",
#     y = "Value",
#     title = "SWC and SWP across selected depths: regclim vs redprec",
#     color = "Climate"
#   )
# 
# print(p_layers_monthly)
# 
# # ============================================================
# # 6 & 7. Difference redprec - regclim across layers
# # ============================================================
# 
# all_layers_diff <- all_layers_monthly %>%
#   pivot_wider(
#     names_from = climate,
#     values_from = value
#   ) %>%
#   mutate(
#     diff_redprec_minus_regclim = redprec - regclim,
#     relative_diff = redprec / regclim
#   )
# 
# p_layers_diff <- all_layers_diff %>%
#   filter(layer %in% layers_to_plot) %>%
#   ggplot(aes(x = month_date, y = diff_redprec_minus_regclim, color = as.factor(layer))) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_line(linewidth = 0.7, alpha = 0.8) +
#   facet_grid(variable ~ scenario_base, scales = "free_y") +
#   theme_bw(base_size = 12) +
#   labs(
#     x = "Date",
#     y = "redprec - regclim",
#     title = "Difference in SWC/SWP (redprec - regclim) by Layer",
#     color = "Layer"
#   )
# 
# print(p_layers_diff)
# 
# # ============================================================
# # 8. Summary Table for ALL layers
# # ============================================================
# 
# # Now you get the statistics for every single layer independently
# all_layers_summary <- all_layers_diff %>%
#   group_by(scenario_base, layer, variable) %>%
#   summarise(
#     mean_regclim = mean(regclim, na.rm = TRUE),
#     mean_redprec = mean(redprec, na.rm = TRUE),
#     mean_difference = mean(diff_redprec_minus_regclim, na.rm = TRUE),
#     min_difference = min(diff_redprec_minus_regclim, na.rm = TRUE),
#     max_difference = max(diff_redprec_minus_regclim, na.rm = TRUE),
#     mean_relative_diff = mean(relative_diff, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(scenario_base, variable, layer)
# 
# print(all_layers_summary, n = Inf, width = Inf)
# # write.csv(all_layers_summary,"/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/infiltration/veg/all_layers_summary.csv")
# 
# 
# library(tidyverse)
# 
# wb_inf2 <- wb_inf %>%
#   mutate(
#     scenario = as.character(scenario),
#     
#     climate = case_when(
#       str_detect(scenario, regex("redprec", ignore_case = TRUE)) ~ "redprec",
#       str_detect(scenario, regex("regclim", ignore_case = TRUE)) ~ "regclim",
#       TRUE ~ NA_character_
#     ),
#     
#     scenario_base = scenario %>%
#       str_replace_all(regex("redprec|regclim", ignore_case = TRUE), "") %>%
#       str_replace_all("__+", "_") %>%
#       str_remove("^_") %>%
#       str_remove("_$")
#   )
# 
# inf_summary <- wb_inf2 %>%
#   group_by(scenario_base, climate) %>%
#   summarise(
#     total_throughfall = sum(throughfall, na.rm = TRUE),
#     total_infiltration = sum(infiltration_real, na.rm = TRUE),
#     total_runoff = sum(runoff, na.rm = TRUE),
#     
#     infil_frac_total = total_infiltration / total_throughfall,
#     runoff_frac_total = total_runoff / total_throughfall,
#     
#     mean_daily_throughfall = mean(throughfall, na.rm = TRUE),
#     mean_daily_infiltration = mean(infiltration_real, na.rm = TRUE),
#     mean_daily_runoff = mean(runoff, na.rm = TRUE),
#     
#     wet_days = sum(throughfall > 0, na.rm = TRUE),
#     runoff_days = sum(runoff > 0 & throughfall > 0, na.rm = TRUE),
#     frac_wet_days_with_runoff = runoff_days / wet_days,
#     
#     .groups = "drop"
#   )
# 
# print(inf_summary, n = Inf, width = Inf)
# 
# wb_budget_summary <- wb_inf2 %>%
#   mutate(
#     transpiration_total =
#       transpiration_0 + transpiration_1 + transpiration_2 +
#       transpiration_3 + transpiration_4 + transpiration_5 +
#       transpiration1016
#   ) %>%
#   group_by(scenario_base, climate) %>%
#   summarise(
#     total_throughfall = sum(throughfall, na.rm = TRUE),
#     total_infiltration = sum(infiltration_real, na.rm = TRUE),
#     total_runoff = sum(runoff, na.rm = TRUE),
#     total_evaporation = sum(evaporation, na.rm = TRUE),
#     total_transpiration = sum(transpiration_total, na.rm = TRUE),
#     total_leak = sum(leak, na.rm = TRUE),
#     
#     mean_SWC_0 = mean(SWC_0, na.rm = TRUE),
#     mean_SWP_0 = mean(SWP_0, na.rm = TRUE),
#     
#     .groups = "drop"
#   )
# 
# print(wb_budget_summary, n = Inf, width = Inf)
