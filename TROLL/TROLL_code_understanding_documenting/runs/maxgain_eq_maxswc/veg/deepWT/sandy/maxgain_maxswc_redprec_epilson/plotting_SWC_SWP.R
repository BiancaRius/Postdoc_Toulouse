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
experiment_table <- tibble(
  family = c("units_vertmov", "units_vertmov"),
  vegetation = c("veg", "veg"),
  water_table = c("deepWT", "deepWT"),
  texture = c("sandy", "sandy"),
  experiment_name = c(
    "debug_hydraulic_locking",
    "test_nan"
  )
)

res3 <- analyze_troll_runs_table(
  project_dir = project_dir,
  experiment_table = experiment_table,
  save_figures = FALSE
)
