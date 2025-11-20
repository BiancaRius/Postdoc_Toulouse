# -------------------------------------------
# Script to plot results from TROLL simulations with fixed WTD values and
# Capillarity implementation
# Author: Bianca Fazio Rius
# -------------------------------------------

library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)
library(tidyr) # Added for data reshaping (pivot_longer)

# ========== 1) Main path and scenario definitions ==========
# Define one main path, then add subfolders for each scenario.
main_path <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/Capillarity_implementation/"


# Example: you can add 2, 3, 5, ... scenarios.
# If you do NOT name them, the script will use folder names as labels.
scenario_paths <- c("notUnified_noWT_fcSWC", "unified_noWT_fcSWC", "notUnified_noWT_capOn_fcSWC", "Unified_noWT_maxSWC",
                    "Unified_maxSWC_shallowWT",  "Unified_maxSWC_deepWT")
#scenario_paths <- c("wtOn_capOn_vegetation_deepWT", "wtOn_capOn_vegetation_shallowWT", "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/WT_implementation/regular_climate/deep_WTD/")
#scenario_paths <- c("wtOn_capOn_vegetation_shallowWT", "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/WT_implementation/regular_climate/shallow_WTD/")

# ========== 2) Variable groups and corresponding input file types ==========
biogeochemical_vars   <- c("npp", "gpp", "agb", "sum1", "sum10", "sum30", "ba", "ba10", "litterfall")
soil_water_content    <- c("SWC_0", "SWC_1", "SWC_2", "SWC_3", "SWC_4")
soil_water_potential  <- c("SWP_0", "SWP_1", "SWP_2", "SWP_3", "SWP_4")
transpiration_layers  <- c("transpiration_0", "transpiration_1", "transpiration_2", "transpiration_3", "transpiration_4")
water_flux_vars       <- c("precipitation", "interception", "throughfall", "runoff", "leak", "evaporation")
water_change_volume   <- c("wcv_0", "wcv_1", "wcv_2", "wcv_3", "wcv_4")
water_upward_volume   <- c("wupv_interface_0_1", "wupv_interface_1_2", "wupv_interface_2_3", "wupv_interface_3_4")

variable_groups <- list(
  biogeochemical       = biogeochemical_vars,
  soil_water_content   = soil_water_content,
  soil_water_potential = soil_water_potential,
  transpiration_layers = transpiration_layers,
  water_fluxes         = water_flux_vars,
  water_change_volume  = water_change_volume,
  water_upward_volume  = water_upward_volume
)

variable_file_map <- c(
  biogeochemical       = "sumstats",
  soil_water_content   = "water_balance",
  soil_water_potential = "water_balance",
  transpiration_layers = "water_balance",
  water_fluxes         = "water_balance",
  water_change_volume  = "vertical_water_flux",
  water_upward_volume  = "vertical_water_flux"
)

# ========== 3) Helper functions ==========
# Check if path is absolute or relative
is_absolute <- function(p) grepl("^(/|~)", p)

# Resolve full path for a given child directory
resolve_path <- function(child) {
  if (is_absolute(child)) {
    return(child)
  } else {
    return(file.path(main_path, child))
  }
}

# Generate scenario labels: use provided names if they exist, otherwise use folder names
scenario_labels <- function(paths_vec) {
  if (!is.null(names(paths_vec)) && any(nzchar(names(paths_vec)))) {
    return(names(paths_vec))
  } else {
    return(basename(normalizePath(sapply(paths_vec, resolve_path), mustWork = FALSE)))
  }
}

# Safe file reader (checks existence and catches errors)
safe_read <- function(full_dir, file_type, label) {
  file_path <- file.path(full_dir, paste0("(null)_0_", file_type, ".txt"))
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }
  df <- tryCatch(
    {
      x <- read.table(file_path, header = TRUE)
      x$scenario <- label
      x
    },
    error = function(e) {
      warning("Failed to read: ", file_path, " -> ", conditionMessage(e))
      NULL
    }
  )
  df
}

# ========== 4) Plot function for multiple scenarios ==========
plot_variable_multi <- function(variable, scenario_paths_vec, file_type, seed = NULL) {
  labs <- scenario_labels(scenario_paths_vec)
  dirs <- sapply(scenario_paths_vec, resolve_path, USE.NAMES = FALSE)
  
  dfs <- Map(function(d, lab) safe_read(d, file_type, lab), dirs, labs)
  all_data <- bind_rows(dfs[!sapply(dfs, is.null)])
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No valid data found for file type: ", file_type)
  }
  if (!variable %in% names(all_data)) {
    stop("Column '", variable, "' does not exist in file type '", file_type, "'.")
  }
  
  # Random color palette (reproducible if seed is provided)
  if (!is.null(seed)) set.seed(seed)
  nsc <- length(unique(all_data$scenario))
  base_cols <- grDevices::hcl.colors(n = max(nsc, 3), palette = "Dark 3")
  scen_levels <- unique(all_data$scenario)
  scen_colors <- setNames(sample(base_cols, nsc), scen_levels)
  
  ggplot(all_data, aes(x = iter, y = .data[[variable]], color = scenario)) +
    geom_line(linewidth = 0.7, alpha = 0.9) +
    scale_color_manual(values = scen_colors, name = NULL) +
    scale_x_continuous(
      name = "Year",
      breaks = seq(0, 10000, by = 365 * 3),
      labels = function(x) floor(x / 365) + 1
    ) +
    labs(y = variable, title = variable) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}


# ========== 5) Main wrapper function ==========
# This function remains unchanged.
plot_results_multi <- function(plot_all = TRUE,
                               selected_variable = NULL,
                               save_pdf = FALSE,
                               pdf_file = "TROLL_multi_scenarios.pdf",
                               seed_colors = NULL) {
  
  if (save_pdf) {
    pdf(pdf_file, width = 12, height = 8)
    on.exit(dev.off(), add = TRUE)
  }
  
  for (group in names(variable_groups)) {
    vars <- variable_groups[[group]]
    file_type <- variable_file_map[[group]]
    message("\n==> Group: ", group, " | file_type: ", file_type)
    
    if (plot_all) {
      for (var in vars) {
        message("Plotting: ", var)
        print(plot_variable_multi(var, scenario_paths, file_type, seed = seed_colors))
      }
    } else if (!is.null(selected_variable) && selected_variable %in% vars) {
      message("Plotting selected variable: ", selected_variable)
      print(plot_variable_multi(selected_variable, scenario_paths, file_type, seed = seed_colors))
    }
  }
  
  if (save_pdf) {
    message("PDF saved at: ", normalizePath(pdf_file))
  }
}


# ========== 6) Grid plot function with SHARED Y-AXIS ==========
# UPDATED FUNCTION
plot_SWC_grid <- function(scenario_paths_vec = scenario_paths,
                          seed_colors = 123,
                          ncol = 2,
                          file_type = "water_balance") {
  
  # --- Step 1: Load data for all scenarios ONCE ---
  labs <- scenario_labels(scenario_paths_vec)
  dirs <- sapply(scenario_paths_vec, resolve_path, USE.NAMES = FALSE)
  dfs <- Map(function(d, lab) safe_read(d, file_type, lab), dirs, labs)
  all_data <- bind_rows(dfs[!sapply(dfs, is.null)])
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No valid data found for file type: ", file_type)
  }
  
  # --- Step 2: Determine the common y-axis range ---
  swc_vars <- paste0("SWC_", 0:4)
  
  # Reshape data to long format to easily find the overall min/max
  long_data <- all_data %>%
    select(iter, scenario, all_of(swc_vars)) %>%
    tidyr::pivot_longer(
      cols = all_of(swc_vars),
      names_to = "variable",
      values_to = "swc_value"
    )
  
  # Calculate the shared y-axis limits across all SWC variables
  y_limits <- range(long_data$swc_value, na.rm = TRUE)
  message("Common y-axis range for SWC plots: ", round(y_limits[1], 4), " to ", round(y_limits[2], 4))
  
  # --- Step 3: Generate colors and individual plots ---
  # Generate reproducible colors
  if (!is.null(seed_colors)) set.seed(seed_colors)
  nsc <- length(unique(long_data$scenario))
  scen_levels <- unique(long_data$scenario)
  scen_colors <- setNames(sample(grDevices::hcl.colors(n = max(nsc, 3), palette = "Dark 3"), nsc), scen_levels)
  
  # Build one ggplot per variable, applying the shared y-axis limit
  plots <- lapply(swc_vars, function(v) {
    plot_data <- filter(long_data, variable == v)
    
    ggplot(plot_data, aes(x = iter, y = swc_value, color = scenario)) +
      geom_line(linewidth = 0.7, alpha = 0.9) +
      scale_color_manual(values = scen_colors, name = NULL) +
      # Apply the shared y-axis limits to each plot
      scale_y_continuous(limits = y_limits) +
      scale_x_continuous(
        name = "Year",
        breaks = seq(0, 10000, by = 365 * 3),
        labels = function(x) floor(x / 365) + 1
      ) +
      labs(title = v, y = "Soil Water Content (m3/m3)") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
  })
  
  # --- Step 4: Assemble the final grid ---
  fig <- patchwork::wrap_plots(plots, ncol = ncol, guides = "collect") &
    theme(legend.position = "bottom")
  
  return(fig)
}


# ========== 7) Examples ==========
# Your original examples will still work for single plots
# plot_results_multi(plot_all = FALSE, selected_variable = "SWC_0", seed_colors = 123)

# EXAMPLE USE for the updated grid plot function:
# This will now generate a grid where all subplots have the same y-axis range.
p <- plot_SWC_grid(scenario_paths, seed_colors = 123, ncol = 2)
print(p)

# If you want to save the grid plot:
# ggsave("SWC_grid_shared_yaxis.pdf", p, width = 12, height = 8)


# ===== Upward fluxes: grid de todas as interfaces (Y compartilhado) =====
plot_upward_flux_grid <- function(scenario_paths_vec = scenario_paths,
                                  seed_colors = 123,
                                  ncol = 2,
                                  file_type = "vertical_water_flux",
                                  ylabel = "Upward water volume (m³)") {
  
  labs <- scenario_labels(scenario_paths_vec)
  dirs <- sapply(scenario_paths_vec, resolve_path, USE.NAMES = FALSE)
  dfs  <- Map(function(d, lab) safe_read(d, file_type, lab), dirs, labs)
  all_data <- dplyr::bind_rows(dfs[!sapply(dfs, is.null)])
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No valid data found for file type: ", file_type)
  }
  
  # interfaces esperadas
  wupv_vars_expected <- c("wupv_interface_0_1","wupv_interface_1_2",
                          "wupv_interface_2_3","wupv_interface_3_4")
  
  wupv_vars <- intersect(wupv_vars_expected, names(all_data))
  if (length(wupv_vars) == 0) {
    stop("Nenhuma coluna wupv_interface_*_* encontrada em '", file_type, "'.")
  } else if (length(wupv_vars) < length(wupv_vars_expected)) {
    warning("Colunas faltando: ",
            paste(setdiff(wupv_vars_expected, wupv_vars), collapse=", "))
  }
  
  long_data <- all_data %>%
    dplyr::select(iter, scenario, dplyr::all_of(wupv_vars)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(wupv_vars),
                        names_to = "variable", values_to = "value")
  
  # eixo Y compartilhado
  y_limits <- range(long_data$value, na.rm = TRUE)
  
  # cores por cenário
  if (!is.null(seed_colors)) set.seed(seed_colors)
  scen_levels <- unique(long_data$scenario)
  scen_colors <- setNames(
    sample(grDevices::hcl.colors(n = max(length(scen_levels), 3), palette = "Dark 3"),
           length(scen_levels)),
    scen_levels
  )
  
  plots <- lapply(wupv_vars, function(v) {
    plot_data <- dplyr::filter(long_data, variable == v)
    ggplot2::ggplot(plot_data, ggplot2::aes(x = iter, y = value, color = scenario)) +
      ggplot2::geom_line(linewidth = 0.7, alpha = 0.9) +
      ggplot2::scale_color_manual(values = scen_colors, name = NULL) +
      ggplot2::scale_y_continuous(limits = y_limits) +
      ggplot2::scale_x_continuous(
        name = "Year",
        breaks = seq(0, 10000, by = 365 * 3),
        labels = function(x) floor(x / 365) + 1
      ) +
      ggplot2::labs(title = v, y = ylabel) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
  })
  
  fig <- patchwork::wrap_plots(plots, ncol = ncol, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  
  fig
}

# Grid com todas as interfaces (0_1, 1_2, 2_3, 3_4)
p_grid <- plot_upward_flux_grid(scenario_paths, seed_colors = 123, ncol = 2)
print(p_grid)
# ggsave("upward_flux_grid.pdf", p_grid, width = 12, height = 8)


# =========================
# 3.b) STABILITY DIAGNOSTICS
# =========================

# --- Utility: Load all data for a given file type across scenarios ---
load_all_data <- function(scenario_paths_vec, file_type) {
  labs <- scenario_labels(scenario_paths_vec)
  dirs <- sapply(scenario_paths_vec, resolve_path, USE.NAMES = FALSE)
  dfs  <- Map(function(d, lab) safe_read(d, file_type, lab), dirs, labs)
  dplyr::bind_rows(dfs[!sapply(dfs, is.null)])
}

# -------- D1: Temporal variability (stability over time) --------
# Computes per-layer and per-scenario metrics:
# - sd_diff: standard deviation of time increments (x[t]-x[t-1])
# - tv: total variation (sum of absolute changes)
# - iqr: interquartile range (robust to outliers)
temporal_variability <- function(df, var_prefix = c("SWC","SWP")) {
  var_prefix <- match.arg(var_prefix)
  vars <- grep(paste0("^", var_prefix, "_[0-9]+$"), names(df), value = TRUE)
  if (length(vars) == 0) stop("No columns found matching ", var_prefix, "_*.")
  
  out <- lapply(vars, function(v) {
    df |>
      dplyr::group_by(scenario) |>
      dplyr::arrange(iter, .by_group = TRUE) |>
      dplyr::summarise(
        variable = v,
        sd_diff = stats::sd(diff(.data[[v]]), na.rm = TRUE),
        tv      = sum(abs(diff(.data[[v]])), na.rm = TRUE),
        iqr     = IQR(.data[[v]], na.rm = TRUE),
        .groups = "drop"
      )
  })
  dplyr::bind_rows(out) |>
    dplyr::mutate(layer = as.integer(gsub(paste0(var_prefix, "_"), "", variable)))
}

# -------- D2: Vertical coherence (SWP monotonicity) --------
# Checks how often the soil water potential profile violates monotonicity.
# Expected (non-frozen soil): SWP becomes less negative with depth.
vertical_coherence_swp <- function(df) {
  swp_vars <- grep("^SWP_[0-9]+$", names(df), value = TRUE)
  if (length(swp_vars) == 0) stop("No SWP_* columns found.")
  swp_vars <- swp_vars[order(as.integer(gsub("SWP_", "", swp_vars)))]
  
  df |>
    dplyr::select(iter, scenario, dplyr::all_of(swp_vars)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(swp_vars),
                        names_to = "layer", values_to = "swp") |>
    dplyr::mutate(layer = as.integer(gsub("SWP_", "", layer))) |>
    dplyr::arrange(scenario, iter, layer) |>
    dplyr::group_by(scenario, iter) |>
    dplyr::summarise(
      breaks = sum(diff(swp) < 0, na.rm = TRUE), # <0 means inversion
      .groups = "drop"
    ) |>
    dplyr::group_by(scenario) |>
    dplyr::summarise(
      frac_timesteps_with_inversion = mean(breaks > 0, na.rm = TRUE),
      avg_inversions_per_timestep = mean(breaks, na.rm = TRUE),
      .groups = "drop"
    )
}

# -------- D3: Water balance check --------
# Compares changes in storage (ΔStorage) with input/output fluxes.
# Uses available columns among: precipitation, interception, throughfall,
# runoff, leak, evaporation, transpiration_*, wcv_*.
water_balance_check <- function(df_wb, df_vf = NULL) {
  has_wcv <- !is.null(df_vf) && any(grepl("^wcv_[0-9]+$", names(df_vf)))
  
  dat <- dplyr::left_join(
    df_wb, if (has_wcv) dplyr::select(df_vf, dplyr::any_of(c("iter","scenario",
                                                             grep("^wcv_[0-9]+$", names(df_vf), value=TRUE)))), 
    by = c("iter","scenario")
  )
  
  # Compute ΔStorage
  wcv_cols <- grep("^wcv_[0-9]+$", names(dat), value = TRUE)
  dat$delta_storage <- if (length(wcv_cols) > 0) {
    rowSums(dat[, wcv_cols], na.rm = TRUE)
  } else {
    # fallback: estimate from SWC time differences
    swc_cols <- grep("^SWC_[0-9]+$", names(dat), value = TRUE)
    if (length(swc_cols) == 0) stop("No wcv_* or SWC_* to estimate storage change.")
    dat <- dat |>
      dplyr::group_by(scenario) |>
      dplyr::arrange(iter, .by_group = TRUE) |>
      dplyr::mutate(!!!rlang::set_names(
        lapply(swc_cols, \(cname) c(NA, diff(.data[[cname]]))),
        paste0("d_", swc_cols))) |>
      dplyr::ungroup()
    rowSums(dat[, grep("^d_SWC_[0-9]+$", names(dat), value = TRUE)], na.rm = TRUE)
  }
  
  # Main fluxes
  get <- function(nm) if (nm %in% names(dat)) dat[[nm]] else 0
  P   <- get("precipitation")
  I   <- get("interception")
  TF  <- if ("throughfall" %in% names(dat)) dat$throughfall else (P - I)
  Ev  <- get("evaporation")
  Tr  <- rowSums(dat[, grep("^transpiration_[0-9]+$", names(dat), value = TRUE), drop = FALSE], na.rm = TRUE)
  Qr  <- get("runoff")
  Lk  <- get("leak")
  
  # Simple column water balance:
  #  TF - (Ev + Tr + Qr + Lk) ≈ ΔStorage
  resid <- TF - (Ev + Tr + Qr + Lk) - dat$delta_storage
  
  tibble::tibble(scenario = dat$scenario, resid = resid) |>
    dplyr::group_by(scenario) |>
    dplyr::summarise(
      mean_error   = mean(resid, na.rm = TRUE),
      mean_abs_error = mean(abs(resid), na.rm = TRUE),
      p95_abs_error = quantile(abs(resid), 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

# -------- Master function: runs all diagnostics --------
run_stability_diagnostics <- function(scenario_paths_vec = scenario_paths) {
  wb  <- load_all_data(scenario_paths_vec, "water_balance")
  vf  <- try(load_all_data(scenario_paths_vec, "vertical_water_flux"), silent = TRUE)
  if (inherits(vf, "try-error")) vf <- NULL
  
  d1_swc <- temporal_variability(wb, var_prefix = "SWC")
  d1_swp <- if (any(grepl("^SWP_[0-9]+$", names(wb)))) temporal_variability(wb, "SWP") else NULL
  d2_swp <- if (any(grepl("^SWP_[0-9]+$", names(wb)))) vertical_coherence_swp(wb) else NULL
  d3_bal <- try(water_balance_check(wb, vf), silent = TRUE)
  if (inherits(d3_bal, "try-error")) d3_bal <- NULL
  
  list(
    D1_temporal_variability_SWC = d1_swc,
    D1_temporal_variability_SWP = d1_swp,
    D2_vertical_coherence_SWP   = d2_swp,
    D3_water_balance            = d3_bal
  )
}

# -------- Quick visualization of D1 (bar plots per layer/scenario) --------
plot_temporal_variability <- function(d1_table, metric = c("sd_diff","tv","iqr")) {
  metric <- match.arg(metric)
  ggplot(d1_table, aes(x = factor(layer), y = .data[[metric]], fill = scenario)) +
    geom_col(position = position_dodge(width = 0.8)) +
    labs(x = "Layer", y = metric, title = paste("Temporal variability (", metric, ")", sep="")) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
}


# === Run diagnostics ===
diag <- run_stability_diagnostics(scenario_paths)

# 1) Temporal stability (SWC): smaller = more stable
print(plot_temporal_variability(diag$D1_temporal_variability_SWC, metric = "sd_diff"))
print(plot_temporal_variability(diag$D1_temporal_variability_SWC, metric = "tv"))

# 2) Vertical coherence (SWP): fewer inversions = better physical consistency
if (!is.null(diag$D2_vertical_coherence_SWP)) {
  print(diag$D2_vertical_coherence_SWP)
}

# 3) Water balance: mean_error ~ 0 and small p95_abs_error = good conservation
if (!is.null(diag$D3_water_balance)) {
  print(diag$D3_water_balance)
}


# ========== SWC grid: one panel per layer x scenario ==========
plot_SWC_grid_facets <- function(scenario_paths_vec = scenario_paths,
                                 file_type = "water_balance",
                                 ncol = NULL) {
  
  # --- 1) Load data for all scenarios ---
  labs <- scenario_labels(scenario_paths_vec)
  dirs <- sapply(scenario_paths_vec, resolve_path, USE.NAMES = FALSE)
  dfs  <- Map(function(d, lab) safe_read(d, file_type, lab), dirs, labs)
  all_data <- dplyr::bind_rows(dfs[!sapply(dfs, is.null)])
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No valid data found for file type: ", file_type)
  }
  
  # --- 2) Select SWC columns and reshape to long format ---
  swc_vars <- paste0("SWC_", 0:4)
  
  long_data <- all_data %>%
    dplyr::select(iter, scenario, dplyr::all_of(swc_vars)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(swc_vars),
      names_to = "variable",
      values_to = "swc_value"
    ) %>%
    dplyr::mutate(
      # Extract numeric layer index (0,1,2,3,4) for nicer facet labels
      layer = as.integer(gsub("SWC_", "", variable))
    )
  
  # --- 3) Common y-axis range across all layers and scenarios ---
  y_limits <- range(long_data$swc_value, na.rm = TRUE)
  message("Common y-axis range for SWC facets: ",
          round(y_limits[1], 4), " to ", round(y_limits[2], 4))
  
  # --- 4) Single ggplot with facet_grid(layer ~ scenario) ---
  p <- ggplot2::ggplot(long_data, ggplot2::aes(x = iter, y = swc_value)) +
    ggplot2::geom_line(linewidth = 0.6) +  # one line per scenario in each facet
    ggplot2::scale_y_continuous(limits = y_limits) +
    ggplot2::scale_x_continuous(
      name = "Year",
      breaks = seq(0, 10000, by = 365 * 3),
      labels = function(x) floor(x / 365) + 1
    ) +
    ggplot2::facet_grid(
      rows   = ggplot2::vars(layer),    # one row per soil layer
      cols   = ggplot2::vars(scenario)  # one column per scenario
    ) +
    ggplot2::labs(
      title = "Soil Water Content by layer and scenario",
      y     = "Soil Water Content (m³/m³)",
      x     = "Year"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}

p_swc_facets <- plot_SWC_grid_facets(scenario_paths)
print(p_swc_facets)

# ========== SWP grid: one panel per layer x scenario ==========


plot_SWP_grid_facets <- function(scenario_paths_vec = scenario_paths,
                                 file_type = "water_balance",
                                 ncol = NULL) {
  
  # --- 1) Load data for all scenarios ---
  labs <- scenario_labels(scenario_paths_vec)
  dirs <- sapply(scenario_paths_vec, resolve_path, USE.NAMES = FALSE)
  dfs  <- Map(function(d, lab) safe_read(d, file_type, lab), dirs, labs)
  all_data <- dplyr::bind_rows(dfs[!sapply(dfs, is.null)])
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No valid data found for file type: ", file_type)
  }
  
  # --- 2) Select SWP columns and reshape to long format ---
  swp_vars <- paste0("SWP_", 0:4)
  
  long_data <- all_data %>%
    dplyr::select(iter, scenario, dplyr::all_of(swp_vars)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(swp_vars),
      names_to = "variable",
      values_to = "swp_value"
    ) %>%
    dplyr::mutate(
      # Extract numeric layer index (0,1,2,3,4) for nicer facet labels
      layer = as.integer(gsub("SWP_", "", variable))
    )
  
  # --- 3) Common y-axis range across all layers and scenarios ---
  y_limits <- range(long_data$swp_value, na.rm = TRUE)
  message("Common y-axis range for SWP facets: ",
          round(y_limits[1], 4), " to ", round(y_limits[2], 4))
  
  # --- 4) Single ggplot with facet_grid(layer ~ scenario) ---
  p <- ggplot2::ggplot(long_data, ggplot2::aes(x = iter, y = swp_value)) +
    ggplot2::geom_line(linewidth = 0.6) +  # one line per scenario in each facet
    ggplot2::scale_y_continuous(limits = y_limits) +
    ggplot2::scale_x_continuous(
      name = "Year",
      breaks = seq(0, 10000, by = 365 * 3),
      labels = function(x) floor(x / 365) + 1
    ) +
    ggplot2::facet_grid(
      rows   = ggplot2::vars(layer),    # one row per soil layer
      cols   = ggplot2::vars(scenario)  # one column per scenario
    ) +
    ggplot2::labs(
      title = "Soil Water Potential by layer and scenario",
      y     = "Soil Water Potential",
      x     = "Year"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}

p_swp_facets <- plot_SWP_grid_facets(scenario_paths)
print(p_swp_facets)

plot_SWP_zoom_layers <- function(scenario_paths_vec = scenario_paths,
                                 layers_to_plot = 1:4,
                                 file_type = "water_balance",
                                 ylim_padding = 0.1) {
  
  # --- 1) Load data ---
  labs <- scenario_labels(scenario_paths_vec)
  dirs <- sapply(scenario_paths_vec, resolve_path, USE.NAMES = FALSE)
  dfs  <- Map(function(d, lab) safe_read(d, file_type, lab), dirs, labs)
  all_data <- dplyr::bind_rows(dfs[!sapply(dfs, is.null)])
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No valid data found for file type: ", file_type)
  }
  
  # --- 2) Select SWP columns and reshape ---
  swp_vars <- paste0("SWP_", layers_to_plot)
  
  long_data <- all_data %>%
    dplyr::select(iter, scenario, dplyr::all_of(swp_vars)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(swp_vars),
                        names_to = "variable",
                        values_to = "swp_value") %>%
    dplyr::mutate(
      layer = as.integer(gsub("SWP_", "", variable))
    )
  
  # --- 3) Determine zoomed y-axis range ---
  y_min <- min(long_data$swp_value, na.rm = TRUE)
  y_max <- max(long_data$swp_value, na.rm = TRUE)
  
  # add a small padding to make the plot nicer
  y_range <- c(y_min - abs(y_min)*ylim_padding,
               y_max + abs(y_max)*ylim_padding)
  
  message("Zoom Y-limits: ", round(y_range[1], 4), " to ", round(y_range[2], 4))
  
  # --- 4) Facet plot ---
  p <- ggplot2::ggplot(long_data,
                       ggplot2::aes(x = iter, y = swp_value, color = scenario)) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.9) +
    ggplot2::scale_y_continuous(limits = y_range) +
    ggplot2::scale_x_continuous(
      name = "Year",
      breaks = seq(0, 10000, by = 365 * 5),
      labels = function(x) floor(x / 365) + 1
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(layer),
      cols = ggplot2::vars(scenario),
      scales = "free_x"
    ) +
    ggplot2::labs(
      title = "SWP (zoomed view) for layers 1–4",
      y = "Soil Water Potential",
      x = "Year"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}
p_zoom <- plot_SWP_zoom_layers(scenario_paths)
print(p_zoom)
