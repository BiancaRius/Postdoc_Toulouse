  suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
    library(arrow)
  })
  
  # -----------------------------
  # Params
  # -----------------------------
  resolution <- 0.1
  dry_months <- 10:11
  wet_months <- 6:7
  start_date <- as.Date("2004-01-01")
  
  # -----------------------------
  # Paths (ONE situation)
  # -----------------------------
  model_path <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting"
  out_prefix <- "(null)"
  
  files_df <- tibble(
    model_name    = "run1",  # qualquer nome (só pra facetar)
    experiment_dir = model_path,
    pedology_path  = file.path(model_path, "Paracou_input_pedology.txt"),
    waterbal_path  = file.path(model_path, paste0(out_prefix, "_0_water_balance.txt"))
  )
  
  stopifnot(file.exists(files_df$pedology_path))
  stopifnot(file.exists(files_df$waterbal_path))
  
  # -----------------------------
  # Robust reader
  # -----------------------------
  read_waterbal_clean <- function(path) {
    tmp <- tempfile(fileext = ".tsv")
    lines <- readLines(path, warn = FALSE)
    lines <- sub("\r$", "", lines)
    lines <- sub("[\t ]+$", "", lines)
    writeLines(lines, tmp, useBytes = TRUE)
    
    readr::read_tsv(tmp, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::select(-dplyr::matches("^\\.\\.\\.[0-9]+$"))
  }
  
  # -----------------------------
  # Build depths from pedology
  # -----------------------------
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
  
  depths <- build_depths_from_pedology(files_df$pedology_path)
  
  # -----------------------------
  # Read + pivot water balance (SWC/SWP by layer)
  # -----------------------------
  wb <- read_waterbal_clean(files_df$waterbal_path)
  
  if (!"iter" %in% names(wb)) stop("Column 'iter' not found in water balance file.")
  
  data_long <- wb %>%
    mutate(
      model_name = files_df$model_name[1],
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
  
  # check
  stopifnot(all(c("depth_min","depth_max") %in% names(data_long)))
  
  # -----------------------------
  # Rasterize
  # -----------------------------
  raster <- data_long %>%
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
      sim_year = days_run / 365
    )
  
  # -----------------------------
  # PLOTS
  # -----------------------------
  
  # SWC raster
  p_swc_raster <- raster %>%
    filter(variable == "SWC") %>%
    ggplot(aes(sim_year, depth, fill = value)) +
    geom_tile(width = 10/365) +
    scale_y_reverse() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(expression("SWC ["~m^3~m^3~"]"), direction = -1) +
    labs(x = "Time [years]", y = "Depth [m]", title = "SWC") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  print(p_swc_raster)
  
  # SWP raster (log scale, NA as grey)
  swp_vals <- raster %>%
    filter(variable == "SWP") %>%
    transmute(val = -value) %>%
    filter(is.finite(val), val > 0)
  
  swp_min <- min(swp_vals$val, na.rm = TRUE)
  swp_max <- max(swp_vals$val, na.rm = TRUE)
  
  p_swp_raster <- raster %>%
    filter(variable == "SWP") %>%
    mutate(val = -value,
           val = ifelse(val <= 0, NA, val)) %>%
    ggplot(aes(sim_year, depth, fill = val)) +
    geom_tile(width = 10/365) +
    scale_y_reverse() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(
      expression("|SWP| ["~MPa~"]"),
      direction = 1,
      trans = "log",
      limits = c(swp_min, swp_max),
      na.value = "grey50"
    ) +
    labs(x = "Time [years]", y = "Depth [m]", title = "SWP") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  print(p_swp_raster)
  
  # SWC profiles wet vs dry
  p_swc_profile <- raster %>%
    filter(variable == "SWC", month(date) %in% c(dry_months, wet_months)) %>%
    mutate(season = ifelse(month(date) %in% dry_months, "dry", "wet")) %>%
    group_by(depth_min, depth_max, season) %>%
    summarise(
      m  = median(value, na.rm = TRUE),
      ll = quantile(value, 0.05, na.rm = TRUE),
      hh = quantile(value, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ggplot(aes(depth_max, m, color = season, group = season)) +
    geom_line(linewidth = 0.6) +
    geom_linerange(aes(ymin = ll, ymax = hh), linewidth = 0.5) +
    geom_point(size = 1) +
    coord_flip() +
    scale_x_reverse("Depth [m]") +
    scale_y_reverse(expression("SWC ["~m^3~m^3~"]")) +
    labs(title = "SWC profiles (wet vs dry)", color = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  print(p_swc_profile)
  
  # SWP profiles wet vs dry (log)
  p_swp_profile <- raster %>%
    filter(variable == "SWP", month(date) %in% c(dry_months, wet_months)) %>%
    mutate(season = ifelse(month(date) %in% dry_months, "dry", "wet"),
           val = -value) %>%
    filter(is.finite(val), val > 0) %>%
    group_by(depth_min, depth_max, season) %>%
    summarise(
      m  = median(val, na.rm = TRUE),
      ll = quantile(val, 0.05, na.rm = TRUE),
      hh = quantile(val, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ggplot(aes(depth_max, m, color = season, group = season)) +
    geom_line(linewidth = 0.6) +
    geom_linerange(aes(ymin = ll, ymax = hh), linewidth = 0.5) +
    geom_point(size = 1) +
    coord_flip() +
    scale_x_reverse("Depth [m]") +
    scale_y_log10(expression("|SWP| ["~MPa~"]")) +
    labs(title = "SWP profiles (wet vs dry)", color = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  print(p_swp_profile)
  
  # (Optional) save figures
  # ggsave(file.path(model_path, "SWC_raster.png"), p_swc_raster, width = 10, height = 6, dpi = 300)
  # ggsave(file.path(model_path, "SWP_raster.png"), p_swp_raster, width = 10, height = 6, dpi = 300)
  # ggsave(file.path(model_path, "SWC_profiles.png"), p_swc_profile, width = 7, height = 5, dpi = 300)
  # ggsave(file.path(model_path, "SWP_profiles.png"), p_swp_profile, width = 7, height = 5, dpi = 300)