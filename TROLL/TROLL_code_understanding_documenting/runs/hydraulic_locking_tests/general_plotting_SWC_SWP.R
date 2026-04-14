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
out_prefix <- "(null)"

# -----------------------------
# Paths
# One row = one experiment
# Add as many rows as you want
# -----------------------------
files_df <- tribble(
  ~model_name,               ~panel_label,                    ~experiment_dir,                                                                                                                                         ~pedology_path,
  "baseline_noveg_redprec",  "baseline_noveg_redprec",        "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/noveg/theta/shallowWT/baseline_noveg_redprec/output", "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/noveg/theta/shallowWT/common_inputs/Paracou_input_pedology.txt",
  "noveg_noWT", "noveg_noWT", "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/no_WT/noveg/noveg_noWT/output", "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/no_WT/noveg/common_inputs/Paracou_input_pedology.txt",
  "theta1e-6_ksfloor1e-12_noveg_redprec", "ks1e-12_noveg_redprec_shallow", "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/noveg/ks/shallowWT/theta1e-6_ksfloor1e-10_noveg_redprec/output", "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/noveg/ks/shallowWT/common_inputs/Paracou_input_pedology.txt",
  "theta1e-6_ksfloor1e-12_veg_redprec_shallowwt", "ks1e-12_veg_redprec_shallow", "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/veg/ks/shallowWT/sandy/theta1e-6_ksfloor1e-10_veg_redprec_shallowwt/output/","~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/veg/ks/shallowWT/sandy/common_inputs/Paracou_input_pedology.txt",
  "theta1e-6_ksfloor1e-12_noveg_redprec_deepwt", "ks1e-12_noveg_redprec_deep", "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/noveg/ks/deepWT/theta1e-6_ksfloor1e-12_noveg_redprec_deepwt/output", "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/noveg/ks/deepWT/common_inputs/Paracou_input_pedology.txt"

  # Exemplo de como adicionar mais runs:
  # ,
  # "run2", "run2_label",
  # "~/caminho/para/run2/output",
  # "~/caminho/para/run2/Paracou_input_pedology.txt"
  # ,
  # "run3", "run3_label",
  # "~/caminho/para/run3/output",
  # "~/caminho/para/run3/Paracou_input_pedology.txt"
) %>%
  mutate(
    experiment_dir = path.expand(experiment_dir),
    pedology_path  = path.expand(pedology_path),
    waterbal_path  = file.path(experiment_dir, paste0(out_prefix, "_0_water_balance.txt"))
  )

stopifnot(all(file.exists(files_df$pedology_path)))
stopifnot(all(file.exists(files_df$waterbal_path)))

# -----------------------------
# Output dir
# -----------------------------
output_dir <- path.expand(
  "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/compare_formulations"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

# -----------------------------
# Read all runs
# -----------------------------
data_long_all <- purrr::map_dfr(seq_len(nrow(files_df)), function(i) {
  
  depths <- build_depths_from_pedology(files_df$pedology_path[i])
  wb <- read_waterbal_clean(files_df$waterbal_path[i])
  
  if (!"iter" %in% names(wb)) {
    stop("Column 'iter' not found in water balance file: ", files_df$waterbal_path[i])
  }
  
  wb %>%
    mutate(
      model_name  = files_df$model_name[i],
      panel_label = files_df$panel_label[i],
      date        = start_date + iter
    ) %>%
    select(model_name, panel_label, date, matches("^(SWC|SWP)_[0-9]+$")) %>%
    pivot_longer(
      cols = -c(model_name, panel_label, date),
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

# -----------------------------
# Rasterize
# -----------------------------
raster_all <- data_long_all %>%
  mutate(
    depth_min = round(depth_min, 1),
    depth_max = round(depth_max, 1)
  ) %>%
  rowwise() %>%
  mutate(depth = list(seq(depth_min, depth_max, by = resolution))) %>%
  unnest(depth) %>%
  ungroup() %>%
  mutate(
    days_run    = as.numeric(date - start_date),
    sim_year    = days_run / 365,
    panel_label = factor(panel_label, levels = unique(files_df$panel_label))
  )

# -----------------------------
# Shared scales
# -----------------------------
swc_vals <- raster_all %>%
  filter(variable == "SWC") %>%
  filter(is.finite(value))

swc_min <- min(swc_vals$value, na.rm = TRUE)
swc_max <- max(swc_vals$value, na.rm = TRUE)

swp_vals <- raster_all %>%
  filter(variable == "SWP") %>%
  transmute(val = -value) %>%
  filter(is.finite(val), val > 0)

swp_min <- min(swp_vals$val, na.rm = TRUE)
swp_max <- max(swp_vals$val, na.rm = TRUE)

# -----------------------------
# Optional: save rasterized data
# -----------------------------
arrow::write_parquet(
  raster_all,
  sink = file.path(output_dir, "raster_all.parquet")
)

# -----------------------------
# Plot sizes adapt to number of panels
# -----------------------------
n_panels <- nrow(files_df)
ncol_facets <- ifelse(n_panels == 1, 1, 2)
nrow_facets <- ceiling(n_panels / ncol_facets)

plot_width  <- ifelse(ncol_facets == 1, 7, 12)
plot_height <- max(4, 3.5 * nrow_facets)

# -----------------------------
# PLOTS: ALL SWC together
# -----------------------------
p_swc_all <- raster_all %>%
  filter(variable == "SWC") %>%
  ggplot(aes(sim_year, depth, fill = value)) +
  geom_tile(width = 10 / 365) +
  scale_y_reverse() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(
    expression("SWC [" ~ m^3 ~ m^-3 ~ "]"),
    direction = -1,
    limits = c(swc_min, swc_max),
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 1.5
    )
  ) +
  facet_wrap(~ panel_label, ncol = ncol_facets) +
  labs(
    x = "Time [years]",
    y = "Depth [m]",
    title = "SWC"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 17),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
  )

print(p_swc_all)

# -----------------------------
# PLOTS: ALL SWP together
# -----------------------------
p_swp_all <- raster_all %>%
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
    na.value = "grey50",
    guide = guide_colorbar(
      barwidth = 20,
      barheight = 1.5
    )
  ) +
  facet_wrap(~ panel_label, ncol = ncol_facets) +
  labs(
    x = "Time [years]",
    y = "Depth [m]",
    title = "SWP"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 17),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

print(p_swp_all)

# limites globais
swc_vals <- raster_all %>%
  filter(variable == "SWC") %>%
  filter(is.finite(value))

swc_min <- min(swc_vals$value, na.rm = TRUE)
swc_max <- max(swc_vals$value, na.rm = TRUE)

swp_vals <- raster_all %>%
  filter(variable == "SWP") %>%
  transmute(val = -value) %>%
  filter(is.finite(val), val > 0)

swp_min <- min(swp_vals$val, na.rm = TRUE)
swp_max <- max(swp_vals$val, na.rm = TRUE)

make_plot_swc <- function(run_name) {
  raster_all %>%
    filter(variable == "SWC", panel_label == run_name) %>%
    ggplot(aes(sim_year, depth, fill = value)) +
    geom_tile(width = 10 / 365) +
    scale_y_reverse() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(
      expression("SWC [" ~ m^3 ~ m^-3 ~ "]"),
      direction = -1,
      limits = c(swc_min, swc_max),
      guide = guide_colorbar(
        barwidth = 20,
        barheight = 1.5
      )
    ) +
    labs(
      x = "Time [years]",
      y = "Depth [m]",
      title = paste("SWC -", run_name)
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(angle = 45, hjust = 1, size = 15),
      panel.grid = element_blank(),
      legend.title = element_text(size = 17),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    )
}

make_plot_swp <- function(run_name) {
  raster_all %>%
    filter(variable == "SWP", panel_label == run_name) %>%
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
      na.value = "grey50",
      guide = guide_colorbar(
        barwidth = 20,
        barheight = 1.5
      )
    ) +
    labs(
      x = "Time [years]",
      y = "Depth [m]",
      title = paste("SWP -", run_name)
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(angle = 45, hjust = 1, size = 15),
      panel.grid = element_blank(),
      legend.title = element_text(size = 17),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    )
}

p1 <- make_plot_swc("ks1e-12_noveg_redprec_shallow")
p2 <- make_plot_swp("ks1e-12_noveg_redprec_shallow")

print(p1)
print(p2)

p3 <- make_plot_swc("ks1e-12_noveg_redprec_shallow")
p4 <- make_plot_swp("ks1e-12_noveg_redprec_shallow")

print(p3)
print(p4)

p5 <- make_plot_swp("ks1e-12_noveg_redprec_deep")
plot(p5)

p6 <- make_plot_swp("ks1e-12_veg_redprec_shallow")
p6
# -----------------------------
# Save figures
# -----------------------------
# ggsave(
#   filename = file.path(output_dir, "SWC_all_configs.png"),
#   plot = p_swc_all,
#   width = plot_width,
#   height = plot_height,
#   dpi = 300
# )
# 
# ggsave(
#   filename = file.path(output_dir, "SWP_all_configs.png"),
#   plot = p_swp_all,
#   width = plot_width,
#   height = plot_height,
#   dpi = 300
# )
# 
# cat("Saved files in:\n", output_dir, "\n")