suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(arrow)
})

# -----------------------------
# Params
# -----------------------------
resolution <- 0.1
start_date <- as.Date("2004-01-01")
out_prefix <- "(null)"

# -----------------------------
# Paths (MULTIPLE situations)
# -----------------------------
base_dir <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/noveg/"

files_df <- tibble(
  model_name = c(
    "baseline_noveg",
    "legacy_theta0_reset1e-3_noveg",
    "theta_thres1e-6_noveg",
    "theta_thres1e-5_noveg",
    "theta_thres1e-4_noveg",
    "theta_thres1e-3_noveg"
  ),
  experiment_dir = file.path(base_dir, c(
    "baseline_noveg/output",
    "legacy_theta0_reset1e-3_noveg/output",
    "theta_thres1e-6_noveg/output",
    "theta_thres1e-5_noveg/output",
    "theta_thres1e-4_noveg/output",
    "theta_thres1e-3_noveg/output"
  )),
  pedology_path = rep(
    file.path(base_dir, "common_inputs/Paracou_input_pedology.txt"),
    6
  )
) %>%
  mutate(
    waterbal_path = file.path(experiment_dir, paste0(out_prefix, "_0_water_balance.txt"))
  )

stopifnot(all(file.exists(files_df$pedology_path)))
stopifnot(all(file.exists(files_df$waterbal_path)))

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
    days_run = as.numeric(date - start_date),
    sim_year = days_run / 365
  )

# -----------------------------
# Shared SWP scale
# -----------------------------
swp_vals <- raster_all %>%
  filter(variable == "SWP") %>%
  transmute(val = -value) %>%
  filter(is.finite(val), val > 0)

swp_min <- min(swp_vals$val, na.rm = TRUE)
swp_max <- max(swp_vals$val, na.rm = TRUE)

# -----------------------------
# PLOTS: ALL SWC together
# -----------------------------
p_swc_all <- raster_all %>%
  filter(variable == "SWC") %>%
  ggplot(aes(sim_year, depth, fill = value)) +
  geom_tile(width = 10/365) +
  scale_y_reverse() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(expression("SWC ["~m^3~m^-3~"]"), direction = -1) +
  facet_wrap(~ model_name, ncol = 2) +
  labs(x = "Time [years]", y = "Depth [m]", title = "SWC - all configurations") +
  theme_bw() +
  theme(legend.position = "bottom")

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
  facet_wrap(~ model_name, ncol = 2) +
  labs(x = "Time [years]", y = "Depth [m]", title = "SWP - all configurations") +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_swp_all)

# -----------------------------
# Optional save
# -----------------------------
ggsave("SWC_all_configs.png", p_swc_all, width = 12, height = 8, dpi = 300)
ggsave("SWP_all_configs.png", p_swp_all, width = 12, height = 8, dpi = 300)
