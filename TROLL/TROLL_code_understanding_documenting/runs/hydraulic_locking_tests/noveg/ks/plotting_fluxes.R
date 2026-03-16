suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------
# Params
# -----------------------------
out_prefix  <- "(null)"
save_plots  <- TRUE
start_date  <- as.Date("2004-01-01")

# -----------------------------
# Base paths
# -----------------------------
base_dir <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/noveg/ks"
plot_dir <- file.path(base_dir, "plots_vertical_flux_faceted")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Runs to compare
# -----------------------------
files_df <- tibble(
  model_name = c(
    "theta1e-6_ksfloor1e-12_noveg",
    "theta1e-6_ksfloor1e-10_noveg",
    "theta1e-6_ksfloor1e-8_noveg"
  )
) %>%
  mutate(
    experiment_dir = file.path(base_dir, model_name, "output"),
    flux_path = file.path(
      experiment_dir,
      paste0(out_prefix, "_0_vertical_water_flux.txt")
    )
  )

print(files_df)

# -----------------------------
# Check files
# -----------------------------
missing_files <- files_df %>%
  filter(!file.exists(flux_path))

if (nrow(missing_files) > 0) {
  stop(
    "These flux files were not found:\n",
    paste(missing_files$flux_path, collapse = "\n")
  )
}

# -----------------------------
# Reader
# -----------------------------
read_flux_file <- function(path) {
  df <- read.table(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    fill = TRUE
  )
  
  empty_idx <- which(names(df) == "")
  if (length(empty_idx) > 0) {
    names(df)[empty_idx] <- paste0("empty_col_", seq_along(empty_idx))
  }
  
  df %>%
    select(-starts_with("empty_col_"))
}

# -----------------------------
# Read all runs
# -----------------------------
flux_all <- purrr::map_dfr(seq_len(nrow(files_df)), function(i) {
  df <- read_flux_file(files_df$flux_path[i])
  
  if (!"iter" %in% names(df)) {
    stop("Column 'iter' not found in: ", files_df$flux_path[i])
  }
  
  df %>%
    mutate(
      iter = as.numeric(iter),
      model_name = files_df$model_name[i],
      date = start_date + iter,
      sim_year = as.numeric(date - start_date) / 365
    )
})

# -----------------------------
# Pretty model labels
# -----------------------------
model_labels <- c(
  "theta1e-6_ksfloor1e-12_noveg" = "Ks floor = 1e-12",
  "theta1e-6_ksfloor1e-10_noveg" = "Ks floor = 1e-10",
  "theta1e-6_ksfloor1e-8_noveg"  = "Ks floor = 1e-8"
)

model_levels <- c(
  "theta1e-6_ksfloor1e-12_noveg",
  "theta1e-6_ksfloor1e-10_noveg",
  "theta1e-6_ksfloor1e-8_noveg"
)

# -----------------------------
# Long format
# -----------------------------
flux_long <- flux_all %>%
  pivot_longer(
    cols = -c(model_name, iter, date, sim_year),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    value = as.numeric(value),
    metric = case_when(
      str_detect(variable, "^mean_flux_layers") ~ "mean_flux",
      str_detect(variable, "^mean_abs_flux_layers") ~ "mean_abs_flux",
      str_detect(variable, "^mean_delta_swp_layers") ~ "mean_delta_swp",
      str_detect(variable, "^mean_ks_harmonic_layers") ~ "mean_ks_harmonic",
      str_detect(variable, "^net_volumetric_change_layers") ~ "net_volumetric_change",
      str_detect(variable, "^gross_volumetric_change_layers") ~ "gross_volumetric_change",
      str_detect(variable, "^mean_swp_layer") ~ "mean_swp_layer",
      str_detect(variable, "^mean_ks_layer") ~ "mean_ks_layer",
      TRUE ~ NA_character_
    ),
    interface = str_extract(variable, "\\d+_\\d+"),
    layer = str_extract(variable, "(?<=layer)\\d+"),
    model_name = factor(model_name, levels = model_levels, labels = model_labels)
  )

# -----------------------------
# Separate interface-level and layer-level variables
# -----------------------------
df_interface <- flux_long %>%
  filter(!is.na(metric), !is.na(interface)) %>%
  mutate(
    interface = factor(interface, levels = unique(interface))
  )

df_layer <- flux_long %>%
  filter(!is.na(metric), !is.na(layer)) %>%
  mutate(
    layer = factor(layer, levels = unique(layer))
  )

# -----------------------------
# Plot functions
# -----------------------------
plot_interface_faceted <- function(data, metric_name, title_txt, y_lab, add_zero = FALSE) {
  p <- ggplot(
    data %>% filter(metric == metric_name),
    aes(x = sim_year, y = value)
  ) +
    geom_line(linewidth = 0.6) +
    facet_grid(model_name ~ interface, scales = "free_y") +
    labs(
      title = title_txt,
      x = "Time [years]",
      y = y_lab
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  if (add_zero) {
    p <- p + geom_hline(yintercept = 0, linetype = 2)
  }
  
  p
}

plot_layer_faceted <- function(data, metric_name, title_txt, y_lab, add_zero = FALSE) {
  p <- ggplot(
    data %>% filter(metric == metric_name),
    aes(x = sim_year, y = value)
  ) +
    geom_line(linewidth = 0.6) +
    facet_grid(model_name ~ layer, scales = "free_y") +
    labs(
      title = title_txt,
      x = "Time [years]",
      y = y_lab
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  if (add_zero) {
    p <- p + geom_hline(yintercept = 0, linetype = 2)
  }
  
  p
}

# -----------------------------
# Interface plots
# -----------------------------
p_flux <- plot_interface_faceted(
  df_interface,
  metric_name = "mean_flux",
  title_txt = "Mean flux by interface",
  y_lab = "Mean flux",
  add_zero = TRUE
)

p_abs_flux <- plot_interface_faceted(
  df_interface,
  metric_name = "mean_abs_flux",
  title_txt = "Mean absolute flux by interface",
  y_lab = "Mean absolute flux"
)

p_dphi <- plot_interface_faceted(
  df_interface,
  metric_name = "mean_delta_swp",
  title_txt = "Mean delta SWP by interface",
  y_lab = "Mean delta SWP",
  add_zero = TRUE
)

p_kharm <- plot_interface_faceted(
  df_interface,
  metric_name = "mean_ks_harmonic",
  title_txt = "Mean harmonic Ks by interface",
  y_lab = "Mean harmonic Ks"
)

p_net <- plot_interface_faceted(
  df_interface,
  metric_name = "net_volumetric_change",
  title_txt = "Net volumetric change by interface",
  y_lab = "Net volumetric change",
  add_zero = TRUE
)

p_gross <- plot_interface_faceted(
  df_interface,
  metric_name = "gross_volumetric_change",
  title_txt = "Gross volumetric change by interface",
  y_lab = "Gross volumetric change"
)

# -----------------------------
# Layer plots
# -----------------------------
p_swp_layer <- plot_layer_faceted(
  df_layer,
  metric_name = "mean_swp_layer",
  title_txt = "Mean SWP by layer",
  y_lab = "Mean SWP"
)

p_ks_layer <- plot_layer_faceted(
  df_layer,
  metric_name = "mean_ks_layer",
  title_txt = "Mean Ks by layer",
  y_lab = "Mean Ks"
)

# -----------------------------
# Combined panel - interface variables
# -----------------------------
df_panel_interface <- df_interface %>%
  filter(metric %in% c(
    "mean_flux",
    "mean_abs_flux",
    "mean_delta_swp",
    "mean_ks_harmonic",
    "net_volumetric_change",
    "gross_volumetric_change"
  )) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "mean_flux",
        "mean_abs_flux",
        "mean_delta_swp",
        "mean_ks_harmonic",
        "net_volumetric_change",
        "gross_volumetric_change"
      ),
      labels = c(
        "Mean flux",
        "Mean absolute flux",
        "Mean delta SWP",
        "Mean harmonic Ks",
        "Net volumetric change",
        "Gross volumetric change"
      )
    )
  )

p_panel_interface <- ggplot(
  df_panel_interface,
  aes(x = sim_year, y = value)
) +
  geom_line(linewidth = 0.4) +
  facet_grid(metric + model_name ~ interface, scales = "free_y") +
  labs(
    title = "Vertical water flux diagnostics by configuration and interface",
    x = "Time [years]",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# Combined panel - layer variables
# -----------------------------
df_panel_layer <- df_layer %>%
  filter(metric %in% c("mean_swp_layer", "mean_ks_layer")) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("mean_swp_layer", "mean_ks_layer"),
      labels = c("Mean SWP", "Mean Ks")
    )
  )

p_panel_layer <- ggplot(
  df_panel_layer,
  aes(x = sim_year, y = value)
) +
  geom_line(linewidth = 0.5) +
  facet_grid(metric + model_name ~ layer, scales = "free_y") +
  labs(
    title = "Layer-level diagnostics by configuration and layer",
    x = "Time [years]",
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# Save plots
# -----------------------------
if (save_plots) {
  ggsave(file.path(plot_dir, "mean_flux_by_interface_faceted.png"), p_flux, width = 14, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "mean_abs_flux_by_interface_faceted.png"), p_abs_flux, width = 14, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "mean_delta_swp_by_interface_faceted.png"), p_dphi, width = 14, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "mean_ks_harmonic_by_interface_faceted.png"), p_kharm, width = 14, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "net_volumetric_change_by_interface_faceted.png"), p_net, width = 14, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "gross_volumetric_change_by_interface_faceted.png"), p_gross, width = 14, height = 8, dpi = 300)
  
  ggsave(file.path(plot_dir, "mean_swp_by_layer_faceted.png"), p_swp_layer, width = 14, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "mean_ks_by_layer_faceted.png"), p_ks_layer, width = 14, height = 8, dpi = 300)
  
  ggsave(file.path(plot_dir, "vertical_water_flux_panel_interface_faceted.png"), p_panel_interface, width = 18, height = 14, dpi = 300)
  ggsave(file.path(plot_dir, "vertical_water_flux_panel_layer_faceted.png"), p_panel_layer, width = 16, height = 8, dpi = 300)
}

# -----------------------------
# Print plots
# -----------------------------
p_flux
p_abs_flux
p_dphi
p_kharm
p_net
p_gross
p_swp_layer
p_ks_layer
p_panel_interface
p_panel_layer

########################
###### Checking ks and ks harmonics = or ~zero ##########
library(dplyr)

df <- df_layer

ks_vars <- c(
  "mean_ks_layer0",
  "mean_ks_layer1",
  "mean_ks_layer2",
  "mean_ks_layer3",
  "mean_ks_layer4",
  "mean_ks_layer5"
)

ks_zero_check <- df %>%
  filter(variable %in% ks_vars) %>%
  group_by(model_name, variable) %>%
  summarise(
    n_total = n(),
    n_non_na = sum(!is.na(value)),
    n_zero_exact = sum(value == 0, na.rm = TRUE),
    pct_zero_exact = 100 * n_zero_exact / n_total,
    n_negative = sum(value < 0, na.rm = TRUE),
    n_na = sum(is.na(value)),
    min_val = if (all(is.na(value))) NA_real_ else min(value, na.rm = TRUE),
    max_val = if (all(is.na(value))) NA_real_ else max(value, na.rm = TRUE),
    min_positive = {
      x <- value[value > 0 & !is.na(value)]
      if (length(x) == 0) NA_real_ else min(x)
    },
    n_lt_1e20 = sum(value > 0 & value < 1e-20, na.rm = TRUE),
    n_lt_1e15 = sum(value > 0 & value < 1e-15, na.rm = TRUE),
    n_lt_1e12 = sum(value > 0 & value < 1e-12, na.rm = TRUE),
    n_lt_1e10 = sum(value > 0 & value < 1e-10, na.rm = TRUE),
    n_lt_1e8  = sum(value > 0 & value < 1e-8,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(variable, model_name)

ks_zero_check

df %>%
  filter(variable %in% ks_vars) %>%
  group_by(model_name) %>%
  summarise(
    min_value_all_layers = min(value, na.rm = TRUE),
    n_zero_exact_all_layers = sum(value == 0, na.rm = TRUE),
    min_positive_all_layers = min(value[value > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(width = Inf)


##### ks harmonic
library(dplyr)

ksh_vars <- c(
  "mean_ks_harmonic_layers0_1",
  "mean_ks_harmonic_layers1_2",
  "mean_ks_harmonic_layers2_3",
  "mean_ks_harmonic_layers3_4",
  "mean_ks_harmonic_layers4_5"
)

ksh_summary_all <- df_interface %>%
  filter(variable %in% ksh_vars) %>%
  group_by(model_name) %>%
  summarise(
    min_value_all_interfaces = min(value, na.rm = TRUE),
    n_zero_exact_all_interfaces = sum(value == 0, na.rm = TRUE),
    min_positive_all_interfaces = min(value[value > 0], na.rm = TRUE),
    .groups = "drop"
  )

ksh_summary_all

ksh_by_interface <- df_interface %>%
  filter(variable %in% ksh_vars) %>%
  mutate(interface_id = str_extract(variable, "\\d_\\d$")) %>%
  group_by(model_name, interface_id, variable) %>%
  summarise(
    n_total = n(),
    n_non_na = sum(!is.na(value)),
    n_zero_exact = sum(value == 0, na.rm = TRUE),
    pct_zero_exact = 100 * n_zero_exact / n_total,
    n_negative = sum(value < 0, na.rm = TRUE),
    n_na = sum(is.na(value)),
    min_val = if (all(is.na(value))) NA_real_ else min(value, na.rm = TRUE),
    max_val = if (all(is.na(value))) NA_real_ else max(value, na.rm = TRUE),
    min_positive = {
      x <- value[value > 0 & !is.na(value)]
      if (length(x) == 0) NA_real_ else min(x)
    },
    n_lt_1e20 = sum(value > 0 & value < 1e-20, na.rm = TRUE),
    n_lt_1e15 = sum(value > 0 & value < 1e-15, na.rm = TRUE),
    n_lt_1e12 = sum(value > 0 & value < 1e-12, na.rm = TRUE),
    n_lt_1e10 = sum(value > 0 & value < 1e-10, na.rm = TRUE),
    n_lt_1e8  = sum(value > 0 & value < 1e-8,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(interface_id, model_name)

ksh_by_interface

df_interface %>%
  filter(variable %in% ksh_vars) %>%
  arrange(value) %>%
  select(model_name, date, variable, interface, value) %>%
  head(40)

######## o que acontece quando Ksh está no mínimo
library(dplyr)
library(stringr)
library(purrr)

interfaces <- c("0_1", "1_2", "2_3", "3_4", "4_5")

diag_floor <- map_dfr(interfaces, function(intf) {
  ksh_var  <- paste0("mean_ks_harmonic_layers", intf)
  flux_var <- paste0("mean_abs_flux_layers", intf)
  swp_var  <- paste0("mean_delta_swp_layers", intf)
  
  ksh <- df_interface %>%
    filter(variable == ksh_var) %>%
    select(model_name, date, interface, ksh = value)
  
  flux <- df_interface %>%
    filter(variable == flux_var) %>%
    select(model_name, date, abs_flux = value)
  
  swp <- df_interface %>%
    filter(variable == swp_var) %>%
    select(model_name, date, delta_swp = value)
  
  joined <- ksh %>%
    left_join(flux, by = c("model_name", "date")) %>%
    left_join(swp, by = c("model_name", "date"))
  
  joined %>%
    group_by(model_name, interface) %>%
    summarise(
      min_ksh = min(ksh, na.rm = TRUE),
      n_at_floor = sum(ksh == min(ksh, na.rm = TRUE), na.rm = TRUE),
      mean_abs_flux_at_floor = mean(abs_flux[ksh == min(ksh, na.rm = TRUE)], na.rm = TRUE),
      median_abs_flux_at_floor = median(abs_flux[ksh == min(ksh, na.rm = TRUE)], na.rm = TRUE),
      mean_delta_swp_at_floor = mean(delta_swp[ksh == min(ksh, na.rm = TRUE)], na.rm = TRUE),
      median_delta_swp_at_floor = median(delta_swp[ksh == min(ksh, na.rm = TRUE)], na.rm = TRUE),
      .groups = "drop"
    )
})

diag_floor

