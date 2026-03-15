suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------
# Params
# -----------------------------
out_prefix <- "(null)"
save_plots <- TRUE   # mude para TRUE se quiser salvar png
start_date <- as.Date("2004-01-01")

# -----------------------------
# Base path
# -----------------------------
base_dir <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/hydraulic_locking_tests/noveg"

# -----------------------------
# Runs to compare
# -----------------------------
files_df <- tibble(
  model_name = c(
    "baseline_noveg",
    "legacy_theta0_reset1e-3_noveg",
    "theta_thres1e-6_noveg",
    "theta_thres1e-5_noveg",
    "theta_thres1e-4_noveg",
    "theta_thres1e-3_noveg"
  )
) %>%
  mutate(
    experiment_dir = file.path(base_dir, model_name, "output"),
    debug_hyd_path = file.path(experiment_dir, paste0(out_prefix, "_0_debug_hydraulics_by_iter.txt"))
  )

print(files_df)

# -----------------------------
# Check files
# -----------------------------
missing_files <- files_df %>%
  filter(!file.exists(debug_hyd_path))

if (nrow(missing_files) > 0) {
  stop(
    "These debug files were not found:\n",
    paste(missing_files$debug_hyd_path, collapse = "\n")
  )
}

# -----------------------------
# Reader for old broken header
# -----------------------------
read_debug_fixed <- function(path) {
  col_names <- c(
    "iter",
    "total_voxels",
    "total_interfaces",
    "n_theta_cap_lt_1e6",
    "n_theta_cap_lt_1e5",
    "n_theta_cap_lt_1e4",
    "n_theta_cap_lt_1e3",
    "n_ks_cap_lt_1e12",
    "n_ks_cap_lt_1e10",
    "n_ks_cap_lt_1e8",
    "n_ks_cap_eq_0",
    "n_ksh_cap_lt_1e12",
    "n_ksh_cap_lt_1e10",
    "n_ksh_cap_lt_1e8",
    "n_ksh_cap_eq_0"
  )
  
  readr::read_tsv(
    path,
    skip = 1,
    col_names = col_names,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    dplyr::select(-dplyr::matches("^X[0-9]+$"))
}

# -----------------------------
# Read all runs
# -----------------------------
debug_all <- purrr::map_dfr(seq_len(nrow(files_df)), function(i) {
  df <- read_debug_fixed(files_df$debug_hyd_path[i])
  
  if (!all(c("iter", "total_voxels", "total_interfaces") %in% names(df))) {
    stop(
      "Columns 'iter', 'total_voxels' and/or 'total_interfaces' not found in: ",
      files_df$debug_hyd_path[i]
    )
  }
  
  df %>%
    mutate(
      model_name = files_df$model_name[i],
      date = start_date + iter
    )
})

# -----------------------------
# Long format
# -----------------------------
debug_long <- debug_all %>%
  pivot_longer(
    cols = -c(model_name, iter, total_voxels, total_interfaces, date),
    names_to = "metric",
    values_to = "n_voxels"
  ) %>%
  mutate(
    group = case_when(
      str_detect(metric, "^n_theta_cap") ~ "theta",
      str_detect(metric, "^n_ks_cap")    ~ "Ks",
      str_detect(metric, "^n_ksh_cap")   ~ "Ksh",
      TRUE ~ "other"
    ),
    pct_voxels = 100 * n_voxels / total_voxels,
    pct_interfaces = 100 * n_voxels / total_interfaces,
    sim_year = as.numeric(date - start_date) / 365
  ) %>%
  filter(group %in% c("theta", "Ks", "Ksh"))

# -----------------------------
# Pretty labels
# -----------------------------
metric_labels <- c(
  "n_theta_cap_lt_1e6" = "theta < 1e-6",
  "n_theta_cap_lt_1e5" = "theta < 1e-5",
  "n_theta_cap_lt_1e4" = "theta < 1e-4",
  "n_theta_cap_lt_1e3" = "theta < 1e-3",
  
  "n_ks_cap_eq_0"      = "Ks = 0",
  "n_ks_cap_lt_1e12"   = "Ks < 1e-12",
  "n_ks_cap_lt_1e10"   = "Ks < 1e-10",
  "n_ks_cap_lt_1e8"    = "Ks < 1e-8",
  
  "n_ksh_cap_eq_0"     = "Ksh = 0",
  "n_ksh_cap_lt_1e12"  = "Ksh < 1e-12",
  "n_ksh_cap_lt_1e10"  = "Ksh < 1e-10",
  "n_ksh_cap_lt_1e8"   = "Ksh < 1e-8"
)

debug_long <- debug_long %>%
  mutate(metric_pretty = recode(metric, !!!metric_labels))

# -----------------------------
# Split data
# -----------------------------
theta_df <- debug_long %>% filter(group == "theta")
ks_df    <- debug_long %>% filter(group == "Ks")
ksh_df   <- debug_long %>% filter(group == "Ksh")

# -----------------------------
# Plot function
# -----------------------------
plot_thresholds <- function(df, title_txt, y_var = c("pct_voxels", "n_voxels")) {
  y_var <- match.arg(y_var)
  
  y_lab <- if (y_var == "pct_voxels") "% of voxels" else "Number of voxels"
  
  ggplot(df, aes(x = sim_year, y = .data[[y_var]], color = metric_pretty)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ model_name, ncol = 2) +
    labs(
      x = "Time [years]",
      y = y_lab,
      color = NULL,
      title = title_txt
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

# -----------------------------
# Main plots (% of voxels)
# -----------------------------
p_theta_pct <- plot_thresholds(
  theta_df,
  title_txt = "Theta thresholds across configurations",
  y_var = "pct_voxels"
)

p_ks_pct <- plot_thresholds(
  ks_df,
  title_txt = "Ks thresholds across configurations",
  y_var = "pct_voxels"
)

p_ksh_pct <- plot_thresholds(
  ksh_df,
  title_txt = "Ksh thresholds across configurations",
  y_var = "pct_voxels"
)

print(p_theta_pct)
print(p_ks_pct)
print(p_ksh_pct)

# -----------------------------
# Optional: absolute counts
# -----------------------------
p_theta_n <- plot_thresholds(
  theta_df,
  title_txt = "Theta thresholds across configurations (absolute counts)",
  y_var = "n_voxels"
)

p_ks_n <- plot_thresholds(
  ks_df,
  title_txt = "Ks thresholds across configurations (absolute counts)",
  y_var = "n_voxels"
)

p_ksh_n <- plot_thresholds(
  ksh_df,
  title_txt = "Ksh thresholds across configurations (absolute counts)",
  y_var = "n_voxels"
)

# Uncomment if you want to see them too
# print(p_theta_n)
# print(p_ks_n)
# print(p_ksh_n)

# -----------------------------
# Optional save
# -----------------------------
if (save_plots) {
  ggsave("theta_thresholds_pct.png", p_theta_pct, width = 13, height = 8, dpi = 300)
  ggsave("ks_thresholds_pct.png",    p_ks_pct,    width = 13, height = 8, dpi = 300)
  ggsave("ksh_thresholds_pct.png",   p_ksh_pct,   width = 13, height = 8, dpi = 300)
  
  ggsave("theta_thresholds_n.png", p_theta_n, width = 13, height = 8, dpi = 300)
  ggsave("ks_thresholds_n.png",    p_ks_n,    width = 13, height = 8, dpi = 300)
  ggsave("ksh_thresholds_n.png",   p_ksh_n,   width = 13, height = 8, dpi = 300)
}