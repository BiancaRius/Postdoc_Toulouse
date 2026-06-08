# ============================================================
# Multitrait syndrome analysis
# Functional trait PCA for final communities
# ============================================================

# Load required packages
library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)

# ------------------------------------------------------------
# 1. Define traits, scenarios and labels
# ------------------------------------------------------------

traits_for_pca <- c(
  "wsg",
  "LMA",
  "Nmass",
  "Pmass",
  "tlp",
  "dbhmax",
  "hmax",
  "leafarea",
  "g1"
)

trait_titles <- c(
  "wsg" = "WD",
  "LMA" = "LMA",
  "Nmass" = "Nmass",
  "Pmass" = "Pmass",
  "tlp" = "TLP",
  "dbhmax" = "DBHmax",
  "hmax" = "Hmax",
  "leafarea" = "Leaf area",
  "g1" = "g1"
)

scenarios_keep <- c(
  "deep_clayey_regclim",
  "shallow_sandy_regclim"
)

scenario_labels <- c(
  "deep_clayey_regclim" = "Deep WT | Clayey soil",
  "shallow_sandy_regclim" = "Shallow WT | Sandy soil"
)

dbh_keep <- "DBH > 10 cm"

# ------------------------------------------------------------
# 2. Prepare wide trait table
# ------------------------------------------------------------
base_dir <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke"

scenario_table <- tibble(
  scenario = c(
    "deep_clayey_regclim",
    "deep_clayey_regclim_ks12",
    "shallow_sandy_regclim",
    "shallow_sandy_regclim_ks12"
  ),
  
  run_dir = c(
    "deepWT_clayeysoil/ducke_deep_clayey_regclim",
    "deepWT_clayeysoil/ducke_deep_clayey_regclim_ks12",
    "shallowWT_sandysoil/ducke_shallow_sandy_regclim",
    "shallowWT_sandysoil/ducke_shallow_sandy_regclim_ks12"
  )
) %>%
  mutate(
    path = file.path(base_dir, run_dir),
    final_pattern_file = file.path(path, "output", "(null)_0_final_pattern.txt")
  )


df_all_raw <- map2_dfr(
  scenario_table$scenario,
  scenario_table$final_pattern_file,
  function(scenario_i, file_i) {
    
    read.table(file_i, header = TRUE) %>%
      as_tibble() %>%
      mutate(
        scenario = scenario_i,
        entity_id = paste(scenario_i, row_number(), sep = "_")
      )
  }
)

# ------------------------------------------------------------
# 2. Add weights and DBH thresholds
# ------------------------------------------------------------

df_all <- df_all_raw %>%
  mutate(
    basal_area = pi * (dbh / 2)^2,
    abundance_weight = 1,
    basal_area_weight = basal_area,
    biomass_weight = AGB
  )

df_threshold <- bind_rows(
  df_all %>%
    filter(dbh > 0.01) %>%
    mutate(dbh_threshold = "DBH > 1 cm"),
  
  df_all %>%
    filter(dbh > 0.10) %>%
    mutate(dbh_threshold = "DBH > 10 cm")
)

# ------------------------------------------------------------
# 4. Prepare PCA input table
# ------------------------------------------------------------

df_pca_input <- df_threshold %>%
  filter(
    dbh_threshold == dbh_keep,
    scenario %in% scenarios_keep
  ) %>%
  select(
    scenario,
    dbh_threshold,
    entity_id,
    abundance_weight,
    basal_area_weight,
    biomass_weight,
    all_of(traits_for_pca)
  ) %>%
  drop_na(all_of(traits_for_pca))

# ------------------------------------------------------------
# 6. Run PCA on standardized traits
# ------------------------------------------------------------

trait_matrix <- df_pca_input %>%
  select(all_of(traits_for_pca)) %>%
  as.matrix()

pca_traits <- prcomp(
  trait_matrix,
  center = TRUE,
  scale. = TRUE
)

pca_variance <- tibble(
  PC = paste0("PC", seq_along(pca_traits$sdev)),
  eigenvalue = pca_traits$sdev^2,
  variance_explained = eigenvalue / sum(eigenvalue),
  cumulative_variance = cumsum(variance_explained)
)

pca_variance

# ------------------------------------------------------------
# 7. Extract PCA scores and loadings
# ------------------------------------------------------------

pca_scores <- df_pca_input %>%
  bind_cols(as_tibble(pca_traits$x)) %>%
  mutate(
    scenario_label = recode(scenario, !!!scenario_labels)
  )

pca_loadings <- as_tibble(
  pca_traits$rotation,
  rownames = "trait"
) %>%
  mutate(
    trait_label = recode(trait, !!!trait_titles)
  )

pca_loadings

# ------------------------------------------------------------
# 8. Inspect PCA loadings
# ------------------------------------------------------------

# Print loadings for the first three PCA axes
pca_loadings %>%
  select(trait, trait_label, PC1, PC2, PC3)

# ------------------------------------------------------------
# 9. Prepare loadings in long format
# ------------------------------------------------------------

pca_loadings_long <- pca_loadings %>%
  select(trait, trait_label, PC1, PC2, PC3) %>%
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "loading"
  )

pca_loadings_long

# ------------------------------------------------------------
# 10. Plot trait loadings for PC1, PC2 and PC3
# ------------------------------------------------------------

ggplot(
  pca_loadings_long,
  aes(x = fct_reorder(trait_label, loading), y = loading)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_col(alpha = 0.75) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Trait loadings on the first PCA axes",
    x = NULL,
    y = "Loading"
  )

# ------------------------------------------------------------
# 11. Convert weights to long format
# ------------------------------------------------------------

pca_scores_weighted <- pca_scores %>%
  pivot_longer(
    cols = c(abundance_weight, basal_area_weight, biomass_weight),
    names_to = "weighting",
    values_to = "weight"
  ) %>%
  mutate(
    weighting = recode(
      weighting,
      abundance_weight = "Abundance-weighted",
      basal_area_weight = "Basal-area-weighted",
      biomass_weight = "Biomass-weighted"
    )
  )

pca_scores_weighted

# ------------------------------------------------------------
# 12. Calculate weighted community centroids in PCA space
# ------------------------------------------------------------

weighted_mean_safe <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  weighted.mean(x[keep], w[keep])
}

cwm_pca <- pca_scores_weighted %>%
  group_by(weighting, scenario, scenario_label) %>%
  summarise(
    CWM_PC1 = weighted_mean_safe(PC1, weight),
    CWM_PC2 = weighted_mean_safe(PC2, weight),
    CWM_PC3 = weighted_mean_safe(PC3, weight),
    .groups = "drop"
  )

cwm_pca

# ------------------------------------------------------------
# 13. Plot multitrait functional space
# ------------------------------------------------------------

ggplot(
  pca_scores_weighted,
  aes(
    x = PC1,
    y = PC2,
    color = scenario_label
  )
) +
  geom_point(
    aes(size = sqrt(weight)),
    alpha = 0.08
  ) +
  geom_point(
    data = cwm_pca,
    aes(
      x = CWM_PC1,
      y = CWM_PC2,
      color = scenario_label
    ),
    size = 4,
    shape = 21,
    fill = "white",
    stroke = 1.2,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ weighting) +
  theme_bw() +
  labs(
    title = "Multitrait functional space - DBH > 10 cm",
    x = paste0(
      "PC1 (",
      round(pca_variance$variance_explained[pca_variance$PC == "PC1"] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance$variance_explained[pca_variance$PC == "PC2"] * 100, 1),
      "%)"
    ),
    color = "WT | Pedology",
    size = "sqrt(weight)"
  )

# ------------------------------------------------------------
# 14. Compare weighted PCA centroids between scenarios
# ------------------------------------------------------------

pca_centroid_difference <- cwm_pca %>%
  select(weighting, scenario, CWM_PC1, CWM_PC2, CWM_PC3) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(CWM_PC1, CWM_PC2, CWM_PC3)
  ) %>%
  mutate(
    diff_PC1_shallow_minus_deep =
      CWM_PC1_shallow_sandy_regclim - CWM_PC1_deep_clayey_regclim,
    diff_PC2_shallow_minus_deep =
      CWM_PC2_shallow_sandy_regclim - CWM_PC2_deep_clayey_regclim,
    diff_PC3_shallow_minus_deep =
      CWM_PC3_shallow_sandy_regclim - CWM_PC3_deep_clayey_regclim
  )

pca_centroid_difference


# ------------------------------------------------------------
# 15. Weighted quantile function
# ------------------------------------------------------------

weighted_quantile <- function(x, w, prob) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  cw <- cumsum(w) / sum(w)
  x[which(cw >= prob)[1]]
}

# ------------------------------------------------------------
# 16. Weighted boxplot summary for PC1
# ------------------------------------------------------------

df_pc1_boxplot <- pca_scores_weighted %>%
  group_by(weighting, scenario, scenario_label) %>%
  summarise(
    ymin = weighted_quantile(PC1, weight, 0.10),
    lower = weighted_quantile(PC1, weight, 0.25),
    middle = weighted_quantile(PC1, weight, 0.50),
    upper = weighted_quantile(PC1, weight, 0.75),
    ymax = weighted_quantile(PC1, weight, 0.90),
    weighted_mean = weighted_mean_safe(PC1, weight),
    .groups = "drop"
  )

ggplot(
  df_pc1_boxplot,
  aes(x = scenario_label, fill = scenario_label)
) +
  geom_boxplot(
    aes(
      ymin = ymin,
      lower = lower,
      middle = middle,
      upper = upper,
      ymax = ymax
    ),
    stat = "identity",
    width = 0.55,
    alpha = 0.65
  ) +
  geom_point(
    aes(y = weighted_mean),
    shape = 21,
    size = 2.7,
    fill = "white"
  ) +
  facet_wrap(~ weighting) +
  theme_bw() +
  labs(
    title = "Weighted distribution of PC1 scores",
    x = NULL,
    y = "PC1 score",
    fill = "WT | Pedology"
  ) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

# ------------------------------------------------------------
# 17. Weighted boxplots for PC1, PC2 and PC3
# ------------------------------------------------------------

df_pca_scores_long <- pca_scores_weighted %>%
  select(
    scenario,
    scenario_label,
    weighting,
    weight,
    PC1,
    PC2,
    PC3
  ) %>%
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "score"
  )

df_pca_boxplot <- df_pca_scores_long %>%
  group_by(weighting, scenario, scenario_label, PC) %>%
  summarise(
    ymin = weighted_quantile(score, weight, 0.10),
    lower = weighted_quantile(score, weight, 0.25),
    middle = weighted_quantile(score, weight, 0.50),
    upper = weighted_quantile(score, weight, 0.75),
    ymax = weighted_quantile(score, weight, 0.90),
    weighted_mean = weighted_mean_safe(score, weight),
    .groups = "drop"
  )

ggplot(
  df_pca_boxplot,
  aes(x = scenario_label, fill = scenario_label)
) +
  geom_boxplot(
    aes(
      ymin = ymin,
      lower = lower,
      middle = middle,
      upper = upper,
      ymax = ymax
    ),
    stat = "identity",
    width = 0.55,
    alpha = 0.65
  ) +
  geom_point(
    aes(y = weighted_mean),
    shape = 21,
    size = 2.4,
    fill = "white"
  ) +
  facet_grid(PC ~ weighting, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Weighted distribution of PCA scores",
    x = NULL,
    y = "PCA score",
    fill = "WT | Pedology"
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# ------------------------------------------------------------
# 18. Helper functions to quantify functional space
# ------------------------------------------------------------

convex_hull_area_2d <- function(df, axes = c("PC1", "PC2")) {
  
  xy <- df %>%
    select(all_of(axes)) %>%
    drop_na() %>%
    distinct()
  
  if (nrow(xy) < 3) {
    return(NA_real_)
  }
  
  h <- chull(xy[[axes[1]]], xy[[axes[2]]])
  
  x <- xy[[axes[1]]][c(h, h[1])]
  y <- xy[[axes[2]]][c(h, h[1])]
  
  area <- 0.5 * abs(sum(x[-1] * y[-length(y)] - x[-length(x)] * y[-1]))
  
  return(area)
}


weighted_centroid <- function(mat, w) {
  
  keep <- complete.cases(mat) & is.finite(w) & w > 0
  mat <- as.matrix(mat[keep, , drop = FALSE])
  w <- w[keep]
  
  if (nrow(mat) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }
  
  w <- w / sum(w)
  colSums(mat * w)
}


weighted_fdis <- function(df, axes = c("PC1", "PC2"), weight_col = "weight") {
  
  mat <- df %>%
    select(all_of(axes)) %>%
    as.matrix()
  
  w <- df[[weight_col]]
  
  keep <- complete.cases(mat) & is.finite(w) & w > 0
  mat <- mat[keep, , drop = FALSE]
  w <- w[keep]
  
  if (nrow(mat) < 2) {
    return(NA_real_)
  }
  
  w <- w / sum(w)
  centroid <- colSums(mat * w)
  
  distances <- sqrt(rowSums(sweep(mat, 2, centroid, "-")^2))
  
  sum(w * distances)
}


weighted_cov_matrix <- function(df, axes = c("PC1", "PC2"), weight_col = "weight") {
  
  mat <- df %>%
    select(all_of(axes)) %>%
    as.matrix()
  
  w <- df[[weight_col]]
  
  keep <- complete.cases(mat) & is.finite(w) & w > 0
  mat <- mat[keep, , drop = FALSE]
  w <- w[keep]
  
  if (nrow(mat) < 2) {
    return(matrix(NA_real_, nrow = length(axes), ncol = length(axes)))
  }
  
  w <- w / sum(w)
  centroid <- colSums(mat * w)
  mat_centered <- sweep(mat, 2, centroid, "-")
  
  cov_mat <- t(mat_centered) %*% (mat_centered * w)
  
  return(cov_mat)
}


weighted_ellipse_area_2d <- function(df, axes = c("PC1", "PC2"), weight_col = "weight", level = 0.95) {
  
  cov_mat <- weighted_cov_matrix(df, axes = axes, weight_col = weight_col)
  
  if (any(!is.finite(cov_mat))) {
    return(NA_real_)
  }
  
  det_cov <- det(cov_mat)
  
  if (!is.finite(det_cov) || det_cov <= 0) {
    return(NA_real_)
  }
  
  chi_val <- qchisq(level, df = 2)
  
  area <- pi * chi_val * sqrt(det_cov)
  
  return(area)
}

# ------------------------------------------------------------
# 19. Quantify functional space in PC1-PC2
# ------------------------------------------------------------

functional_space_2d <- pca_scores_weighted %>%
  group_by(weighting, scenario, scenario_label) %>%
  group_modify(~ {
    
    tibble(
      n_entities = nrow(.x),
      convex_hull_area_PC1_PC2 = convex_hull_area_2d(.x, axes = c("PC1", "PC2")),
      weighted_FDis_PC1_PC2 = weighted_fdis(.x, axes = c("PC1", "PC2"), weight_col = "weight"),
      weighted_ellipse_area95_PC1_PC2 = weighted_ellipse_area_2d(.x, axes = c("PC1", "PC2"), weight_col = "weight", level = 0.95)
    )
  }) %>%
  ungroup()

functional_space_2d

# ------------------------------------------------------------
# 20. Compare functional space between scenarios
# ------------------------------------------------------------

functional_space_comparison <- functional_space_2d %>%
  select(
    weighting,
    scenario,
    convex_hull_area_PC1_PC2,
    weighted_FDis_PC1_PC2,
    weighted_ellipse_area95_PC1_PC2
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(
      convex_hull_area_PC1_PC2,
      weighted_FDis_PC1_PC2,
      weighted_ellipse_area95_PC1_PC2
    )
  ) %>%
  mutate(
    diff_hull_shallow_minus_deep =
      convex_hull_area_PC1_PC2_shallow_sandy_regclim -
      convex_hull_area_PC1_PC2_deep_clayey_regclim,
    
    diff_FDis_shallow_minus_deep =
      weighted_FDis_PC1_PC2_shallow_sandy_regclim -
      weighted_FDis_PC1_PC2_deep_clayey_regclim,
    
    diff_ellipse_shallow_minus_deep =
      weighted_ellipse_area95_PC1_PC2_shallow_sandy_regclim -
      weighted_ellipse_area95_PC1_PC2_deep_clayey_regclim,
    
    ratio_hull_shallow_deep =
      convex_hull_area_PC1_PC2_shallow_sandy_regclim /
      convex_hull_area_PC1_PC2_deep_clayey_regclim,
    
    ratio_FDis_shallow_deep =
      weighted_FDis_PC1_PC2_shallow_sandy_regclim /
      weighted_FDis_PC1_PC2_deep_clayey_regclim,
    
    ratio_ellipse_shallow_deep =
      weighted_ellipse_area95_PC1_PC2_shallow_sandy_regclim /
      weighted_ellipse_area95_PC1_PC2_deep_clayey_regclim
  )

functional_space_comparison

# ------------------------------------------------------------
# 21. Plot functional space metrics
# ------------------------------------------------------------

functional_space_plot <- functional_space_2d %>%
  pivot_longer(
    cols = c(
      convex_hull_area_PC1_PC2,
      weighted_FDis_PC1_PC2,
      weighted_ellipse_area95_PC1_PC2
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      convex_hull_area_PC1_PC2 = "Convex hull area",
      weighted_FDis_PC1_PC2 = "Weighted FDis",
      weighted_ellipse_area95_PC1_PC2 = "Weighted 95% ellipse area"
    )
  )

ggplot(
  functional_space_plot,
  aes(x = scenario_label, y = value, fill = scenario_label)
) +
  geom_col(alpha = 0.75, width = 0.65) +
  facet_grid(metric ~ weighting, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Functional space size in PC1-PC2",
    x = NULL,
    y = "Functional space metric",
    fill = "WT | Pedology"
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )
