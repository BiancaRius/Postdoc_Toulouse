

## TRAITS DISTRIBUTION ##

library(dplyr)
library(ggplot2)
library(tidyverse)

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


df_all_raw <- scenario_table %>%
  mutate(data = map(final_pattern_file, ~ read.table(.x, header = TRUE))) %>%
  select(scenario, data) %>%
  unnest(data)

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

trait_vars <- c(
  "LMA",
  "wsg",
  "Nmass",
  "Pmass",
  "tlp",
  "dbhmax",
  "hmax",
  "leafarea",
  "g1"
)

trait_vars <- trait_vars[trait_vars %in% names(df_threshold)]

cwm_traits <- df_threshold %>%
  mutate(tree_id = row_number()) %>%
  pivot_longer(
    cols = all_of(trait_vars),
    names_to = "trait",
    values_to = "trait_value"
  ) %>%
  group_by(scenario, dbh_threshold, trait) %>%
  summarise(
    n_trees = n_distinct(tree_id[!is.na(trait_value)]),
    
    CWM_abundance = weighted.mean(trait_value, abundance_weight),
    CWM_basal_area = weighted.mean(trait_value, basal_area_weight),
    CWM_biomass = weighted.mean(trait_value, biomass_weight),
    
    .groups = "drop"
  )

cwm_traits_long <- cwm_traits %>%
  pivot_longer(
    cols = starts_with("CWM_"),
    names_to = "weighting",
    values_to = "CWM"
  ) %>%
  mutate(
    weighting = recode(
      weighting,
      "CWM_abundance" = "Abundance",
      "CWM_basal_area" = "Basal area",
      "CWM_biomass" = "Biomass"
    )
  )

df_trait_dist <- df_threshold %>%
  select(
    scenario,
    dbh_threshold,
    all_of(trait_vars),
    abundance_weight,
    basal_area_weight,
    biomass_weight
  ) %>%
  pivot_longer(
    cols = all_of(trait_vars),
    names_to = "trait",
    values_to = "trait_value"
  ) %>%
  pivot_longer(
    cols = c(abundance_weight, basal_area_weight, biomass_weight),
    names_to = "weighting",
    values_to = "weight"
  ) %>%
  mutate(
    weighting = recode(
      weighting,
      "abundance_weight" = "Abundance",
      "basal_area_weight" = "Basal area",
      "biomass_weight" = "Biomass"
    )
  ) %>%
  filter(
    !is.na(trait_value),
    !is.na(weight),
    weight > 0
  )

# select dbh > 10
cwm_dbh_gt_10 <- df_trait_dist %>% filter(dbh_threshold == "DBH > 10 cm")
cwm_dbh_gt_10_ks14 <- cwm_dbh_gt_10 %>% filter(scenario == "deep_clayey_regclim" | scenario == "shallow_sandy_regclim")

# select dbh > 1
cwm_dbh_gt_1 <- df_trait_dist %>% filter(dbh_threshold == "DBH > 1 cm")
cwm_dbh_gt_1_ks14 <- cwm_dbh_gt_1 %>% filter(scenario == "deep_clayey_regclim" | scenario == "shallow_sandy_regclim")

scenario_colors <- c(
  "deep_clayey_regclim" = "#2F4858",     # dark blue/grey
  "shallow_sandy_regclim" = "#A67C52"    # muted brown/sandy
)

scenario_labels <- c(
  "deep_clayey_regclim" = "Deep WT | Clayey soil",
  "shallow_sandy_regclim" = "Shallow WT | Sandy soil"
)

plot_trait_distribution <- function(data,
                                    trait_name,
                                    dbh_threshold = 10,
                                    plot_type = c("density", "histogram"),
                                    binwidth = NULL,
                                    save_plot = FALSE,
                                    output_dir = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke/trait_plots",
                                    scenario_keep = c(
                                      "deep_clayey_regclim",
                                      "shallow_sandy_regclim"
                                    ),
                                    colors = scenario_colors,
                                    labels = scenario_labels,
                                    width = 8,
                                    height = 5,
                                    dpi = 300) {
  
  plot_type <- match.arg(plot_type)
  
  # transforma 10 em "DBH > 10 cm" e 1 em "DBH > 1 cm"
  dbh_label <- if (is.numeric(dbh_threshold)) {
    paste0("DBH > ", dbh_threshold, " cm")
  } else {
    dbh_threshold
  }
  
  # nomes bonitos para título e eixo x
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
  
  trait_xlabs <- c(
    "wsg" = "Wood density (g cm-3)",
    "LMA" = "LMA",
    "Nmass" = "Nmass",
    "Pmass" = "Pmass",
    "tlp" = "TLP",
    "dbhmax" = "DBHmax",
    "hmax" = "Hmax",
    "leafarea" = "Leaf area",
    "g1" = "g1"
  )
  
  trait_title <- if (trait_name %in% names(trait_titles)) {
    unname(trait_titles[trait_name])
  } else {
    trait_name
  }
  
  x_label <- if (trait_name %in% names(trait_xlabs)) {
    unname(trait_xlabs[trait_name])
  } else {
    trait_name
  }
  
  df_plot <- data %>%
    filter(
      trait == trait_name,
      dbh_threshold == dbh_label,
      scenario %in% scenario_keep
    )
  
  if (nrow(df_plot) == 0) {
    stop("No data after filtering. Check trait_name, dbh_threshold, and scenario names.")
  }
  
  if (plot_type == "density") {
    
    p <- df_plot %>%
      ggplot(aes(
        x = trait_value,
        fill = scenario,
        color = scenario,
        weight = weight
      )) +
      geom_density(alpha = 0.25, linewidth = 0.8) +
      facet_grid(weighting ~ ., scales = "free") +
      theme_bw() +
      scale_color_manual(
        values = scenario_colors,
        labels = scenario_labels,
        name = "WT | Pedology"
      ) +
      scale_fill_manual(
        values = scenario_colors,
        labels = scenario_labels,
        name = "WT | Pedology"
      ) +
      labs(
        title = paste0("Final ", trait_title, " distribution - ", dbh_label),
        x = x_label,
        y = "Weighted density"
      )
    
  } else if (plot_type == "histogram") {
    
    if (is.null(binwidth)) {
      range_trait <- range(df_plot$trait_value, na.rm = TRUE)
      binwidth <- diff(range_trait) / 30
      
      if (!is.finite(binwidth) || binwidth == 0) {
        binwidth <- 0.01
      }
    }
    
    df_bins <- df_plot %>%
      mutate(
        bin = floor(trait_value / binwidth) * binwidth,
        bin_mid = bin + binwidth / 2
      ) %>%
      group_by(weighting, scenario, bin_mid) %>%
      summarise(
        bin_weight = sum(weight, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      group_by(weighting, scenario) %>%
      mutate(
        prop_weight = bin_weight / sum(bin_weight, na.rm = TRUE)
      ) %>%
      ungroup()
    
    p <- ggplot(df_bins, aes(
      x = bin_mid,
      y = prop_weight,
      fill = scenario
    )) +
      geom_col(
        position = "identity",
        alpha = 0.55,
        width = binwidth
      ) +
      facet_grid(weighting ~ .) +
      theme_bw() +
      scale_fill_manual(
        values = scenario_colors,
        labels = scenario_labels,
        name = "WT | Pedology"
      ) +
      labs(
        title = paste0("Final ", trait_title, " distribution - ", dbh_label),
        x = x_label,
        y = "Weighted proportion"
      )
  }
  
  if (save_plot) {
    
    if (is.null(output_dir)) {
      stop("You need to provide output_dir when save_plot = TRUE.")
    }
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    dbh_tag <- gsub("[^A-Za-z0-9]+", "_", dbh_label)
    
    filename <- file.path(
      output_dir,
      paste0(plot_type, "_", trait_name, "_", dbh_tag, ".png")
    )
    
    ggsave(
      filename = filename,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  return(p)
}

wsg_density_gt10 <- plot_trait_distribution(
  data = df_trait_dist,
  trait_name = "wsg",
  dbh_threshold = 10,
  plot_type = "density",
  save_plot = FALSE
)

wsg_density_gt10

wsg_density_gt1 <- plot_trait_distribution(
  data = df_trait_dist,
  trait_name = "wsg",
  dbh_threshold = 1,
  plot_type = "density",
  save_plot = FALSE
)

wsg_density_gt1

# wsg_hist_gt10 <- plot_trait_distribution(
#   data = df_trait_dist,
#   trait_name = "wsg",
#   dbh_threshold = 10,
#   plot_type = "histogram",
#   save_plot = TRUE
# )
# 
# wsg_hist_gt10

# wsg_hist_gt1 <- plot_trait_distribution(
#   data = df_trait_dist,
#   trait_name = "wsg",
#   dbh_threshold = 1,
#   plot_type = "histogram",
#   save_plot = TRUE
# )
# 
# wsg_hist_gt1

# 
lma_density_gt10 <- plot_trait_distribution(
  data = df_trait_dist,
  trait_name = "LMA",
  dbh_threshold = 10,
  plot_type = "density",
  save_plot = FALSE
)

lma_density_gt10
# 
# lma_density_gt1 <- plot_trait_distribution(
#   data = df_trait_dist,
#   trait_name = "LMA",
#   dbh_threshold = 1,
#   plot_type = "density",
#   save_plot = TRUE
# )
# 
# lma_density_gt1

# lma_hist_gt10 <- plot_trait_distribution(
#   data = df_trait_dist,
#   trait_name = "LMA",
#   dbh_threshold = 10,
#   plot_type = "histogram",
#   save_plot = TRUE
# )
# 
# lma_hist_gt10

# lma_hist_gt1 <- plot_trait_distribution(
#   data = df_trait_dist,
#   trait_name = "LMA",
#   dbh_threshold = 1,
#   plot_type = "histogram",
#   save_plot = TRUE
# )
# lma_hist_gt1

# Calculating the quantiles
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

# Boxplot
plot_trait_weighted_boxplot <- function(data,
                                        trait_name,
                                        dbh_threshold = 10,
                                        save_plot = FALSE,
                                        output_dir = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke/trait_plots",
                                        scenario_keep = c(
                                          "deep_clayey_regclim",
                                          "shallow_sandy_regclim"
                                        ),
                                        colors = NULL,
                                        labels = NULL,
                                        width = 8,
                                        height = 5,
                                        dpi = 300) {
  
  # Convert 10 into "DBH > 10 cm" and 1 into "DBH > 1 cm"
  dbh_label <- if (is.numeric(dbh_threshold)) {
    paste0("DBH > ", dbh_threshold, " cm")
  } else {
    dbh_threshold
  }
  
  # Default colors if no color vector is provided
  if (is.null(colors)) {
    colors <- c(
      "deep_clayey_regclim" = "#2F4858",
      "shallow_sandy_regclim" = "#A67C52"
    )
  }
  
  # Default labels if no label vector is provided
  if (is.null(labels)) {
    labels <- c(
      "deep_clayey_regclim" = "Deep WT | Clayey soil",
      "shallow_sandy_regclim" = "Shallow WT | Sandy soil"
    )
  }
  
  # Define nicer trait names for plot titles
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
  
  # Define y-axis labels for each trait
  trait_ylabs <- c(
    "wsg" = "Wood density (g cm-3)",
    "LMA" = "LMA",
    "Nmass" = "Nmass",
    "Pmass" = "Pmass",
    "tlp" = "TLP",
    "dbhmax" = "DBHmax",
    "hmax" = "Hmax",
    "leafarea" = "Leaf area",
    "g1" = "g1"
  )
  
  # Use formatted trait name if available
  trait_title <- if (trait_name %in% names(trait_titles)) {
    unname(trait_titles[trait_name])
  } else {
    trait_name
  }
  
  # Use formatted y-axis label if available
  y_label <- if (trait_name %in% names(trait_ylabs)) {
    unname(trait_ylabs[trait_name])
  } else {
    trait_name
  }
  
  # Weighted quantile function
  weighted_quantile <- function(x, w, prob) {
    keep <- is.finite(x) & is.finite(w) & w > 0
    x <- x[keep]
    w <- w[keep]
    
    if (length(x) == 0) {
      return(NA_real_)
    }
    
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    
    cw <- cumsum(w) / sum(w)
    x[which(cw >= prob)[1]]
  }
  
  # Filter data
  df_plot <- data %>%
    filter(
      trait == trait_name,
      dbh_threshold == dbh_label,
      scenario %in% scenario_keep
    )
  
  # Stop if there is no data after filtering
  if (nrow(df_plot) == 0) {
    stop("No data after filtering. Check trait_name, dbh_threshold, and scenario names.")
  }
  
  # Calculate weighted boxplot statistics
  df_boxplot <- df_plot %>%
    group_by(weighting, scenario) %>%
    summarise(
      ymin = weighted_quantile(trait_value, weight, 0.10),
      lower = weighted_quantile(trait_value, weight, 0.25),
      middle = weighted_quantile(trait_value, weight, 0.50),
      upper = weighted_quantile(trait_value, weight, 0.75),
      ymax = weighted_quantile(trait_value, weight, 0.90),
      weighted_mean = weighted.mean(trait_value, weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create weighted boxplot
  p <- ggplot(df_boxplot, aes(
    x = scenario,
    fill = scenario,
    color = scenario
  )) +
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
      alpha = 0.6
    ) +
    geom_point(
      aes(y = weighted_mean),
      shape = 21,
      size = 2.5,
      fill = "white"
    ) +
    facet_wrap(~ weighting) +
    theme_bw() +
    scale_fill_manual(
      values = colors,
      labels = labels,
      name = "WT | Pedology"
    ) +
    scale_color_manual(
      values = colors,
      labels = labels,
      name = "WT | Pedology"
    ) +
    scale_x_discrete(
      labels = labels
    ) +
    labs(
      title = paste0("Weighted ", trait_title, " distribution - ", dbh_label),
      x = NULL,
      y = y_label
    )
  
  # Save plot if requested
  if (save_plot) {
    
    if (is.null(output_dir)) {
      stop("You need to provide output_dir when save_plot = TRUE.")
    }
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    dbh_tag <- gsub("[^A-Za-z0-9]+", "_", dbh_label)
    
    filename <- file.path(
      output_dir,
      paste0("boxplot_", trait_name, "_", dbh_tag, ".png")
    )
    
    ggsave(
      filename = filename,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  return(p)
}

# wsg_boxplot_gt10 <- plot_trait_weighted_boxplot(
#   data = df_trait_dist,
#   trait_name = "wsg",
#   dbh_threshold = 10,
#   save_plot = FALSE
# )
# 
# wsg_boxplot_gt10

  # wsg_boxplot_gt1 <- plot_trait_weighted_boxplot(
#   data = df_trait_dist,
#   trait_name = "wsg",
#   dbh_threshold = 1,
#   save_plot = TRUE
# )
# 
# wsg_boxplot_gt1

# LMA_boxplot_gt10 <- plot_trait_weighted_boxplot(
#   data = df_trait_dist,
#   trait_name = "LMA",
#   dbh_threshold = 10,
#   save_plot = FALSE
# )
# # 
# LMA_boxplot_gt10
# 
# LMA_boxplot_gt1 <- plot_trait_weighted_boxplot(
#   data = df_trait_dist,
#   trait_name = "LMA",
#   dbh_threshold = 1,
#   save_plot = TRUE
# )
# 
# LMA_boxplot_gt1


# wsg_summary_gt10 <- df_trait_dist %>%
#   filter(
#     trait == "wsg",
#     dbh_threshold == "DBH > 10 cm",
#     scenario %in% c("deep_clayey_regclim", "shallow_sandy_regclim")
#   ) %>%
#   group_by(weighting, scenario) %>%
#   summarise(
#     weighted_mean = weighted.mean(trait_value, weight, na.rm = TRUE),
#     q10 = weighted_quantile(trait_value, weight, 0.10),
#     median = weighted_quantile(trait_value, weight, 0.50),
#     q90 = weighted_quantile(trait_value, weight, 0.90),
#     prop_wsg_below_060 = sum(weight[trait_value < 0.60], na.rm = TRUE) / sum(weight, na.rm = TRUE),
#     prop_wsg_below_070 = sum(weight[trait_value < 0.70], na.rm = TRUE) / sum(weight, na.rm = TRUE),
#     .groups = "drop"
#   )
# wsg_summary_gt10
# write_csv(wsg_summary_gt10,"~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke/trait_plots/wsg_summary_gt10.csv")
# 
# LMA_summary_gt10 <- df_trait_dist %>%
#   filter(
#     trait == "LMA",
#     dbh_threshold == "DBH > 10 cm",
#     scenario %in% c("deep_clayey_regclim", "shallow_sandy_regclim")
#   ) %>%
#   group_by(weighting, scenario) %>%
#   summarise(
#     weighted_mean = weighted.mean(trait_value, weight, na.rm = TRUE),
#     q10 = weighted_quantile(trait_value, weight, 0.10),
#     median = weighted_quantile(trait_value, weight, 0.50),
#     q90 = weighted_quantile(trait_value, weight, 0.90),
#     prop_LMA_below_060 = sum(weight[trait_value < 0.60], na.rm = TRUE) / sum(weight, na.rm = TRUE),
#     prop_LMA_below_070 = sum(weight[trait_value < 0.70], na.rm = TRUE) / sum(weight, na.rm = TRUE),
#     .groups = "drop"
#   )
# LMA_summary_gt10
# write_csv(LMA_summary_gt10,"~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke/trait_plots/LMA_summary_gt10.csv")

for (trait in trait_vars){
  plot_trait_distribution(
    data = df_trait_dist,
    trait_name = trait,
    dbh_threshold = 10,
    plot_type = "density",
    save_plot = TRUE
  )
}

for (trait in trait_vars){
  plot_trait_distribution(
    data = df_trait_dist,
    trait_name = trait,
    dbh_threshold = 10,
    plot_type = "histogram",
    save_plot = TRUE
  )
}

for (trait in trait_vars){
  plot_trait_weighted_boxplot(
      data = df_trait_dist,
      trait_name = trait,
      dbh_threshold = 10,
      save_plot = TRUE
    )
}
