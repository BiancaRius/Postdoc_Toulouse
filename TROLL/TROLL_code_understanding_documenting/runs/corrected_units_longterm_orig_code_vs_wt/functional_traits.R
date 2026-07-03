## TRAITS DISTRIBUTION ##

library(dplyr)
library(ggplot2)
library(tidyverse)

base_dir <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/corrected_units_longterm_orig_code_vs_wt/"

scenario_table <- tibble(
  scenario = c(
    "longterm_origcode_sandy_regclim",
    "longterm_deep_sandy_regclim",
    "longterm_shallow_sandy_regclim"
  ),
  run_dir = c(
    "longterm_origcode_sandy_regclim",
    "longterm_deep_sandy_regclim",
    "longterm_shallow_sandy_regclim"  )
) %>%
  mutate(
    path = file.path(base_dir, run_dir),
    final_pattern_file = file.path(path, "output", "(null)_0_final_pattern.txt")
  )


safe_read <- possibly(
  ~ read.table(.x, header = TRUE),
  otherwise = NULL
)

df_all_raw <- scenario_table %>%
  mutate(data = map(final_pattern_file, safe_read)) %>%
  filter(!map_lgl(data, is.null)) %>%
  select(scenario, data) %>%
  unnest(data)

df_all <- df_all_raw %>%
  mutate(
    basal_area        = pi * (dbh / 2)^2,
    abundance_weight  = 1,
    basal_area_weight = basal_area,
    biomass_weight    = AGB
  )

df_threshold <- bind_rows(
  df_all %>%
    filter(dbh > 0.01) %>%
    mutate(dbh_threshold = "DBH > 1 cm"),
  df_all %>%
    filter(dbh > 0.10) %>%
    mutate(dbh_threshold = "DBH > 10 cm")
)

trait_vars <- c("LMA", "wsg", "Nmass", "Pmass", "tlp", "dbhmax", "hmax", "leafarea", "g1")
trait_vars <- trait_vars[trait_vars %in% names(df_threshold)]

cwm_traits <- df_threshold %>%
  mutate(tree_id = row_number()) %>%
  pivot_longer(
    cols      = all_of(trait_vars),
    names_to  = "trait",
    values_to = "trait_value"
  ) %>%
  group_by(scenario, dbh_threshold, trait) %>%
  summarise(
    n_trees        = n_distinct(tree_id[!is.na(trait_value)]),
    CWM_abundance  = weighted.mean(trait_value, abundance_weight),
    CWM_basal_area = weighted.mean(trait_value, basal_area_weight),
    CWM_biomass    = weighted.mean(trait_value, biomass_weight),
    .groups = "drop"
  )

cwm_traits_long <- cwm_traits %>%
  pivot_longer(
    cols      = starts_with("CWM_"),
    names_to  = "weighting",
    values_to = "CWM"
  ) %>%
  mutate(
    weighting = recode(
      weighting,
      "CWM_abundance"  = "Abundance",
      "CWM_basal_area" = "Basal area",
      "CWM_biomass"    = "Biomass"
    )
  )

df_trait_dist <- df_threshold %>%
  select(
    scenario, dbh_threshold,
    all_of(trait_vars),
    abundance_weight, basal_area_weight, biomass_weight
  ) %>%
  pivot_longer(
    cols      = all_of(trait_vars),
    names_to  = "trait",
    values_to = "trait_value"
  ) %>%
  pivot_longer(
    cols      = c(abundance_weight, basal_area_weight, biomass_weight),
    names_to  = "weighting",
    values_to = "weight"
  ) %>%
  mutate(
    weighting = recode(
      weighting,
      "abundance_weight"  = "Abundance",
      "basal_area_weight" = "Basal area",
      "biomass_weight"    = "Biomass"
    )
  ) %>%
  filter(!is.na(trait_value), !is.na(weight), weight > 0)

scenario_colors <- c(
  "longterm_origcode_sandy_regclim" = "black",
  "longterm_deep_sandy_regclim"     = "#1B6CA8",
  "longterm_shallow_sandy_regclim" = "#CC5500"
  
)

scenario_labels <- c(
  "longterm_origcode_sandy_regclim" = "Original code",# | Sandy soil | Regclim",
  "longterm_deep_sandy_regclim"     = "Deep WT", # | Sandy soil | Regclim",
  "longterm_shallow_sandy_regclim" = "Shallow WT" #| Sandy soil | Regclim"
)

# -----------------------------------------------------------------------------
# FUNCTION: plot_trait_distribution
# Plots weighted density or histogram of a trait, with dashed vertical line
# at the weighted mode per scenario and faceted by weighting scheme.
# MODIFICATION 2: dashed vertical line at weighted mode, shown in legend
# -----------------------------------------------------------------------------
plot_trait_distribution <- function(data,
                                    trait_name,
                                    dbh_threshold = 10,
                                    plot_type     = c("density", "histogram"),
                                    binwidth      = NULL,
                                    save_plot     = FALSE,
                                    output_dir    = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/corrected_units_longterm_orig_code_vs_wt/plot_traits",
                                    scenario_keep = c(
                                      "longterm_origcode_sandy_regclim",
                                      "longterm_deep_sandy_regclim",
                                      "longterm_shallow_sandy_regclim"
                                    ),
                                    colors = scenario_colors,
                                    labels = scenario_labels,
                                    width  = 8,
                                    height = 5,
                                    dpi    = 300) {
  
  plot_type <- match.arg(plot_type)
  
  # Convert numeric threshold to label string
  dbh_label <- if (is.numeric(dbh_threshold)) {
    paste0("DBH > ", dbh_threshold, " cm")
  } else {
    dbh_threshold
  }
  
  # Formatted trait names for title and x-axis
  trait_titles <- c(
    "wsg"      = "WD",       "LMA"      = "LMA",
    "Nmass"    = "Nmass",    "Pmass"    = "Pmass",
    "tlp"      = "TLP",      "dbhmax"   = "DBHmax",
    "hmax"     = "Hmax",     "leafarea" = "Leaf area",
    "g1"       = "g1"
  )
  
  trait_xlabs <- c(
    "wsg"      = "Wood density (g cm-3)", "LMA"      = "LMA",
    "Nmass"    = "Nmass",                 "Pmass"    = "Pmass",
    "tlp"      = "TLP",                   "dbhmax"   = "DBHmax",
    "hmax"     = "Hmax",                  "leafarea" = "Leaf area",
    "g1"       = "g1"
  )
  
  trait_title <- if (trait_name %in% names(trait_titles)) unname(trait_titles[trait_name]) else trait_name
  x_label     <- if (trait_name %in% names(trait_xlabs))  unname(trait_xlabs[trait_name])  else trait_name
  
  # Filter to the requested trait, DBH threshold, and scenarios
  df_plot <- data %>%
    filter(
      trait         == trait_name,
      dbh_threshold == dbh_label,
      scenario      %in% scenario_keep
    )
  
  if (nrow(df_plot) == 0) {
    stop("No data after filtering. Check trait_name, dbh_threshold, and scenario names.")
  }
  
  if (plot_type == "density") {
    
    # --- MODIFICATION 2: compute weighted mode (density peak) per group ---
    df_mode <- df_plot %>%
      group_by(scenario, weighting) %>%
      summarise(
        mode_val = {
          
          # Keep only valid trait values and valid weights
          keep <- is.finite(trait_value) & is.finite(weight) & weight > 0
          x <- trait_value[keep]
          w <- weight[keep]
          
          # If there are too few values, return NA
          if (length(x) < 2 || sum(w, na.rm = TRUE) <= 0) {
            NA_real_
          } else {
            
            # Normalize weights so they sum to 1
            w <- w / sum(w, na.rm = TRUE)
            
            # Compute a numeric bandwidth.
            # This avoids the warning caused by bw = "nrd0".
            bw_value <- bw.nrd0(x)
            
            # If the bandwidth is invalid, use a small fallback value
            if (!is.finite(bw_value) || bw_value <= 0) {
              bw_value <- diff(range(x, na.rm = TRUE)) / 30
            }
            
            if (!is.finite(bw_value) || bw_value <= 0) {
              bw_value <- 0.01
            }
            
            # Estimate weighted density using a numeric bandwidth
            d <- density(
              x,
              weights = w,
              bw = bw_value,
              na.rm = TRUE
            )
            
            # Return the x value where the weighted density is highest
            d$x[which.max(d$y)]
          }
        },
        .groups = "drop"
      )
    
    p <- df_plot %>%
      ggplot(aes(
        x      = trait_value,
        fill   = scenario,
        color  = scenario,
        weight = weight
      )) +
      geom_density(alpha = 0.25, linewidth = 0.8) +
      # Dashed vertical line at the weighted mode for each scenario
      geom_vline(
        data        = df_mode,
        aes(xintercept = mode_val, color = scenario),
        linetype    = "dashed",
        linewidth   = 0.6
      ) +
      facet_grid(weighting ~ ., scales = "free") +
      theme_bw() +
      scale_color_manual(
        values = colors,
        labels = labels,
        name   = "WT | Pedology"
      ) +
      scale_fill_manual(
        values = colors,
        labels = labels,
        name   = "WT | Pedology"
      ) +
      # Show dashed line style in the legend key
      guides(
        color = guide_legend(
          override.aes = list(
            linetype  = "dashed",
            linewidth = 0.8,
            fill      = NA
          )
        ),
        fill = guide_legend(
          override.aes = list(alpha = 0.25)
        )
      ) +
      labs(
        title    = paste0("Final ", trait_title, " distribution - ", dbh_label),
        x        = x_label,
        y        = "Weighted density",
        # Explain dashed line meaning in the caption
        caption  = "Dashed vertical line = weighted mode"
      )
    
  } else if (plot_type == "histogram") {
    
    if (is.null(binwidth)) {
      range_trait <- range(df_plot$trait_value, na.rm = TRUE)
      binwidth    <- diff(range_trait) / 30
      if (!is.finite(binwidth) || binwidth == 0) binwidth <- 0.01
    }
    
    df_bins <- df_plot %>%
      mutate(
        bin     = floor(trait_value / binwidth) * binwidth,
        bin_mid = bin + binwidth / 2
      ) %>%
      group_by(weighting, scenario, bin_mid) %>%
      summarise(bin_weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
      group_by(weighting, scenario) %>%
      mutate(prop_weight = bin_weight / sum(bin_weight, na.rm = TRUE)) %>%
      ungroup()
    
    # --- MODIFICATION 2: compute weighted mode per group for histogram ---
    df_mode <- df_plot %>%
      group_by(scenario, weighting) %>%
      summarise(
        mode_val = {
          
          # Keep only valid trait values and valid weights
          keep <- is.finite(trait_value) & is.finite(weight) & weight > 0
          x <- trait_value[keep]
          w <- weight[keep]
          
          # If there are too few values, return NA
          if (length(x) < 2 || sum(w, na.rm = TRUE) <= 0) {
            NA_real_
          } else {
            
            # Normalize weights so they sum to 1
            w <- w / sum(w, na.rm = TRUE)
            
            # Compute a numeric bandwidth.
            # This avoids the warning caused by bw = "nrd0".
            bw_value <- bw.nrd0(x)
            
            # If the bandwidth is invalid, use a small fallback value
            if (!is.finite(bw_value) || bw_value <= 0) {
              bw_value <- diff(range(x, na.rm = TRUE)) / 30
            }
            
            if (!is.finite(bw_value) || bw_value <= 0) {
              bw_value <- 0.01
            }
            
            # Estimate weighted density using a numeric bandwidth
            d <- density(
              x,
              weights = w,
              bw = bw_value,
              na.rm = TRUE
            )
            
            # Return the x value where the weighted density is highest
            d$x[which.max(d$y)]
          }
        },
        .groups = "drop"
      )
    
    p <- ggplot(df_bins, aes(x = bin_mid, y = prop_weight, fill = scenario)) +
      geom_col(position = "identity", alpha = 0.55, width = binwidth) +
      # Dashed vertical line at the weighted mode
      geom_vline(
        data      = df_mode,
        aes(xintercept = mode_val, color = scenario),
        linetype  = "dashed",
        linewidth = 0.6
      ) +
      facet_grid(weighting ~ .) +
      theme_bw() +
      scale_fill_manual(
        values = colors,
        labels = labels,
        name   = "WT | Pedology"
      ) +
      scale_color_manual(
        values = colors,
        labels = labels,
        name   = "WT | Pedology"
      ) +
      labs(
        title   = paste0("Final ", trait_title, " distribution - ", dbh_label),
        x       = x_label,
        y       = "Weighted proportion",
        caption = "Dashed vertical line = weighted mode"
      )
  }
  
  if (save_plot) {
    if (is.null(output_dir)) stop("You need to provide output_dir when save_plot = TRUE.")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    dbh_tag  <- gsub("[^A-Za-z0-9]+", "_", dbh_label)
    filename <- file.path(output_dir, paste0(plot_type, "_", trait_name, "_", dbh_tag, ".png"))
    ggsave(filename = filename, plot = p, width = width, height = height, dpi = dpi)
  }
  
  return(p)
}


# -----------------------------------------------------------------------------
# HELPER: weighted_quantile
# Returns the weighted quantile at probability `prob`
# -----------------------------------------------------------------------------
weighted_quantile <- function(x, w, prob) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x    <- x[keep]
  w    <- w[keep]
  if (length(x) == 0) return(NA_real_)
  ord  <- order(x)
  x    <- x[ord]
  w    <- w[ord]
  cw   <- cumsum(w) / sum(w)
  x[which(cw >= prob)[1]]
}


# -----------------------------------------------------------------------------
# FUNCTION: plot_trait_weighted_boxplot
# Plots weighted boxplots (10th–90th percentile) with weighted mean as a
# white dot, faceted by weighting scheme.
# MODIFICATION 1: x-axis labels rotated 45 degrees
# MODIFICATION 3: white dot (weighted mean) explained in legend caption
# -----------------------------------------------------------------------------
plot_trait_weighted_boxplot <- function(data,
                                        trait_name,
                                        dbh_threshold = 10,
                                        save_plot     = FALSE,
                                        output_dir    = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/corrected_units_longterm_orig_code_vs_wt/plot_traits",
                                        scenario_keep = c(
                                          "longterm_origcode_sandy_regclim",
                                          "longterm_deep_sandy_regclim",
                                          "longterm_shallow_sandy_regclim"
                                        ),
                                        colors = NULL,
                                        labels = NULL,
                                        width  = 8,
                                        height = 5,
                                        dpi    = 300) {
  
  # Convert numeric threshold to label string
  dbh_label <- if (is.numeric(dbh_threshold)) {
    paste0("DBH > ", dbh_threshold, " cm")
  } else {
    dbh_threshold
  }
  
  # Default color palette
  if (is.null(colors)) {
    colors <- c(
      "longterm_origcode_sandy_regclim" = "black",
      "longterm_deep_sandy_regclim"     = "#1B6CA8",
      "longterm_shallow_sandy_regclim" = "#CC5500"
      
      
    )
  }
  
  # Default scenario labels
  if (is.null(labels)) {
    labels <- c(
      "longterm_origcode_sandy_regclim" = "Original code",# | Sandy soil | Regclim",
      "longterm_deep_sandy_regclim"     = "Deep WT", # | Sandy soil | Regclim",
      "longterm_shallow_sandy_regclim" = "Shallow WT" #| Sandy soil | Regclim"
    )
  }
  
  # Formatted trait names for title and y-axis
  trait_titles <- c(
    "wsg"      = "WD",       "LMA"      = "LMA",
    "Nmass"    = "Nmass",    "Pmass"    = "Pmass",
    "tlp"      = "TLP",      "dbhmax"   = "DBHmax",
    "hmax"     = "Hmax",     "leafarea" = "Leaf area",
    "g1"       = "g1"
  )
  
  trait_ylabs <- c(
    "wsg"      = "Wood density (g cm-3)", "LMA"      = "LMA",
    "Nmass"    = "Nmass",                 "Pmass"    = "Pmass",
    "tlp"      = "TLP",                   "dbhmax"   = "DBHmax",
    "hmax"     = "Hmax",                  "leafarea" = "Leaf area",
    "g1"       = "g1"
  )
  
  trait_title <- if (trait_name %in% names(trait_titles)) unname(trait_titles[trait_name]) else trait_name
  y_label     <- if (trait_name %in% names(trait_ylabs))  unname(trait_ylabs[trait_name])  else trait_name
  
  # Weighted quantile (defined locally to keep function self-contained)
  weighted_quantile <- function(x, w, prob) {
    keep <- is.finite(x) & is.finite(w) & w > 0
    x    <- x[keep]
    w    <- w[keep]
    if (length(x) == 0) return(NA_real_)
    ord  <- order(x)
    x    <- x[ord]
    w    <- w[ord]
    cw   <- cumsum(w) / sum(w)
    x[which(cw >= prob)[1]]
  }
  
  # Filter to the requested trait, DBH threshold, and scenarios
  df_plot <- data %>%
    filter(
      trait         == trait_name,
      dbh_threshold == dbh_label,
      scenario      %in% scenario_keep
    )
  
  if (nrow(df_plot) == 0) {
    stop("No data after filtering. Check trait_name, dbh_threshold, and scenario names.")
  }
  
  # Compute weighted boxplot statistics (10th–90th percentile range)
  df_boxplot <- df_plot %>%
    group_by(weighting, scenario) %>%
    summarise(
      ymin          = weighted_quantile(trait_value, weight, 0.10),
      lower         = weighted_quantile(trait_value, weight, 0.25),
      middle        = weighted_quantile(trait_value, weight, 0.50),
      upper         = weighted_quantile(trait_value, weight, 0.75),
      ymax          = weighted_quantile(trait_value, weight, 0.90),
      weighted_mean = weighted.mean(trait_value, weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(df_boxplot, aes(x = scenario, fill = scenario, color = scenario)) +
    geom_boxplot(
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat  = "identity",
      width = 0.55,
      alpha = 0.6
    ) +
    # White dot = weighted mean
    geom_point(
      aes(y = weighted_mean),
      shape = 21,
      size  = 2.5,
      fill  = "white",
      color = "black"
    ) +
    facet_wrap(~ weighting) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(
      values = colors,
      labels = labels,
      name   = "WT | Pedology"
    ) +
    scale_color_manual(
      values = colors,
      labels = labels,
      name   = "WT | Pedology"
    ) +
    scale_x_discrete(labels = labels) +
    # --- MODIFICATION 3: explain white dot in caption ---
    labs(
      title   = paste0("Weighted ", trait_title, " distribution - ", dbh_label),
      x       = NULL,
      y       = y_label,
      caption = "Boxes: 25th–75th percentile | Whiskers: 10th–90th percentile | ● Weighted mean"
    )
  
  if (save_plot) {
    if (is.null(output_dir)) stop("You need to provide output_dir when save_plot = TRUE.")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    dbh_tag  <- gsub("[^A-Za-z0-9]+", "_", dbh_label)
    filename <- file.path(output_dir, paste0("boxplot_", trait_name, "_", dbh_tag, ".png"))
    ggsave(filename = filename, plot = p, width = width, height = height, dpi = dpi)
  }
  
  return(p)
}


# # --- Example calls ---
# 
# wsg_density_gt10 <- plot_trait_distribution(
#   data          = df_trait_dist,
#   trait_name    = "wsg",
#   dbh_threshold = 10,
#   plot_type     = "density",
#   save_plot     = FALSE
# )
# wsg_density_gt10
# 
# wsg_boxplot_gt10 <- plot_trait_weighted_boxplot(
#   data          = df_trait_dist,
#   trait_name    = "wsg",
#   dbh_threshold = 10,
#   save_plot     = FALSE
# )
# wsg_boxplot_gt10
# 
# 
# # --- Batch export all traits ---
# 
# for (trait in trait_vars) {
#   plot_trait_distribution(
#     data          = df_trait_dist,
#     trait_name    = trait,
#     dbh_threshold = 10,
#     plot_type     = "density",
#     save_plot     = TRUE
#   )
# }
# 
# for (trait in trait_vars) {
#   plot_trait_distribution(
#     data          = df_trait_dist,
#     trait_name    = trait,
#     dbh_threshold = 10,
#     plot_type     = "histogram",
#     save_plot     = TRUE
#   )
# }
# 
# for (trait in trait_vars) {
#   plot_trait_weighted_boxplot(
#     data          = df_trait_dist,
#     trait_name    = trait,
#     dbh_threshold = 10,
#     save_plot     = TRUE
#   )
# }

# -----------------------------------------------------------------------------
# HELPER: weighted mode based on weighted density
# -----------------------------------------------------------------------------
weighted_mode_density <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  
  if (length(x) < 2 || sum(w, na.rm = TRUE) <= 0) {
    return(NA_real_)
  }
  
  w <- w / sum(w, na.rm = TRUE)
  
  bw_value <- bw.nrd0(x)
  
  if (!is.finite(bw_value) || bw_value <= 0) {
    bw_value <- diff(range(x, na.rm = TRUE)) / 30
  }
  
  if (!is.finite(bw_value) || bw_value <= 0) {
    bw_value <- 0.01
  }
  
  d <- density(
    x,
    weights = w,
    bw = bw_value,
    na.rm = TRUE
  )
  
  d$x[which.max(d$y)]
}


# -----------------------------------------------------------------------------
# FUNCTION: plot all trait weighted densities in one figure
# -----------------------------------------------------------------------------
plot_all_traits_distribution <- function(data,
                                         dbh_threshold = 10,
                                         save_plot = FALSE,
                                         output_dir = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/corrected_units_longterm_orig_code_vs_wt/plot_traits",
                                         scenario_keep = c(
                                           "longterm_origcode_sandy_regclim",
                                           "longterm_deep_sandy_regclim",
                                           "longterm_shallow_sandy_regclim"
                                         ),
                                         colors = scenario_colors,
                                         labels = scenario_labels,
                                         width = 18,
                                         height = 8,
                                         dpi = 300) {
  
  dbh_label <- if (is.numeric(dbh_threshold)) {
    paste0("DBH > ", dbh_threshold, " cm")
  } else {
    dbh_threshold
  }
  
  trait_labels <- c(
    "wsg"      = "WD",
    "LMA"      = "LMA",
    "Nmass"    = "Nmass",
    "Pmass"    = "Pmass",
    "tlp"      = "TLP",
    "dbhmax"   = "DBHmax",
    "hmax"     = "Hmax",
    "leafarea" = "Leaf area",
    "g1"       = "g1"
  )
  
  df_plot <- data %>%
    filter(
      dbh_threshold == dbh_label,
      scenario %in% scenario_keep
    ) %>%
    mutate(
      trait_label = recode(trait, !!!trait_labels),
      trait_label = factor(trait_label, levels = trait_labels[trait_vars])
    )
  
  if (nrow(df_plot) == 0) {
    stop("No data after filtering. Check dbh_threshold and scenario names.")
  }
  
  df_mode <- df_plot %>%
    group_by(trait, trait_label, scenario, weighting) %>%
    summarise(
      mode_val = weighted_mode_density(trait_value, weight),
      .groups = "drop"
    )
  
  p <- ggplot(
    df_plot,
    aes(
      x = trait_value,
      fill = scenario,
      color = scenario,
      weight = weight
    )
  ) +
    geom_density(alpha = 0.25, linewidth = 0.7, na.rm = TRUE) +
    geom_vline(
      data = df_mode,
      aes(xintercept = mode_val, color = scenario),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    facet_wrap(
      ~ weighting + trait_label,
      scales = "free",
      ncol = length(unique(df_plot$trait_label))
    ) +
    scale_color_manual(
      values = colors,
      labels = labels,
      name = "Scenario"
    ) +
    scale_fill_manual(
      values = colors,
      labels = labels,
      name = "Scenario"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.x = element_blank()
    ) +
    labs(
      title = paste0("Final trait distributions - ", dbh_label),
      y = "Weighted density",
      caption = "Dashed vertical line = weighted mode"
    )
  
  if (save_plot) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    dbh_tag <- gsub("[^A-Za-z0-9]+", "_", dbh_label)
    
    ggsave(
      filename = file.path(output_dir, paste0("density_all_traits_", dbh_tag, ".png")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  return(p)
}

all_traits_density_gt10 <- plot_all_traits_distribution(
  data = df_trait_dist,
  dbh_threshold = 10,
  save_plot = TRUE
)

all_traits_density_gt10

all_traits_density_gt1 <- plot_all_traits_distribution(
  data = df_trait_dist,
  dbh_threshold = 1,
  save_plot = TRUE
)

all_traits_density_gt1

 # -----------------------------------------------------------------------------
# FUNCTION: plot all trait weighted boxplots in one figure
# -----------------------------------------------------------------------------
plot_all_traits_weighted_boxplot <- function(data,
                                             dbh_threshold = 10,
                                             save_plot = FALSE,
                                             output_dir = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/corrected_units_longterm_orig_code_vs_wt/plot_traits",
                                             scenario_keep = c(
                                               "longterm_origcode_sandy_regclim",
                                               "longterm_deep_sandy_regclim",
                                               "longterm_shallow_sandy_regclim"
                                             ),
                                             colors = scenario_colors,
                                             labels = scenario_labels,
                                             width = 12,
                                             height = 14,
                                             dpi = 300) {
  
  dbh_label <- if (is.numeric(dbh_threshold)) {
    paste0("DBH > ", dbh_threshold, " cm")
  } else {
    dbh_threshold
  }
  
  trait_labels <- c(
    "wsg"      = "WD",
    "LMA"      = "LMA",
    "Nmass"    = "Nmass",
    "Pmass"    = "Pmass",
    "tlp"      = "TLP",
    "dbhmax"   = "DBHmax",
    "hmax"     = "Hmax",
    "leafarea" = "Leaf area",
    "g1"       = "g1"
  )
  
  df_plot <- data %>%
    filter(
      dbh_threshold == dbh_label,
      scenario %in% scenario_keep
    ) %>%
    mutate(
      trait_label = recode(trait, !!!trait_labels),
      trait_label = factor(trait_label, levels = trait_labels[trait_vars])
    )
  
  if (nrow(df_plot) == 0) {
    stop("No data after filtering. Check dbh_threshold and scenario names.")
  }
  
  df_boxplot <- df_plot %>%
    group_by(trait, trait_label, weighting, scenario) %>%
    summarise(
      ymin          = weighted_quantile(trait_value, weight, 0.10),
      lower         = weighted_quantile(trait_value, weight, 0.25),
      middle        = weighted_quantile(trait_value, weight, 0.50),
      upper         = weighted_quantile(trait_value, weight, 0.75),
      ymax          = weighted_quantile(trait_value, weight, 0.90),
      weighted_mean = weighted.mean(trait_value, weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(df_boxplot, aes(x = scenario, fill = scenario, color = scenario)) +
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
      size = 2,
      fill = "white",
      color = "black"
    ) +
    facet_grid(trait_label ~ weighting, scales = "free_y") +
    scale_fill_manual(
      values = colors,
      labels = labels,
      name = "Scenario"
    ) +
    scale_color_manual(
      values = colors,
      labels = labels,
      name = "Scenario"
    ) +
    scale_x_discrete(labels = labels) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      axis.text.y = element_text(size = 7)
    ) +
    labs(
      title = paste0("Weighted trait distributions - ", dbh_label),
      x = NULL,
      y = "Trait value",
      caption = "Boxes: 25th–75th percentile | Whiskers: 10th–90th percentile | white dot = weighted mean"
    )
  
  if (save_plot) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    dbh_tag <- gsub("[^A-Za-z0-9]+", "_", dbh_label)
    
    ggsave(
      filename = file.path(output_dir, paste0("boxplot_all_traits_", dbh_tag, ".png")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  
  return(p)
}


all_traits_boxplot_gt10 <- plot_all_traits_weighted_boxplot(
  data = df_trait_dist,
  dbh_threshold = 10,
  save_plot = TRUE
)

all_traits_boxplot_gt10