## TRAITS DISTRIBUTION ##

library(dplyr)
library(ggplot2)
library(tidyverse)

base_dir <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke"

scenario_table <- tibble(
  scenario = c(
    "ducke_deep_clayey_regclim_ks14",
    "ducke_deep_clayey_redprec30_ks14",
    "ducke_deep_sandy_redprec30_ks14",
    "ducke_deep_sandy_regclim_ks14",
    "ducke_shallow_sandy_regclim_ks14",
    "ducke_shallow_sandy_redprec30_ks14",
    "ducke_shallow_clayey_redprec30_ks14",
    "ducke_shallow_clayey_regclim_ks14"
  ),
  run_dir = c(
    "deepWT_clayeysoil/ducke_deep_clayey_regclim_ks14",
    "deepWT_clayeysoil/ducke_deep_clayey_redprec30_ks14",
    "deepWT_sandysoil/ducke_deep_sandy_redprec30_ks14",
    "deepWT_sandysoil/ducke_deep_sandy_regclim_ks14",
    "shallowWT_sandysoil/ducke_shallow_sandy_regclim_ks14",
    "shallowWT_sandysoil/ducke_shallow_sandy_redprec30_ks14",
    "shallowWT_clayeysoil/ducke_shallow_clayey_redprec30_ks14",
    "shallowWT_clayeysoil/ducke_shallow_clayey_regclim_ks14"
  )
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
  "ducke_deep_clayey_regclim_ks14"     = "#1B6CA8",
  "ducke_deep_clayey_redprec30_ks14"   = "#56B4E9",
  "ducke_deep_sandy_redprec30_ks14"   = "green4",
  "ducke_deep_sandy_regclim_ks14"   = "purple",
  "ducke_shallow_sandy_regclim_ks14"   = "#CC5500",
  "ducke_shallow_sandy_redprec30_ks14" = "#E8A268",
  "ducke_shallow_clayey_redprec30_ks14" = "chocolate4",
  "ducke_shallow_clayey_regclim_ks14" = "red"
  
)

scenario_labels <- c(
  "ducke_deep_clayey_regclim_ks14"     = "Deep WT | Clayey soil | Regclim",
  "ducke_deep_clayey_redprec30_ks14"   = "Deep WT | Clayey soil | Redprec",
  "ducke_deep_sandy_redprec30_ks14"   = "Deep WT | Sandy soil | Redprec",
  "ducke_deep_sandy_regclim_ks14"   = "Deep WT | Sandy soil | Regclim",
  "ducke_shallow_sandy_regclim_ks14"   = "Shallow WT | Sandy soil | Regclim",
  "ducke_shallow_sandy_redprec30_ks14" = "Shallow WT | Sandy soil | Redprec",
  "ducke_shallow_clayey_redprec30_ks14" = "Shallow WT | Clayey soil | Redprec",
  "ducke_shallow_clayey_regclim_ks14" = "Shallow WT | Clayey soil | Regclim"
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
                                    output_dir    = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke/trait_plots",
                                    scenario_keep = c(
                                      "ducke_deep_clayey_regclim_ks14",
                                      "ducke_deep_clayey_redprec30_ks14",
                                      "ducke_deep_sandy_redprec30_ks14",
                                      "ducke_deep_sandy_regclim_ks14",
                                      "ducke_shallow_sandy_regclim_ks14",
                                      "ducke_shallow_sandy_redprec30_ks14",
                                      "ducke_shallow_clayey_redprec30_ks14",
                                      "ducke_shallow_clayey_regclim_ks14"
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
                                        output_dir    = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/reserva_Ducke/trait_plots",
                                        scenario_keep = c(
                                          "ducke_deep_clayey_regclim_ks14",
                                          "ducke_deep_clayey_redprec30_ks14",
                                          "ducke_deep_sandy_redprec30_ks14",
                                          "ducke_deep_sandy_regclim_ks14",
                                          "ducke_shallow_sandy_regclim_ks14",
                                          "ducke_shallow_sandy_redprec30_ks14",
                                          "ducke_shallow_clayey_redprec30_ks14",
                                          "ducke_shallow_clayey_regclim_ks14"
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
      "ducke_deep_clayey_regclim_ks14"     = "#1B6CA8",
      "ducke_deep_clayey_redprec30_ks14"   = "#56B4E9",
      "ducke_deep_sandy_redprec30_ks14"  = "green4",
      "ducke_deep_sandy_regclim_ks14"  = "purple",
      "ducke_shallow_sandy_regclim_ks14"   = "#CC5500",
      "ducke_shallow_sandy_redprec30_ks14" = "#E8A268",
      "ducke_shallow_clayey_redprec30_ks14" = "chocolate4",
      "ducke_shallow_clayey_regclim_ks14" = "red"
      
    )
  }
  
  # Default scenario labels
  if (is.null(labels)) {
    labels <- c(
      "ducke_deep_clayey_regclim_ks14"     = "Deep WT | Clayey soil | Regclim",
      "ducke_deep_clayey_redprec30_ks14"   = "Deep WT | Clayey soil | Redprec",
      "ducke_deep_sandy_redprec30_ks14"   = "Deep WT | Sandy soil | Redprec",
      "ducke_deep_sandy_regclim_ks14"   = "Deep WT | Sandy soil | Regclim",
      "ducke_shallow_sandy_regclim_ks14"   = "Shallow WT | Sandy soil | Regclim",
      "ducke_shallow_sandy_redprec30_ks14" = "Shallow WT | Sandy soil | Redprec",
      "ducke_shallow_clayey_redprec30_ks14" = "Shallow WT | Clayey soil | Redprec",
      "ducke_shallow_clayey_regclim_ks14" = "Shallow WT | Clayey soil | Regclim"
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


# --- Example calls ---

wsg_density_gt10 <- plot_trait_distribution(
  data          = df_trait_dist,
  trait_name    = "wsg",
  dbh_threshold = 10,
  plot_type     = "density",
  save_plot     = FALSE
)
wsg_density_gt10

wsg_boxplot_gt10 <- plot_trait_weighted_boxplot(
  data          = df_trait_dist,
  trait_name    = "wsg",
  dbh_threshold = 10,
  save_plot     = FALSE
)
wsg_boxplot_gt10


# --- Batch export all traits ---

for (trait in trait_vars) {
  plot_trait_distribution(
    data          = df_trait_dist,
    trait_name    = trait,
    dbh_threshold = 10,
    plot_type     = "density",
    save_plot     = TRUE
  )
}

for (trait in trait_vars) {
  plot_trait_distribution(
    data          = df_trait_dist,
    trait_name    = trait,
    dbh_threshold = 10,
    plot_type     = "histogram",
    save_plot     = TRUE
  )
}

for (trait in trait_vars) {
  plot_trait_weighted_boxplot(
    data          = df_trait_dist,
    trait_name    = trait,
    dbh_threshold = 10,
    save_plot     = TRUE
  )
}
