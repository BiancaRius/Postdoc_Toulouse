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
scenario_paths <- c("wtOn_capOff")

# ========== 2) Variable groups and corresponding input file types ==========
biogeochemical_vars   <- c("npp", "gpp", "agb", "sum1", "sum10", "sum30", "ba", "ba10", "litterfall")
soil_water_content    <- c("SWC_0", "SWC_1", "SWC_2", "SWC_3", "SWC_4")
soil_water_potential  <- c("SWP_0", "SWP_1", "SWP_2", "SWP_3", "SWP_4")
transpiration_layers  <- c("transpiration_0", "transpiration_1", "transpiration_2", "transpiration_3", "transpiration_4")
water_flux_vars       <- c("precipitation", "interception", "throughfall", "runoff", "leak", "evaporation")

variable_groups <- list(
  biogeochemical       = biogeochemical_vars,
  soil_water_content   = soil_water_content,
  soil_water_potential = soil_water_potential,
  transpiration_layers = transpiration_layers,
  water_fluxes         = water_flux_vars
)

variable_file_map <- c(
  biogeochemical       = "sumstats",
  soil_water_content   = "water_balance",
  soil_water_potential = "water_balance",
  transpiration_layers = "water_balance",
  water_fluxes         = "water_balance"
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
      labs(title = v, y = "Soil Water Content") +
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