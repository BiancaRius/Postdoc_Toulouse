---
title: "Bianca"
author: "S. Schmitt"
date: last-modified
format: 
  html:
    embed-resources: true
    toc: true
---

# Setup

```{r libs}
#| message: false
#| warning: false
library(tidyverse)
library(lubridate)
library(arrow)
library(dplyr)
library(ggplot2)
library(units)
library(scales)
library(patchwork)

```

The parameters are the vertical resolution of the soil depth (here 10 cm),
the layer depths used in the model,
and the dry and wet season month used for inter annual summary.

> Note that I used end of the rain season where precipiation is not the highest but soil water is.

```{r pars}
#| message: false
#| warning: false

# -----------------------------
# Global parameters
# -----------------------------
resolution <- 0.1
dry_months <- 10:11
wet_months <- 6:7

# # -----------------------------
# # Experiment definitions
# # -----------------------------
model_path <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/sensitivity_tests"
model_name <- c("shallow_variable_sandy", "deep_variable_sandy",
"shallow_variable_clayey", "deep_variable_clayey",
"shallow_intermediate_sandy", "deep_intermediate_sandy",
"shallow_intermediate_clayey", "deep_intermediate_clayey",
"shallow_thin_sandy", "deep_thin_sandy",
"shallow_thin_clayey", "deep_thin_clayey")

files_df <- tibble(model_name = model_name) %>%
  mutate(
    experiment_dir = file.path(model_path, model_name),
    pedology_path  = file.path(experiment_dir, "Paracou_input_pedology.txt"),
    waterbal_name  = paste0(model_name, "_0_water_balance.txt"),
    waterbal_path  = file.path(experiment_dir, "outputs", waterbal_name)
)

 
# -----------------------------
# Robust reader: remove trailing tabs/spaces before parsing
# -----------------------------
read_waterbal_clean <- function(path) {
  tmp <- tempfile(fileext = ".tsv")
  lines <- readLines(path, warn = FALSE)
  lines <- sub("\r$", "", lines)        # remove Windows CR if present
  lines <- sub("[\t ]+$", "", lines)    # remove trailing tabs/spaces
  writeLines(lines, tmp, useBytes = TRUE)

  readr::read_tsv(tmp, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::select(-dplyr::matches("^\\.\\.\\.[0-9]+$")) # drop any phantom cols if still created
}

```

# Data (uncomment it if you need to reset data)

First we need to compute layers depth structure with accumulated minimum and maximum values.
The pedology files provide only layer thickness. Here we convert thicknesses into cumulativedepth bounds for each layer: depth_max is the lower boundary (cumulative sum), and depth_min is the upper boundary (previous depth_max, with 0 at the surface). These depth intervals let
us map model outputs to real depths and compare experiments with different layer discretizations.

Then we open the TROLL output of water balance.
And we extract the soil water content and potentials (SWC & SWP) and add the corresponding layer depths.
```{r depths}
#| message: false
#| warning: false

# # -----------------------------
# # Build depth table from a pedology file
# # -----------------------------
# build_depths_from_pedology <- function(pedology_path) {
#   ped <- readr::read_tsv(pedology_path, show_col_types = FALSE, progress = FALSE)
# 
#   layer_thickness <- as.numeric(ped$layer_thickness)
# 
#   # (optional safety)
#   if (any(!is.finite(layer_thickness)) || any(layer_thickness <= 0, na.rm = TRUE)) {
#     stop("Non-finite or non-positive layer_thickness detected in: ", pedology_path)
#   }
# 
#   depth_max <- cumsum(layer_thickness)
#   depth_min <- c(0, head(depth_max, -1))
# 
#   tibble(
#     layer     = 0:(length(layer_thickness) - 1),
#     depth_min = depth_min,
#     depth_max = depth_max
#   )
# }
# 
# # -----------------------------
# # Read + pivot each experiment separately, then bind
# # -----------------------------
# data <- files_df %>%
#   mutate(dat = purrr::pmap(
#     list(model_name, waterbal_path, pedology_path),
#     function(model_name, waterbal_path, pedology_path) {
# 
#       depths <- build_depths_from_pedology(pedology_path)
# 
# # -----------------------------
# # Read all water balance files and attach depths by model_name + layer
# # -----------------------------
#       read_waterbal_clean(waterbal_path) %>%
#         mutate(
#           model_name = model_name,
#           date = as.Date("2004-01-01") + iter
#         ) %>%
#         select(model_name, date, matches("^(SWC|SWP)_[0-9]+$")) %>%
#         pivot_longer(
#           cols = -c(model_name, date),
#           names_to = "var_layer",
#           values_to = "value"
#         ) %>%
#         separate(var_layer, into = c("variable", "layer"), convert = TRUE) %>%
#         left_join(depths, by = "layer") %>%
#         filter(!is.na(depth_min), !is.na(depth_max)) %>%     # keep only real layers for that pedology
#         mutate(simulation = model_name) %>%
#         select(model_name, simulation, date, layer, depth_min, depth_max, variable, value)
#     }
#   )) %>%
#   select(dat) %>%
#   unnest(dat)
# 
# depths_df <- files_df %>%
#   mutate(depths = purrr::map(pedology_path, build_depths_from_pedology)) %>%
#   select(model_name, depths) %>%
#   tidyr::unnest(depths)
```

Resolution (uncomment it if you need to reset data)
To compute figures at a finer vertical resolution,
we first need to round up the depth to the resolution we will use,
and then we simply expand observed values to sublayers.
"Rasterize" each soil layer to a common fine depth grid (e.g., 0.1 m):
we round layer boundaries to avoid floating-point misalignment, then for each (model_name, date, layer, variable) we create a sequence of depths between depth_min and depth_max and replicate the same layer-mean value across all those sub-depths (piecewise-constant within the layer). This makes depth plots comparable across experiments with different layer thicknesses.


```{r raster}
#| message: false
#| warning: false
# raster <- data %>%
#   mutate(depth_min = round(depth_min, 1),
#          depth_max = round(depth_max, 1)) %>%
#   rowwise() %>%
#   mutate(depth = list(seq(depth_min,
#                           depth_max,
#                           by = resolution))) %>%
#   unnest(depth)
```

Save raster file and depths_df (uncomment it if you need to reset data)
```{r raster}

 
# write_parquet(raster,paste0(model_path,"/raster_data.parquet"))
# 
# write_parquet(depths_df,paste0(model_path,"/depths_metadata.parquet"))

```

Import and open the complete rasterized dataset:

```{r raster}
# Import raster data
raster <- read_parquet(paste0(model_path,"/raster_data.parquet"))
depths_df <- read_parquet(paste0(model_path,"/depths_metadata.parquet"))
```


# Figures

And here are the codes to represent SWC and SWP as raster with depth and time or
mean dry and wet season profiles.

Plot pairs without saving it. Uncomment if you want.
```{r }
#| message: false
#| warning: false
#| fig-cap: "Soil water content (SWC) across time (x-axis) and depth (y-axis), plotted pairwise (shallow vs deep) for each pedology setup."

# # -----------------------------
# # Define the pairs (shallow vs deep)
# # -----------------------------
# pairs <- list(
#   list(key = "variable_sandy",        label = "Variable — Sandy",        models = c("shallow_variable_sandy",        "deep_variable_sandy")),
#   list(key = "variable_clayey",       label = "Variable — Clayey",       models = c("shallow_variable_clayey",       "deep_variable_clayey")),
#   list(key = "intermediate_sandy",    label = "Intermediate — Sandy",    models = c("shallow_intermediate_sandy",    "deep_intermediate_sandy")),
#   list(key = "intermediate_clayey",   label = "Intermediate — Clayey",   models = c("shallow_intermediate_clayey",   "deep_intermediate_clayey")),
#   list(key = "thin_sandy",            label = "Thin — Sandy",            models = c("shallow_thin_sandy",            "deep_thin_sandy")),
#   list(key = "thin_clayey",           label = "Thin — Clayey",           models = c("shallow_thin_clayey",           "deep_thin_clayey"))
# )
# 
# # -----------------------------
# # Open Arrow dataset once
# # -----------------------------
# ds_raster <- open_dataset(file.path(model_path, "raster_data.parquet"))
# 
# # -----------------------------
# # Loop over pairs and plot
# # -----------------------------
# for (pp in pairs) {
# 
#   pair_models <- pp$models
# 
#   # Layer boundaries for this pair
#   layer_bounds_pair <- depths_df %>%
#     filter(model_name %in% pair_models) %>%
#     distinct(model_name, depth_max) %>%
#     mutate(model_name = factor(model_name, levels = pair_models))
# 
#   # Collect ONLY this pair (much lighter)
#   raster_pair <- ds_raster %>%
#     filter(variable == "SWC", model_name %in% pair_models) %>%
#     select(model_name, date, depth, value) %>%
#     collect() %>%
#     mutate(model_name = factor(model_name, levels = pair_models))
# 
#   # Plot (2 panels)
#   p <- ggplot(raster_pair, aes(date, depth, fill = value)) +
#     geom_raster() +
#     scale_y_reverse() +
#     scale_fill_viridis_c(
#       expression("SWC ["~m^3~m^3~"]"),
#       direction = -1,
#       limits = c(0.1, 0.6)
#     ) +
#     labs(
#       title = pp$label,
#       x = NULL,
#       y = "Depth [ m ]"
#     ) +
#     facet_wrap(~model_name, ncol = 2) +
#     theme_bw() +
#     theme(
#       strip.text = element_text(size = 8),
#       legend.position = "bottom",
#       panel.spacing = unit(0.5, "lines"),
#       plot.title = element_text(size = 10, face = "bold")
#     )
# 
#   print(p)
# 
#   # Free memory between pairs
#   rm(raster_pair)
#   gc()
# }

```


Plot and save SWC by pairs
```{r}

# Define the pairs in the desired order
pairs <- list(
  list(label = "Variable — Sandy",      models = c("shallow_variable_sandy",      "deep_variable_sandy")),
  list(label = "Variable — Clayey",     models = c("shallow_variable_clayey",     "deep_variable_clayey")),
  list(label = "Intermediate — Sandy",  models = c("shallow_intermediate_sandy",  "deep_intermediate_sandy")),
  list(label = "Intermediate — Clayey", models = c("shallow_intermediate_clayey", "deep_intermediate_clayey")),
  list(label = "Thin — Sandy",          models = c("shallow_thin_sandy",          "deep_thin_sandy")),
  list(label = "Thin — Clayey",         models = c("shallow_thin_clayey",         "deep_thin_clayey"))
)

# List to store plots
plot_list <- list()

for (i in seq_along(pairs)) {
  
  pp <- pairs[[i]]
  pair_models <- pp$models
  
  # 1. Prepare layer boundaries (only for this pair)
  layer_bounds_pair <- depths_df %>%
    filter(model_name %in% pair_models) %>%
    distinct(model_name, depth_max) %>%
    mutate(model_name = factor(model_name, levels = pair_models))

  # 2. Load data
  raster_pair <- raster %>%
  filter(variable == "SWC", model_name %in% pair_models) %>%
  select(model_name, date, depth, value) %>%
  collect() %>%
  mutate(
    model_name = factor(model_name, levels = pair_models),
    days_run = as.numeric(date - start_date),
    sim_year = days_run / 365
  )

  # 3. Plotting with fixes
  p <- ggplot(raster_pair, aes(sim_year, depth, fill = value)) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(0, ceiling(max(raster_pair$sim_year, na.rm = TRUE)), by = 5)
  ) +
  geom_tile(width = 10/365) +  # 10 dias em "anos"
  scale_y_reverse() +
  scale_fill_viridis_c(expression("SWC ["~m^3~m^3~"]"), direction = -1, limits = c(0.1, 0.6)) +
  labs(title = pp$label, x = NULL, y = "Depth [m]") +
  facet_wrap(~model_name, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

  
  # If this is the last plot (i == 6), bring the x-axis back
  if (i == length(pairs)) {
    p <- p + theme(axis.text.x = element_text(), axis.ticks.x = element_line())
  }
  
  plot_list[[i]] <- p
  
  # Memory cleanup
  rm(raster_pair)
  gc()
}

# COMBINE AND SAVE

# Patchwork to stack vertically (ncol = 1)
final_figure <- wrap_plots(plot_list, ncol = 1) + 
  plot_layout(guides = "collect") & # Collect the legend
  theme(legend.position = "bottom")

# Save
# ggsave("Figure_SWC_allpairs.png", final_figure, width = 10, height = 18, dpi = 300)
# message("Figure saved successfully!")
```


Plot SWP by pairs - viridis gradient
```{r}

pairs <- list(
  list(label = "",      models = c("shallow_variable_sandy",      "shallow_thin_sandy")), list(label = "",     models = c("deep_variable_sandy",     "deep_thin_sandy")))
# pairs <- list(
#   list(label = "Variable — Sandy",      models = c("shallow_variable_sandy",      "deep_variable_sandy")),
#   list(label = "Variable — Clayey",     models = c("shallow_variable_clayey",     "deep_variable_clayey")),
#   list(label = "Intermediate — Sandy",  models = c("shallow_intermediate_sandy",  "deep_intermediate_sandy")),
#   list(label = "Intermediate — Clayey", models = c("shallow_intermediate_clayey", "deep_intermediate_clayey")),
#   list(label = "Thin — Sandy",          models = c("shallow_thin_sandy",          "deep_thin_sandy")),
#   list(label = "Thin — Clayey",         models = c("shallow_thin_clayey",         "deep_thin_clayey"))
# )

swp_vals <- raster %>%
  filter(variable == "SWP") %>%
  select(value) %>%          # ONLY this column (light)
  collect()

val_plot <- -swp_vals$value  # EXACTLY like your simple plot: fill = -value
val_plot <- val_plot[is.finite(val_plot) & val_plot > 0]  # keep valid for log scale

swp_min <- min(val_plot, na.rm = TRUE)
swp_max <- max(val_plot, na.rm = TRUE)

swp_scale <- scale_fill_viridis_c(
  expression("|SWP| ["~MPa~"]"),
  direction = 1,
  trans = "log",
  limits = c(swp_min, swp_max),
  na.value = "grey" # <--- color the water table
)


rm(swp_vals, val_plot); gc()

start_date <- as.Date("2004-01-01")
  
plot_list <- list()

for (i in seq_along(pairs)) {

  pp <- pairs[[i]]
  pair_models <- pp$models

  raster_pair <- raster %>%
    filter(variable == "SWP", model_name %in% pair_models) %>%
    select(model_name, date, depth, value) %>%
    collect() %>%
    mutate(
      model_name = factor(model_name, levels = pair_models),
      
      # --- TRANSFORMATION STEPS ---
      # 1. Recover 'iter' (days since start)
      days_run = as.numeric(date - start_date),
      
      # 2. Convert to Years (using 365 days)
      sim_year = days_run / 365,
      
      # 3. Handle SWP sign and Saturation (blue color)
      value_plot = -value,
      value_plot = ifelse(value_plot <= 0, NA, value_plot) 
    )

  # Use sim_year on x-axis
  p <- ggplot(raster_pair, aes(sim_year, depth, fill = value_plot)) +
    geom_tile() + 
    theme_bw() +
    scale_y_reverse() +
    swp_scale + 
    xlab("") + # Axis label hidden for all except the last one
    ylab("Depth [ m ]") +
    facet_wrap(~model_name, ncol = 2, scales = "free_y") +
    ggtitle(pp$label) +
    
    # Force the x-axis to start exactly at 0 if desired
    scale_x_continuous(expand = c(0, 0),
      breaks = seq(0,    
      ceiling(max(raster_pair$sim_year)),   by = 5)
  ) + 
    
    theme(
      strip.text = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold"),
      panel.spacing = unit(0.5, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  # Only add the "Years" label and ticks to the bottom-most plot
  if (i == length(pairs)) {
    p <- p + 
      xlab("Time [ years ]") + 
      theme(axis.text.x = element_text(), axis.ticks.x = element_line())
  }

  plot_list[[i]] <- p

  rm(raster_pair)
  gc()
}

final_figure_swp <- wrap_plots(plot_list, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text  = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.title = element_text(angle = 0)  # opcional
  )


plot(final_figure_swp)

ggsave("Figure_SWP_variable_thin.png",
        final_figure_swp, width = 10, height = 10, dpi = 300)

# ggsave("Figure_SWP_deep_variable.png",
#         final_figure_swp, width = 10, height = 8,
# dpi = 300)

```


```{r swp_rast}
#| message: false
#| warning: false
#| fig-cap: "Soil water potential across time (x-axis) and depth (y-axis)."


raster %>% 
  filter(variable == "SWP") %>% 
  ggplot(aes(date, depth, fill = -value)) +
  geom_raster() +
  theme_bw() +
  scale_y_reverse() +
  scale_fill_viridis_c(expression("|SWP| ["~MPa~"]"),
                       direction = 1, trans = "log") +
  xlab("") +
  ylab("Depth [ m ]")+facet_wrap(~model_name, nrow = 6)
```

```{r}
library(ggplot2)
library(dplyr)
library(stringr) # Para limpar o nome do arquivo

# Loop para iterar sobre cada item da lista 'pairs'
for (pair in pairs) {
  
  # 1. Extrair informações do par atual
  plot_label <- pair$label
  model_names <- pair$models
  
  # 2. Filtrar os dados apenas para os modelos desse par
  plot_data <- raster %>% 
    filter(variable == "SWP") %>% 
    filter(model_name %in% model_names)
  
  # 3. Forçar a ordem dos fatores para aparecer na ordem da lista (Esq -> Dir)
  # Isso garante que o primeiro modelo da lista fique na esquerda
  plot_data$model_name <- factor(plot_data$model_name, levels = model_names)
  
  # 4. Criar o gráfico
  p <- ggplot(plot_data, aes(date, depth, fill = -value)) +
    geom_raster() +
    theme_bw() +
    scale_y_reverse() +
    # Ajustei a escala de cores baseada no seu código anterior
    scale_fill_viridis_c(expression("|SWP| ["~MPa~"]"),
                         direction = 1, 
                         trans = "log",
                         # Sugestão: adicione limits se quiser fixar a escala igual para todos
                         # limits = c(0.001, 1000), 
                         na.value = "grey50") + 
    xlab("") +
    ylab("Depth [ m ]") +
    facet_wrap(~model_name, nrow = 1, scales = "free_y") +
    ggtitle(plot_label) # Adiciona o título do par
  
  # 5. Gerar um nome de arquivo seguro (sem espaços ou caracteres especiais)
  # Exemplo: "Variable_Sandy.png"
  file_name <- paste0("SWP2_", str_replace_all(plot_label, "[^[:alnum:]]", "_"), ".png")
  
  # 6. Salvar
  # Ajuste width e height conforme necessário
  # ggsave(filename = file_name, plot = p, width = 10, height = 5, dpi = 300)
  # 
  # print(paste("Salvo:", file_name))
}
```

## Profiles

> Profiles are made with daily means across years, but maybe we should use first a season mean within year and then summarise it across years.

```{r swc_profile}
#| message: false
#| warning: false
#| fig-cap: "Mean soil water content for wet and dry season across depth (y-axis). Point represent the median across days and bars the 90% interval."
raster %>% 
  filter(variable == "SWC") %>% 
  filter(month(date) %in% c(dry_months, wet_months)) %>% 
  mutate(season = ifelse(month(date) %in% dry_months, "dry", "wet")) %>% 
  group_by(simulation,variable, depth_min, depth_max, season) %>% 
  summarise(m = median(value), 
            ll = quantile(value, .05), 
            hh = quantile(value, .95)) %>% 
  mutate(depth = (depth_min + depth_max)/2) %>% 
  ggplot(aes(depth_max, m, col = paste(simulation, season))) +
  geom_line() +
  geom_linerange(aes(ymin = ll, ymax = hh)) +
  geom_point() +
  theme_bw() +
  coord_flip() +
  xlab("Depth [ m ]") +
  scale_y_reverse(expression("SWC ["~m^3~m^3~"]"))
```


```{r}
library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)

pairs <- list(
  list(label = "Variable — Sandy",      models = c("shallow_variable_sandy",      "deep_variable_sandy")),
  list(label = "Variable — Clayey",     models = c("shallow_variable_clayey",     "deep_variable_clayey")),
  list(label = "Intermediate — Sandy",  models = c("shallow_intermediate_sandy",  "deep_intermediate_sandy")),
  list(label = "Intermediate — Clayey", models = c("shallow_intermediate_clayey", "deep_intermediate_clayey")),
  list(label = "Thin — Sandy",          models = c("shallow_thin_sandy",          "deep_thin_sandy")),
  list(label = "Thin — Clayey",         models = c("shallow_thin_clayey",         "deep_thin_clayey"))
)

plot_list <- list()

for (i in seq_along(pairs)) {

  pp <- pairs[[i]]
  pair_models <- pp$models

  df_pair <- raster %>%
    filter(variable == "SWC") %>%
    filter(model_name %in% pair_models) %>%
    filter(month(date) %in% c(dry_months, wet_months)) %>%
    mutate(
      season = ifelse(month(date) %in% dry_months, "dry", "wet"),
      simulation = factor(model_name, levels = pair_models)
    ) %>%
    group_by(simulation, depth_min, depth_max, season) %>%
    summarise(
      m  = median(value, na.rm = TRUE),
      ll = quantile(value, 0.05, na.rm = TRUE),
      hh = quantile(value, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(depth = (depth_min + depth_max) / 2)

  p <- ggplot(df_pair, aes(depth_max, m, color = season, group = season)) +
    geom_line(linewidth = 0.6) +
    geom_linerange(aes(ymin = ll, ymax = hh), linewidth = 0.5) +
    geom_point(size = 1) +
    coord_flip() +
    scale_x_reverse("Depth [ m ]") +
    scale_y_reverse(expression("SWC ["~m^3~m^3~"]"))+
    labs(title = pp$label, color = NULL) +
    facet_wrap(~simulation, ncol = 2) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 8)
    )

  plot_list[[i]] <- p
}

final_fig <- wrap_plots(plot_list, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("SWC_profiles_pairs_wet_dry.png", final_fig, width = 10, height = 18, dpi = 300)
# message("Saved: SWC_profiles_pairs_wet_dry_freeY.png")


```

```{r swp_profile}
#| message: false
#| warning: false
#| fig-cap: "Mean soil water potential for wet and dry season across depth (y-axis). Point represent the median across days and bars the 90% interval."
raster %>% 
  filter(variable == "SWP") %>% 
  filter(month(date) %in% c(dry_months, wet_months)) %>% 
  mutate(season = ifelse(month(date) %in% dry_months, "dry", "wet")) %>% 
  group_by(simulation,variable, depth_min, depth_max, season) %>% 
  mutate(value = -value) %>% 
  summarise(m = median(value), 
            ll = quantile(value, .05), 
            hh = quantile(value, .95)) %>% 
  mutate(depth = (depth_min + depth_max)/2) %>% 
  ggplot(aes(depth, m, col = paste(simulation,season))) +
  geom_line() +
  geom_linerange(aes(ymin = ll, ymax = hh)) +
  geom_point() +
  theme_bw() +
  coord_flip() +
  scale_x_reverse("Depth [ m ]") +
  scale_y_log10(expression("|SWP| ["~MPa~"]"))
```

```{r}
library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)

pairs <- list(
  list(label = "Variable — Sandy",      models = c("shallow_variable_sandy",      "deep_variable_sandy")), list(label = "Variable — Clayey",     models = c("shallow_variable_clayey",     "deep_variable_clayey")))

  # pairs <- list(
#   list(label = "Variable — Sandy",      models = c("shallow_variable_sandy",      "deep_variable_sandy")),
#   list(label = "Variable — Clayey",     models = c("shallow_variable_clayey",     "deep_variable_clayey")),
#   list(label = "Intermediate — Sandy",  models = c("shallow_intermediate_sandy",  "deep_intermediate_sandy")),
#   list(label = "Intermediate — Clayey", models = c("shallow_intermediate_clayey", "deep_intermediate_clayey")),
#   list(label = "Thin — Sandy",          models = c("shallow_thin_sandy",          "deep_thin_sandy")),
#   list(label = "Thin — Clayey",         models = c("shallow_thin_clayey",         "deep_thin_clayey"))
# )

plot_list <- list()

for (i in seq_along(pairs)) {

  pp <- pairs[[i]]
  pair_models <- pp$models

  df_pair <- raster %>%
    filter(variable == "SWP") %>%
    filter(model_name %in% pair_models) %>%   # only this pair
    filter(month(date) %in% c(dry_months, wet_months)) %>%
    mutate(
      season = ifelse(month(date) %in% dry_months, "dry", "wet"),
      simulation = factor(model_name, levels = pair_models),
      value = -value
    ) %>%
    # IMPORTANT for log: remove non-positive after flipping sign
    filter(is.finite(value), value > 0) %>%
    group_by(simulation, depth_min, depth_max, season) %>%
    summarise(
      m  = median(value, na.rm = TRUE),
      ll = quantile(value, 0.05, na.rm = TRUE),
      hh = quantile(value, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(depth = (depth_min + depth_max) / 2)

  p <- ggplot(df_pair, aes(depth_max, m, color = season, group = season)) +
    geom_line(linewidth = 0.6) +
    geom_linerange(aes(ymin = ll, ymax = hh), linewidth = 0.5) +
    geom_point(size = 1) +
    theme_bw() +
    coord_flip() +
    scale_x_reverse("Depth [ m ]") +
    scale_y_log10(expression("|SWP| ["~MPa~"]")) +
    labs(title = pp$label, color = NULL) +
    facet_wrap(~simulation, ncol = 2) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 8)
    )

  plot_list[[i]] <- p
}

final_fig <- wrap_plots(plot_list, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("SWP_profiles_pairs_clayey_sandy_variable.png", final_fig, width = 10, height = 8, dpi = 300)
# ggsave("SWP_profiles_pairs_wet_dry_log10.png", final_fig, width = 10, height = 18, dpi = 300)
# message("Saved: SWP_profiles_pairs_wet_dry_log10_freeY.png")

```