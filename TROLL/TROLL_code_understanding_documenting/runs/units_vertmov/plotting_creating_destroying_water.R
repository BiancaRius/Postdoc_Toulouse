library(tidyverse)
library(stringr)

project_dir <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/"

experiment_table <- tibble(
  family = c("units_vertmov", "safety_margin_capacities", "units_vertmov", "safety_margin_capacities", "units_vertmov", "safety_margin_capacities", "units_vertmov", "safety_margin_capacities"),
  vegetation = c("veg", "veg", "veg", "veg", "veg", "veg", "veg", "veg"),
  water_table = c("shallowWT", "shallowWT", "deepWT", "deepWT", "shallowWT", "shallowWT", "deepWT", "deepWT"),
  texture = c("sandy", "sandy","sandy", "sandy", "clayey","clayey","clayey","clayey"),
  experiment_name = c(
    "shallowWT_sandy_ksgeom",
    "safetymargin_shallowWT_sandy",
    "deepWT_sandy_ksgeom",
    "safetymargin_deepWT_sandy",
    "shallowWT_clayey_ksgeom",
    "safetymargin_shallowWT_clayey",
    "deepWT_clayey_ksgeom",
    "safetymargin_deepWT_clayey"
    
  )
)

experiment_table <- experiment_table %>%
  mutate(
    log_file = file.path(
      project_dir,
      family,
      vegetation,
      water_table,
      texture,
      experiment_name,
      "run.log"
    )
  )

experiment_table %>%
  mutate(file_exists = file.exists(log_file))

parse_troll_clamp_log <- function(log_file) {
  
  lines <- readLines(log_file, warn = FALSE)
  
  water_lines <- lines[str_detect(lines, "annual clamp-created water")]
  
  tibble(raw = water_lines) %>%
    mutate(
      year = as.numeric(str_extract(raw, "(?<=Year )\\d+")),
      created_mm = as.numeric(str_extract(raw, "(?<=created water \\(mm\\) = )[-0-9\\.eE]+")),
      destroyed_mm = as.numeric(str_extract(raw, "(?<=destroyed water \\(mm\\) = )[-0-9\\.eE]+")),
      net_mm = as.numeric(str_extract(raw, "(?<=net water change \\(mm\\) = )[-0-9\\.eE]+"))
    ) %>%
    select(year, created_mm, destroyed_mm, net_mm)
}

df_clamps <- experiment_table %>%
  mutate(data = map(log_file, parse_troll_clamp_log)) %>%
  unnest(data)

df_clamps

summary_clamps <- df_clamps %>%
  group_by(family, vegetation, water_table, texture, experiment_name) %>%
  summarise(
    total_created_mm = sum(created_mm, na.rm = TRUE),
    total_destroyed_mm = sum(destroyed_mm, na.rm = TRUE),
    total_net_mm = sum(net_mm, na.rm = TRUE),
    mean_created_mm_year = mean(created_mm, na.rm = TRUE),
    mean_destroyed_mm_year = mean(destroyed_mm, na.rm = TRUE),
    mean_net_mm_year = mean(net_mm, na.rm = TRUE),
    max_created_mm_year = max(created_mm, na.rm = TRUE),
    max_destroyed_mm_year = max(destroyed_mm, na.rm = TRUE),
    max_abs_net_mm_year = max(abs(net_mm), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(total_net_mm)))

summary_clamps

ggplot(df_clamps, aes(x = year, y = net_mm, color = experiment_name)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_grid(texture ~ water_table, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Simulation year",
    y = "Net artificial water change by clamp (mm)",
    color = "Experiment"
  )

df_clamps_long <- df_clamps %>%
  select(family, vegetation, water_table, texture, experiment_name,
         year, created_mm, destroyed_mm, net_mm) %>%
  pivot_longer(
    cols = c(created_mm, destroyed_mm, net_mm),
    names_to = "variable",
    values_to = "water_mm"
  )

ggplot(df_clamps_long, aes(x = year, y = water_mm, color = variable)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ experiment_name, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Simulation year",
    y = "Artificial water change by clamp (mm)",
    color = NULL
  )

df_clamps_long <- df_clamps_long %>%
  mutate(
    variable = recode(
      variable,
      created_mm = "Created water",
      destroyed_mm = "Destroyed water",
      net_mm = "Net change"
    )
  )

summary_clamps