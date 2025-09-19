# -------------------------------------------
# Script to plot results from TROLL simulations with fixed WTD values and
# Capillarity implementation
# Author: Bianca Fazio Rius
# -------------------------------------------

library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)


# ---- 1) Caminhos dos dois cenários ----
#scenario_paths <- list(
  #maximum of water in the layer = max_SWC
 # max_maxswc = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/Capillarity_implementation/deep_WTD_capON/SWC_updated/max_eq_max_swc/",
  #maximum of water in the layer = field_capacity
  #max_fc = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/Capillarity_implementation/deep_WTD_capON/SWC_updated/max_eq_field_capacity/"
#)

#scenario_paths <- list(
  #maximum of water in the layer = max_SWC
#max_maxswc = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/Capillarity_implementation/deep_WTD_capON/SWC_updated/max_eq_max_swc/",
#maximum of water in the layer = field_capacity, donor capacity = infinity for WT, cseedrain = 1
#max_fc = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/Capillarity_implementation/deep_WTD_capON/SWC_updated/max_eq_field_capacity_don_infinity/"
#)

#scenario_paths <- list(
  ##maximum of water in the layer = field_capacity, donor capacity = not infinity for WT, cseedrain = 1
 # max_fc_notinf = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/Capillarity_implementation/deep_WTD_capON/SWC_updated/max_eq_field_capacity_don_notinfinity",
  #maximum of water in the layer = field_capacity, donor capacity = infinity for WT, cseedrain = 1
  #max_fc = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/Capillarity_implementation/deep_WTD_capON/SWC_updated/max_eq_field_capacity_don_infinity/"
#)

scenario_paths <- list(
#maximum of water in the layer = field_capacity, donor capacity = not infinity for WT, cseedrain = 1
  max_fc_notinf = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/Capillarity_implementation/deep_WTD_capON/SWC_updated/max_eq_field_capacity_don_notinfinity",
#maximum of water in the layer = field_capacity, donor capacity = not infinity for WT, cseedrain = 1, wt = max_SWC commented
  max_fc_notinf_wt = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/Capillarity_implementation/deep_WTD_capON/SWC_updated/max_eq_field_capacity_don_notinfinity_wt/"
)

# ---- 2) Grupos de variáveis e mapeamento p/ arquivo ----
biogeochemical_vars   <- c("npp", "gpp", "agb", "sum1", "sum10", "sum30", "ba", "ba10", "litterfall")
soil_water_content    <- c("SWC_0", "SWC_1", "SWC_2", "SWC_3", "SWC_4")
soil_water_potential  <- c("SWP_0", "SWP_1", "SWP_2", "SWP_3", "SWP_4")
# Corrija aqui se o nome das colunas no arquivo tem typo:
transpiration_layers  <- c("transpiration_0", "transpiration_1", "transpiration_2", "transpiration_3", "transpiration_4")
water_flux_vars       <- c("precipitation", "interception", "throughfall", "runoff", "leak", "evaporation")

variable_groups <- list(
  biogeochemical        = biogeochemical_vars,
  soil_water_content    = soil_water_content,
  soil_water_potential  = soil_water_potential,
  transpiration_layers  = transpiration_layers,
  water_fluxes          = water_flux_vars
)

variable_file_map <- c(
  biogeochemical        = "sumstats",
  soil_water_content    = "water_balance",
  soil_water_potential  = "water_balance",
  transpiration_layers  = "water_balance",
  water_fluxes          = "water_balance"
)

# ---- 3) Função robusta de leitura por cenário ----
read_scenario_data <- function(path, file_type, scenario_label) {
  file_path <- file.path(path, paste0("(null)_0_", file_type, ".txt"))
  if (!file.exists(file_path)) {
    stop("Arquivo não encontrado: ", file_path)
  }
  # read.table é OK; readr::read_tsv/read_csv pode ser mais rápido se os separadores permitirem
  df <- read.table(file_path, header = TRUE)
  df$scenario <- scenario_label
  df
}

# ---- 4) Função de plot para comparar cenários ----
plot_variable_two_scenarios <- function(variable, paths, file_type) {
  scenarios <- names(paths)
  
  all_data <- lapply(scenarios, function(sc) {
    read_scenario_data(paths[[sc]], file_type, sc)
  }) %>% bind_rows()
  
  # Checa existência da coluna solicitada
  if (!variable %in% names(all_data)) {
    stop("A coluna '", variable, "' não existe no arquivo '", file_type, "'.")
  }
  
  # Paleta fixa por cenário
  scen_labels <- c(max_fc_notinf = "max = fc not inf", max_fc_notinf_wt = "max = fc not inf wt")
  scen_colors <- c("max = fc not inf" = "white", "max = fc not inf wt " = "purple")
  
  all_data <- all_data %>%
    mutate(scenario = factor(scenario, levels = names(scen_labels), labels = unname(scen_labels)))
  
  p <- ggplot(all_data, aes(x = iter, y = .data[[variable]], color = scenario)) +
    geom_line(size = 0.7, alpha = 0.7) +
    scale_color_manual(values = scen_colors, name = NULL) +
    scale_x_continuous(
      name = "Year",
      breaks = seq(0, 10000, by = 365 * 3),
      labels = function(x) floor(x / 365) + 1
    ) +
    labs(y = variable, title = variable) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  # Exemplo de ajuste de limites (se quiser aplicar só para SWC_4, como no script original)
  # if (variable == "SWC_4") {
  #   y_limits <- range(all_data$SWC_4, na.rm = TRUE)
  #   p <- p + ylim(y_limits)
  # }
  
  p
}

# ---- 5) Função principal ----
plot_results_two_scenarios <- function(plot_all = TRUE, selected_variable = NULL, save_pdf = FALSE, pdf_file = "TROLL_2scenarios_plots.pdf") {
  
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
        print(plot_variable_two_scenarios(var, scenario_paths, file_type))
      }
    } else if (!is.null(selected_variable) && selected_variable %in% vars) {
      message("Plotting selected variable: ", selected_variable)
      print(plot_variable_two_scenarios(selected_variable, scenario_paths, file_type))
    }
  }
  
  if (save_pdf) {
    message("PDF salvo em: ", normalizePath(pdf_file))
  }
}

# ---- 6) Exemplos de chamada ----
# Um único gráfico (na tela)
plot_results_two_scenarios(plot_all = FALSE, selected_variable = "npp")
plot_results_two_scenarios(plot_all = FALSE, selected_variable = "SWC_0")
plot_results_two_scenarios(plot_all = FALSE, selected_variable = "SWC_1")
plot_results_two_scenarios(plot_all = FALSE, selected_variable = "SWC_2")
plot_results_two_scenarios(plot_all = FALSE, selected_variable = "SWC_3")
plot_results_two_scenarios(plot_all = FALSE, selected_variable = "SWC_4")




# Todos os gráficos, salvando em PDF
# plot_results_two_scenarios(plot_all = TRUE, save_pdf = TRUE, pdf_file = "TROLL_2scenarios_allvars.pdf")
