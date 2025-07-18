# -------------------------------------------
# Script to plot results from TROLL simulations with fixed WTD values
# Comparing 4 climate conditions: regular climate, 30% red.prec. (reduced precipitation),
# 50% red. prec. and 70% red. prec.
# Author: Bianca Fazio Rius
# -------------------------------------------

library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)

# # ---- 1. Define paths to simulation outputs ----

climate_paths = list(
  regular = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/",
  redprec30 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_30/",
  redprec50 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_50/",
  redprec70 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_70/"
)


# ---- 2. Define variable groups and types ----
biogeochemical_vars <- c("npp", "gpp", "agb", "sum1", "sum10", "sum30", "ba", "ba10", "litterfall")
soil_water_content  <- c("SWC_0", "SWC_1", "SWC_2", "SWC_3", "SWC_4")
soil_water_potential <- c("SWP_0", "SWP_1", "SWP_2", "SWP_3", "SWP_4")
transpiration_layers <- c("transpitation_0", "transpitation_1", "transpitation_2", "transpitation_3", "transpitation_4")
water_flux_vars     <- c("precipitation", "interception", "throughfall", "runoff", "leak", "evaporation")

# Variable group to metadata mapping
variable_types <- c("biogeochemical", "soil water content", "soil water potential", "transpiration layers", "water fluxes")

get_variable_info <- function(variable_type) {
  switch(variable_type,
         "biogeochemical"     = list(vars = biogeochemical_vars, file = "sumstats"),
         "soil water content" = list(vars = soil_water_content, file = "water_balance"),
         "soil water potential" = list(vars = soil_water_potential, file = "water_balance"),
         "transpiration layers" = list(vars = transpiration_layers, file = "water_balance"),
         "water fluxes"         = list(vars = water_flux_vars, file = "water_balance"),
         stop("Unknown variable type.")
  )
}

# ---- 3. Function to read data for a given WTD condition ----
read_simulation_data <- function(path, wtd, file_type, climate) {
  file_path <- file.path(path, wtd, paste0("(null)_0_", file_type, ".txt"))
  df <- read.table(file_path, header = TRUE)
  df$WTD <- wtd
  df$climate <- climate
  return(df)
}


# ---- 4. Function to create plots ----
plot_variable_facet <- function(variable, paths, file_type, by_WTD = FALSE) {
  climate_keys <- names(paths)
  climate_labels <- c("Reg. Clim.", "Red. Prec. 30%", "Red. Prec. 50%", "Red. Prec. 70%")
  wtd_levels <- c("no_WTD", "shallow_WTD", "deep_WTD")
  wtd_labels <- c("No WTD", "Shallow WTD", "Deep WTD")
  
  # Definir cores SEMPRE pela WTD
  color_var <- "WTD"
  color_values <- c(
    "No WTD" = "darkblue",
    "Shallow WTD" = "darkgreen",
    "Deep WTD" = "darkorange"
  )
  
  # Ler todos os dados
  all_data <- lapply(seq_along(climate_keys), function(i) {
    climate <- climate_keys[i]
    path <- paths[[climate]]
    
    bind_rows(
      read_simulation_data(path, "no_WTD", file_type, climate),
      read_simulation_data(path, "shallow_WTD", file_type, climate),
      read_simulation_data(path, "deep_WTD", file_type, climate)
    )
  }) %>% bind_rows()
  
  # Formatar fatores
  all_data <- all_data %>%
    mutate(
      Climate = factor(climate, levels = climate_keys, labels = climate_labels),
      WTD = factor(WTD, levels = wtd_levels, labels = wtd_labels)
    )
  
  # Ajustar limites de y para SWC_4
  if (variable %in% soil_water_content) {
    swc4_vals <- all_data$SWC_4
    y_limits <- range(swc4_vals, na.rm = TRUE)
  } else {
    y_limits <- NULL
  }
  
  if (!by_WTD) {
    # Modo original: comparando WTD dentro de cada clima
    p <- ggplot(all_data, aes(x = iter, y = .data[[variable]], color = .data[[color_var]])) +
      geom_line(size = 0.7, alpha = 0.5) +
      facet_wrap(~ Climate, nrow = 1) +
      scale_color_manual(values = color_values) +
      scale_x_continuous(
        name = "Year",
        breaks = seq(0, 10000, by = 365 * 3),
        labels = function(x) floor(x / 365) + 1
      ) +
      labs(y = variable, color = "Water Table Depth") +
      theme_minimal() +
      theme(
        # Remova ou comente a linha abaixo para mostrar os rótulos de clima
        # strip.text = element_blank(), 
        legend.position = 'bottom'
      )
    
    if (!is.null(y_limits)) {
      p <- p + ylim(y_limits)
    }
    
    return(p)
    
  } else {
    # Novo modo: um gráfico para cada WTD
    plot_list <- lapply(seq_along(levels(all_data$WTD)), function(i) {
      wtd_level <- levels(all_data$WTD)[i]
      subset_data <- all_data %>% filter(WTD == wtd_level)
      
      p <- ggplot(subset_data, aes(x = iter, y = .data[[variable]], color = .data[[color_var]])) +
        geom_line(size = 0.7, alpha = 0.7) +
        facet_wrap(~ Climate, nrow = 1) +
        scale_color_manual(values = color_values) +
        scale_x_continuous(
          name = "Year",
          breaks = seq(0, 10000, by = 365 * 3),
          labels = function(x) floor(x / 365) + 1
        ) +
        labs(
          title = paste(variable, "-", wtd_level), # Título para cada WTD
          y = variable,
          color = "Water Table Depth"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          # Remova ou comente a linha abaixo para mostrar os rótulos de clima
          # strip.text = element_blank(), 
          legend.position = "none"
        )
      
      if (!is.null(y_limits)) {
        p <- p + ylim(y_limits)
      }
      
      return(p)
    })
    
    return(plot_list)
  }
}
# Definir grupos de variáveis
variable_groups <- list(
  biogeochemical = biogeochemical_vars,
  soil_water_content = soil_water_content,
  soil_water_potential = soil_water_potential,
  transpiration_layers = transpiration_layers,
  water_fluxes = water_flux_vars
)

# Mapeamento de arquivo de cada grupo de variáveis
variable_file_map <- c(
  biogeochemical = "sumstats",
  soil_water_content = "water_balance",
  soil_water_potential = "water_balance",
  transpiration_layers = "water_balance",
  water_fluxes = "water_balance"
)

# Agora, o restante do código pode ser executado normalmente

# ---- 5. Main function ----
plot_results <- function(plot_all = TRUE, selected_variable = NULL, by_WTD = FALSE) {
  for (group in names(variable_groups)) {
    vars <- variable_groups[[group]]
    file_type <- variable_file_map[[group]]
    
    message("\n==> Processing group: ", group)
    
    if (plot_all) {
      for (var in vars) {
        message("Plotting: ", var)
        plots <- plot_variable_facet(var, climate_paths, file_type, by_WTD)
        if (by_WTD) { # Se by_WTD for TRUE, então 'plots' é uma lista de gráficos
          print(wrap_plots(plots, ncol = 1, axis = "tblr") + plot_layout(guides = "collect"))
        } else { # Caso contrário, 'plots' é um único objeto ggplot
          print(plots)
        }
      }
    } else if (!is.null(selected_variable) && selected_variable %in% vars) {
      message("Plotting selected variable: ", selected_variable)
      plots <- plot_variable_facet(selected_variable, climate_paths, file_type, by_WTD)
      if (by_WTD) { # Se by_WTD for TRUE, então 'plots' é uma lista de gráficos
        print(wrap_plots(plots, ncol = 1, axis = "tblr") + plot_layout(guides = "collect"))
      } else { # Caso contrário, 'plots' é um único objeto ggplot
        print(plots)
      }
    }
  }
}


# ---- 6. Example: call plotting ----

# Define o nome do arquivo PDF de saída
output_pdf_file <- "TROLL_simulation_plots.pdf"

# Inicia a gravação no arquivo PDF
pdf(output_pdf_file, width = 12, height = 8) # Ajuste a largura (width) e altura (height) conforme necessário

# Gera os plots
plot_results(plot_all = FALSE, selected_variable = "npp", by_WTD = FALSE)
plot_results(plot_all = FALSE, selected_variable = "litterfall", by_WTD = FALSE)

plot_results(plot_all = FALSE, selected_variable = "npp")

plot_results(plot_all = FALSE, selected_variable = "gpp")

plot_results(plot_all = FALSE, selected_variable = "agb")

plot_results(plot_all = FALSE, selected_variable = "sum1")

plot_results(plot_all = FALSE, selected_variable = "sum10")

plot_results(plot_all = FALSE, selected_variable = "sum30")

plot_results(plot_all = FALSE, selected_variable = "ba")

plot_results(plot_all = FALSE, selected_variable = "ba10")

plot_results(plot_all = FALSE, selected_variable = "litterfall")

plot_results(plot_all = FALSE, selected_variable = "SWC_0", by_WTD = FALSE)

plot_results(plot_all = FALSE, selected_variable = "SWC_1", by_WTD = FALSE)

plot_results(plot_all = FALSE, selected_variable = "SWC_2", by_WTD = FALSE)

plot_results(plot_all = FALSE, selected_variable = "SWC_3", by_WTD = FALSE)

plot_results(plot_all = FALSE, selected_variable = "SWC_4", by_WTD = FALSE)

plot_results(plot_all = FALSE, selected_variable = "transpiration_0", by_WTD = FALSE)

plot_results(plot_all = FALSE, selected_variable = "transpiration_1", by_WTD = FALSE)

plot_results(plot_all = FALSE, selected_variable = "transpiration_2", by_WTD = FALSE)

plot_results(plot_all = FALSE, selected_variable = "transpiration_3", by_WTD = FALSE)


# Fecha o dispositivo gráfico e salva o PDF
dev.off()

cat("Todos os gráficos foram salvos em:", output_pdf_file, "\n")


# # -------------------------------------------
# # Script to plot results from TROLL simulations with fixed WTD values
# # Comparing 4 climate conditions: regular climate, 30% red.prec. (reduced precipitation),
# # 50% red. prec. and 70% red. prec.
# # Author: Bianca Fazio Rius
# # -------------------------------------------
# 
# library(ggplot2)
# library(gridExtra)
# library(dplyr)
# library(patchwork)
# 
# # # ---- 1. Define paths to simulation outputs ----
# # climate_paths = list(
# #   regular = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/",
# #   redprec30 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_30/",
# #   redprec50 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_50/",
# #   redprec70 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_70/"
# # )
# 
# climate_paths = list(
#   regular = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/regular_climate/",
#   redprec30 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_30/",
#   redprec50 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_50/",
#   redprec70 = "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/prec_red_70/"
# )
# 
# 
# # ---- 2. Define variable groups and types ----
# biogeochemical_vars <- c("npp", "gpp", "agb", "sum1", "sum10", "sum30", "ba", "ba10", "litterfall")
# soil_water_content  <- c("SWC_0", "SWC_1", "SWC_2", "SWC_3", "SWC_4")
# soil_water_potential <- c("SWP_0", "SWP_1", "SWP_2", "SWP_3", "SWP_4")
# transpiration_layers <- c("transpitation_0", "transpitation_1", "transpitation_2", "transpitation_3", "transpitation_4")
# water_flux_vars     <- c("precipitation", "interception", "throughfall", "runoff", "leak", "evaporation")
# 
# # Variable group to metadata mapping
# variable_types <- c("biogeochemical", "soil water content", "soil water potential", "transpiration layers", "water fluxes")
# 
# get_variable_info <- function(variable_type) {
#   switch(variable_type,
#          "biogeochemical"       = list(vars = biogeochemical_vars, file = "sumstats"),
#          "soil water content"   = list(vars = soil_water_content, file = "water_balance"),
#          "soil water potential" = list(vars = soil_water_potential, file = "water_balance"),
#          "transpiration layers" = list(vars = transpiration_layers, file = "water_balance"),
#          "water fluxes"         = list(vars = water_flux_vars, file = "water_balance"),
#          stop("Unknown variable type.")
#   )
# }
# 
# # ---- 3. Function to read data for a given WTD condition ----
# read_simulation_data <- function(path, wtd, file_type, climate) {
#   file_path <- file.path(path, wtd, paste0("(null)_0_", file_type, ".txt"))
#   df <- read.table(file_path, header = TRUE)
#   df$WTD <- wtd
#   df$climate <- climate
#   return(df)
# }
# 
# 
# # ---- 4. Function to create plots ----
# plot_variable_facet <- function(variable, paths, file_type, by_WTD = FALSE) {
#   climate_keys <- names(paths)
#   climate_labels <- c("Reg. Clim.", "Red. Prec. 30%", "Red. Prec. 50%", "Red. Prec. 70%")
#   wtd_levels <- c("no_WTD", "shallow_WTD", "deep_WTD")
#   wtd_labels <- c("No WTD", "Shallow WTD", "Deep WTD")
#   
#   # Definir cores SEMPRE pela WTD
#   color_var <- "WTD"
#   color_values <- c(
#     "No WTD" = "darkblue",
#     "Shallow WTD" = "darkgreen",
#     "Deep WTD" = "darkorange"
#   )
#   
#   # Ler todos os dados
#   all_data <- lapply(seq_along(climate_keys), function(i) {
#     climate <- climate_keys[i]
#     path <- paths[[climate]]
#     
#     bind_rows(
#       read_simulation_data(path, "no_WTD", file_type, climate),
#       read_simulation_data(path, "shallow_WTD", file_type, climate),
#       read_simulation_data(path, "deep_WTD", file_type, climate)
#     )
#   }) %>% bind_rows()
#   
#   # Formatar fatores
#   all_data <- all_data %>%
#     mutate(
#       Climate = factor(climate, levels = climate_keys, labels = climate_labels),
#       WTD = factor(WTD, levels = wtd_levels, labels = wtd_labels)
#     )
#   
#   # Ajustar limites de y para SWC_4
#   if (variable %in% soil_water_content) {
#     swc4_vals <- all_data$SWC_4
#     y_limits <- range(swc4_vals, na.rm = TRUE)
#   } else {
#     y_limits <- NULL
#   }
#   
#   if (!by_WTD) {
#     # Modo original: comparando WTD dentro de cada clima
#     p <- ggplot(all_data, aes(x = iter, y = .data[[variable]], color = .data[[color_var]])) +
#       geom_line(size = 0.7, alpha = 0.5) +
#       facet_wrap(~ Climate, nrow = 1) +
#       scale_color_manual(values = color_values) +
#       scale_x_continuous(
#         name = "Year",
#         breaks = seq(0, 10000, by = 365 * 3),
#         labels = function(x) floor(x / 365) + 1
#       ) +
#       labs(y = variable, color = "Water Table Depth") +
#       theme_minimal() +
#       theme(
#         strip.text = element_blank(),
#         legend.position = 'none'
#       )
#     
#     if (!is.null(y_limits)) {
#       p <- p + ylim(y_limits)
#     }
#     
#     return(p)
#     
#   } else {
#     # Novo modo: um gráfico para cada WTD
#     plot_list <- lapply(seq_along(levels(all_data$WTD)), function(i) {
#       wtd_level <- levels(all_data$WTD)[i]
#       subset_data <- all_data %>% filter(WTD == wtd_level)
#       
#       p <- ggplot(subset_data, aes(x = iter, y = .data[[variable]], color = .data[[color_var]])) +
#         geom_line(size = 0.7, alpha = 0.7) +
#         facet_wrap(~ Climate, nrow = 1) +
#         scale_color_manual(values = color_values) +
#         scale_x_continuous(
#           name = "Year",
#           breaks = seq(0, 10000, by = 365 * 3),
#           labels = function(x) floor(x / 365) + 1
#         ) +
#         labs(
#           title = if (i == 1) paste(variable, "-", wtd_level) else NULL,  # Título só para o primeiro
#           y = variable,
#           color = "Water Table Depth"
#         ) +
#         theme_minimal() +
#         theme(
#           plot.title = element_text(hjust = 0.5, face = "bold"),
#           strip.text = element_blank(),
#           legend.position = "none"
#         )
#       
#       if (!is.null(y_limits)) {
#         p <- p + ylim(y_limits)
#       }
#       
#       return(p)
#     })
#     
#     return(plot_list)
#   }
# }
# 
# # Definir grupos de variáveis
# variable_groups <- list(
#   biogeochemical = biogeochemical_vars,
#   soil_water_content = soil_water_content,
#   soil_water_potential = soil_water_potential,
#   transpiration_layers = transpiration_layers,
#   water_fluxes = water_flux_vars
# )
# 
# # Mapeamento de arquivo de cada grupo de variáveis
# variable_file_map <- c(
#   biogeochemical = "sumstats",
#   soil_water_content = "water_balance",
#   soil_water_potential = "water_balance",
#   transpiration_layers = "water_balance",
#   water_fluxes = "water_balance"
# )
# 
# # Agora, o restante do código pode ser executado normalmente
# 
# # ---- 5. Main function ----
# plot_results <- function(plot_all = TRUE, selected_variable = NULL, by_WTD = FALSE) {
#   for (group in names(variable_groups)) {
#     vars <- variable_groups[[group]]
#     file_type <- variable_file_map[[group]]
#     
#     message("\n==> Processing group: ", group)
#     
#     if (plot_all) {
#       for (var in vars) {
#         message("Plotting: ", var)
#         plots <- plot_variable_facet(var, climate_paths, file_type, by_WTD)
#         if (by_WTD) { # Se by_WTD for TRUE, então 'plots' é uma lista de gráficos
#           print(wrap_plots(plots, ncol = 1, axis = "tblr") + plot_layout(guides = "collect"))
#         } else { # Caso contrário, 'plots' é um único objeto ggplot
#           print(plots)
#         }
#       }
#     } else if (!is.null(selected_variable) && selected_variable %in% vars) {
#       message("Plotting selected variable: ", selected_variable)
#       plots <- plot_variable_facet(selected_variable, climate_paths, file_type, by_WTD)
#       if (by_WTD) { # Se by_WTD for TRUE, então 'plots' é uma lista de gráficos
#         print(wrap_plots(plots, ncol = 1, axis = "tblr") + plot_layout(guides = "collect"))
#       } else { # Caso contrário, 'plots' é um único objeto ggplot
#         print(plots)
#       }
#     }
#   }
# }
# 
# 
# # ---- 6. Example: call plotting ----
# 
# plot_results(plot_all = FALSE, selected_variable = "npp", by_WTD = FALSE)
# plot_results(plot_all = FALSE, selected_variable = "litterfall", by_WTD = TRUE)
# 
# plot_results(plot_all = FALSE, selected_variable = "npp")
# 
# plot_results(plot_all = FALSE, selected_variable = "gpp")
# 
# plot_results(plot_all = FALSE, selected_variable = "agb")
# 
# plot_results(plot_all = FALSE, selected_variable = "sum1")
# 
# plot_results(plot_all = FALSE, selected_variable = "sum10")
# 
# plot_results(plot_all = FALSE, selected_variable = "sum30")
# 
# plot_results(plot_all = FALSE, selected_variable = "ba")
# 
# plot_results(plot_all = FALSE, selected_variable = "ba10")
# 
# plot_results(plot_all = FALSE, selected_variable = "litterfall")
# 
# plot_results(plot_all = FALSE, selected_variable = "SWC_0", by_WTD = TRUE)
# 
# plot_results(plot_all = FALSE, selected_variable = "SWC_1", by_WTD = FALSE)
# 
# plot_results(plot_all = FALSE, selected_variable = "SWC_2", by_WTD = FALSE)
# 
# plot_results(plot_all = FALSE, selected_variable = "SWC_3", by_WTD = FALSE)
# 
# plot_results(plot_all = FALSE, selected_variable = "SWC_4", by_WTD = FALSE)
# 
# plot_results(plot_all = FALSE, selected_variable = "transpiration_0", by_WTD = FALSE)
# 
# plot_results(plot_all = FALSE, selected_variable = "transpiration_1", by_WTD = FALSE)
# 
# plot_results(plot_all = FALSE, selected_variable = "transpiration_2", by_WTD = FALSE)
# 
# plot_results(plot_all = FALSE, selected_variable = "transpiration_3", by_WTD = FALSE)
# 
# 
# 
# 
