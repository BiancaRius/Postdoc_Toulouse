#code to plot results from tests regarding WTD
library(ggplot2)
library(gridExtra)
library(dplyr)


# Define the plotting function
plot_time_series <- function(data, time_col, value_col, line_color = "blue") {
  ggplot(data, aes_string(x = time_col, y = value_col)) +
    geom_line(color = line_color, size = 1) +
    theme_minimal() +
    labs(title = paste(value_col, "Over Time"),
         x = "Iteration",
         y = value_col) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Function to read data and generate plots
plot_variables <- function(file_path, variables, ncol, nrow, time_col = "iter", line_color = "darkblue") {
  data <- read.table(file_path, header = TRUE)
  plots <- lapply(variables, function(col) plot_time_series(data, time_col, col, line_color))
  grid.arrange(grobs = plots, ncol = ncol, nrow = nrow)
}

# File paths
biogem_file <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/transfer_8110311_files_8a5bba64/(null)_0_sumstats.txt"
water_file <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/transfer_8110311_files_8a5bba64/(null)_0_water_balance.txt"

# Biogeochemical variables
plot_variables(biogem_file, 
               c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                 'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall'),
               ncol = 3, nrow = 5)

# Water balance variables
plot_variables(water_file, 
               c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4'),
               ncol = 2, nrow = 3)

#########
biogem_file_nowtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/no_WTD/(null)_0_sumstats.txt"
water_file_nowtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/no_WTD/(null)_0_water_balance.txt"

biogem_file_shallowwtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/shallow_WTD/(null)_0_sumstats.txt"
water_file_shallowwtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/shallow_WTD/(null)_0_water_balance.txt"

biogem_file_deepwtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/deep_WTD/(null)_0_sumstats.txt"
water_file_deepwtd <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/deep_WTD/(null)_0_water_balance.txt"


# Biogeochemical variables
plot_variables(biogem_file_shallowwtd, 
               c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                 'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall'),
               ncol = 3, nrow = 5)

# Water balance variables
plot_variables(water_file_shallowwtd, 
               c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4'),
               ncol = 2, nrow = 3)

plot_variables(biogem_file_deepwtd, 
               c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                 'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall'),
               ncol = 3, nrow = 5)

# Water balance variables
plot_variables(water_file_deepwtd, 
               c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4'),
               ncol = 2, nrow = 3)

plot_variables(biogem_file_nowtd, 
               c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                 'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall'),
               ncol = 3, nrow = 5)

# Water balance variables
plot_variables(water_file_nowtd, 
               c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4'),
               ncol = 2, nrow = 3)







library(ggplot2)
library(gridExtra)

# Função para criar gráficos com cores diferenciadas
plot_time_series <- function(data, time_col, value_col, line_color, title, x_label) {
  ggplot(data, aes_string(x = time_col, y = value_col)) +
    geom_line(color = line_color, alpha = 0.6, size = 1) +  # Ajuste de transparência
    ggtitle(title) +
    xlab(x_label) +
    ylab(value_col) +
    theme_minimal()
}

# Lendo os dados
table_shallow <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/shallow_WTD/(null)_0_water_balance.txt", header = TRUE)
table_deep <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/deep_WTD/(null)_0_water_balance.txt", header = TRUE)

columns_to_plot <- c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4')

plots <- list()

# Criando gráficos para cada coluna
for (col in columns_to_plot) {
  combined_plot <- ggplot() +
    geom_line(data = table_shallow, aes_string(x = "iter", y = col), color = "blue", alpha = 0.6, size = 1) +  # Azul para shallow WTD
    geom_line(data = table_deep, aes_string(x = "iter", y = col), color = "red", alpha = 0.6, size = 1) +  # Vermelho para deep WTD
    ggtitle(paste(col, "- Shallow (Blue) vs Deep (Red) WTD")) +
    xlab("Iteration") +
    ylab(col) +
    theme_minimal()
  
  plots[[col]] <- combined_plot
}

# Arranjando os gráficos em uma grade
grid.arrange(grobs = plots, ncol = 2, nrow = 3)


library(ggplot2)
library(gridExtra)
library(dplyr)

# Função para criar gráficos combinados
plot_time_series <- function(shallow, deep, nowtd, time_col, value_col) {
  ggplot() +
    geom_line(data = shallow, aes_string(x = time_col, y = value_col), color = "green", alpha = 0.6, size = 1) +  # Verde para shallow WTD
    geom_line(data = deep, aes_string(x = time_col, y = value_col), color = "red", alpha = 0.6, size = 1) +  # Vermelho para deep WTD
    geom_line(data = nowtd, aes_string(x = time_col, y = value_col), color = "blue", alpha = 0.6, size = 1) +  # Azul para no WTD
    ggtitle(paste(value_col, "- No WTD (Blue), Shallow (Green), Deep (Red)")) +
    xlab("Iteration") +
    ylab(value_col) +
    theme_minimal()
}

# Lendo os dados
table_shallow <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/shallow_WTD/(null)_0_water_balance.txt", header = TRUE)
table_deep <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/deep_WTD/(null)_0_water_balance.txt", header = TRUE)
table_nowtd <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/(null)_0_water_balance.txt", header = TRUE)

# Filtrar as 5000 primeiras observações para no WTD
table_nowtd <- table_nowtd %>% filter(iter <= 5000)

columns_to_plot <- c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4')

plots <- list()

# Criando gráficos para cada coluna
for (col in columns_to_plot) {
  combined_plot <- plot_time_series(table_shallow, table_deep, table_nowtd, "iter", col)
  plots[[col]] <- combined_plot
}

# Arranjando os gráficos em uma grade
grid.arrange(grobs = plots, ncol = 2, nrow = 3)

library(ggplot2)
library(gridExtra)
library(dplyr)

# Função para criar gráficos combinados
plot_time_series <- function(shallow, nowtd, time_col, value_col) {
  ggplot() +
    geom_line(data = shallow, aes_string(x = time_col, y = value_col), color = "black", alpha = 0.8, size = 1) +  # Verde para shallow WTD
    geom_line(data = nowtd, aes_string(x = time_col, y = value_col), color = "blue", alpha = 0.5, size = 1) +  # Azul para no WTD
    ggtitle(paste(value_col, "- No WTD (Blue) & Shallow WTD (Green)")) +
    xlab("Iteration") +
    ylab(value_col) +
    theme_minimal()
}

# Lendo os dados
table_shallow <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/shallow_WTD/(null)_0_water_balance.txt", header = TRUE)
table_nowtd <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/(null)_0_water_balance.txt", header = TRUE)

# Filtrar as 5000 primeiras observações para no WTD
table_nowtd <- table_nowtd %>% filter(iter <= 5000)

columns_to_plot <- c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4')

plots <- list()

# Criando gráficos para cada coluna
for (col in columns_to_plot) {
  combined_plot <- plot_time_series(table_shallow, table_nowtd, "iter", col)
  plots[[col]] <- combined_plot
}

# Arranjando os gráficos em uma grade
grid.arrange(grobs = plots, ncol = 2, nrow = 3)

library(ggplot2)
library(gridExtra)
library(dplyr)

plot_time_series <- function(shallow, deep, nowtd, time_col, value_col) {
  ggplot() +
    geom_line(data = shallow, aes_string(x = time_col, y = value_col), color = "green", alpha = 0.6, size = 1) +  # Verde para shallow WTD
    geom_line(data = deep, aes_string(x = time_col, y = value_col), color = "red", alpha = 0.6, size = 1) +  # Vermelho para deep WTD
    geom_line(data = nowtd, aes_string(x = time_col, y = value_col), color = "blue", alpha = 0.6, size = 1) +  # Azul para no WTD
    ggtitle(paste(value_col, "- No WTD (Blue), Shallow (Green), Deep (Red)")) +
    xlab("Iteration") +
    ylab(value_col) +
    theme_minimal()
}

# Lendo os dados
table_shallow <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/shallow_WTD/(null)_0_sumstats.txt", header = TRUE)
table_deep <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/deep_WTD/(null)_0_sumstats.txt", header = TRUE)
table_nowtd <- read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/(null)_0_sumstats.txt", header = TRUE)

# Filtrar as 5000 primeiras observações para no WTD
table_nowtd <- table_nowtd %>% filter(iter <= 5000)

# Variáveis para plotar
columns_to_plot <- c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                     'gpp', 'rday', 'rnight', 'rstem', 'litterfall')

plots <- list()

for (col in columns_to_plot) {
  combined_plot <- plot_time_series(table_shallow, table_nowtd,table_deep, "iter", col)
  plots[[col]] <- combined_plot
}

grid.arrange(grobs = plots, ncol = 3, nrow = 5)

### current run
table_Paracou_current = read.table("/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/(null)_0_water_balance.txt", header = TRUE) # nolint


columns_to_plot = c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4')

plots = list()

# Loop through the column names and create plots
for (col in columns_to_plot) {
  plot <- plot_time_series(table_Paracou_current, time_col = "iter", value_col = col,
                           line_color = "darkblue",
                           title = paste(col, "Over Time"),
                           x_label = "Iteration")
  plots[[col]] <- plot  # Store the plot in the list
}

# Arrange the plots side by side
grid.arrange(grobs = plots, ncol = 2, nrow = 3)



library(ggplot2)
library(gridExtra)
library(dplyr)

# Function to read data
read_data <- function(base_path, condition, file_type) {
  file_path <- file.path(base_path, condition, paste0("(null)_0_", file_type, ".txt"))
  read.table(file_path, header = TRUE)
}

plot_time_series <- function(shallow, deep, nowtd, time_col, value_col) {
  ggplot() +
    geom_line(data = nowtd, aes_string(x = time_col, y = value_col), 
              color = "darkblue", linetype = "solid", size = 0.7, alpha = 0.9) +  # Sólida
    geom_line(data = shallow, aes_string(x = time_col, y = value_col), 
              color = "darkgreen", linetype = "solid", size = 0.7,alpha = 0.6) +  # Tracejada
    geom_line(data = deep, aes_string(x = time_col, y = value_col), 
              color = "darkorange", linetype = "solid", size = 0.7,alpha = 0.4) +  # Pontilhada
    
    ggtitle(paste(value_col, "- No WTD (Blue, Solid), Shallow (Green, Dashed), Deep (Red, Dotted)")) +
    xlab("Iteration") +
    ylab(value_col) +
    theme_minimal()
}

# Define paths
base_path <- "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs"

# Variable lists
biogeochem_vars <- c('npp', 'sum1', 'sum10', 'sum30', 'ba', 'ba10', 'agb', 
                     'gpp', 'npp', 'rday', 'rnight', 'rstem', 'litterfall')

water_soil_content <- c('SWC_0', 'SWC_1', 'SWC_2', 'SWC_3', 'SWC_4')

water_soil_potential <- c('SWP_0', 'SWP_1', 'SWP_2', 'SWP_3', 'SWP_4')

transpiration_layer <- c('transpitation_0', 'transpitation_1', 'transpitation_2', 'transpitation_3', 'transpitation_4')

waterflux_vars <- c('precipitation', 'interception', 'throughfall',	'runoff',	'leak',	'evaporation')


# Function to plot selected variables
plot_selected_variables <- function(variable_type) {
  if (variable_type == "biogeochemical") {
    vars <- biogeochem_vars
    file_type <- "sumstats"
  } else if (variable_type == "water soil content") {
    vars <- water_soil_content
    file_type <- "water_balance"
  } else if (variable_type == "water soil potential") {
    vars <- water_soil_potential
    file_type <- "water_balance"
  } else if (variable_type == "transpiration layer") {
    vars <- transpiration_layer
    file_type <- "water_balance"
  } else if (variable_type == "water fluxes") {
    vars <- waterflux_vars
    file_type <- "water_balance"
  } else {
    stop("Invalid variable type. Choose 'biogeochemical' or 'water'.")
  }
  
  # Read data
  shallow <- read_data(base_path, "shallow_WTD", file_type)
  deep <- read_data(base_path, "deep_WTD", file_type)
  nowtd <- read_data(base_path, "no_WTD", file_type)
  
  # Generate and arrange plots
  plots <- lapply(vars, function(var) plot_time_series(shallow, deep, nowtd, "iter", var))
  grid.arrange(grobs = plots, ncol = 3, nrow = ceiling(length(vars) / 3))
}

# Example: Call the function with "biogeochemical" or "water"
plot_selected_variables("biogeochemical")  
plot_selected_variables("water soil content")         
plot_selected_variables("water soil potential")        
plot_selected_variables("transpiration layer")         
plot_selected_variables("water fluxes")        


