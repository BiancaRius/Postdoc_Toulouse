# File to create the different input_pedology needed for the sensitivity analysis
# All the analysis combination can be found at: 
#https://docs.google.com/spreadsheets/d/16yrdHod3HMiNctPirBy3NFTWqjIWo7QIL1rO2uqs6vo/edit?usp=sharing

# Read original file
df_pedo <- read.table("~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/transfer_8110311_files_8a5bba64/Paracou_input_pedology.txt", header = TRUE)

template_path <- "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/transfer_8110311_files_8a5bba64/Paracou_input_pedology.txt"


make_pedology_from_template <- function(template_path,
                                        new_thickness,
                                        output_path) {
  
  # Stop execution if the output file already exists
  if (file.exists(output_path)) {
    stop(
      paste(
        "ERROR: Pedology file already exists at:",
        output_path,
        "\nDelete or rename the file before running this script again."
      ),
      call. = FALSE
    )
  }
  
  # Read the pedology template file
  template <- read.table(
    template_path,
    header = TRUE,
    sep = "",
    stringsAsFactors = FALSE
  )
  
  # Number of layers in the new configuration
  n_layers <- length(new_thickness)
  
  # Create a new pedology data.frame by replicating the first row
  # This preserves all soil properties (SOC, DBD, pH, CEC, etc.)
  new_df <- template[rep(1, n_layers), ]
  
  # Replace only the layer thickness column
  new_df$layer_thickness <- new_thickness
  
  # Reset row names to avoid duplicated indices
  rownames(new_df) <- NULL
  
  return(new_df)
}


## ------------------------ ##
## shallow_variable_sandy ##

# new_thickness <- c(0.1, 0.2, 0.4, 0.7, 1.1, 2.5)
# 
# 
# df_name = "shallow_variable_sandy"
# 
# out_path <- file.path(
#   "~",
#   "Desktop/Postdoc_Toulouse/Postdoc_Toulouse",
#   "TROLL/TROLL_code_understanding_documenting/runs/sensitivity_tests",
#   df_name,
#   "Paracou_input_pedology.txt"
# )
# 
# new_df <- make_pedology_from_template(
#   template_path = template_path,
#   new_thickness = new_thickness,
#   output_path   = out_path
# )
# 
# write.table(
#   new_df,
#   out_path,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )

# Check total soil depth
# total_depth <- sum(new_thickness)
# 
# message(
#   "Scenario: shallow_variable_sandy | ",
#   "Number of layers: ", length(new_thickness),
#   " | Total soil depth: ", total_depth, " m"
# )

## ------------------------ ##
## deep_variable_sandy ##
# 
# new_thickness <- c(0.1, 0.5, 1.0, 2.4, 6.0, 10.0)
# 
# 
# df_name = "deep_variable_sandy"
# 
# out_path <- file.path(
#   "~",
#   "Desktop/Postdoc_Toulouse/Postdoc_Toulouse",
#   "TROLL/TROLL_code_understanding_documenting/runs/sensitivity_tests",
#   df_name,
#   "Paracou_input_pedology.txt"
# )
# 
# new_df <- make_pedology_from_template(
#   template_path = template_path,
#   new_thickness = new_thickness,
#   output_path   = out_path
# )
# 
# write.table(
#   new_df,
#   out_path,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )
# 
# # Check total soil depth
# total_depth <- sum(new_thickness)
# 
# message(
#   "Scenario: deep_variable_sandy | ",
#   "Number of layers: ", length(new_thickness),
#   " | Total soil depth: ", total_depth, " m"
# )

## ------------------------ ##
## shallow_thin_sandy ##
# new_thickness <- c(
#   rep(0.1, 25),  # 25 thin layers of 0.1 m
#   2.5            # last deep layer
# )

# df_name = "shallow_thin_sandy"
# 
# out_path <- file.path(
#   "~",
#   "Desktop/Postdoc_Toulouse/Postdoc_Toulouse",
#   "TROLL/TROLL_code_understanding_documenting/runs/sensitivity_tests",
#   df_name,
#   "Paracou_input_pedology.txt"
# )
# 
# new_df <- make_pedology_from_template(
#   template_path = template_path,
#   new_thickness = new_thickness,
#   output_path   = out_path
# )
# 
# write.table(
#   new_df,
#   out_path,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )
# 
# # Check total soil depth
# total_depth <- sum(new_thickness)
# 
# message(
#   "Scenario: deep_variable_sandy | ",
#   "Number of layers: ", length(new_thickness),
#   " | Total soil depth: ", total_depth, " m"
# )


## ------------------------ ##
## deep_thin_sandy ##
# new_thickness <- c(
#   rep(0.1, 100),  # 100 thin layers of 0.1 m
#   10.0            # last deep layer
# )
# 
# df_name = "deep_thin_sandy"
# 
# out_path <- file.path(
#   "~",
#   "Desktop/Postdoc_Toulouse/Postdoc_Toulouse",
#   "TROLL/TROLL_code_understanding_documenting/runs/sensitivity_tests",
#   df_name,
#   "Paracou_input_pedology.txt"
# )
# 
# new_df <- make_pedology_from_template(
#   template_path = template_path,
#   new_thickness = new_thickness,
#   output_path   = out_path
# )
# 
# write.table(
#   new_df,
#   out_path,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )
# 
# # Check total soil depth
# total_depth <- sum(new_thickness)
# 
# message(
#   "Scenario: deep_thin_sandy | ",
#   "Number of layers: ", length(new_thickness),
#   " | Total soil depth: ", total_depth, " m"
# )


## ------------------------ ##
# shallow_intermediate_sandy
# new_thickness <- c(0.1, 0.15,
#   rep(0.25, 9),  # 9 thin layers of 0.25
#   2.5            # last deep layer
# )
# 
# df_name = "shallow_intermediate_sandy"
# 
# out_path <- file.path(
#   "~",
#   "Desktop/Postdoc_Toulouse/Postdoc_Toulouse",
#   "TROLL/TROLL_code_understanding_documenting/runs/sensitivity_tests",
#   df_name,
#   "Paracou_input_pedology.txt"
# )
# 
# new_df <- make_pedology_from_template(
#   template_path = template_path,
#   new_thickness = new_thickness,
#   output_path   = out_path
# )
# 
# write.table(
#   new_df,
#   out_path,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = "\t"
# )
# 
# # Check total soil depth
# total_depth <- sum(new_thickness)
# 
# message(
#   "Scenario: ", df_name,
#   "Number of layers: ", length(new_thickness),
#   " | Total soil depth: ", total_depth, " m"
# )

## ------------------------ ##
# deep_intermediate_sandy
new_thickness <- c(0.1, 0.15,
                   rep(0.25, 39),  # 39 thin layers of 0.25
                   10.0            # last deep layer
)

df_name = "deep_intermediate_sandy"

out_path <- file.path(
  "~",
  "Desktop/Postdoc_Toulouse/Postdoc_Toulouse",
  "TROLL/TROLL_code_understanding_documenting/runs/sensitivity_tests",
  df_name,
  "Paracou_input_pedology.txt"
)

new_df <- make_pedology_from_template(
  template_path = template_path,
  new_thickness = new_thickness,
  output_path   = out_path
)

write.table(
  new_df,
  out_path,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# Check total soil depth
total_depth <- sum(new_thickness)

message(
  "Scenario: ", df_name,
  "Number of layers: ", length(new_thickness),
  " | Total soil depth: ", total_depth, " m"
)

# shallow_variable_clayey

# shallow_intermediate_clayey

# shallow_thin_clayey
# deep_variable_clayey

# deep_intermediate_clayey

# deep_thin_clayey

# Pega a segunda linha como referência (0.1 m), pois todas as colunas são iguais
ref <- df[1, ]

# Cria 25 layers de 0.1 m
layers_0.1 <- ref[rep(1, 25), ]
layers_0.1$layer_thickness <- 0.1

# Cria última camada de 2.5 m
layer_2.5 <- ref
layer_2.5$layer_thickness <- 2.5

# Junta tudo
new_df <- rbind(layers_0.1, layer_2.5)

# Salva o arquivo final
write.table(new_df, "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/Paracou_input_pedology.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")
