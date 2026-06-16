#######################################
# Main simulation results from "Joint Spectral Clustering in
# Multilayer Degree Corrected Blockmodles"
#######################################

# Load all methods for simulations
source("Experiments/run_all_methods.R")
library(dplyr)
# Note: the  code excludes the method graph-tool by default.
# To run graph-tool, install the Python package and uncomment
# the corresponding lines in "R/run_all_methods.R"

#######################################
# Simulation settings
#######################################

# num_replications <- 100
# num_layers <- list( 5, 10, 15, 20, 25,30,35, 40, 45, 50)
num_replications <- 100
num_layers <- list(1,2,3,5,7,10,15,20,30,40,50)

parameters_list <- num_layers
param_iter = parameters_list


#######################################
# Run different scenarios
#######################################

# core-periphery
results_simulation8 <- iterate_parameters(sim_setting = simulation8, parameters_list, param_iter, num_replications)

resume <- results_simulation8 %>%
  group_by(parameter) %>%
  summarise(
    across(
      .cols = c(FROST_MF, FROST_US, US, MF, OLMF, DC_MASE,graph.tool,Sum.A,S.A.2.Bias.adj,MASE),
      .fns = list(
        moyenne = ~ mean(.x, na.rm = TRUE),
        ecart_type = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )
resume[-1] <- lapply(resume[-1], round, digits = 4)
write.table(
  resume,
  file = "simulation8.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)





