
# Note: the  code excludes the method graph-tool by default.
# To run graph-tool, install the Python package and uncomment
# the corresponding lines in "R/run_all_methods.R"

#######################################
# Simulation settings
#######################################
source("Experiments/run_all_methods.R")
library(dplyr)
num_replications <- 100
# num_layers <- list(1, 2, 3, 5, 7, 10,15,20,30, 40, 50) #list(1, 2)
num_layers <- list(1,2,3,5,7,10,15,20,30,40,50) #list(1, 2)
parameters_list <- num_layers
param_iter = parameters_list

results_simulationC1 <- iterate_parameters(sim_setting = simulation_identifiability_m, parameters_list, param_iter, num_replications)
resume <- results_simulationC1 %>%
  group_by(parameter) %>%
  summarise_at(
    vars(
      FROST_MF,
      FROST_US,
      US,
      MF,
      OLMF,
      DC_MASE,
      graph.tool,
      Sum.A,
      S.A.2.Bias.adj,
      MASE
    ),
    list(
      moyenne = ~mean(.x, na.rm = TRUE),
      ecart_type = ~sd(.x, na.rm = TRUE)
    )
  )
resume[-1] <- lapply(resume[-1], round, digits=4)

write.table(
  resume,
  "simulationC1.txt",
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)

save(results_simulationC1, file = "./Added-Simulations/Results-addC1-rep100-miscerror.RData")


#######################################
# Plot simulation results
#######################################
load("./Added-Simulations/Results-addC1-rep100-miscerror.RData")
png("./Added-Simulations/Figures/Simulation-C1-rep100-flipped-modified.png", width = 1200, height = 1400, res = 200)
make_ggplot_single(results_simulationC1, "Number of graphs", xbreaks = c(1, 2, 3, 5, 7, 10,15,20,30, 40, 50),
                       methodnames = c("FROST_MF", "FROST_US","US","MF","OLMF","DC_MASE","graph.tool","Sum.A","S.A.2.Bias.adj","MASE"))
dev.off()