#######################################
# Main simulation results from "Joint Spectral Clustering in
# Multilayer Degree Corrected Blockmodles"
#######################################

# Load all methods for simulations
source("./Experiments/run_all_methods.R")

# Note: the  code excludes the method graph-tool by default.
# To run graph-tool, install the Python package and uncomment
# the corresponding lines in "R/run_all_methods.R"

#######################################
# Simulation settings
#######################################

num_replications <- 100
delta <- list(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
parameters_list <- delta
param_iter = parameters_list

results_simulationC2 <- iterate_parameters(sim_setting = simulation_identifiability_delta, parameters_list, param_iter, num_replications)

resume <- results_simulationC2 %>%
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
  "simulationC2.txt",
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)

save(results_simulationC2, file = "./Added-Simulations/Results-addC2-rep100-miscerror.RData")


#######################################
# Plot simulation results
#######################################
load("./Added-Simu-Eva/Results-addC2-rep100-miscerror.RData")
png("./Added-Simu-Eva/Figures/Simulation-C2-rep100-flipped-modified.png", width = 1200, height = 1400, res = 200)
make_ggplot_single(results_simulationC2, "Proportion of A-type collapse layers", xbreaks = c(1, seq(0, 1, 0.1)),
                       methodnames = c("DC-MASE", "Sum of adj. matrices", "Bias-adjusted SoS",  "MASE", "OLMF", "graph-tool"))
dev.off()