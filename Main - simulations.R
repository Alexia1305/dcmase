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

# Same B same theta
results_simulation1 <- iterate_parameters(sim_setting = simulation1, parameters_list, param_iter, num_replications)

resume <- results_simulation1 %>%
  group_by(parameter) %>%
  summarise(
    across(
      .cols = c(FROST_MF, FROST_US,FROST_DCMASE, US, MF, OLMF, DC_MASE,graph.tool,Sum.A,S.A.2.Bias.adj,MASE),
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
  file = "simulation1.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

browser()

# Different B same theta
results_simulation2 <- iterate_parameters(sim_setting = simulation2, parameters_list, param_iter, num_replications)
resume <- results_simulation2 %>%
  group_by(parameter) %>%
  summarise(
    across(
      .cols = c(FROST_MF, FROST_US,FROST_DCMASE, US, MF, OLMF, DC_MASE,graph.tool,Sum.A,S.A.2.Bias.adj,MASE),
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
  file = "simulation2.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# Diff B diff theta
results_simulation3 <- iterate_parameters(simulation3, parameters_list, param_iter, num_replications)
resume <- results_simulation3 %>%
  group_by(parameter) %>%
  summarise(
    across(
      .cols = c(FROST_MF, FROST_US,FROST_DCMASE, US, MF, OLMF, DC_MASE,graph.tool,Sum.A,S.A.2.Bias.adj,MASE),
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
  file = "simulation3.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
# Same B different theta
results_simulation4 <- iterate_parameters(simulation4, parameters_list, param_iter, num_replications)
resume <- results_simulation4 %>%
  group_by(parameter) %>%
  summarise(
    across(
      .cols = c(FROST_MF, FROST_US,FROST_DCMASE, US, MF, OLMF, DC_MASE,graph.tool,Sum.A,S.A.2.Bias.adj,MASE),
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
  file = "simulation4.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
# # Core-periphery B same theta
# results_simulation8 <- iterate_parameters(sim_setting = simulation8, parameters_list, param_iter, num_replications)
# 
# 
# # Core)periphery B different theta
# results_simulation9 <- iterate_parameters(sim_setting = simulation9, parameters_list, param_iter, num_replications)
# # 


parameters_list <- lapply(param_iter,
                          function(x) c(x, 150, 3))
# Same B alternating theta
results_simulation6 <- iterate_parameters(simulation6, parameters_list, param_iter, num_replications)
resume <- results_simulation6 %>%
  group_by(parameter) %>%
  summarise(
    across(
      .cols = c(FROST_MF, FROST_US,FROST_DCMASE, US, MF, OLMF, DC_MASE,graph.tool,Sum.A,S.A.2.Bias.adj,MASE),
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
  file = "simulation6.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Different B alternating theta
results_simulation7 <- iterate_parameters(simulation7, parameters_list, param_iter, num_replications)
resume <- results_simulation7 %>%
  group_by(parameter) %>%
  summarise(
    across(
      .cols = c(FROST_MF, FROST_US,FROST_DCMASE, US, MF, OLMF, DC_MASE,graph.tool,Sum.A,S.A.2.Bias.adj,MASE),
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
  file = "simulation7.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


#######################################
# Combine results
#######################################
results_simulation1$scenarioB <- "Same B"
results_simulation2$scenarioB <- "Different B"
results_simulation3$scenarioB <- "Different B"
results_simulation4$scenarioB <- "Same B"
results_simulation6$scenarioB <- "Same B"
results_simulation7$scenarioB <- "Different B"

results_simulation1$scenarioT <- "Same \u0398"
results_simulation2$scenarioT <- "Same \u0398"
results_simulation3$scenarioT <- "Different \u0398"
results_simulation4$scenarioT <- "Different \u0398"
results_simulation6$scenarioT <- "Alternating \u0398"
results_simulation7$scenarioT <- "Alternating \u0398"

different_scenarios <- rbind(results_simulation1, results_simulation2, results_simulation4, results_simulation3, 
                             results_simulation6, results_simulation7)

different_scenarios$scenarioB <- factor(different_scenarios$scenarioB,
                                       levels = c("Same B", "Different B"))
different_scenarios$scenarioT <- factor(different_scenarios$scenarioT,
                                        levels = c("Same \u0398", "Different \u0398", "Alternating \u0398"))
save(different_scenarios, file = "Results-testcomplete.RData")

source("R/make_ggplot.R")
#######################################
# Plot simulation results
#######################################
#load("Results-testcomplete.RData")
png("Simulation-rep100-6scenarios-flipped.png", width = 1200, height = 1500, res = 200)
different_scenarios[] <- lapply(different_scenarios, function(col) {
  if (is.list(col)) {
    as.numeric(unlist(col))
  } else {
    col
  }
})

p <- make_ggplot_multipleBT2(different_scenarios, "Number of graphs", xbreaks = c(1, seq(10, 50, 10)),methodnames = c("FROST_MF", "FROST_US","FROST_DCMASE", "US", "MF", "OLMF","DC_MASE","graph.tool","Sum.A","S.A.2.Bias.adj","MASE"))#, "graph-tool"))
ggsave(
  filename = "figure_paper.png",
  plot = p,              # ton objet ggplot
  width = 10.5,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)
#make_ggplot_multipleBT2(different_scenarios, "Number of graphs", xbreaks = c(1, seq(10, 50, 10)),
                        #methodnames = c("DC-MASE", "Sum of adj. matrices", "Bias-adjusted SoS",  "MASE", "OLMF"))#, "graph-tool"))

#dev.off()


