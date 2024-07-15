#----------------------------
# GALLERIA IMMUNE EXPERIMENTS
#----------------------------


# setup:
source("galleria_immune_params.R")
source("galleria_immune_model.R")
source("galleria_immune_simulations.R")


#----------------------------------------------------------------------------------
# simulation figures

# inocula tested and time range:
inocs <- c(10**c(1:4), 10**seq(4.5,4.8, 0.1), 10**c(6:8))
times <- seq(0, 24, length = 240) 


#----------------------------------------------------------------------------------
# figure dynamics for different inocula
png("dynamics_fig.png", width = 18, height = 18, unit = "cm", res = 300)
dynamics_fig_data <- int_inoc_layers(inocs) 
dev.off()
write.csv(dynamics_fig_data, "dynamics_fig_data.csv")


#----------------------------------------------------------------------------------
# figure end state comparison for different inocula
png("inocs_vs_endstate.png", width = 12, height = 12, unit = "cm", res = 300)
inocs_vs_endstate_fig_data <- plot_inoc_vs_endstate(inocs, model = integrated_handling, params = params, y_0 = y, solver_method = "lsoda")
dev.off()
inocs_vs_endstate_fig_data$time <- NULL
write.csv(inocs_vs_endstate_fig_data, "inocs_vs_endstate_fig_data.csv")

