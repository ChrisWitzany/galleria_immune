#------------------------------------------
# GALLERIA IMMUNE SENSITIVITY SIMULATIONS
#------------------------------------------

# generate parameters and run simulations

#------------------------------------------
# setup 
set.seed(42) # for reproducibility
require(tidyverse)
require(data.table)
require(lhs)
require(tictoc)
require(grid)
require(gridExtra)
require(matrixStats)
require(RColorBrewer)
## for first time install of SAFE-R:
#require(devtools)
#install.packages("http://www.maths.bris.ac.uk/~mazjcr/calibrater_0.51.tar.gz", repos = NULL, type = "source")
#install_github("SAFEtoolbox/SAFE-R")
require(SAFER)


#--------------------------------------
# which parameters should be sampled?
# (dont just sample them all; this would be computationally expensive)
params_to_sample = c("K_P", "f", "b", "h_1", "h_2", "d", "a", "g", "s") 


# generate and scale lhs - NOTE this can take a while!
#lhs = generate_lhs(N = 10000, params_to_sample, params)
#scaled_lhs <- scale_lhs(lhs, params_space)
#saveRDS(scaled_lhs, file = "scaled_lhs.rds") # save LHS

# read in already generated LHS
scaled_lhs <- readRDS("scaled_lhs.rds")

#--------------------------------------
# # quick visual check whether sampling looks good
# plot_list <- list()
# 
# for (i in 1:ncol(scaled_lhs)){
#   
#   p <- ggplot(data = data.frame(x = log10(scaled_lhs[, i])), aes(x = x))+
#     geom_histogram(binwidth = 0.1, fill = "blue", color = "black")+
#     ggtitle(colnames(scaled_lhs)[i])+
#     theme_minimal()
#   
#   plot_list[[i]] <- p
# 
#   }
# 
# do.call("grid.arrange", c(plot_list, ncol = ceiling(sqrt(length(plot_list))), nrow = ceiling(length(plot_list) / ceiling(sqrt(length(plot_list))))))


#--------------------------------------
# run simulations

# time frame, and inocs tested
times <- seq(0, 72, length = 201) 
inocs <- floor(10**seq(1, 7.5, 0.1)) 

# run simulations for all sampled parameters and save it in a list
all_sim_res <- c()
tic("now running simulations for sensitivity analysis...")
for(i in 1:nrow(scaled_lhs)){

  all_sim_res[i] <- get_threshold(inocs, model = integrated_handling, params = scaled_lhs[i,], y_0 = y, solver_method = "lsoda")
  
  if(i %% 1001 == 0){cat("\014")} # this clears the console after 10001 iterations - to avoid flooding the console and causing a crash
  
}
toc()

# combine simulation results with scaled_LHS
all_sim_res <- as.data.frame(all_sim_res)
sim_res_w_paras <- cbind(all_sim_res, scaled_lhs)

# save results
saveRDS(sim_res_w_paras, file = "simulation_results.rds") 




