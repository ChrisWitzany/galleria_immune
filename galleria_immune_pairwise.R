#------------------------------------------------
# GALLERIA IMMUNE PAIRWISE PARAMETER VARIATION 
#------------------------------------------------
# (i.e. two parameters are varied, all others stay the same)

# setup:
source("galleria_immune_params.R")
source("galleria_immune_model.R")
source("galleria_immune_simulations.R")

require(zoo)
require(viridisLite)
require(scales)


# time frame, and inocs tested
times <- seq(0, 72, length = 201) 
inocs <- floor(10**seq(1, 7.5, 0.1)) 

# function for running pairwise simulations for para_1 (x) and para_2 (y) - where all other params are kept constant
pairwise_simulations <- function(para_1, para_2, by = 0.1, params = params){
  
  # create ranges
  x_range <- 10**seq(log10(as.numeric(params_space[para_1, 1])), 
                     log10(as.numeric(params_space[para_1, 2])), 
                     by = by)

  para_2 = "h_2"
  y_range <- 10**seq(log10(as.numeric(params_space[para_2, 1])), 
                     log10(as.numeric(params_space[para_2, 2])), 
                     by = by)
  
  # create pairwise combinations
  combinations <- expand.grid(x = x_range, y = y_range)
  names(combinations) <- c(para_1, para_2)
  
  # keep all parameters the same for these pairwise combinations
  pair_params <- data.frame(t(params))
  pair_params <- pair_params[rep(1:nrow(pair_params), each = nrow(combinations)), ]
  rownames(pair_params) <- NULL

  pair_params[,para_1] <- combinations[,para_1]
  pair_params[,para_2] <- combinations[,para_2]
  
  # run simulations for all pairs and save it in a list
  all_sim_res <- c()
  
  tic("now running simulations for sensitivity analysis...")
  for(i in 1:nrow(pair_params)){
  
    result <- get_threshold(inocs, model = integrated_handling, params = pair_params[i,], y_0 = y, solver_method = "lsoda")
  
    if(!is.na(result) & result == "terminated"){ # switch to stiff solver
      
      print("lsoda terminated - trying vode instead:")
      result <- get_threshold(inocs, model = integrated_handling, params = pair_params[i,], y_0 = y, solver_method = "vode")
      
    }
    
    all_sim_res[i] <- result
    
  if(i %% 1001 == 0){cat("\014")} # this clears the console after 10001 iterations - to avoid flooding the console and causing a crash
  
  }
  toc()

  # combine simulation results with pair_params
  sim_res <- cbind(as.data.frame(all_sim_res), pair_params)
  
  return(sim_res)
}



# function for plotting a nice looking heatmap
pairwise_heatmap <- function(sim_res, para_1, para_2){
  
  # some data prepping
  sim_res$res <- ifelse(is.na(sim_res$all_sim_res), 
                        10**7.5, 
                        sim_res$all_sim_res)
  
  sim_res$res <- round(log10(as.numeric(sim_res$res)), 2)
  
  # create nice looking plot
  gg_pairwise <- ggplot(sim_res, aes(log10(h_1), log10(h_2), fill = (res)))+ 
    geom_tile()+
    scale_x_continuous(labels = function(x) parse(text = paste0("10^", x)), expand = c(0, 0))+ 
    scale_y_continuous(labels = function(x) parse(text = paste0("10^", x)), expand = c(0, 0))+ 
    scale_fill_gradientn(colours = rev(c("blue", "dodgerblue", "skyblue", "lightblue", 
                                         "lightpink", "pink", "salmon", "red", "darkred")),
                         values = rescale(c(1, 2, 3, 4, 5, 6, 7.5)),
                         na.value = "darkgrey",
                         guide = "colorbar",
                         trans = "reverse",
                         breaks = c(1, 2, 3, 4, 5, 6, 7.5),
                         labels = parse(text = c(bquote(10^.(1)),
                                                 bquote(10^.(2)),
                                                 bquote(10^.(3)),
                                                 bquote(10^.(4)),
                                                 bquote(10^.(5)),
                                                 bquote(10^.(6)),
                                                 bquote(">"*10^.(7.5)))))+
    xlab(paste0(para_1))+
    ylab(paste0(para_2))+
    labs(fill = "threshold")+
    theme_bw()+
    theme(axis.title = element_text(size = 12),       
          axis.text = element_text(size = 12),      
          legend.title = element_text(size = 10),   
          legend.text = element_text(size = 8),
          legend.key.size = unit(1, "lines"))
  
  print(gg_pairwise)
  return(gg_pairwise)
  
}


#--------------------------------------------------------------
# run pairwise simulations (NOTE: that this takes a while!)
#sim_res <- pairwise_simulations(para_1 = "h_1", para_2 = "h_2", by = 0.1, params = params)

# save results
#saveRDS(sim_res, file = "pairwise_simulation_results.rds")


#-------------------------------------------------------------
# read in results
#sim_res <- readRDS("pairwise_simulation_results.rds")

# create pairwise plot
#gg_pairwise <- pairwise_heatmap(sim_res = sim_res, para_1 = "h_1", para_2 = "h_2")
#ggsave("h_1_h_2_pairwise.png", gg_pairwise, height = 9, width = 11, unit = "cm", dpi = 300 )




