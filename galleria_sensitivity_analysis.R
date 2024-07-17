#--------------------------------------
# GALLERIA IMMUNE SENSITIVITY ANALYSIS
#--------------------------------------
# identify which parameters are influential

#--------------------------------------
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
# read in the simulation results
sim_res_w_paras <- readRDS("simulation_results_new.rds")

# which parameters were sampled? (see sensitivity_simulations file)
params_to_sample = c("K_P", "f", "b", "h_1", "h_2", "d", "a", "g", "s") 

# subset only found thresholds, i.e., no NAs or "termination"s or "threshold but drops after"
summary(as.factor(sim_res_w_paras$all_sim_res))
sim_res_w_paras2 <- subset(sim_res_w_paras, !is.na(sim_res_w_paras$all_sim_res) & !sim_res_w_paras$all_sim_res  == "terminated" & !sim_res_w_paras$all_sim_res  == "threshold but drops after" )
sim_res_w_paras2$all_sim_res <- as.numeric(sim_res_w_paras2$all_sim_res)

#----------------------------------------------
# set up model parameters (X) and output (Y)
X <- as.matrix(sim_res_w_paras2[,params_to_sample])
Y <- as.matrix((sim_res_w_paras2$all_sim_res))


#----------------------------------------------
# PAWN parameters, and plotting
n <- 10 # number of conditioning intervals 
N <- nrow(Y) # number of samples - rule of thumb: N/n>80 
NN <- seq(N/(N/100*2), N, by = N/(N/100*2)) # subsample size for convergence - Note: NN - floor(NN) < 10^-6 needs to be met!
Nboot <- 10000 # bootstraps for CIs
alfa <- 0.05 # significance level for CIs estimated by bootstrapping 
x_labels <- params_to_sample # for plotting
x_labels_dummy <- c(params_to_sample, "dummy")

#----------------------------------------------
# generate and plot cdf
pawn_cdf <- pawn_plot_CDF(log10(X), log10(Y), n = n, n_col = 3, y_label = "threshold", labelinput = x_labels)

pawn_ks <- pawn_plot_ks(pawn_cdf$YF, pawn_cdf$FU, pawn_cdf$FC, pawn_cdf$xc, n_col = 3, x_labels = x_labels)


#----------------------------------------------
# calc PAWN indexes with CIs, and bootstrapping

# identification of influential and non-influential inputs
# only KS_max for screening purposes (as recommended) 

# NOTE: this might take a while as CDF computations is costly
pawn_ind <- pawn_indices(X, Y, n, Nboot, dummy = TRUE)

KS_max <- pawn_ind$KS_max # KS_max has dim (Nboot, M)
KS_dummy <- pawn_ind$KS_dummy # KS_dummy has dim (Nboot, 1)

# compute means and Cis across the bootstraps
KS_max_m <- colMeans(KS_max)
KS_max_lb <- colQuantiles(KS_max, probs = alfa/2) 
KS_max_ub <- colQuantiles(KS_max, probs = 1-alfa/2)

KS_dummy_m <- mean(KS_dummy)
KS_dummy_lb <-  quantile(KS_dummy, probs = alfa/2)
KS_dummy_ub <- quantile(KS_dummy, probs = 1-alfa/2) 

# combine KS max for all inputs and for the dummy variable for plot function
KS_max_d_m <- c(KS_max_m, KS_dummy_m) 
KS_max_d_lb <- c(KS_max_lb, KS_dummy_lb)
KS_max_d_ub <- c(KS_max_ub, KS_dummy_ub)

# create plot
gg_pawn_indices <- boxplot1_dummy(mu = KS_max_d_m, lb = KS_max_d_lb, ub = KS_max_d_ub, prnam = x_labels)
gg_pawn_indices


# prettify plot (adapt original function from SAFE-R package boxplot1_dummy function)
boxplot1_dummy_touchup <- function(mu, lb = NULL, ub = NULL, prnam = NULL){
  
  mu1 <- mu[1:(length(mu)-1)]
  lb1 <- lb[1:(length(lb)-1)]
  ub1 <- ub[1:(length(ub)-1)]
  mu2 <- tail(mu,n=1) 
  lb2 <- tail(lb,n=1)
  ub2 <- tail(ub,n=1)
  
  dat <- data.frame(x = factor(prnam, levels = prnam ), mu = mu1)
  
  .pl <- ggplot(data = dat, mapping = aes(x = reorder(x,-mu1), y = mu1, color = reorder(x,-mu1)))+ # added reorder to plot by decreasing mu1
            #geom_point(size = 5) + 
            #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lb2, ymax = ub2), 
            #  fill = "grey49", alpha = 0.01)+ # semitransparent segment
            geom_col(fill = NA, size = 1.3, width = 0.8)+ # changed to barplot
            scale_color_brewer(palette = "Set1")+# nicer colors
            geom_errorbar(mapping = aes(ymin = lb1, ymax = ub1), size = 0.5, width = 0.5, color = "black")+ # draw errors on top
            xlab(NULL)+ 
            ylab("PAWN sensitivity index")+
            scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1), expand = c(0, 0)) + # expand to put bars directly on x-axis
            geom_hline(yintercept=mu2, linewidth = 0.5, color = "grey49") +
            geom_hline(yintercept=lb2, linewidth = 0.5, linetype = "dashed", color = "grey49") +
            geom_hline(yintercept=ub2, linewidth = 0.5, linetype = "dashed", color = "grey49") +
            #annotate("text", x = length(dat$x)/2, y = 0, label = "threshold for non-influential input factors",
            #size = 20/.pt, color = "grey49")+
            theme_classic()+
            theme(legend.position = "none", 
                  text = element_text(size = 12))

  return( .pl )
  
}

gg_pawn_indices <- boxplot1_dummy_touchup(mu = KS_max_d_m, lb = KS_max_d_lb, ub = KS_max_d_ub, prnam = x_labels)
gg_pawn_indices
#save this one
ggsave("pawn_sensitivity.png", gg_pawn_indices, height = 12, width = 12, unit = "cm", dpi = 300 )



#----------------------------------------------
# assess convergence using bootstrapped CIs

# NOTE: this might take a while as CDF computations is costly (~overnight)
pawn_conv <- pawn_convergence(X, Y, n, NN, Nboot, dummy = TRUE)

KS_median_c <- pawn_conv$KS_median
KS_mean_c <- pawn_conv$KS_mean
KS_max_c <- pawn_conv$KS_max
KS_dummy_c <- pawn_conv$KS_dummy

# calculate statistics across bootstrap resamples
KS_stat <- KS_max_c
KS_max_c_m <- t(sapply(KS_stat, colMeans)) 
KS_max_c_lb <-  t(sapply(KS_stat, colQuantiles, probs = alfa/2))
KS_max_c_ub <- t(sapply(KS_stat, colQuantiles, probs = 1-alfa/2)) 

KS_stat <- KS_dummy_c
KS_dummy_c_m <- sapply(KS_stat, mean)
KS_dummy_c_lb <-  sapply(KS_dummy_c, quantile, alfa/2) 
KS_dummy_c_ub <- sapply(KS_dummy_c, quantile, 1-alfa/2) 

# combine KS max for all inputs and for dummy variable to plot
KS_max_d_c_m <- unname(cbind(KS_max_c_m, KS_dummy_c_m))
KS_max_d_c_lb <- unname(cbind(KS_max_c_lb, KS_dummy_c_lb))
KS_max_d_c_ub <- unname(cbind(KS_max_c_ub, KS_dummy_c_ub))

# plot convergence results:
plot_convergence(NN, KS_max_d_c_m, KS_max_d_c_lb, KS_max_d_c_ub, xlab = "no of model executions", ylab = "KS", labels = x_labels_dummy, panel.first = grid())


