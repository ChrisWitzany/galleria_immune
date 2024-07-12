# working on the sensitivity analysis

#setup 
library(tidyverse)
library(data.table)


# read in the data 
setwd("/Users/witzany/polybox/zzz_misc/galleria_immune_model/galleria_immune")
sens_data <- readRDS("sensitivity_results_integrated_handling.rds") # save thsi data!

str(sens_data)
summary((log10(sens_data$end)))
hist(log10(sens_data$end))
max(sens_data$end)

# clean sensitivity data
# drop prematurely terminated simulations 
result_agg <- aggregate(time ~ n_sample, data=sens_data, FUN = function(x) any(x != max(x)), na.action = NULL)
names(result_agg)[2] <- "maxstep_timeout" # rename to count of maxstep timeout 
sens_data <- merge(sens_data, result_agg, by=c("n_sample")) # add count of NAs back to original data
sens_data <-  subset(sens_data, sens_data$maxstep_timeout==FALSE)

# ok this does not seem as easy as I thought.
threshold <- sens_data %>%
  group_by(n_sample) %>% # for each parameter set 
  mutate(jump_over_1e6 = end > 1e7, # total pop over 10 to the power of 6
         fall_below_1e6_after_jump = if_else(jump_over_1e6, lead(end < 1e6, default = FALSE), FALSE)) %>%  # and correct this if it ever goes down again
  summarise(
    first_jump_value = ifelse(any(jump_over_1e6), first(inoc[jump_over_1e6]), NA_real_), # get first inoculum for which that happens
    falls_below_after_jump = any(fall_below_1e6_after_jump)
  ) %>%
  ungroup()

View(threshold)
summary(threshold$first_jump_value)
summary(threshold$falls_below_after_jump)


# join the information back to the original data
all_result <- sens_data  %>%
  left_join(threshold, by = "n_sample")
View(result)

first_match_result <- all_result %>%
  group_by(n_sample) %>%
  slice(1) %>%
  ungroup()

result <- first_match_result[, c("K_P", "f", "b", "h_1", "h_2", "d", "a", "g", "s", "first_jump_value")]

View(result)
length(result$first_jump_value)
# ok now I have result which should be what I want :P


hist(result$first_jump_value)

# these do not work here out of the box...
#library(sensobol) # hmmm this package is better but cannot use pre-run data :/

#library(sensitivity)
#sobol_result <- sobol(model = NULL, X1 = X, y = Y, nboot = 100, conf = 0.95)
#?sobol

# also does not work well with already generated data - lets do 
library(devtools)
#install.packages("http://www.maths.bris.ac.uk/~mazjcr/calibrater_0.51.tar.gz", repos = NULL, type = "source")
#install_github("SAFEtoolbox/SAFE-R")
library(SAFER)


result <- as.data.table(result)

X <- as.matrix(result[,.(K_P, f, b, h_1,h_2, d, a, g,s)])
Y <- as.matrix(result$first_jump_value)
x_labels <- c("K_P", "f","b", "h_1","h_2", "d", "a", "g", "s")
pawn_cdf <- pawn_plot_CDF(X, Y, n=10, n_col=3, y_label='output y', labelinput=x_labels)

YF <- pawn_cdf$YF
FU <- pawn_cdf$FU
FC <- pawn_cdf$FC
xc <- pawn_cdf$xc

# Compute and plot KS statistics for each conditioning interval:
KS <- pawn_ks(YF, FU, FC)
dev.new()
KS_all <- pawn_plot_ks(YF, FU, FC, xc, n_col=3, x_labels = x_labels)

# Compute PAWN sensitivity indices:
pawn_ind <- pawn_indices(X, Y, n=10)

KS_median <- pawn_ind$KS_median
KS_mean <- pawn_ind$KS_mean
KS_max <- pawn_ind$KS_max

# Plot results:
library(grid)
library(gridExtra)
dev.new()
p1 <- boxplot1(as.vector(KS_median), prnam = x_labels) + ylab("KS (median)") +
  ggtitle("threshold") + theme(plot.title = element_text(hjust = 0.5))
p2 <- boxplot1(as.vector(KS_mean), prnam = x_labels) + ylab("KS (mean)") +
  ggtitle("threshold") + theme(plot.title = element_text(hjust = 0.5))
p3 <- boxplot1(as.vector(KS_max), prnam = x_labels) + ylab("KS (max)") +
  ggtitle("threshold") + theme(plot.title = element_text(hjust = 0.5))
grid.arrange(grobs = list(p1, p2, p3), ncol = 3)

# Use bootstrapping to derive confidence bounds:
Nboot <- 10000

# Compute sensitivity indices for Nboot bootstrap resamples
# (Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
pawn_ind <- pawn_indices(X, Y, n=10, Nboot)
KS_median <- pawn_ind$KS_median
KS_mean <- pawn_ind$KS_mean
KS_max <- pawn_ind$KS_max

# KS_median and KS_mean and KS_max have shape (Nboot, M)
# Compute mean and confidence intervals of the sensitivity indices across the
# bootstrap resamples:
alfa <- 0.05 # Significance level for the confidence intervals estimated by bootstrapping 

require(matrixStats)
KS_stat <- KS_median
KS_median_m <- colMeans(KS_stat) # mean
KS_median_lb <-  colQuantiles(KS_stat,probs=alfa/2) # Lower bound
KS_median_ub <- colQuantiles(KS_stat,probs=1-alfa/2) # Upper bound

KS_stat <- KS_mean
KS_mean_m <- colMeans(KS_stat) # mean
KS_mean_lb <-  colQuantiles(KS_stat, probs = alfa/2) # Lower bound
KS_mean_ub <- colQuantiles(KS_stat, probs = 1 - alfa/2) # Upper bound

KS_stat <- KS_max
KS_max_m <- colMeans(KS_stat) # mean
KS_max_lb <- colQuantiles(KS_stat, probs = alfa/2) # Lower bound
KS_max_ub <- colQuantiles(KS_stat, probs = 1-alfa/2) # Upper bound

# Plot bootstrapping results:
dev.new()
p1 <- boxplot1(mu = KS_median_m, lb = KS_median_lb, ub = KS_median_ub, prnam = x_labels) + ylab("KS (median)") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- boxplot1(mu = KS_mean_m, lb = KS_mean_lb, ub = KS_mean_ub, prnam = x_labels) + ylab("KS (mean)") +
  theme(plot.title = element_text(hjust = 0.5))
p3 <- boxplot1(mu = KS_max_m, lb = KS_max_lb, ub = KS_max_ub, prnam = x_labels) + ylab("KS (max)") +
  theme(plot.title = element_text(hjust = 0.5))
grid.arrange(grobs = list(p1, p2, p3), ncol = 3)


############

# Analyze convergence of sensitivity indices:
N <- nrow(result) # Number of samples
NN <- seq(N / 5, N, by = N / 5 )
n <- 10
#( Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
pawn_conv <- pawn_convergence(X, Y, n, NN)
KS_median_c <- pawn_conv$KS_median
KS_mean_c <- pawn_conv$KS_mean
KS_max_c <- pawn_conv$KS_max

KS_median_c <- do.call("rbind",KS_median_c)
KS_mean_c <- do.call("rbind",KS_mean_c)
KS_max_c <- do.call("rbind",KS_max_c)

# Plot convergence
dev.new()
par(mfrow=c(3,1))
plot_convergence(NN, KS_median_c, xlab = "no of model executions", ylab = "KS (median)", labels = x_labels, panel.first = grid())
plot_convergence(NN, KS_mean_c, xlab = "no of model executions", ylab = "KS (mean)", labels = x_labels, panel.first = grid())
plot_convergence(NN, KS_max_c, xlab = "no of model executions", ylab = "KS (max)", labels = x_labels, panel.first = grid())

# Analyze convergence using bootstrapping to derive confidence intervals
#( Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
pawn_conv <- pawn_convergence(X, Y, n, NN, Nboot)

KS_median_c <- pawn_conv$KS_median
KS_mean_c <- pawn_conv$KS_mean
KS_max_c <- pawn_conv$KS_max

# Compute mean and confidence intervals of the sensitivity indices across the
# bootstrap resamples:
alfa <- 0.05 # Significance level for the confidence intervals estimated by bootstrapping 

KS_stat <- KS_median_c
KS_median_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_median_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_median_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

KS_stat <- KS_mean_c
KS_mean_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_mean_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_mean_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

KS_stat <- KS_max_c
KS_max_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_max_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_max_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

# Plot convergence results:
dev.new()
par(mfrow=c(3,1))
plot_convergence(NN, KS_median_c_m, KS_median_c_lb, KS_median_c_ub, xlab = "no of model executions", ylab = "KS (median)", labels = x_labels, panel.first = grid())
plot_convergence(NN, KS_mean_c_m, KS_mean_c_lb, KS_mean_c_ub, xlab = "no of model executions", ylab = "KS (mean)", labels = x_labels, panel.first = grid())
plot_convergence(NN, KS_max_c_m, KS_max_c_lb, KS_max_c_ub, xlab = "no of model executions", ylab = "KS (max)", labels = x_labels, panel.first = grid())


#### Step 6 (identification of influential and non-influential inputs) ####
# This is done by adding an artificial 'dummy' input to the list of the model inputs. 
# The sensitivity indices for the dummy parameter estimate the approximation error of the
# sensitivity indices. For reference and more details, see help of the function pawn_indices

# Sensitivity indices using bootstrapping for the model inputs and the dummy input:
# Use bootstrapping to derive confidence bounds:
#Nboot <- 100
# Compute sensitivity indices for Nboot bootstrap resamples. We analyse KS_max
# only (and not KS_median and KS_mean) for screening purposes.
# (Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
pawn_ind <- pawn_indices(X, Y, n, Nboot, dummy = TRUE)

KS_max <- pawn_ind$KS_max # KS_max has dim (Nboot, M)
KS_dummy <- pawn_ind$KS_dummy # KS_dummy has dim (Nboot, 1)

# Compute mean and confidence intervals of the sensitivity indices across the
# bootstrap resamples:
alfa <- 0.05 # Significance level for the confidence intervals estimated by bootstrapping 

KS_stat <- KS_max
KS_max_m <- colMeans(KS_stat) # mean
KS_max_lb <- colQuantiles(KS_stat,probs=alfa/2) # Lower bound
KS_max_ub <- colQuantiles(KS_stat,probs=1-alfa/2) # Upper bound

KS_stat <- KS_dummy
KS_dummy_m <- mean(KS_stat) # mean
KS_dummy_lb <-  quantile(KS_dummy,alfa/2) # Lower bound
KS_dummy_ub <- quantile(KS_dummy,1-alfa/2) # Upper bound

# Combine KS max for all inputs and for dummy to plot
KS_max_d_m <- c(KS_max_m,KS_dummy_m) 
KS_max_d_lb <- c(KS_max_lb,KS_dummy_lb)
KS_max_d_ub <- c(KS_max_ub,KS_dummy_ub)

# Plot bootstrapping results:
dev.new()
boxplot1_dummy(mu = KS_max_d_m, lb = KS_max_d_lb, ub = KS_max_d_ub, prnam = x_labels) + ylab("KS") +
  theme(plot.title = element_text(hjust = 0.5))


# saving this one
reso <- 150
length <- 3.25*reso/72
png("sensitivity_plot_1000_boot.png",units="in",res=reso,height=length,width=length)
sens_plot <- boxplot1_dummy(mu = KS_max_d_m, lb = KS_max_d_lb, ub = KS_max_d_ub, prnam = x_labels) + ylab("KS") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
dev.off()  
ggsave("gg_sens_10000_boots.png",sens_plot, dpi = 150)
# can I heat map the biggest factors?
library(cowplot)

shoddy_heatmap <- ggplot(result, aes(log10(h_1), log10(h_2), color= as.factor(round(log10(first_jump_value),1)))) + 
  geom_point(size = 5, alpha = 0.5)+
  ylab("h_2 (log10)")+
  xlab("h_1 (log10)")+
  scale_color_discrete("threshold (log10)") +
  theme_bw()+
  theme(text = element_text(size = 12))

shoddy_heatmap
ggsave("quick_heatmap.png",shoddy_heatmap, dpi = 150)
ggplot(result, aes(log10(h_1), log10(h_2), fill= as.factor(round(log10(first_jump_value),1)))) + 
 # geom_density_2d_filled(show.legend = FALSE) +
#  coord_cartesian(expand = FALSE) 
  geom_tile(width=0.1,height=0.1)
  


result <- as.data.table(result)

X <- as.matrix(result[,.(K_P, f, b, h_1,h_2, d, a, g,s)])
Y <- as.matrix(result$first_jump_value)

# Analyze convergence using bootstrapping to derive confidence intervals
#( Warning: the following line may take some time to run, as the computation of
# CDFs is costly):
NN <- seq(N / 5, N, by = N / 5 )
pawn_conv <- pawn_convergence(X, Y, n, NN, Nboot, dummy = TRUE)

KS_median_c <- pawn_conv$KS_median
KS_mean_c <- pawn_conv$KS_mean
KS_max_c <- pawn_conv$KS_max
KS_dummy_c <- pawn_conv$KS_dummy

# Calculate statistics across bootstrap resamples (mean, lower and upper bounds of sensitivity indices):
alfa <- 0.05 # Significance level for the confidence intervals estimated by bootstrapping 

KS_stat <- KS_max_c
KS_max_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_max_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_max_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

KS_stat <- KS_dummy_c
KS_dummy_c_m <- sapply(KS_stat,mean) # mean
KS_dummy_c_lb <-  sapply(KS_dummy_c,quantile,alfa/2) # Lower bound
KS_dummy_c_ub <- sapply(KS_dummy_c,quantile,1-alfa/2) # Upper bound

# Combine KS max for all inputs and for dummy to plot
KS_max_d_c_m <- unname(cbind(KS_max_c_m,KS_dummy_c_m))
KS_max_d_c_lb <- unname(cbind(KS_max_c_lb,KS_dummy_c_lb))
KS_max_d_c_ub <- unname(cbind(KS_max_c_ub,KS_dummy_c_ub))

x_labels_dummy <- c("K_P", "f","b","h_1","h_2", "d", "a", "g", "s", "dummy")

# Plot convergence results:
dev.new()
plot_convergence(NN, KS_max_d_c_m, KS_max_d_c_lb, KS_max_d_c_ub, xlab = "no of model executions", ylab = "KS", labels = x_labels_dummy, panel.first = grid())




