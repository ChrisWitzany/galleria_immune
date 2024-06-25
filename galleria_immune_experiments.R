#----------------------------
# GALLERIA IMMUNE EXPERIMENTS
#----------------------------

set.seed(42) # for reproducibility

# setup:
source("galleria_immune_params.R")
source("galleria_immune_model.R")
source("galleria_immune_simulations.R")


# figure making:

# inocula tested and time range:
inocs <- c(10**c(1:4), 10**seq(4.5,4.8, 0.1), 10**c(6:8))
times <- seq(0, 24, length = 240) 

# figure 1 - example dynamics for different inocula
png("example_dynamics.png", width = 18, height = 18, unit = "cm", res = 300)
fig_1_data_for_brandon <- int_inoc_layers(inocs) # does this function take times?
dev.off()

# sanity check 
library(ggplot2)
ggplot(fig_1_data_for_brandon, aes(x = time, y = log10(P+1), group = inoculum))+
  geom_line()

# lets loose A for brandon
fig_1_data_for_brandon$A <- NULL 
write.csv(fig_1_data_for_brandon, "fig_1_data_for_brandon.csv")

# figure 2 - endstate comparison for different inocula
png("example_inocs_vs_endstate.png", width = 12, height = 12, unit = "cm", res = 300)
fig_2_data_for_brandon<- plot_inoc_vs_endstate(inocs, model = integrated_handling, params = params, y_0 = y, solver_method = "lsoda")
dev.off()

fig_2_data_for_brandon$time <- NULL
write.csv(fig_2_data_for_brandon, "fig_2_data_for_brandon.csv")

#-------
# make the code below prettier!
# why is there an inocs and an inoculum column in the sensitivity output data?
# make psi(A) bug save by simply setting A to zero at some point? 
#-------

# figure 3 - sensitivity analysis - this one is a bit more complicated!

# set separte parameters - bigger inoc range and more time!
inocs <- floor(10**seq(1, log10(params[["K_U"]]), 0.5))
times <- seq(0, 72, length = 201) 

# this will take a while! I run this overnight - set n_samples lower for test runs
#saving_df = sensitivity_analysis(n_samples = 10000, inocs = inocs, tspan = times, model = simple_seperate_handling, solver_method = "lsoda")
#saveRDS(saving_df, file = "sensitivity_results.rds") # save thsi data!

# this will take a while! I run this overnight - set n_samples lower for test runs
saving_df = sensitivity_analysis(n_samples = 10000, inocs = inocs, tspan = times, model = integrated_handling, solver_method = "lsoda")
saveRDS(saving_df, file = "sensitivity_results_integrated_handling.rds") # save thsi data!

# it is possible that some simulations terminate early because max step number is reached
summary(saving_df$time == max(times)) # how often does that happen?

# drop all the simulations, i.e., all inocula values where that happened:
# count NAs by name-surname combos (na.action arg is important!)
result_agg <- aggregate(time ~ n_sample, data=saving_df, FUN = function(x) any(x != max(times)), na.action = NULL)
names(result_agg)[2] <- "maxstep_timeout" # rename to count of maxstep timeout 
saving_df <- merge(saving_df, result_agg, by=c("n_sample")) # add count of NAs back to original data

# only plot the ones where that doesnt happen:
succ_sims = subset(saving_df, saving_df$maxstep_timeout == FALSE)
failed_sims = subset(saving_df, saving_df$maxstep_timeout == TRUE)

# making plot pretty
font_import() # this will take a while! 
fonts()

library(scales)

# Manually creating labels for x and y axis
x_breaks <- log10(seq(10, 100000000, by = 10^1))  # adjust by as needed
x_labels <- sapply(10^x_breaks, function(x) bquote(10^.(log10(x))))

y_breaks <- seq(1,8,1) #log10(seq(10, 100000000, by = 10^1))  # adjust by as needed
y_labels <- sapply(10^y_breaks, function(x) bquote(10^.(log10(x))))

# plot this!
sens_visu <- ggplot(succ_sims, aes(x = log10(inoc), y = log10(end), group = n_sample))+
  geom_hline(aes(yintercept = log10(params[["K_U"]])), col = "grey", linetype = "solid", size = 2)+
  geom_line(alpha = 0.01)+
  geom_point(alpha = 0.01)+
  ylab("Population size at endstate")+
  xlab("Initial inoculum")+
  # setting "nicer" axis labels 
  scale_y_continuous(breaks = seq(1, 8, 1), 
                     labels = math_format(10^.x))+
  scale_x_continuous(breaks = seq(1, 8, 1), 
                     labels = math_format(10^.x))+
  annotate(geom = "text", family = "Arial", x = 7, y = log10(params[["K_U"]])+0.4, label = "K[U]", parse = T, size = 7, col = "grey")+
  #ylim(0, 8.5)+
  #xlim(1.5, 8.5)+
  theme_half_open()+
  theme(text = element_text(family = "Arial", size = 12))

sens_visu 
ggsave("sensitivity_analysis_arial.png", sens_visu, dpi = 300, width = 12, unit = "cm", height = 12)
print("hi")




# ok one of those looks interesting lets see what is happening there

look <- subset(saving_df, saving_df$inoc > 10**4)
look <- subset(look, look$inoc < 10**5)
look <- subset(look, look$end < 10**3)
View(look)

interst <- subset(saving_df, saving_df$n_sample ==89)

# plot this!
ggplot(interst, aes(x = log10(inoc), y = log10(end), group = n_sample))+
  geom_hline(aes(yintercept = log10(params[["K_U"]])), col = "grey", linetype = "solid", size = 2)+
  geom_line(alpha = 0.1)+
  geom_point(alpha = 0.1)+
  ylab("population size at endstate (log10)")+
  xlab("initial inoculum (log10)")+
  annotate(geom = "text", x = 7, y = log10(params[["K_U"]])+0.4, label = "K[U]", parse = T,size = 7, col = "grey")+
  #ylim(0, 8.5)+
  #xlim(2, 8.5)+
  theme_half_open()

interst

# well this one looks interesting!
# how do the parameters look like?
sub = tail(interst[, c(5:23)], 1)
as.list(sub)
sub["Amax"]
params = sub

library(tibble)
library(dplyr)
deframe(sub) %>%
  as.list


# turn them into params for the sims!
paras_of_interest = split(x=, f=names(params))

# these are the params of interest:
structure(list(inoculum = 0, Amax = 0, k = 1, psimax = 1, psimin = -5, 
               kappa = 1.5, MIC = 1, K_U = 1e+08, K_P = 21.6307672110914, 
               f = 8.32207217523862e-07, b = 6.36528080865859e-07, E_0 = 94340.8540647187, 
               E_tot = 1e+05, h_1 = 4.70596409477393e-05, h_2 = 1.93660867635099e-05, 
               d = 0.0720759004499978, a = 5.17776803070887e-05, g = 0.390540995741218, 
               s = 6.82334109156235e-06), row.names = 1335L, class = "data.frame")


################################################
################################################
# TO DO 
# why is there an inoculum in the readout?
# remove antibiotics (A)

