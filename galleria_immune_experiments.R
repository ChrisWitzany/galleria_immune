#----------------------------
# GALLERIA IMMUNE EXPERIMENTS
#----------------------------

set.seed(42) # for reproducibility

# setup:
source("galleria_immune_params.R")
source("galleria_immune_model.R")
source("galleria_immune_simulations.R")


inocs <- floor(10**seq(1,log10(params[["K_U"]]), 0.1))
times <- seq(0, 48, length = 201) 
times <- seq(0, 24 ,0.001)


params["h_1"] = 0.0004 # killing rate for bacteria 
params["h_2"] = 0.0001001 # rate by which they end up in enganged
params["d"] = 0.0001 # return rate to resting state
params["a"] = 0.01 #0.001, # background activation rate
params["g"] = 0.0049 #handling time = 1/g
params["s"] = 0.000035 # per meeting rate/activation rate

int_inoc_layers(inocs) # does this funciton take times?


# integrated handling:
plot_inoc_vs_endstate(inocs = inocs, model = integrated_handling, params, y)

# seperete handling:
plot_inoc_vs_endstate(inocs, model = integrated_handling, params = params, y_0 = y, solver_method = "lsoda")



#---------------------------
# the max timeframe might actually matter quite a bit here!
saving_df = sensitivity_analysis(n_samples = 100, inocs = 10**seq(3,9, 1), tspan = seq(0, 48, by = 0.1), model = simple_seperate_handling, solver_method = "lsoda")

# plot this!
with_NAs <- ggplot(saving_df, aes(x = log10(inoc+1), y = log10(end+1), group = n_sample))+
  geom_line(alpha = 0.1)+
  geom_point(alpha = 0.1)+
  ylab("endstate")+
  xlab("inoculum")+
  #ylim(0, 10)+
  #xlim(3, 9)+
  theme_half_open()

with_NAs 












# the max timeframe might actually matter quite a bit here!
saving_df = sensitivity_analysis(n_samples = 100, inocs = 10**seq(3, 9, 1), tspan = seq(0, 24, by = 0.1), model = integrated_handling, solver_method = "lsoda")

# plot this!
with_NAs <- ggplot(saving_df, aes(x = log10(inoc+1), y = log10(end+1), group = n_sample))+
  geom_line(alpha = 0.1)+
  geom_point(alpha = 0.1)+
  ylab("endstate")+
  xlab("inoculum")+
  #ylim(0, 10)+
  #xlim(3, 9)+
  theme_half_open()

with_NAs 









# lets drop all the NA values? 

# count NAs by name-surname combos (na.action arg is important!)
result_agg <- aggregate(end ~ n_sample, data=saving_df, FUN=function(x) sum(is.na(x)), na.action=NULL)

# rename is count of NAs column
names(result_agg)[2] <- "number_of_na"

#add count of NAs back to original data
df <- merge(saving_df, result_agg, by=c("n_sample"))

# subset the original data
result <- df[df$number_of_na == 0, ]

# plot wo NAs!
wo_NAs <- ggplot(result, aes(x = log10(inoc), y = log10(end), group = n_sample))+
  geom_line(alpha = 0.1)+
  geom_point(alpha = 0.1)+
  ylab("endstate")+
  xlab("inoculum")+
  #ylim(0, 10)+
  #xlim(3, 9)+
  theme_half_open()

length(result[,1])
length(saving_df[,1])

#side by side
with_NAs + wo_NAs


# the max timeframe might actually matter quite a bit here!
saving_df = sensitivity_analysis(n_samples = 10, inocs = 10**seq(3,9, 1), tspan = seq(0, 24, length = 201), model = integrated_handling)

# plot this!
with_NAs <- ggplot(saving_df, aes(x = log10(inoc), y = log10(end+1), group = n_sample))+
  geom_line(alpha = 0.1)+
  geom_point(alpha = 0.1)+
  ylab("endstate")+
  xlab("inoculum")+
  #ylim(0, 10)+
  #xlim(3, 9)+
  theme_half_open()

with_NAs 

# lets drop all the NA values? 

# count NAs by name-surname combos (na.action arg is important!)
result_agg <- aggregate(end ~ n_sample, data=saving_df, FUN=function(x) sum(is.na(x)), na.action=NULL)

# rename is count of NAs column
names(result_agg)[2] <- "number_of_na"

#add count of NAs back to original data
df <- merge(saving_df, result_agg, by=c("n_sample"))

# subset the original data
result <- df[df$number_of_na == 0, ]

# plot wo NAs!
wo_NAs <- ggplot(result, aes(x = log10(inoc), y = log10(end), group = n_sample))+
  geom_line(alpha = 0.1)+
  geom_point(alpha = 0.1)+
  ylab("endstate")+
  xlab("inoculum")+
  #ylim(0, 10)+
  #xlim(3, 9)+
  theme_half_open()

length(result[,1])
length(saving_df[,1])

#side by side
with_NAs + wo_NAs



