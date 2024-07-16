#---------------------------------
# GALLERIA IMMUNE SIMULATIONS
#---------------------------------
# simulation and plotting functions


#--------------------------------
# run and plot multiple dynamics  

int_inoc_layers <- function(inocs, ylimit = log10(params[["K_U"]]), model=integrated_handling, y_0 = y){
  
  tic("for simulating the inoc layers plot")
  
  colfunc_red <- colorRampPalette(c("blue", "red"))
  
  for(i in 1:length(inocs)){
    
    y_0["U"] = inocs[i]
    out = ode(times = times, y = y_0, func = model, parms = params, atol = 10**-8, rtol = 10**-8, maxsteps = 10**4)
    out <- data.frame(out)
    
    if(i == 1){
      par(mfrow = c(2, 2))
      plot(log10(out$U+out$P+1)~out$time, type = "l", col = colfunc_red(length(inocs))[i], lwd = 3, ylim = c(0, ylimit), main = "Total Bacteria\n(Unprotected + Protected)", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inocs))[i], lwd = 3, ylim = c(0, ylimit), main = "Effectors", ylab = "", xlab = "time")
      plot(log10(out$U+1)~out$time, type = "l", col = colfunc_red(length(inocs))[i], lwd = 3, ylim = c(0, ylimit), main = "Unprotected", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$P+1)~out$time, type = "l", col = colfunc_red(length(inocs))[i], lwd = 3, ylim = c(0, ylimit), main = "Protected", ylab = "", xlab = "time")
    }
    par(mfg = c(1, 1)) # pick which plot is being drawn onto
    points(log10(out$U+out$P+1)~out$time, type = "l", col = colfunc_red(length(inocs))[i], lwd = 3)
    par(mfg = c(1, 2))
    points(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inocs))[i], lwd = 3)
    par(mfg = c(2, 1))
    points(log10(out$U+1)~out$time, type = "l", col = colfunc_red(length(inocs))[i], lwd = 3)
    par(mfg = c(2, 2))
    points(log10(out$P+1)~out$time, type = "l", col = colfunc_red(length(inocs))[i], lwd = 3)
    
    # if first run create saving dataframe - otherwise add to it 
    if(i == 1){
      saving_df <- out # create saving df 
    }else{
      saving_df <- rbind(saving_df, out)
    }
  }
  par(mfrow = c(1, 1)) # set window management back to default 
  
  toc()
  return(saving_df)
  
}


#----------------------------------
# generate inoculum vs end state population sizes data 
inoc_vs_endstate <- function(inocs, tspan = times, ylimit = log10(params[["K_U"]]), model=integrated_handling, y_0 = y, params, solver_method = "lsoda"){
  
  tic("simulating the inoc vs. endstate ")
  
  df <- data.frame("inoc" = double(), "end" = double()) 
  
  for(i in 1:length(inocs)){
    
    y_0["U"] = inocs[i]
    out = ode(times = tspan, y = y_0, func = model, parms = params, method = solver_method, maxsteps = 10**4) #, atol = 10**-6, rtol = 10**-6,  
    out <- data.frame(out)
    
    df[i, "inoc"] = inocs[i]
    df[i, "end"] = tail(out, 1)$P+tail(out, 1)$U
    df[i, "time"] = tail(out, 1)$time
    
  }
  
  toc()
  return(df)
}


# a wrapper function to simulate and plot inoculum vs endstate  
plot_inoc_vs_endstate <- function(inocs, model, tspan = times, params, y_0, solver_method = "lsoda", thres = 10**7){
  
  data = inoc_vs_endstate(inocs, tspan = tspan, model = model, params = params, y_0 = y_0, solver_method = solver_method)
  
  thres_val <- data$inoc[which(data$end > thres)[1]]
  #print(thres_val) # debugging
  
  plot(log10(data$end+1)~log10(data$inoc), 
       type = "l", 
       main = paste("final population ~ inoculum"), 
       ylab = "final population size (log10)", 
       xlab = "inoculum size (log10)", 
       ylim = c(0, max(log10(data+1)+1, na.rm = T)))
  abline(v = log10(thres_val), lty = 2, col = "red", lwd = 2)
  abline(h=log10(params[["K_U"]]), col = "grey", lwd = 1.5, lty = 2)
  text(x = 2, y =log10(params[["K_U"]])+0.4, label = expression('K'[U]), col = "grey")
  abline(h=log10(params[["K_P"]]), col = "grey", lwd = 1.5, lty = 2)
  text(x = log10(max(inocs))-1, y =log10(params[["K_P"]])+0.4, label = expression('K'[P]), col = "grey", lwd = 2)
  points(log10(data$end)~log10(data$inoc), pch = 18)
  points(x = log10(thres_val), y = log10(data$end[data$inoc==thres_val]), col = "red", lwd = 2, cex = 2)
  points(log10(data$end)~log10(data$inoc), pch = 18, type = "l", lwd = 2)
  
  return(data)
}


# a wrapper function that returns only the threshold value when given a set of parameters
get_threshold <- function(inocs, model, tspan = times, params, y_0, solver_method = "lsoda", cutoff = 1e+07){
  
  # building this function
  data = inoc_vs_endstate(inocs, tspan = tspan, model = model, params = params, y_0 = y_0, solver_method = solver_method)
  #data = inoc_vs_endstate(inocs, model = integrated_handling, params = params, y_0 = y, solver_method = "lsoda")
  
  # check for terminations
  if(any(data$time < max(data$time))){
    return("terminated")
  }

  threshold_inoc <- data%>%
    mutate(
            value = end >= cutoff, # total pop over cutoff (default 10**7)
            fall_below_after = if_else(value, lead(end < cutoff, default = FALSE), FALSE)
          )%>%  # and correct this if it ever goes down again
    summarise(
            first_jump_value = ifelse(any(value), first(inoc[value]), NA_real_), # get first inoculum for which that happens
            falls_below_after_jump = any(fall_below_after)
             )
  
  if(threshold_inoc$falls_below_after_jump == FALSE){
    
    return(threshold_inoc$first_jump_value)
  
  }else{
    
    return("threshold but drops after")
  
  }

}


