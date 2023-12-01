#----------------------------
# GALLERIA IMMUNE SIMULATIONS
#----------------------------
# simulation and plotting functions


#----------------------------
# run and plot multiple dynamics  

int_inoc_layers <- function(inoculum, ylimit = log10(params[["K_U"]]), model=integrated_handling, y_0 = y){
  
  tic("for simulating the inoc layers plot")
  
  colfunc_red <- colorRampPalette(c("blue", "red"))
  
  for(i in 1:length(inoculum)){
    
    y_0["U"] = inoculum[i]
    out = ode(times = times, y = y_0, func = model, parms = params, atol = 10**-8, rtol = 10**-8, maxsteps = 10**4)
    out <- data.frame(out)
    
    if(i == 1){
      par(mfrow = c(2, 2))
      plot(log10(out$U+out$P+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim = c(0, ylimit), main = "Total Bacteria\n(Unprotected + Protected)", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim = c(0, ylimit), main = "Effectors", ylab = "", xlab = "time")
      plot(log10(out$U+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim = c(0, ylimit), main = "Unprotected", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$P+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim = c(0, ylimit), main = "Protected", ylab = "", xlab = "time")
    }
    par(mfg = c(1, 1)) # pick which plot is being drawn onto
    points(log10(out$U+out$P+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(1, 2))
    points(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(2, 1))
    points(log10(out$U+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(2, 2))
    points(log10(out$P+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    
  }
  par(mfrow = c(1, 1)) # set window management back to default 
  toc()

}


#----------------------------------
# generate inoculum vs end state population sizes data 
inoc_vs_endstate <- function(inoculum, tspan, ylimit = log10(params[["K_U"]]), model=integrated_handling, y_0 = y, params, solver_method = "lsoda"){
  
  tic("simulating the inoc vs. endstate ")
  
  df <- data.frame("inoc"=double(), "end"=double()) 
  
  for(i in 1:length(inoculum)){
    
    y_0["U"] = inoculum[i]
    out = ode(times = tspan, y = y_0, func = model, parms = params, method = solver_method, maxsteps = 10**5) #, atol = 10**-6, rtol = 10**-6,  
    out <- data.frame(out)
    
    df[i, "inoc"] = inoculum[i]
    df[i, "end"] = tail(out,1)$P+tail(out,1)$U
    
  }
  
  toc()
  return(df)
}

# a wrapper function to simulate and plot inoculum vs endstate  
plot_inoc_vs_endstate <- function(inocs, model, tspan = times, params, y_0, solver_method = "lsoda"){
  
  data = inoc_vs_endstate(inocs, tspan = tspan, model = model, params = params, y_0 = y_0, solver_method = solver_method)
  
  plot(log10(data$end+1)~log10(data$inoc), 
       type = "l", 
       main = paste("final population ~ inoculum"), 
       ylab = "final population size (log10)", 
       xlab = "inoculum size (log10)", 
       ylim = c(0, max(log10(data+1)+1, na.rm = T)))
  abline(h=log10(params[["K_U"]]), col = "grey", lwd = 1.5, lty = 2)
  text(x = 1, y =log10(params[["K_U"]])+0.3, label = expression('K'[U]), col = "grey")
  abline(h=log10(params[["K_P"]]), col = "grey", lwd = 1.5, lty = 2)
  text(x = log10(max(inocs))-1, y =log10(params[["K_P"]])+0.3, label = expression('K'[P]), col = "grey", lwd = 2)
  points(log10(data$end)~log10(data$inoc), pch = 18)
  points(log10(data$end)~log10(data$inoc), pch = 18, type = "l")
  
  return(data)
}


# function to generate parameters for sensitivity analysis
generate_parameters <- function(n_samples, lhs_type = "random"){
  
  tic("parameter generation:")
  
  if(lhs_type == "random"){
    
    # quicker for testing:
    lhs <- randomLHS(n = n_samples,          # n = row = number of samples 
                     k = length(params))  # k = col =  number of parameters
    
  } else if(lhs_type == "genetic"){
    
    # generate random parameter set here
    lhs <- geneticLHS(n = n_samples,          # n = row = number of samples 
                      k = length(params))  # k = col =  number of parameters
    
  }else{ stop("lhs_type has to be random or genetic!")}
  
  print("-----------------------")
  print("generating LHS...")
  lhs <- matrix(lhs, nrow = nrow(lhs), ncol = ncol(lhs))  # each row will consist of one sample of the parameter space.
  colnames(lhs) <- names(params) # easier indexing
  
  # transform lhs via params_space to actual parameter values
  for (i in 1:nrow(params_space)){
    
    if(params_space[i,3]=="unif"){
      
      lhs[,i] <- qunif(lhs[,i], 
                       min = as.numeric(params_space[i,1]),
                       max = as.numeric(params_space[i,2])) 
      
    } else if(params_space[i,3]=="log"){
      
      lhs[,i] <- qlunif(lhs[,i],
                        min = as.numeric(params_space[i,1]), 
                        max = as.numeric(params_space[i,2]))
    }
  }
  
  print("...LHS generated!")
  
  toc()
  
  return(lhs)
  
}

# run sensitivity analysis 
sensitivity_analysis <- function(n_samples, inocs = 10**seq(3, 9, 0.25), tspan = seq(0, 48, length = 201), model = simple_seperate_handling, lhs_type = "random", solver_method = "lsoda"){
  
  lhs = generate_parameters(n_samples = n_samples, lhs_type = lhs_type)
  
  print("-----------------------")
  print("running simulations...")
  tic("all senstitivity simulations:")
  
  for(i in 1:n_samples){
    
    params = lhs[i, ]   # change params to values in lhs 
    
    if(i == 1){
      
      saving_df = inoc_vs_endstate(inocs, tspan = tspan, params = params, model = simple_seperate_handling, y_0 = y2, solver_method = solver_method) # if no saving_df create it...
      saving_df[, names(params)] = data.frame(t(params))
      saving_df[, "n_sample"] = i # add an identifier for plotting 
      
    }else{ 
      
      to_append = inoc_vs_endstate(inocs, tspan = tspan, params = params, model = simple_seperate_handling, y_0 = y2, solver_method = solver_method)
      to_append[ , names(params)] = data.frame(t(params))
      to_append[, "n_sample"] = i
      saving_df = dplyr::bind_rows(saving_df, to_append)  # ...otherwise append to it
      
    }
    
  }
  
  print("...simulations r done!")
  print("-----------------------")
  toc()
  
  return(saving_df)
  
}
