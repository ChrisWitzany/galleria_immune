#----------------------------
# IMMUNO MODEL GALLERIA
#----------------------------

# libraries
library(tidyverse)
library(deSolve)


params <- c(
  
  inoculum = 10**3, # starting population size of bacteria
  K = 10**7, # carrying capacity
  
  # pharmacokinetics
  Amax = 0, # antibiotic dose
  k = 1, # decay/excretion rate of the antibiotic
  
  # pharmacodynamics
  psimax = 1,
  psimin = -5,
  kappa = 1.5,
  MIC = 1, 
  
  # bacteria refugee
  f = 1*10**-2,#0.0001,
  b = 1*0.1,
  K_R = 10**3, # carrying capacity inside refugee
  
  # immune system
  BL = 1, # baseline level of immune effectors
  reservoir = 10*10**4, # total number of cells in the body
  
  delta = 0.000035, # induction rate 
  beta = 0.0001,  # max immune killing rate
  gamma = 0.01, # effector expenditure (per killing)
  zeta = 0.005, #immune handling time    
  h = 10**-6,

  # trying Pilyugin and Antia model kind of stuff
  h_1 = 0.00051, # killing rate for bacteria 
  h_2 = 0.00014, # rate by which they end up in eganged
  d = 0.0005,# return rate to resting state
  a = 0.01, #0.001, # background activation rate
  g = 0.005, #handling time = 1/g
  s = 0.000035 # per meeting rate/activation rate
)


# populations 
y <- c( 
  
  A = params[["Amax"]], # antibiotics
  B = params[["inoculum"]], # bacteria
  E = params[["BL"]], # active immune effectors
  H = 0,
  R = 0 # refugee compartment 
)


# pharmacodynamic function
psi <- function(A){
  with(as.list(params),{
    
    psimax-(psimax-psimin)*((A/MIC)^kappa)/(((A/MIC)^kappa)-psimin/psimax)
    
  })
}

#----------------------------
# setting up the models
#----------------------------

# induction and handling time seperated equations 
sep_handling_pop <- function (t, y, params) {
  
  with(as.list(c(y, params)), {
    
    # antibiotic
    dA <- -k*A
    
    # bacteria 
    dB <- psi(A)*B*(1-B/K)  - f*(1-R/K_R)*B + b*(1-B/K)*R - beta*B*E #1/(1+h*B)*beta*B*E
    
    # immmune effectors (representative for hemocytes, AMPs, NETs, etc.)    
    dE <- delta*B*(reservoir-E-H)  + (1-gamma)*zeta*H - beta*E*B#1/(1+h*B)*beta*E*B 

    # immmune effectors handling/active (representative for hemocytes, AMPs, NETs, etc.)   
    dH <- - zeta*H + beta*E*B#1/(1+h*B)*beta*E*B# h2 is the engagement rate # g is when they kill smth
    
    # refugee
    dR <- 1*psi(0)*R*(1-R/K_R) + f*(1-R/K_R)*B - b*(1-B/K)*R
    
    return(list(c(dA, dB, dE, dH, dR)))
    
  })
}

sep_handling_pop_w_dampening <- function (t, y, params) {
  
  with(as.list(c(y, params)), {
    
    # antibiotic
    dA <- -k*A
    
    # bacteria 
    dB <- psi(A)*B*(1-B/K)  - f*(1-R/K_R)*B + b*(1-B/K)*R - 1/(1+h*B)*beta*B*E
    
    # immmune effectors (representative for hemocytes, AMPs, NETs, etc.)    
    dE <- delta*B*(reservoir-E-H)  + (1-gamma)*zeta*H - 1/(1+h*B)*beta*E*B 
    
    # immmune effectors handling/active (representative for hemocytes, AMPs, NETs, etc.)   
    dH <- - zeta*H + 1/(1+h*B)*beta*E*B# h2 is the engagement rate # g is when they kill smth
    
    # refugee
    dR <- 1*psi(0)*R*(1-R/K_R) + f*(1-R/K_R)*B - b*(1-B/K)*R
    
    return(list(c(dA, dB, dE, dH, dR)))
    
  })
}



# model following Pilyugin and Antia with an integrated handling time 
integrated_handling <- function (t, y, params) {
    
    with(as.list(c(y, params)), {
      
      # antibiotic
      dA <- -k*A
      
      # bacteria 
      dB <- psi(A)*B*(1-B/K) - h_1*B*E - f*(1-R/K_R)*B + b*(1-B/K)*R
      
      # immmune effectors (handling time included here)
      dE <- (a + s*B)*(reservoir-E-(h_2/g)*E*B) - d*E  #gamma is the percentage that gets used up during killing?
      
      # refugee
      dR <- 1*psi(0)*R*(1-R/K_R) + f*(1-R/K_R)*B - b*(1-B/K)*R
      
      return(list(c(dA, dB, dE, dR)))
      
    })
  }


#---------------------------------------
# plotting function 

# color setup (I want to move this honestly)

colfunc_blue <- colorRampPalette(c("lightblue", "orange"))
colfunc_cyan <- colorRampPalette(c("black", "hotpink"))

sep_inoc_layers <- function(inoculum, ylimit = log10(params[["K"]]), model=sep_handling_pop){
  colfunc_red <- colorRampPalette(c("blue", "red"))
  times <- seq(0, 72, length = 201)
  inocs <- floor(10**seq(1,6,0.5))
    for(i in 1:length(inoculum)){
      y["B"] = inoculum[i]
      out = ode(times = times, y = y, func = model, parms = params, atol = 10**-8, rtol = 10**-8)
      out <- data.frame(out)
    if(i == 1){
      par(mfrow = c(1, 4))
      plot(log10(out$B+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Bacteria", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Effectors Free", ylab = "", xlab = "time")
      plot(log10(out$H+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Effectors Handling", ylab = "", xlab = "time")
      plot(log10(out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Refuge", ylab = "", xlab = "time")
    }
      par(mfg = c(1, 1)) # pick which plot is being drawn onto
      points(log10(out$B+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
      par(mfg = c(1, 2))
      points(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
      par(mfg = c(1, 3))
      points(log10(out$H+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
      par(mfg = c(1, 4))
      points(log10(out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    }
  parmfrow=c(1,1) # set window management back to default 
}


int_inoc_layers <- function(inoculum, ylimit = log10(params[["K"]]), model=integrated_handling){
  
  colfunc_red <- colorRampPalette(c("blue", "red"))

    y <- c( 
    A = params[["Amax"]], # antibiotics
    B = params[["inoculum"]], # bacteria
    E = params[["BL"]], # active immune effectors
    R = 0 )# refugee compartment 
  
  for(i in 1:length(inoculum)){
    y["B"] = inoculum[i]
    out = ode(times = times, y = y, func = model, parms = params, atol = 10**-8, rtol = 10**-8)
    out <- data.frame(out)
    if(i == 1){
      par(mfrow = c(1, 3))
      plot(log10(out$B+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Bacteria", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Effectors", ylab = "", xlab = "time")
      plot(log10(out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Refuge", ylab = "", xlab = "time")
    }
    par(mfg = c(1, 1)) # pick which plot is being drawn onto
    points(log10(out$B+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(1, 2))
    points(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(1, 3))
    points(log10(out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
  }
  parmfrow=c(1,1) # set window management back to default 
}


int_inoc_layers_overlay <- function(inoculum, ylimit = log10(params[["K"]]), model=integrated_handling){
  
  colfunc_red <- colorRampPalette(c("blue", "red"))
  
  y <- c( 
    A = params[["Amax"]], # antibiotics
    B = params[["inoculum"]], # bacteria
    E = params[["BL"]], # active immune effectors
    R = 0 )# refugee compartment 
  
  for(i in 1:length(inoculum)){
    y["B"] = inoculum[i]
    out = ode(times = times, y = y, func = model, parms = params, atol = 10**-8, rtol = 10**-8)
    out <- data.frame(out)
    if(i == 1){
      par(mfrow = c(1, 2))
      plot(log10(out$B+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Bacteria", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Effectors", ylab = "", xlab = "time")
      #plot(log10(out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Refuge", ylab = "", xlab = "time")
    }
    par(mfg = c(1, 1)) # pick which plot is being drawn onto
    points(log10(out$B+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(1, 2))
    points(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(1, 1))
    points(log10(out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, lty = 2)
  }
  parmfrow=c(1,1) # set window management back to default 
}





params["BL"] = 10**4

y["E"] = 5688
sep_inoc_layers(inocs)
sep_inoc_layers(inoculum = 28000)
params["h"] = 10**-20
sep_inoc_layers(inoculum = 28000)


times <- seq(0, 240, length = 201)
inocs <- floor(10**seq(3,3.5,0.5))

y["E"] = 1
int_inoc_layers(inocs)
sep_inoc_layers((inocs))
#params["h"] = 10**-5
sep_inoc_layers((inocs), model = sep_handling_pop_w_dampening)


times <- seq(0, 24, length = 201)
# trying Pilyugin and Antia model kind of stuff
params["h_1"] = 0.0004 # killing rate for bacteria 
params["h_2"] = 0.0001001 # rate by which they end up in enganged
params["d"] = 0.0001 # return rate to resting state
params["a"] = 0.01 #0.001, # background activation rate
params["g"] = 0.0049 #handling time = 1/g
params["s"] = 0.000035 # per meeting rate/activation rate


inocs <- floor(10**seq(2,6,0.25))

params["h_2"] = 0.0001001 # rate by which they end up in eganged
int_inoc_layers(inocs)

params["h_2"] = 0.000101 # rate by which they end up in eganged
int_inoc_layers(inocs)




params["h_2"] = 0.0001 # rate by which they end up in eganged
#params["f"] = 0.0
#params["b"] = 0.0
int_inoc_layers(inocs)


