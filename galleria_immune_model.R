#----------------------------
# IMMUNO MODEL GALLERIA
#----------------------------

#----------------------------
# setup
shhh <- suppressPackageStartupMessages # It's a library, so shhh!
shhh(library(tidyverse))
library(deSolve)


#----------------------------
# parameters and initial conditions 

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
  f = 1*10**-3,
  b = 1*0.1,
  K_R = 0.5*10**3, # carrying capacity inside refugee
  
  # immune system
  BL = 1, # baseline level of immune effectors

  loss = 0, # this is a parameter I introduced for playing around
  
  # Pilyugin and Antia model kind of stuff
  reservoir = 1*10**5, # total number of cells in the body
  h_1 = 0.00051, # killing rate for bacteria 
  h_2 = 0.00014, # rate by which they end up in eganged
  d = 0.0005,# return rate to resting state
  a = 0.01, # background activation rate
  g = 0.005, # handling time = 1/g
  s = 0.000035 # per meeting rate/activation rate
)


# populations 
y <- c( 
  
        A = params[["Amax"]], # antibiotics
        B = params[["inoculum"]], # bacteria
        E = params[["BL"]], # active immune effectors
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

# model following Pilyugin and Antia with an integrated handling time 
integrated_handling <- function (t, y, params) {
    
    with(as.list(c(y, params)), {
      
      reservoir = max(c(reservoir-loss*h_1*B*E, 0)) # loss here handles degradation of effectors upon killing

      # antibiotic
      dA <- -k*A
      
      # bacteria 
      dB <- psi(A)*B*(1-B/K) - h_1*B*E - f*(1-R/K_R)*B + b*(1-B/K)*R
      
      # immmune effectors (handling time included here)
      dE <- (a + s*B)*(reservoir-E-((1-loss)*h_2/g)*E*B) - d*E 
      
      # refugee
      dR <- 1*psi(0)*R*(1-R/K_R) + f*(1-R/K_R)*B - b*(1-B/K)*R
      
      return(list(c(dA, dB, dE, dR)))
      
    })
  }


#---------------------------------------
# plotting function 

int_inoc_layers <- function(inoculum, ylimit = log10(params[["K"]]), model=integrated_handling){
  
  colfunc_red <- colorRampPalette(c("blue", "red"))

  for(i in 1:length(inoculum)){
    
    y["B"] = inoculum[i]
    out = ode(times = times, y = y, func = model, parms = params, atol = 10**-8, rtol = 10**-8)
    out <- data.frame(out)
    
    if(i == 1){
      par(mfrow = c(2, 2))
      plot(log10(out$B+out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Total Bacteria\n(Unprotected + Protected)", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Effectors", ylab = "", xlab = "time")
      plot(log10(out$B+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Unprotected", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Protected", ylab = "", xlab = "time")
    }
    par(mfg = c(1, 1)) # pick which plot is being drawn onto
    points(log10(out$B+out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(1, 2))
    points(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(2, 1))
    points(log10(out$B+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
    par(mfg = c(2, 2))
    points(log10(out$R+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3)
  
  }
  parmfrow=c(1, 1) # set window management back to default 
}


#----------------------------
# generating plots 

# default parameters 
times <- seq(0, 25, length = 201) # time span
inocs <- floor(10**seq(3.49, 3.5, 0.001)) # inoc values 
inocs <- floor(c(10**2.5,10**3,10**seq(3.49, 3.5, 0.001),10**4,10**4.5,10**5)) # inoc values 

y["E"] = 10**4 # number of efectors present at time point zero
params["BL"] = 10**4 # this would do the same thing still if the function would recreate y

params["h_1"] = 0.0004 # killing rate for bacteria 
params["h_2"] = 0.000105 # rate by which they end up in enganged
params["d"] = 0.0001 # return rate to resting state
params["a"] = 0.001 #0.001, # background activation rate
params["g"] = 0.0049 #handling time = 1/g
params["s"] = 0.000035 # per meeting rate/activation rate


# creating (and saving) plot
#png("general_visualisation.png", width = 5, height = 5, units = "in", res = 300)
int_inoc_layers(inocs)
#dev.off()







