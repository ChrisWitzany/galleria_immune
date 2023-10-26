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
  K_U = 10**7, # carrying capacity
  
  # pharmacokinetics
  Amax = 0, # antibiotic dose
  k = 1, # decay/excretion rate of the antibiotic
  
  # pharmacodynamics
  psimax = 1,
  psimin = -5,
  kappa = 1.5,
  MIC = 1, 
  
  # protected site
  f = 1*10**-3,
  b = 1*0.1,
  K_P = 0.5*10**3, # carrying capacity inside the protected site
  
  # immune system
  BL = 10**4, # baseline level of immune effectors
  loss = 0, # this is a parameter I introduced for playing around
  
  # Pilyugin and Antia model kind of stuff
  reservoir = 1*10**5, # total number of cells in the body
  h_1 = 0.0004, # killing rate for bacteria 
  h_2 = 0.000105, # rate by which they end up in eganged
  d = 0.0001,# return rate to resting state
  a = 0.001, # background activation rate
  g = 0.0049, # handling time = 1/g
  s = 0.000035 # per meeting rate/activation rate
)


# populations 
y <- c( 
  
        A = params[["Amax"]], # antibiotics
        U = params[["inoculum"]], # bacteria
        E = params[["BL"]], # active immune effectors
        P = 0 # refugee compartment 

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
      
      reservoir = max(c(reservoir-loss*h_1*U*E, 0)) # loss here handles degradation of effectors upon killing

      # antibiotic
      dA <- -k*A
      
      # bacteria 
      dU <- psi(A)*U*(1-U/K_U) - h_1*U*E - f*(1-P/K_P)*U + b*(1-U/K_U)*P
      
      # immmune effectors (handling time included here)
      dE <- (a + s*U)*(reservoir-E-((1-loss)*h_2/g)*E*U) - d*E 
      
      # baprotected site 
      dP <- 1*psi(0)*P*(1-P/K_P) + f*(1-P/K_P)*U - b*(1-U/K_U)*P
      
      return(list(c(dA, dU, dE, dP)))
      
    })
  }


#---------------------------------------
# plotting function 

int_inoc_layers <- function(inoculum, ylimit = log10(params[["K_U"]]), model=integrated_handling){
  
  colfunc_red <- colorRampPalette(c("blue", "red"))

  for(i in 1:length(inoculum)){
    
    y["U"] = inoculum[i]
    out = ode(times = times, y = y, func = model, parms = params, atol = 10**-8, rtol = 10**-8)
    out <- data.frame(out)
    
    if(i == 1){
      par(mfrow = c(2, 2))
      plot(log10(out$U+out$P+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Total Bacteria\n(Unprotected + Protected)", ylab = "pop size (log10)", xlab = "time")
      abline(h=log10(3908.083), col = "black", lwd = 10)
      plot(log10(out$E+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Effectors", ylab = "", xlab = "time")
      plot(log10(out$U+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Unprotected", ylab = "pop size (log10)", xlab = "time")
      plot(log10(out$P+1)~out$time, type = "l", col = colfunc_red(length(inoculum))[i], lwd = 3, ylim=c(0,ylimit), main = "Protected", ylab = "", xlab = "time")
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
  parmfrow=c(1, 1) # set window management back to default 
}


#----------------------------
# generating plots 

# some exemplary parameters 
times <- seq(0, 25, length = 201) 
inocs <- floor(c(10**2.5,10**3,10**seq(3.49, 3.5, 0.001),10**4,10**4.5,10**5)) 

# creating (and saving) plot
png("with_protected_site.png", width = 5, height = 5, units = "in", res = 300)
int_inoc_layers(inocs)
dev.off()



# lets approach this from a bit more of a theoretical side

# first lets turn off the protected site
params[["f"]] = 0 
params[["b"]] = 0


# next we utilize the analytical solution of Pilugyin i.e. 
params[["reservoir"]] = 10**4
params[["s"]] = 3.5e-03
params[["d"]] = 0.1

with(as.list(c(y, params)), {
  
  Phat = sqrt(d/(s*(h_2/g)))-(a/s)
  FatP = (a+s*Phat)*reservoir/(d+(a+s*Phat)*(1+(h_2/g)*Phat))
  
  print(FatP)
  print(psimax/h_1)

  print(paste0("a*reservoir/(d+a)-psimax/h_1 < 0 is ", a*reservoir/(d+a)-psimax/h_1 < 0))
  print(paste0("F(Phat) < r/h_1 is ", FatP < psimax/h_1))
  
  # some logic to help interpretation:
  if(FatP < psimax/h_1 && a*reservoir/(d+a)-psimax/h_1 < 0){print(" --> no control possible - only escape or eradicaiton of pathogens")} # this also meany no steady states for pathogens > 0 because growth is essentially "uncontrolled"
  if(FatP > psimax/h_1 && a*reservoir/(d+a)-psimax/h_1 < 0){print(" --> control possible - and stable")}
  
})

# visualize that control is always possible, but might take some time 
times <- seq(0, 7*24, length = 201) 
params["K_U"] = 10**7

png("without_protected_site.png", width = 5, height = 5, units = "in", res = 300)
int_inoc_layers(floor(c(1**10, 10**seq(2.41, 2.43, 0.001))) )
dev.off()
# extinction is impossible


# note that protected site kills this behaviour at least for these parameters
params[["f"]] = 1*10**-3
params[["b"]] = 0.1
png("same_paras_but_with_protected_site.png", width = 5, height = 5, units = "in", res = 300)
int_inoc_layers(floor(c(1**10, 10**seq(2.41, 2.43, 0.001))) )
dev.off()
