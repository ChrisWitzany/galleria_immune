#----------------------------
# IMMUNO MODEL GALLERIA
#----------------------------

#----------------------------
# setup
shhh <- suppressPackageStartupMessages # it's a library, so shhh!
shhh(library(tidyverse))
library(deSolve)
library(lhs)        # as the name suggests, needed for lhs generation
library(KScorrect)  # needed to generate uniform log distributions via the dlunif() function
library(cowplot)
library(tictoc)
library(patchwork)
#--------------------------------------


# pharmacodynamic function
psi <- function(A){
  with(as.list(params),{
    
    psimax-(psimax-psimin)*((A/MIC)^kappa)/(((A/MIC)^kappa)-psimin/psimax)
    
  })
}


#----------------------------
# setting up the model 
# inspired by Pilyugin and Antia with an integrated handling time 

simple_seperate_handling <- function (t, y, params) {
  
  with(as.list(c(y, params)), {

    # bacteria 
    dU <- psimax*U*(1-U/K_U) - h_1*U*E - f*(1-P/K_P)*U + b*(1-U/K_U)*P
    
    # immune effectors
    dE <- (a + s*U)*(E_tot-E-S) - h_2*E*U + g*S - d*E 
    
    # immune effectors engaged
    dS <- h_2*E*U - g*S
    
    # protected site 
    dP <- psimax*P*(1-P/K_P) + f*(1-P/K_P)*U - b*(1-U/K_U)*P
    
    # debugging
    if(any(is.na(c(dU, dE, dS, dP)))){ browser() } 
    
    return(list(c(dU, dE, dS, dP)))
    
  })
}


integrated_handling <- function (t, y, params) {
    
    with(as.list(c(y, params)), {
      
      # antibiotic
      dA <- -k*A
      
      # bacteria 
      dU <- psi(A)*U*(1-U/K_U) - h_1*U*E - f*(1-P/K_P)*U + b*(1-U/K_U)*P
      
      # immune effectors
      dE <- (a + s*U)*(E_tot-E-(h_2/g)*E*U) - d*E 
      
      # protected site 
      dP <- 1*psi(0)*P*(1-P/K_P) + f*(1-P/K_P)*U - b*(1-U/K_U)*P
      
      return(list(c(dA, dU, dE, dP)))
      
    })
  }