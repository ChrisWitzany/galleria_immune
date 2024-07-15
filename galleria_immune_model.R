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
shhh(library(cowplot))
library(tictoc)
shhh(library(patchwork))
shhh(library(extrafont))


#----------------------------
# setting up the model 
# inspired by Pilyugin and Antia with an integrated handling time 

integrated_handling <- function (t, y, params) {
    
    with(as.list(c(y, params)), {

      # bacteria 
      dU <- r*U*(1-U/K_U) - h_1*U*E - f*(1-P/K_P)*U + b*(1-U/K_U)*P
      
      # immune effectors
      dE <- (a + s*U)*(E_tot-E-(h_2/g)*E*U) - d*E 
      
      # protected site 
      dP <- r*P*(1-P/K_P) + f*(1-P/K_P)*U - b*(1-U/K_U)*P
      
      # debugging
      #if(any(is.na(c(dU, dE, dP)))){ browser() } 
      
      return(list(c(dU, dE, dP)))
      
    })
}

