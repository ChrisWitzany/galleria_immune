#----------------------------
# GALLERIA IMMUNE PARAMETERS
#----------------------------

params <- c(
  
  inoculum = 10**3, # starting population size of bacteria
  
  # pharmacokinetics
  Amax = 0, # antibiotic dose
  k = 1, # decay/excretion rate of the antibiotic
  
  # pharmacodynamics
  psimax = 1,
  psimin = -5,
  kappa = 1.5,
  MIC = 1, 
  
  # migration un/protected site
  K_U = 10**7, # carrying capacity in unprotected site
  K_P = 1.5*10**3, # carrying capacity inside the protected site
  f = 1*10**-3,
  b = 1*0.1,
  
  # immune system - Pilyugin and Antia model kind of stuff
  E_0 = 10**4, # baseline level of immune effectors
  E_tot = 1*10**5, # total number of cells in the body
  h_1 = 0.0004, # killing rate for bacteria 
  h_2 = 0.000105, # rate by which they end up in enganged
  d = 0.0001, # return rate to resting state
  a = 0.01, # background activation rate
  g = 0.00049, # handling time = 1/g
  s = 0.000035 # per meeting rate/activation rate
)


params_space <- c(
  
  inoculum_min = 0, inoculum_max = 0, inoculum_samp = "unif", # starting population size of bacteria
  
  # pharamcokinetics
  Amax_min = 0, Amax_max = 0, A_max_samp = "unif", # antibiotic dose
  k_min = 1, k_max = 1, k_samp = "unif", # decay/excretion rate of the antibiotic
  
  # pharmacodynamics
  psimax_min = 1, psimax_max = 1, psimax_samp = "unif",
  psimin_min = -5, psimin_min = -5 , psimin_samp = "unif",
  kappa_min = 1.5, kappa_max = 1.5 , kappa_samp = "unif",
  MIC_min = 1, MIC_max = 1, MIC_samp = "unif", 
  
  # migration
  K_U_min = 10**9 , K_U_max = 10**9 , K_U_samp = "log", # carrying capacity in unprotected site
  K_P_min = 1*10**3, K_P_max = 1*10**3, K_P_samp = "log",
  f_min = 1*10**-7, f_max = 1*0.1 , f_samp = "log",
  b_min = 1*10**-7, b_max = 1*0.1 , b_samp = "log",
  
  # immune system
  E_0_min = 10**4, E_0_max = 10**4, E_0_samp = "log",  # baseline level of immune effectors
  E_tot_min = 10**6, E_tot_max = 10**6, E_tot_samp = "log", # total number of cells in the body
  h_1_min = 10**-7, h_1_max = 0.1, h_1_samp = "log", # killing rate for bacteria 
  h_2_min = 10**-7, h_2_max = 0.1, h_2_samp = "log", # rate by which they end up in enganged
  d_min = 10**-7, d_max = 1, d_samp = "log", # return rate to resting state
  a_min = 10**-7, a_max = 0.5, a_samp = "log", # background activation rate
  g_min = 10**-7, g_max = 4, g_samp = "log", # handling time = 1/g
  s_min = 10**-7, s_max = 0.1, s_samp = "log" # per meeting rate/activation rate
  
)


# NOTE THAT THE ORDER MATTERS!
# hence some sanity checks 
if(length(params_space)/3 != length(params) || 
   all(gsub('_min', '', names(params_space)[seq(1, length(params_space), 3)]) != names(params))){
  stop("params_space and params dont match")
}

params_space <- matrix(params_space, ncol = 3, byrow = T) # easier to index
if(nrow(params_space) != length(params)){
  stop("something went wrong with the param_space")
} # another sanity check


# populations vector
y <- c( 
  
  A = params[["Amax"]], # antibiotics
  U = params[["inoculum"]], # bacteria
  E = params[["E_0"]], # active immune effectors
  P = 0 # protected site 
  
)

# populations vector
y2 <- c( 
  
  U = params[["inoculum"]], # bacteria
  E = params[["E_0"]], # active immune effectors
  S = 0, 
  P = 0 # protected site 
  
)
