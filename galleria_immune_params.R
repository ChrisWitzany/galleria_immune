#----------------------------
# GALLERIA IMMUNE PARAMETERS
#----------------------------


params <- c(
  
  inoculum = 10**3, # starting population size of bacteria
  
  r = 1, # bacterial growth rate 
  
  # migration un/protected site
  K_U = 10**8, # carrying capacity in unprotected site
  K_P = 1.5*10**3, # carrying capacity inside the protected site
  f = 0.001,
  b = 0.1,
  
  # immune system - Pilyugin and Antia model kind of stuff
  E_0 = 10**3, # baseline level of immune effectors
  E_tot = 1*10**5, # total number of cells in the body
  h_1 = 0.001, # killing rate for bacteria 
  h_2 = 0.001, # rate by which they end up in enganged
  d = 0.01, # return rate to resting state
  a = 0.01, # background activation rate
  g = 0.5, # handling time = 1/g
  s = 0.001 # per meeting rate/activation rate
)


# populations vector
y <- c( 
  
  U = params[["inoculum"]], # bacteria
  E = params[["E_0"]], # active immune effectors
  P = 0 # protected site 
  
)




params_space <- c(
  
  inoculum_min = 0, inoculum_max = 0, inoculum_samp = "unif", # starting population size of bacteria
  
  r_min = 1, r_max = 1, r_samp = "unif", # growth rate 
  
  # migration
  K_U_min = 10**8 , K_U_max = 10**8 , K_U_samp = "log", # carrying capacity in unprotected site
  K_P_min = 1*10**2, K_P_max = 1*10**4, K_P_samp = "log",
  f_min = 1*10**-7, f_max = 1*0.1 , f_samp = "log",
  b_min = 1*10**-7, b_max = 1*0.1 , b_samp = "log",
  
  # immune system
  E_0_min = 1, E_0_max = 10**5, E_0_samp = "log",  # baseline level of immune effectors
  E_tot_min = 10**5, E_tot_max = 10**5, E_tot_samp = "unif", # total number of cells in the body
  h_1_min = 10**-7, h_1_max = 0.1, h_1_samp = "log", # killing rate for bacteria 
  h_2_min = 10**-7, h_2_max = 0.1, h_2_samp = "log", # rate by which they end up in enganged
  d_min = 10**-4, d_max = 1, d_samp = "log", # return rate to resting state
  a_min = 10**-4, a_max = 1, a_samp = "log", # background activation rate
  g_min = 10**-3, g_max = 4, g_samp = "log", # handling time = 1/g
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

# make naming intuitive
rownames(params_space) = names(params)


# function to create LHS with optimization via a genetic algorithm (better coverage)
generate_lhs <- function(N, params_to_sample, params){
  
  k <- length(params_to_sample) # how many parameters?
  warning("Note that params_to_sample need to make sense with params_space - or unexpected behaviour could happen!")
  
  print("creating genetic LHS (this might take a while...)")
  tic("during genetic LHS generation")
  lhs <- geneticLHS(n = N, # n = row = number of samples 
                    k = k, # k = column = number of parameters sampled
                    pop = 50, # how many lhs are generated
                    gen = 5, # how many generations
                    pMut = 0.25, # how likely for a row to change/mutate
                    criterium = "Maximin", # which method for optimizing cooverage
                    verbose = T # give progress
  ) 
  
  # set names for convenience
  colnames(lhs) <- params_to_sample
  
  # add the non-sampled parameters 
  # these will be constant over all simulations and given the default values in the scale function
  constant_params <- names(params)[!rownames(params_space) %in% params_to_sample]
  
  lhs <- as.data.frame(lhs)
  for(n in constant_params){ lhs[, n] <- 1 }
  lhs <- as.matrix(lhs)
  
  toc()
  
  return(lhs)
  
}


# function to transform parameters according to designated param_space 
scale_lhs <- function(lhs, params_space = params_space){
  
  # make sure the order is correct
  lhs <- lhs[, names(params)]
  
  # transform lhs via params_space to actual parameter values
  for (i in 1:nrow(params_space)){
    
    if(params_space[i,3] == "unif"){
      
      lhs[,i] <- qunif(lhs[,i], 
                       min = as.numeric(params_space[i,1]),
                       max = as.numeric(params_space[i,2])) 
      
    } else if(params_space[i,3] == "log"){
      
      lhs[,i] <- qlunif(lhs[,i],
                        min = as.numeric(params_space[i,1]), 
                        max = as.numeric(params_space[i,2]))
    }
  }
  
  return(lhs)
  
}


