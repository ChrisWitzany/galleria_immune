#------------------------------------------
# generating some maybe interesting plots 

# some exemplary parameters 
times <- seq(0, 25, length = 201) 
inocs <- floor(c(10**2.5,10**3,10**seq(3.49, 3.5, 0.001),10**4,10**4.5,10**5)) 

# creating (and saving) plot
#png("with_protected_site.png", width = 5, height = 5, units = "in", res = 300)
int_inoc_layers(inocs)
#dev.off()



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

#png("without_protected_site.png", width = 5, height = 5, units = "in", res = 300)
int_inoc_layers(floor(c(1**10, 10**seq(2.41, 2.43, 0.001))) )
#dev.off()
# extinction is impossible


# note that protected site kills this behaviour at least for these parameters
params[["f"]] = 1*10**-3
params[["b"]] = 0.1
#png("same_paras_but_with_protected_site.png", width = 5, height = 5, units = "in", res = 300)
int_inoc_layers(floor(c(1**10, 10**seq(2.41, 2.43, 0.001))) )
#dev.off()





int_inoc_layers(floor(seq(10**1:9)) )


y["U"] = 10**10
params[["K_U"]] = 10**9
out = ode(times = times, y = y, func = model, parms = params, atol = 10**-8, rtol = 10**-8)
out <- data.frame(out)
out

