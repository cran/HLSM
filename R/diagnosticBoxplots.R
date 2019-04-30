# Function for plotting the diagnostic plots which auto detects the type of HLSM object and plots accordingly

HLSMcovplots<-function(model,burnin=0,thin=1)
{
  model=list(model)
  length(model)
  if(model[[1]]$call[[1]]=="HLSMfixedEF"){
    HLSMfixedboxplots(model[[1]],burnin=burnin,thin=thin)
    
  }else if(model[[1]]$call[[1]]=="HLSMrandomEF"){
   HLSMrandomboxplots(model[[1]],burnin=burnin,thin=thin)
  }
}

# Function that plots the diagnostic boxplots for fixed effects model
HLSMfixedboxplots<-function(model,burnin=burnin,thin=thin){
  if(length(model$draws$ZZ)<=4000){
    t=paste("Looks like your chain length is less than the minimum required value length of 3746, try running a chain length greater than 4000 to see the boxplots")
    return(t)
  }
  else{
  i=getIntercept(model,burnin = burnin,thin = thin)
  betas=getBeta(model,burnin = burnin,thin = thin)
  beta_table=as.data.frame(betas)
  beta_mcmc=as.mcmc(beta_table)
  intercept_mcmc=as.mcmc(i)
  ri=raftery.diag(intercept_mcmc)
  #getting the burnin and the optimum chain length based on the betas
  rb=raftery.diag(beta_table)
  rdf=data.frame(rb[[2]])
  chain_length=max(sapply(rdf[,2],max))
  burnin_value =rdf[which.max(rdf[,2]),1]
  thinnin_value=rdf[which.max(rdf[,2]),4]
  #plotting the boxplots 
  g= plotHLSM.fixed.fit(model,parameter="Beta",burnin=burnin_value,thin=thinnin_value)
  return (c(g))
  }
}

# Function that plots the diagnostic bocplots for random effects model
HLSMrandomboxplots<-function(model,burnin=burnin,thin=thin){
  if(length(model$draws$ZZ)<=4000){
    t=paste("Looks like your chain length is less than the minimum required value length of 3746, try running a chain length greater than 4000 to see the boxplots")
    return(t)
  }
  else{
  intercept=getIntercept(model,burnin=burnin,thin=thin)
  beta_list <- list()
  for(i in 1:ncol(intercept))
  {
    beta_list[[i]]=model$draws$Beta[,,i]
  }
  
  beta_list_df=do.call(cbind.data.frame, beta_list)
  beta_mcmc=as.mcmc(beta_list_df)
  #getting the burnin and the optimum chain length based on the betas
  rb=raftery.diag(beta_list_df)
  rdf=data.frame(rb[[2]])
  chain_length=max(sapply(rdf[,2],max))
  burnin_value =rdf[which.max(rdf[,2]),1]
  thinnin_value=rdf[which.max(rdf[,2]),4]
  #plotting the boxplots 
  g= plotHLSM.random.fit(model,parameter="Beta",burnin=burnin_value,thin=thinnin_value)
  return (c(g))
  }
}