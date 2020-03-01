# Function for plotting the diagnostic plots which auto detects the type of HLSM object and plots accordingly

HLSMcovplots<-function(fitted.model,burnin=0,thin=1)
{
  fitted.model=list(fitted.model)
  length(fitted.model)
  if(fitted.model[[1]]$call[[1]]=="HLSMfixedEF"){
    HLSMfixed.covplots(m=fitted.model[[1]],burninvalue=burnin,thinvalue=thin)
    
  }else if(fitted.model[[1]]$call[[1]]=="HLSMrandomEF"){
   HLSMrandom.covplots(m=fitted.model[[1]],burninvalue=burnin,thinvalue=thin)
  }
}

# Function that plots the diagnostic boxplots for fixed effects model
HLSMfixed.covplots<-function(m,burninvalue,thinvalue){
  if(length(m$draws$ZZ)<=4000){
    t=paste("Try running a chain length greater than 4000 to see parameter plots")
    return(t)
  }
  else{
  i=getIntercept(m,burnin = burninvalue,thin = thinvalue)
  betas=getBeta(m,burnin = burninvalue,thin = thinvalue)
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
  g= plotHLSM.fixed.fit(m,parameter="Beta",burnin=burnin_value,thin=thinnin_value)
  return (c(g))
  }
}

# Function that plots the diagnostic bocplots for random effects model
HLSMrandom.covplots<-function(m,burninvalue,thinvalue){
  if(length(m$draws$ZZ)<=4000){
    t=paste("Try running a chain length greater than 4000 to see the parameter plots")
    return(t)
  }
  else{
  intercept=getIntercept(m,burnin=burninvalue,thin=thinvalue)
  beta_list <- list()
  for(i in 1:ncol(intercept))
  {
    beta_list[[i]]=m$draws$Beta[,,i]
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
  g= plotHLSM.random.fit(m,parameter="Beta",burnin=burnin_value,thin=thinnin_value)
  return (c(g))
  }
}
