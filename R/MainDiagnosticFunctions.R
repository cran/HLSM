
#Main diagnostic function that detects the type of HLSM object passed and provides the correct diagnostic
HLSMDiag<-function(...,burnin=0,thin=1,lag=500,type="0",varnum=1){
 
 model=list(...)
 l=length(model)
 
  if(model[[1]]$call[[1]]=="HLSMfixedEF"){
    HLSMFixedEffectsDiagnostics(model=model,burnin=burnin,thin=thin,lag=lag,type=type)
 
}else if(model[[1]]$call[[1]]=="HLSMrandomEF"){
  HLSMrandomEffectsDiagnostics(model=model,burnin=burnin,thin=thin,lag=lag,type=type,varnum=varnum)
}
  
}
## Helper function for the below Fixed effects diagnostics function
## This function can inform the user of the optimum thinning , burnin and chain length values when one fixed effects object is passed.
## This function will be called from the main fixed effects diagnostic function and hence we do not need to use this function explicitly 

HLSMfixedDiagnostics<-function(model,burnin=burnin,thin=thin,lag=lag,type=type){
  if(length(model$draws$ZZ)<=4000){
    t=paste("Looks like your chain length is less than the minimum required value length of 3746, try running a chain length greater than 4000")
    return(t)
  }else{
    i=getIntercept(model,burnin = burnin,thin = thin)
    
    betas=getBeta(model,burnin = burnin,thin = thin)
    beta_table=as.data.frame(betas)
    beta_mcmc=as.mcmc(beta_table)
    intercept_mcmc=as.mcmc(i)
    ri=raftery.diag(intercept_mcmc)
    rb=raftery.diag(beta_table)
    rdf=data.frame(rb[[2]])
    chain_length=max(sapply(rdf[,2],max))
    burnin_value =rdf[which.max(rdf[,2]),1]
    thinnin_value=rdf[which.max(rdf[,2]),4]
    t=paste("Based on your betas, the optimum chain length is:",chain_length,".The optimum burnin value is:",burnin_value,".The optimum thinning value is :",thinnin_value)
    
    if (type=="traceplot"){
      g=traceplot(beta_mcmc,col=1:ncol(beta_table))
      return (c(t,g))
    }
    else if (type=="autocorr"){
      par(mar=c(1,1,1,1))
      g=acf(beta_mcmc,lag.max =lag)$type
      return (c(g,t))
    }else if(type=="0"){
      return (c(t))
    }
  }
}

#### Fixed effects diagnostic function that calculates the optimal chain length, burnin and thinning value if only one fixed effects model is passed and can also compute the convergence if you pass two or more fixed effect objects
HLSMFixedEffectsDiagnostics<-function(model,burnin=burnin,thin=thin,lag=lag,type=type){
  
  l=length(model)
  
  if(l==1){
    HLSMfixedDiagnostics(model[[1]],burnin=burnin,thin=thin,lag=lag,type=type)
  }
  else if(l==2){
    chain1=model[[1]]$draws$Beta
    chain2=model[[2]]$draws$Beta
    combinedchains=mcmc.list(as.mcmc(as.data.frame(chain1)),as.mcmc(as.data.frame(chain2)))
    par(mar=c(1,1,1,1))
    gelman=plot(combinedchains)
    gelmanDiag=gelman.diag(combinedchains)
    if(gelmanDiag[[1]][[1]]<1.05){
      t=paste("The chains converge with less variance between chains and within chains as the scale reduction factor is less than 1.05, the value is :",round(gelmanDiag[[1]][[1]],2))
    }
    if(gelmanDiag[[1]][[1]]>1.05){
      t=paste("The chains have not converged as the scale reduction factor is greater than 1.5, the value is:",round(gelmanDiag[[1]][[1]],2))
    }
    
    return(c(gelman,t))
  }else if (l==3){
    chain1=model[[1]]$draws$Beta
    chain2=model[[2]]$draws$Beta
    chain3=model[[3]]$draws$Beta
    combinedchains=mcmc.list(as.mcmc(as.data.frame(chain1)),as.mcmc(as.data.frame(chain2)),as.mcmc(as.data.frame(chain3)))
    par(mar=c(1,1,1,1))
    gelman=plot(combinedchains)
    gelmanDiag=gelman.diag(combinedchains)
    if(gelmanDiag[[1]][[1]]<1.05){
      t=paste("The chains converge with less variance between chains and within chains as the scale reduction factor is less than 1.05, the value is :",round(gelmanDiag[[1]][[1]],2))
    }
    if(gelmanDiag[[1]][[1]]>1.05){
      t=paste("The chains have not converged as the scale reduction factor is greater than 1.5, the value is:",round(gelmanDiag[[1]][[1]],2))
    }
    return(c(gelman,t))
  }
  else if (l==4){
    chain1=model[[1]]$draws$Beta
    chain2=model[[2]]$draws$Beta
    chain3=model[[3]]$draws$Beta
    chain4=model[[4]]$draws$Beta
    combinedchains=mcmc.list(as.mcmc(as.data.frame(chain1)),as.mcmc(as.data.frame(chain2)),as.mcmc(as.data.frame(chain3)),as.mcmc(as.data.frame(chain4)))
    par(mar=c(1,1,1,1))
    gelman=plot(combinedchains)
    gelmanDiag=gelman.diag(combinedchains)
    if(gelmanDiag[[1]][[1]]<1.05){
      t=paste("The chains converge with less variance between chains and within chains as the scale reduction factor is less than 1.05, the value is :",round(gelmanDiag[[1]][[1]],2))
    }
    if(gelmanDiag[[1]][[1]]>1.05){
      t=paste("The chains have not converged as the scale reduction factor is greater than 1.5, the value is:",round(gelmanDiag[[1]][[1]],2))
    }
    return(c(gelman,t))
  }else if (l==5){
    chain1=model[[1]]$draws$Beta
    chain2=model[[2]]$draws$Beta
    chain3=model[[3]]$draws$Beta
    chain4=model[[4]]$draws$Beta
    chain5=model[[5]]$draws$Beta
    combinedchains=mcmc.list(as.mcmc(as.data.frame(chain1)),as.mcmc(as.data.frame(chain2)),as.mcmc(as.data.frame(chain3)),as.mcmc(as.data.frame(chain4)),as.mcmc(as.data.frame(chain5)))
    par(mar=c(1,1,1,1))
    gelman=plot(combinedchains)
    gelmanDiag=gelman.diag(combinedchains)
    if(gelmanDiag[[1]][[1]]<1.05){
      t=paste("The chains converge with less variance between chains and within chains as the scale reduction factor is less than 1.05, the value is :",round(gelmanDiag[[1]][[1]],2))
    }
    if(gelmanDiag[[1]][[1]]>1.05){
      t=paste("The chains have not converged as the scale reduction factor is greater than 1.5, the value is:",round(gelmanDiag[[1]][[1]],2))
    }
    return(c(gelman,t))
  }
  
  
}


## Helper function for the below Random effects diagnostics function
## This function can inform the user of the optimum thinning , burnin and chain length values when one random effects object is passed.
## This function will be called from the main random effects diagnostic function and hence we do not need to use this function explicitly 

HLSMrandomDiagnostics<-function(model,burnin=burnin,thin=thin,lag=lag,type=type,varnum=varnum){
  if(length(model$draws$ZZ)<=4000){
    t=paste("Looks like your chain length is less than the minimum required value length of 3746, try running a chain length greater than 4000")
    return(t)
  }
  else{
    intercept=getIntercept(model,burnin = burnin,thin = thin)
    
    
    beta_list <- list()
    for(i in 1:ncol(intercept))
    {
      beta_list[[i]]=model$draws$Beta[,varnum,i]
    }
    
    beta_list_df=do.call(cbind.data.frame, beta_list)
    beta_mcmc=as.mcmc(beta_list_df)
    intercept_mcmc=as.mcmc(intercept)
    ri=raftery.diag(intercept_mcmc)
    rb=raftery.diag( beta_list_df)
    rdf=data.frame(rb[[2]])
    chain_length=max(sapply(rdf[,2],max))
    burnin_value =rdf[which.max(rdf[,2]),1]
    thinnin_value=rdf[which.max(rdf[,2]),4]
    
    
    t=paste("Based on betas across all the networks,the optimum chain length is:",chain_length,"The optimum burnin value is:",burnin_value,".The optimum thinning value is :",thinnin_value)
    if(type=="traceplot"){
      par(mar=c(1,1,1,1))
      g=for(i in 1:ncol(intercept)){plot(model$draws$Beta[,varnum,i],type = "l")}
      
      return(c(g,t))
    }
    else if (type=="autocorr"){
      par(mar=c(1,1,1,1))
      
      a=for(i in 1:ncol(intercept)){acf(model$draws$Beta[,varnum,i], lag.max = lag)}
      return (c(a,t))
    }
    
    else if(type=="0"){
      return (c(t))}
  }
}


#### Random effects diagnostic function that calculates the optimal chain length, burnin and thinning value if only one random effects model is passed and can also compute the convergence if you pass two or more random effect objects

HLSMrandomEffectsDiagnostics<-function(model,burnin=burnin,thin=thin,lag=lag,type=type,varnum=varnum){
  
  l=length(model)
  
  if(l==1){
    HLSMrandomDiagnostics(model[[1]],burnin=burnin,thin=thin,lag=lag,type=type,varnum=varnum)
  }
  
  else if(l==2){
    beta_list1 <- list()
    beta_list2<-list()
    intercept=getIntercept(model[[1]])
    for(i in 1:ncol(intercept))
    {
      beta_list1[[i]]=model[[1]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list2[[i]]=model[[2]]$draws$Beta[,,i]
    }
    
    chain1=do.call(cbind.data.frame, beta_list1)
    chain2=do.call(cbind.data.frame, beta_list2)
    combinedchains=mcmc.list(as.mcmc(chain1),as.mcmc(chain2))
    gelmanDiag=gelman.diag(combinedchains)[[1]]
    g_df=data.frame(gelmanDiag)
    n=nrow(g_df)/15
    n_networks<-seq(from=1,to=15,by=1)
    g_df$network<- rep(n_networks,each=n)
    g_df$Point.est.=round(g_df$Point.est.,2)
    g_higher_text<-paste("The following betas of the following networks have PSRF(Potential scale reduction values greater than 1.05")
    g_higher=g_df[g_df$Point.est.>1.05,]
    num_betas=paste("Your random effect models have",n,"betas")
    g_df$betas<-rep(seq(1,n,1), times = n)
    point_estimate=setNames(aggregate(g_df$Point.est.~g_df$betas,FUN=function(x){
      round(mean(x),2)}),c("Beta Number","Mean PSRF"))
    
    CI=setNames(aggregate(g_df$Upper.C.I.~g_df$betas,FUN=function(x){
      round(mean(x),2)}),c("Beta Number","Mean Upper CI"))
    
    combined_table=cbind.data.frame(point_estimate,CI$`Mean Upper CI`)
    colnames(combined_table)=c("Beta Number","Mean PSRF","Mean Upper CI")
    if(max(point_estimate[,2])<1.05){
      t=paste("The chains converge with less variance between chains and within chains as the average scale reduction factor across all beta chains are less than 1.05")
    }else if(max(point_estimate[,2])>1.05)
    {
      t=paste("The chains do not converge as atleast one beta chain has an average scale reduction greater than 1.05")
    }
    return(list(num_betas,combined_table,t,g_higher_text,g_higher))
  }
  else if(l==3){
    beta_list1 <- list()
    beta_list2<-list()
    beta_list3<-list()
    intercept=getIntercept(model[[1]])
    for(i in 1:ncol(intercept))
    {
      beta_list1[[i]]=model[[1]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list2[[i]]=model[[2]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list3[[i]]=model[[3]]$draws$Beta[,,i]
    }
    
    chain1=do.call(cbind.data.frame, beta_list1)
    chain2=do.call(cbind.data.frame, beta_list2)
    chain3=do.call(cbind.data.frame, beta_list3)
    combinedchains=mcmc.list(as.mcmc(chain1),as.mcmc(chain2),as.mcmc(chain3))
    gelmanDiag=gelman.diag(combinedchains)[[1]]
    g_df=data.frame(gelmanDiag)
    n=nrow(g_df)/15
    n_networks<-seq(from=1,to=15,by=1)
    g_df$network<- rep(n_networks,each=n)
    g_df$Point.est.=round(g_df$Point.est.,2)
    g_higher_text<-paste("The following betas of the following networks have PSRF(Potential scale reduction values greater than 1.05")
    
    g_higher=g_df[g_df$Point.est.>1.05,]
    num_betas=paste("Your random effect models have",n,"betas")
    g_df$betas<-rep(seq(1,n,1), times = n)
    point_estimate=setNames(aggregate(g_df$Point.est.~g_df$betas,FUN=function(x){
      round(mean(x),2)}),c("Beta Number","Mean PSRF"))
    
    CI=setNames(aggregate(g_df$Upper.C.I.~g_df$betas,FUN=function(x){
      round(mean(x),2)}),c("Beta Number","Mean Upper CI"))
    
    combined_table=cbind.data.frame(point_estimate,CI$`Mean Upper CI`)
    colnames(combined_table)=c("Beta Number","Mean PSRF","Mean Upper CI")
    if(max(point_estimate[,2])<1.05){
      t=paste("The chains converge with less variance between chains and within chains as the average scale reduction factor across all beta chains are less than 1.05")
    }else if(max(point_estimate[,2])>1.05)
    {
      t=paste("The chains do not converge as atleast one beta chain has an average scale reduction greater than 1.05")
    }
    return(list(num_betas,combined_table,t,g_higher_text,g_higher))
  }
  else if(l==4){
    beta_list1 <- list()
    beta_list2<-list()
    beta_list3<-list()
    beta_list4<-list()
    intercept=getIntercept(model[[1]])
    for(i in 1:ncol(intercept))
    {
      beta_list1[[i]]=model[[1]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list2[[i]]=model[[2]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list3[[i]]=model[[3]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list4[[i]]=model[[4]]$draws$Beta[,,i]
    }
    
    chain1=do.call(cbind.data.frame, beta_list1)
    chain2=do.call(cbind.data.frame, beta_list2)
    chain3=do.call(cbind.data.frame, beta_list3)
    chain4=do.call(cbind.data.frame, beta_list4)
    combinedchains=mcmc.list(as.mcmc(chain1),as.mcmc(chain2),as.mcmc(chain3),as.mcmc(chain4))
    gelmanDiag=gelman.diag(combinedchains)[[1]]
    g_df=data.frame(gelmanDiag)
    n=nrow(g_df)/15
    n_networks<-seq(from=1,to=15,by=1)
    g_df$network<- rep(n_networks,each=n)
    g_df$Point.est.=round(g_df$Point.est.,2)
    g_higher_text<-paste("The following betas of the following networks have PSRF(Potential scale reduction values greater than 1.05")
    
    g_higher=g_df[g_df$Point.est.>1.05,]
    num_betas=paste("Your random effect models have",n,"betas")
    g_df$betas<-rep(seq(1,n,1), times = n)
    point_estimate=setNames(aggregate(g_df$Point.est.~g_df$betas,FUN=function(x){
      round(mean(x),2)}),c("Beta Number","Mean PSRF"))
    
    CI=setNames(aggregate(g_df$Upper.C.I.~g_df$betas,FUN=function(x){
      round(mean(x),2)}),c("Beta Number","Mean Upper CI"))
    
    combined_table=cbind.data.frame(point_estimate,CI$`Mean Upper CI`)
    colnames(combined_table)=c("Beta Number","Mean PSRF","Mean Upper CI")
    if(max(point_estimate[,2])<1.05){
      t=paste("The chains converge with less variance between chains and within chains as the average scale reduction factor across all beta chains are less than 1.05")
    }else if(max(point_estimate[,2])>1.05)
    {
      t=paste("The chains do not converge as atleast one beta chain has an average scale reduction greater than 1.05")
    }
    return(list(num_betas,combined_table,t,g_higher_text,g_higher))
  }
  else if(l==5){
    beta_list1 <- list()
    beta_list2<-list()
    beta_list3<-list()
    beta_list4<-list()
    beta_list5<-list()
    intercept=getIntercept(model[[1]])
    for(i in 1:ncol(intercept))
    {
      beta_list1[[i]]=model[[1]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list2[[i]]=model[[2]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list3[[i]]=model[[3]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list4[[i]]=model[[4]]$draws$Beta[,,i]
    }
    for(i in 1:ncol(intercept))
    {
      beta_list5[[i]]=model[[5]]$draws$Beta[,,i]
    }
    
    chain1=do.call(cbind.data.frame, beta_list1)
    chain2=do.call(cbind.data.frame, beta_list2)
    chain3=do.call(cbind.data.frame, beta_list3)
    chain4=do.call(cbind.data.frame, beta_list4)
    chain5=do.call(cbind.data.frame, beta_list5)
    combinedchains=mcmc.list(as.mcmc(chain1),as.mcmc(chain2),as.mcmc(chain3),as.mcmc(chain4),as.mcmc(chain5))
    gelmanDiag=gelman.diag(combinedchains)[[1]]
    g_df=data.frame(gelmanDiag)
    n=nrow(g_df)/15
    n_networks<-seq(from=1,to=15,by=1)
    g_df$network<- rep(n_networks,each=n)
    g_df$Point.est.=round(g_df$Point.est.,2)
    g_higher_text<-paste("The following betas of the following networks have PSRF(Potential scale reduction values greater than 1.05")
    g_higher=g_df[g_df$Point.est.>1.05,]
    num_betas=paste("Your random effect models have",n,"betas")
    g_df$betas<-rep(seq(1,n,1), times = n)
    point_estimate=setNames(aggregate(g_df$Point.est.~g_df$betas,FUN=function(x){
      round(mean(x),2)}),c("Beta Number","Mean PSRF"))
    CI=setNames(aggregate(g_df$Upper.C.I.~g_df$betas,FUN=function(x){
      round(mean(x),2)}),c("Beta Number","Mean Upper CI"))
    
    combined_table=cbind.data.frame(point_estimate,CI$`Mean Upper CI`)
    colnames(combined_table)=c("Beta Number","Mean PSRF","Mean Upper CI")
    if(max(point_estimate[,2])<1.05){
      t=paste("The chains converge with less variance between chains and within chains as the average scale reduction factor across all beta chains are less than 1.05")
    }else if(max(point_estimate[,2])>1.05)
    {
      t=paste("The chains do not converge as atleast one beta chain has an average scale reduction greater than 1.05")
    }
    return(list(num_betas,combined_table,t,g_higher_text,g_higher))
  }
  
}


