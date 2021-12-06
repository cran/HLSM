#######################################
##Inputs are all in the form of arrays.
MCMCfunction = function(nn,PP,KK,dd,XX,YY,ZZ,beta,intercept,MuBeta,SigmaBeta,MuZ,VarZ,tuneBetaAll,tuneInt,tuneZAll,niter,PriorA,PriorB){
    accBetaAll = rep(0,KK*PP) 
    accIntAll = rep(0,KK)
    accZAll = rep(0,sum(nn))
    betaFinal = rep(0,niter*length(beta))
    ZZFinal = rep(0,niter*length(ZZ))
    InterceptFinal = rep(0, niter*length(intercept))
    Zvar1 = Zvar2= rep(0, niter)
    postVar = postMu = rep(0,niter*(PP+1))
    likelihood = rep(0,niter)

    out = .C('sampleRandomIntervention',as.integer(niter),
             as.double(XX),as.double(YY),as.double(ZZ),
             as.integer(nn),as.integer(PP),as.integer(dd),as.integer(KK),
             as.double(beta),as.double(intercept),
             as.double(MuBeta),as.double(SigmaBeta), as.double(MuZ),as.double(VarZ), as.double(tuneBetaAll),as.double(tuneInt),as.double(tuneZAll),
             as.double(accBetaAll),as.double(accIntAll),as.double(accZAll), as.double(betaFinal),as.double(ZZFinal),as.double(InterceptFinal), as.double(Zvar1),as.double(Zvar2), as.double(postVar),as.double(postMu), as.double(likelihood),as.double(PriorA), as.double(PriorB))

    ## 1, 2-4, 5-8
    ## 9-10, 11-14
    ## 15-17, 18-20
    ## 21-23, 24-27, 28-30

    betaFinal = array(out[[21]],dim = c(niter,PP,KK ))
    ZZFinal = list()
    accZ = list()
    ZZx = array(out[[22]],dim = c((sum(nn)*dd),niter))
    for(ni in 1:niter ){
	Zsm = list()
        n00 = n0 = 1	
    for(kk in 1:KK ){
        n1 = sum(nn[1:kk])
        Zsm[[kk]] = array(ZZx[n0:(2*n1),ni],dim = c(nn[kk],dd))
        accZ[[kk]] = out[[20]][n00:(n1)]/(niter)
        n00 = (n1) + 1
	n0 = (2*n1) + 1
}
	ZZFinal[[ni]] = Zsm
}


#    n00 = n0 = 1
#    ZZx = array(out[[29]],dim = c(sum(nn),dd,niter))
#    for(kk in 1:KK){
#        n1 = sum(nn[1:kk])
#        ZZFinal[[kk]] = array(ZZx[n00:n1,,],dim = c(nn[kk],dd,niter))
#        accZ[[kk]] = out[[26]][n00:(n1)]/(niter)
#        n00 = (n1) + 1
#}
    InterceptFinal = array(out[[23]],dim=c(niter,KK) )
    Zvar1 = out[[24]]
    Zvar2 = out[[25]]
    postVarBeta = array(out[[26]],dim = c(niter,PP+1))
    postMuBeta = array(out[[27]],dim = c(niter,PP+1))
    likelihood = out[[28]]
    accbeta = array(out[[18]],dim = c(PP,KK))/(niter)
    accint = out[[19]]/(niter)
    Zvar = data.frame(Zvar1,Zvar2)
    draws = list(Intercept = InterceptFinal, Beta= betaFinal, ZZ = ZZFinal, Zvar = Zvar,BetaMu = postMuBeta, VarMu = postVarBeta, likelihood = likelihood)
    accrate = list(intercept = accint, beta = accbeta,Z=accZ)

    return(list(draws = draws, acc = accrate)) 
  #  return(out)  
    rm(out)

}

    
 
