
MCMCfixedEFfixedIntercept = function(nn,PP,KK,dd,XX,YY,ZZ,
		beta,intercept,MuBeta,SigmaBeta,MuZ,VarZ,
		tuneBetaAll,tuneInt,tuneZAll,
		niter,PriorA,PriorB)
{
    accBetaAll = rep(0,PP) 
    accIntAll = 0
    accZAll = rep(0,sum(nn))
    betaFinal = rep(0,niter*length(beta))
    ZZFinal = rep(0,niter*length(ZZ))
    InterceptFinal = rep(0, niter*length(intercept))
    Zvar1 = Zvar2= rep(0, niter)
    likelihood = rep(0,niter)

    out = .C('sampleFixedEFFixedIntercept',as.integer(niter), 
             as.double(XX), as.double(YY), as.double(ZZ), 
             as.integer(nn), as.integer(PP), as.integer(dd), as.integer(KK),
             as.double(beta), as.double(intercept),
             as.double(MuBeta),as.double(SigmaBeta),
 			 as.double(MuZ), as.double(VarZ), as.double(tuneBetaAll), 
             as.double(tuneZAll), as.double(accBetaAll), 
             as.double(accZAll), as.double(betaFinal), as.double(ZZFinal), 
             as.double(Zvar1), as.double(Zvar2),as.double(likelihood),
             as.double(PriorA),as.double(PriorB))

    ## 1, 
    ## 2-4, XX, YY, ZZ
    ## 5-8, nn, PP, dd, KK
    ## 9-10, beta, intercept
    ## 11-12, muBeta, SigmaBeta
    ## 13-15, 16-17,  muZ, VarZ, tuneBetaAll,tuneZall, accBetaAll
    ## 18-20 accZAll, betaFinal, ZZFinal
    ## 21 - 23 Zvar1, Zvar2, likelihood
    ## 24,25 PriorA, PriorB


    betaFinal = array(out[[19]],dim = c(niter,PP,1))  ##betafinal
    ##InterceptFinal = array(out[[23]],dim=c(niter,1) )

    ZZFinal = list()
    accZ = list()
    ZZx = array(out[[20]],dim = c((sum(nn)*dd),niter)) ##ZZFinal
    for(ni in 1:niter ){
	Zsm = list()
        n00 = n0 = 1	
    for(kk in 1:KK ){
        n1 = sum(nn[1:kk] )
        Zsm[[kk]] = array(ZZx[n0:(2*n1),ni],dim = c(nn[kk],dd))
        accZ[[kk]] = out[[18]][n00:(n1)]/(niter) ##accZAll
        n00 = (n1) + 1
	n0 = (2*n1) + 1
}
	ZZFinal[[ni]] = Zsm
}


    Zvar1 = out[[21]] ## Zvar1
    Zvar2 = out[[22]] ##Zvar2
    Zvar = data.frame(Zvar1,Zvar2)
    likelihood = out[[23]] ##likelihood

    accbeta = out[[17]]/(niter) ##accBetaAll
   
    
    draws = list( Beta= betaFinal, ZZ = ZZFinal, Zvar = Zvar,likelihood = likelihood)

    accrate = list(beta = accbeta, Z=accZ)

    return(list(draws = draws, acc = accrate)) 
    rm(out)
 }   
    