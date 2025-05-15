
####Function to run the sampler#####
##SPECIFYING PRIORS######
#########################
##if priors = NULL => uses randomly generated priors
##else: priors is a list of following objects:
        ##MuBeta:= prior mean for betas & intercepts; 
        ##SigmaBeta:= prior variance for betas & intercept;
        ##MuZ; VarZ;
        ##PriorA; PriorB 
##TUNING PARAMETERS#######
##########################
##if tune = NULL => uses auto tuning
##else: tune is a list of following objects:
       ##tuneBeta = vec, PP
       ##tuneInt = vec, len = 1
       ##tuneZ =  list( vec(len = nn[x]])) length of list = KK
##############################################################          
#library(MASS)
LSM= function(Y,edgeCov = NULL, receiverCov = NULL,senderCov =NULL,
     FullX = NULL, initialVals = NULL, priors = NULL, tune = NULL,
        tuneIn = TRUE,dd=2,estimate.intercept=FALSE, niter, verbose=TRUE)
{

    ## Y should be matrix or a dataframe
    if(is(Y, 'matrix')){
        KK = 1
        diag(Y)=0
	if(dim(Y)[1] == dim(Y)[2]){
            Y=list(Y)
            nn =sapply(1:length(Y),function(x) nrow(Y[[x]]))
            nodenames=lapply(1:length(Y), function(x) rownames(Y[[x]]))
            if(sum(sapply(nodenames, is.null))>=1){nodenames=lapply(1:length(Y), function(x) 1:dim(Y[[x]])[1])}}else{stop('Y must be a nxn matrix or dataframe')}}

	
    if(is(Y, 'dataframe')){ 
	if(dim(Y)[1] != dim(Y)[2] & dim(Y)[2] == 3){
            nn = length(unique(c(Y$Receiver,Y$Sender)))
            Y$id=rep(1, nn)
            Y=list(Y)
	nodenames = list(unique(c(Y$Receiver,Y$Sender)))
	}	}


##prepare covariates#####
#########################
	noCOV = FALSE
	if(!is.null(FullX) & !is.null(edgeCov) &!is.null(receiverCov) & !is.null(senderCov))(stop('FullX cannot be used when nodal or edge covariates are provided'))

	if(is.null(FullX) & is.null(edgeCov) & is.null(receiverCov) & is.null(senderCov)){
		X = lapply(1:KK,function(x) array(0, dim = c(nn[x],nn[x],1)))
		noCOV = TRUE
	}

	if(is.null(FullX)){
	if(!is.null(edgeCov) | !is.null(senderCov)| !is.null(receiverCov)){
            if(!is.null(edgeCov)){
		if(is(edgeCov,'data.frame')==FALSE){
                    stop('edgeCov must be of class data.frame')}
                temp.n=dim(edgeCov)[1]
                edgeCov$id=rep(1, temp.n)
		X1 = getEdgeCov(edgeCov, nn,nodenames)
}else(X1 =NULL)
  	  if(!is.null(senderCov)){
		if(is(senderCov, 'data.frame')==FALSE){
                    stop('senderCov must be of class data.frame')}
                temp.n=dim(senderCov)[1]
                senderCov$id=rep(1, temp.n)
		X2 = getSenderCov(senderCov, nn,nodenames)
}else(X2 = NULL)


	  if(!is.null(receiverCov)){
		if(is(receiverCov,'data.frame')==FALSE){
                    stop('receiverCov must be of class data.frame')}
                temp.n=dim(receiverCov)[1]
                receiverCov$id=rep(1, temp.n)
		X3 = getReceiverCov(receiverCov, nn,nodenames)
}else(X3 = NULL)	

	X = lapply(1:KK, function(x){if(!is.null(X1)&!is.null(X2)&!is.null(X3)){
		ncov = dim(X1[[x]])[3]+dim(X2[[x]])[3]+dim(X3[[x]])[3];
		df = array(0, dim = c(nn[x],nn[x],ncov));
		df[,,1:dim(X1[[x]])[3]] = X1[[x]];
		df[,,(dim(X1[[x]])[3]+1):(dim(X1[[x]])[3]+dim(X2[[x]])[3])] = X2[[x]];
		df[,,(dim(X1[[x]])[3]+dim(X2[[x]])[3]+1):(dim(X1[[x]])[3]+dim(X2[[x]])[3]+dim(X3[[x]])[3])] = X3[[x]] };
		if(!is.null(X1)&!is.null(X2) & is.null(X3)){
			ncov = dim(X1[[x]])[3]+dim(X2[[x]])[3];
			df = array(0, dim = c(nn[x],nn[x],ncov));
			df[,,1:dim(X1[[x]])[3]] = X1[[x]];
			df[,,(dim(X1[[x]])[3]+1):(dim(X1[[x]])[3]+dim(X2[[x]])[3])] = X2[[x]]};
		if(!is.null(X1)&!is.null(X3)&is.null(X2)){
			ncov = dim(X1[[x]])[3]+dim(X3[[x]])[3];
			df = array(0, dim = c(nn[x],nn[x],ncov));
			df[,,1:dim(X1[[x]])[3]] = X1[[x]];
			df[,,(dim(X1[[x]])[3]+1):(dim(X1[[x]])[3]+dim(X3[[x]])[3])] = X3[[x]]};
	if(!is.null(X2)&!is.null(X3)&is.null(X1)){
			ncov = dim(X2[[x]])[3]+dim(X3[[x]])[3];
			df = array(0, dim = c(nn[x],nn[x],ncov));
			df[,,1:dim(X2[[x]])[3]] = X2[[x]];
			df[,,(dim(X2[[x]])[3]+1):(dim(X2[[x]])[3]+dim(X3[[x]])[3])] = X3[[x]]};
	if(!is.null(X1)& is.null(X2)& is.null(X3)){
			df = X1[[x]] };
	if(is.null(X1)& !is.null(X2)& is.null(X3)){
			df = X2[[x]] };
	if(is.null(X1)& is.null(X2)& !is.null(X3)){
			df = X3[[x]] };
	return(df) } )
}
}
	if(!is.null(FullX)) X = FullX

if(is(X, 'array')){X=list(X)}
#    nn = sapply(1:length(X),function(x) nrow(X[[x]]))
#    KK = length(X)
    PP = dim(X[[1]])[3]
    XX = unlist(X)
    YY = unlist(Y)
    YY[which(is.na(YY))] = 0  ## Note that missing data is imputed to be 0
    XX[which(is.na(XX))] = 0 ## Note thatmissing data is imputed to be 0
    #Priors

    if(is.null(priors)){
	MuBeta= rep(0,(PP+1)) 
	VarBeta = rep(1,(PP+1)) 
        MuZ = c(0,0)
        VarZ = c(20,20)
        PriorA = 5
        PriorB = 200
     }else{
	if(is(priors,'list')==FALSE)(stop("priors must be of class list, if not NULL"))
	MuBeta = priors$MuBeta
	VarBeta = priors$VarBeta
	MuZ = priors$MuZ
	VarZ = priors$VarZ
	PriorA = priors$PriorA
	PriorB = priors$PriorB
  }
##starting values
##
    ##Procrustean transformation of latent positions
    C = lapply(1:KK,function(tt){
        diag(nn[tt]) - (1/nn[tt]) * array(1, dim = c(nn[tt],nn[tt]))})
    Z0 = lapply(1:KK,function(tt){
        g = graph_from_adjacency_matrix(Y[[tt]]);
        ss = distances(g);
        ss[ss > 4] = 4;
        Z0 = cmdscale(ss, k = dd);
        dimnames(Z0)[[1]] = dimnames(YY[[tt]])[[1]];
        return(Z0)})	
    Z00 = lapply(1:KK,function(tt)C[[tt]]%*%Z0[[tt]])
    

    
    
     if(is.null(initialVals)){
        Z0 = unlist(Z00)
         beta0 = rnorm(PP,0,1)
         if(estimate.intercept==FALSE){
                 intercept0=0}
         else if(is(estimate.intercept, 'logical')==FALSE){intercept0=estimate.intercept}
         else{intercept0 = rnorm(1, 0, 1)}
             
         
        if(verbose){message("Starting Values Set")}
    }else{
	if(is(initialVals, 'list')==FALSE)(stop("initialVals must be of class list, if not NULL"))
	Z0 = initialVals$ZZ
	beta0 = initialVals$beta
	intercept0 = initialVals$intercept
	}

###tuning parameters#####
    if(is.null(tune)){
            a.number = 5
            tuneBeta = rep(1,PP)
            tuneInt = 0.2
            tuneZ =  lapply(1:KK,function(x) rep(1.2,nn[x]))          
            } else{
         	if(is(tune, 'list')==FALSE)(stop("tune must be of class list, if not NULL"))
                 a.number = 1
                 tuneBeta = tune$tuneBeta
                 if(estimate.intercept==TRUE){tuneInt = tune$tuneInt}
                 tuneZ = tune$tuneZ
          }       

    
###Tuning and Running the Sampler####
#### Estimating Intercept = TRUE!########

if(estimate.intercept==TRUE){
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
    while(do.again ==1){
        if(verbose){message('Tuning the Sampler')}
        for(counter in 1:a.number)
{ 
            rslt = MCMCfixedEF(nn=nn,PP=PP,KK=KK,dd=dd,XX = XX,YY = YY,ZZ = Z0,
		beta = beta0 ,intercept = intercept0,
		MuBeta = MuBeta,SigmaBeta = VarBeta,MuZ = MuZ,VarZ = VarZ,tuneBetaAll = tuneBeta, tuneInt = tuneInt,tuneZAll = unlist(tuneZ),niter = 200,PriorA = PriorA, PriorB = PriorB)
     	    tuneZ = lapply(1:KK,function(x)adjust.my.tune(tuneZ[[x]], rslt$acc$Z[[x]], 2))
            tuneBeta = adjust.my.tune(tuneBeta, rslt$acc$beta,1)
            tuneInt =  adjust.my.tune(tuneInt,rslt$acc$intercept,1)
            tuneX = tuneX+1
    }
    extreme = lapply(1:KK,function(x)which.suck(rslt$acc$Z[[x]],2))
    do.again = max(sapply(extreme, length)) > max(nn)

}
    if(verbose){message("Tuning is finished")}  
}

    rslt = MCMCfixedEF(nn=nn,PP=PP,KK=KK,dd=dd,XX = XX,YY = YY,ZZ = Z0,
		beta = beta0 ,intercept = intercept0, MuBeta = MuBeta,SigmaBeta = VarBeta,MuZ = MuZ,VarZ = VarZ,tuneBetaAll = tuneBeta, tuneInt = tuneInt, tuneZAll = unlist(tuneZ),niter = niter,PriorA = PriorA, PriorB = PriorB)
}



#### Estimating Intercept = FALSE!########
##default intercept is 0##
if(estimate.intercept==FALSE){
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
    while(do.again ==1){
        if(verbose){message('Tuning the Sampler')}
        for(counter in 1:a.number)
{ 
            rslt = MCMCfixedIntercept(nn=nn,PP=PP,KK=KK,dd=dd,XX = XX,YY = YY,ZZ = Z0,
		beta = beta0 ,intercept = intercept0,
		MuBeta = MuBeta,SigmaBeta = VarBeta,MuZ = MuZ,VarZ = VarZ,tuneBetaAll = tuneBeta, tuneInt = tuneInt,tuneZAll = unlist(tuneZ),niter = 200,PriorA = PriorA, PriorB = PriorB)
     	    tuneZ = lapply(1:KK,function(x)adjust.my.tune(tuneZ[[x]], rslt$acc$Z[[x]], 2))
            tuneBeta = adjust.my.tune(tuneBeta, rslt$acc$beta,1)
            tuneX = tuneX+1
    }
    extreme = lapply(1:KK,function(x)which.suck(rslt$acc$Z[[x]],2))
    do.again = max(sapply(extreme, length)) > max(nn)

}
    if(verbose){message("Tuning is finished") } 
}

    rslt = MCMCfixedIntercept(nn=nn,PP=PP,KK=KK,dd=dd,XX = XX,YY = YY,ZZ = Z0,
		beta = beta0 ,intercept = intercept0, MuBeta = MuBeta,SigmaBeta = VarBeta,MuZ = MuZ,VarZ = VarZ,tuneBetaAll = tuneBeta, tuneInt = tuneInt, tuneZAll = unlist(tuneZ),niter = niter,PriorA = PriorA, PriorB = PriorB)
}


if(is(estimate.intercept,'numeric')){
	intercept0=estimate.intercept
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
    while(do.again ==1){
        if(verbose){message('Tuning the Sampler')}
        for(counter in 1:a.number)
{ 
            rslt = MCMCfixedIntercept(nn=nn,PP=PP,KK=KK,dd=dd,XX = XX,YY = YY,ZZ = Z0,
		beta = beta0 ,intercept = intercept0,
		MuBeta = MuBeta,SigmaBeta = VarBeta,MuZ = MuZ,VarZ = VarZ,tuneBetaAll = tuneBeta, tuneInt = tuneInt,tuneZAll = unlist(tuneZ),niter = 200,PriorA = PriorA, PriorB = PriorB)
     	    tuneZ = lapply(1:KK,function(x)adjust.my.tune(tuneZ[[x]], rslt$acc$Z[[x]], 2))
            tuneBeta = adjust.my.tune(tuneBeta, rslt$acc$beta,1)
            tuneX = tuneX+1
    }
    extreme = lapply(1:KK,function(x)which.suck(rslt$acc$Z[[x]],2))
    do.again = max(sapply(extreme, length)) > max(nn)

}
    if(verbose){message("Tuning is finished") } 
}

    rslt = MCMCfixedIntercept(nn=nn,PP=PP,KK=KK,dd=dd,XX = XX,YY = YY,ZZ = Z0,
		beta = beta0 ,intercept = intercept0, MuBeta = MuBeta,SigmaBeta = VarBeta,MuZ = MuZ,VarZ = VarZ,tuneBetaAll = tuneBeta, tuneInt = tuneInt, tuneZAll = unlist(tuneZ),niter = niter,PriorA = PriorA, PriorB = PriorB)
}
    

    
    
    ##Procrutes transformation on the final draws of the latent positions
    ##
    Ztransformed = lapply(1:niter, function(ii) {lapply(1:KK,
                                                        function(tt){z= rslt$draws$ZZ[[ii]][[tt]];
                                                        z = C[[tt]]%*%z;
                                                        pr = t(Z00[[tt]])%*% z;
                                                        ssZ = svd(pr)
                                                        tx = ssZ$v%*%t(ssZ$u)
                                                        zfinal = z%*%tx
                                                        return(zfinal)})})    
    rslt$draws$ZZ = Ztransformed
    

    rslt$call = match.call()
    if(noCOV == TRUE ){
		rslt$tune = list(tuneZ = tuneZ, tuneInt = tuneInt)	
		rslt$draws$Beta = NA
    }

    if(noCOV == FALSE){
	    rslt$tune = list(tuneBeta = tuneBeta, tuneZ = tuneZ,tuneInt = tuneInt)
	
}
	
    class(rslt) = 'HLSM'
    rslt
}


    



