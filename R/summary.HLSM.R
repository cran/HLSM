 summary.HLSM = function(object,...){
	if(is(object, 'HLSM')==FALSE){stop('input must be of class HLSM')}
warning("Run HLSMdiag to assess convergence and thinning before reporting")
	call = object$call
type= get_HLSM_type(list(object))

	Betas = getBeta(object, ...)
	if(!is.null(Betas)){
	    if(length(dim(Betas)) == 3){
		est.slopes = lapply(1:dim(Betas)[[3]],
			function(x){ t(sapply(1:ncol(Betas[,,x]),
			function(y) data.frame(min = round(min(Betas[,y,x]),3),
				max = round(max(Betas[,y,x]),3),
				est.mean = round(mean(Betas[,y,x]),3), 
				sd = round(sd(Betas[,y,x]),3), 
				q.25 = round(quantile(Betas[,y,x],0.025),3),
				q.975 = round(quantile(Betas[,y,x],0.975),3))))})
	}else{
		Betas = as.matrix(Betas)
		est.slopes = t(sapply(1:ncol(Betas),function(y) data.frame(min = round(min(Betas[,y]),3), max = round(max(Betas[,y]),3),est.mean = round(mean(Betas[,y]),3), sd = round(sd(Betas[,y]),3), q.25 = round(quantile(Betas[,y],0.025),3), q.975 = round(quantile(Betas[,y],0.975),3)) ) )
	}
	}else(est.slopes = NA)

if(type!='LSM.estInt'){
	Intercept = getIntercept(object,...)
	if(is(Intercept[1], 'matrix')){
		est.intercept = t(sapply(1:ncol(Intercept), function(y) data.frame(min = round(min(Intercept[,y]),3), max = round(max(Intercept[,y]),3), est.mean = round(mean(Intercept[,y]),3), sd = round(sd(Intercept[,y]),3), q.25 = round(quantile(Intercept[,y],0.025),3), q.975 = round(quantile(Intercept[,y],0.975),3) )) )
	}else{
		est.intercept = data.frame(min = round(min(Intercept),3), max = round(max(Intercept),3), est.mean = round(mean(Intercept),3), sd = round(sd(Intercept),3), q.025 = round(quantile(Intercept,0.025),3), q.975 = round(quantile(Intercept, 0.975),3) )
      }
	
	res =  list(call = object$call,est.intercept = est.intercept, est.slopes = est.slopes)
	}else{res =  list(call = object$call,est.slopes = est.slopes)}
		

    class(res) = 'summary.HLSM'
    res
}

