###Plots####

######"BOX" PLOTS########
###chain.list is a list MCMC output, either a vector for one parameter or matrix for several parameters###
get_HLSM_type <- function(object_list) {
  calls <- lapply(object_list, getCall)
  funcs <- sapply(lapply(calls, `[[`, 1), as.character)
  type <- unique(gsub('HLSM(.*)EF', '\\1', funcs))
  if (length(type) > 1) {
    stop("HLSM list must be all of the same type.")
  } else if (!(type %in% c('fixed', 'random', 'LSM'))) {
    stop("Unknown HLSM type found in object.")
  }
  if(type=="LSM"){
  	test=deparse(calls[[1]]) #creates string of the call
  	est.int=grep("estimate.intercept = TRUE", test)
  	if(length(est.int)>0){type="LSM.estInt"}else{type="LSM.fixedInt"}
  }
  return(type)
}


HLSMcovplots=function(fitted.model, burnin=0, thin=1){
	if(class(fitted.model)!="HLSM"){stop("Fitted Model must be of class HLSM")}
	type=get_HLSM_type(list(fitted.model))
	if(type=='random'){
		HLSMplot.random.fit(fitted.model, burnin, thin)
		}else{HLSMplot.fixed.fit(fitted.model,  burnin, thin)}
}

HLSMplot.random.fit=function(fitted.model, burnin = 0, thin = 1){
mycolors=c("navy", "orange", "darkmagenta")
	
	myInt=getIntercept(fitted.model, burnin=burnin, thin=thin)
	myBetas=getBeta(fitted.model,burnin=burnin,thin=thin)

	
	if(is.null(myInt)==TRUE & is.null(myBetas)==FALSE){
		draws=myBetas
		p=dim(draws)[2]
		myTitles=paste('X', 1:p, sep="")
		n.nets=dim(draws)[3]
		print("Covariates are in order of Edge, Sender, Receiver")}else if(is.null(myInt)==FALSE & is.null(myBetas)==TRUE){
	     p=1
	     temp=dim(myInt)
	     draws=array(0, dim=c(temp[1],1,temp[2]))
	     draws[,1,]=myInt	
	     myTitles="Intercept"
	     n.nets=dim(draws)[3]}else if(is.null(myInt)==FALSE & is.null(myBetas)==FALSE){
		temp=dim(myBetas)
		draws=array(data=0, dim=c(temp[1], temp[2]+1, temp[3]))
		draws[,1,]=myInt
		draws[,2:(temp[2]+1),]=myBetas
		p=dim(draws)[2]
		myTitles=c("Int", paste('X', 1:(p-1), sep=""))	
		n.nets=dim(draws)[3]
		print("Covariates are in order of Edge, Sender, Receiver")
	}else{stop("There are not any covariates to plot")}

x=1:n.nets
for(k in 1:length(myTitles)){ ### boxplots across all schools for each parameter ###
	int.quantiles=apply(draws[,k,], 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))   
plot(x, int.quantiles[3,], pch=".", 
			xlim=c(0.75, n.nets+0.25), ylim=c(min(int.quantiles), max(int.quantiles)),
			col="gray", xlab="Network", ylab=expression(beta), main=myTitles[k])	 
			          		
		rect((x-0.05), int.quantiles[1,], (x+0.05), int.quantiles[5,], col=mycolors[3], border=NA)
		rect((x-0.15), int.quantiles[2,], (x+0.15), int.quantiles[4,], col=mycolors[2], border=NA)
		segments((x-0.15), int.quantiles[3,], (x+0.15), int.quantiles[3,], lwd=2, col=mycolors[1])
}
 }               


###############

HLSMplot.fixed.fit=function(fitted.model, burnin =0, thin = 1){
	mycolors=c("navy", "orange", "darkmagenta")
	
	myInt=getIntercept(fitted.model, burnin=burnin, thin=thin)
	myBetas=getBeta(fitted.model,burnin=burnin,thin=thin)
	
	if(is.null(myInt)==TRUE & is.null(myBetas)==TRUE){stop("There are not any covariates to plot")}
	
	if(is.null(myInt)==TRUE & is.null(myBetas)==FALSE){
		draws=myBetas
		p=dim(draws)[2]
		myXlabels=paste('X', 1:p, sep="")
		mytitle="Regression Coefficients"	
		print("Covariates are in order of Edge, Sender, Receiver")
	}
	if(is.null(myInt)==FALSE & is.null(myBetas)==TRUE){
	     p=1
	     draws=myInt	
	     myXlabels=c("")
	     mytitle=c("Intercept")
	}
	if(is.null(myInt)==FALSE && is.null(myBetas)==FALSE){
		draws=cbind(myInt,myBetas)
		p=dim(draws)[2]
		myXlabels=c("Int", paste('X', 1:(p-1), sep=""))	
		mytitle="Regression Coefficients"
		print("Covariates are in order of Edge, Sender, Receiver")
	}

x=1:p
	
if(is.null(myInt)==FALSE & is.null(myBetas)==TRUE){int.quantiles=quantile(draws, c(0.025, 0.25, 0.5, 0.75, 0.975))
plot(x, int.quantiles[3], pch=".", xlim=c(0.75, p+0.25), ylim=c(min(int.quantiles), max(int.quantiles)),
			col="white", ylab="", xlab="", xaxt="n", main=mytitle)
	axis(1, at=x, labels=myXlabels)		           		
		rect((x-0.05), int.quantiles[1], (x+0.05), int.quantiles[5], col=mycolors[3], border=NA)
		rect((x-0.15), int.quantiles[2], (x+0.15), int.quantiles[4], col=mycolors[2], border=NA)
		segments((x-0.15), int.quantiles[3], (x+0.15), int.quantiles[3], lwd=2, col=mycolors[1])
}else{
int.quantiles=apply(draws, 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))   

plot(x, int.quantiles[3,], pch=".", xlim=c(0.75, p+0.25), ylim=c(min(int.quantiles), max(int.quantiles)),
			col="white", ylab=expression(beta), xlab="", xaxt="n", main=mytitle)
	axis(1, at=x, labels=myXlabels)		           		
		rect((x-0.05), int.quantiles[1,], (x+0.05), int.quantiles[5,], col=mycolors[3], border=NA)
		rect((x-0.15), int.quantiles[2,], (x+0.15), int.quantiles[4,], col=mycolors[2], border=NA)
		segments((x-0.15), int.quantiles[3,], (x+0.15), int.quantiles[3,], lwd=2, col=mycolors[1])
}
}





######LATENT SPACE POSTIONS########


HLSMplotLS = function(LS,xx,fitted.model, node.name = FALSE,nodenames = NULL){
	xcor = apply(LS[[xx]][,1,],1,mean)
	ycor = apply(LS[[xx]][,2,],1,mean)
	plot(xcor,ycor,main = paste('Network',xx),cex = 0.5,pch = 19,cex.lab = 1.5,col = 'red')
	if(node.name == TRUE){
		if(!is.null(nodenames)){
			text(xcor,ycor,lab = nodenames) }
	}
}

HLSMplot.fit.LS = function(fitted.model,pdfname = NULL,burnin=0,thin=1,...){
	if(class(fitted.model)!='HLSM'){
		stop('fitted.model must be of class HLSM')}
#	niter = length(fitted.model$draws$Alpha)
#	if(is.null(burnin)){burnin=0.1*niter}
#	if(is.null(thin)){thin = round(0.9*niter/200)}
	LS = getLS(fitted.model,burnin=burnin, thin=thin)
	if(!is.null(pdfname)){pdf(file=paste(pdfname,'.pdf',sep=''))}else(dev.new(height = 10,width = 10))
	if(length(LS) > 5){
		mat = matrix(1:length(LS),ceiling(length(LS)/5),5,byrow = TRUE)
		layout(mat,widths = rep.int(2.5,ncol(mat)), heights = rep.int(1,nrow(mat)))
	}else{
		x1 = ceiling(length(LS)/2)
		par(mfrow = c(2,x1))
	}
	lapply(1:length(LS),function(x) HLSMplotLS(LS,x,fitted.model,...))
	if(!is.null(pdfname))dev.off()
}







        




    
