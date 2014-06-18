\name{plotDiagnostic}
\alias{plotDiagnostic}
\alias{plotLikelihood}

\title{built-in plot functions for HLSM object}

\description{ 
	Functions to plot the likelihood and diagnostic-plot of posterior draws of the parameters from MCMC sample. 
}

\usage{
	plotLikelihood(object,burnin = 0, thin = 1)
	plotDiagnostic(chain)
}


\arguments{

	\item{object}{object of class 'HLSM' obtained as an output from \code{HLSM()} or \code{HLSMfixedEF()}
}

	\item{burnin}{burnin value if want to plot burned chain}

	\item{thin}{a numeric thinning value}

	\item{chain}{a numeric vector of posterior draws of parameter of interest.}
}


\value{
 returns plot objects.
}


\author{Sam Adhikari}

\examples{
#using advice seeking network of teachers in 15 schools
#to fit the data

#Random effect model#
priors = NULL
tune = NULL
initialVals = NULL
niter = 10

random.fit = HLSM(X = ps.edge.vars.mat,Y = ps.advice.mat,
	initialVals = initialVals,priors = priors,
	tune = tune,tuneIn = FALSE,dd = 2,niter = niter,
	intervention = 0)

plotLikelihood(random.fit)

intercept = getIntercept(random.fit)
dim(intercept) ##is an array of dimension niter by 15
plotDiagnostic(intercept[,1])

##look at the diagnostic plot of intercept for the first school


}
