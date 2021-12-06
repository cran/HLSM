\name{HLSMcovplots}
\alias{HLSMplotLikelihood}
\alias{HLSMplot.random.fit}
\alias{HLSMplot.fixed.fit}
\alias{HLSMplot.fit.LS}
\alias{HLSMcovplots}


\title{Plotting functions for HLSM objects}

\description{ 
 Functions for plotting HLSM/LSM model fits of class 'HLSM'. HSLMcovplots is the most recent function to plot posterior distribution summaries. HLSMplotLikelihood( ) plots the likelihood, HLSMcovplots( ) summarizes posterior draws of the parameters from MCMC sample, and HLSMplot.fit.LS( ) is for plotting the mean latent position estimates.
}

\usage{
	HLSMplotLikelihood(object, burnin = 0, thin = 1)
	HLSMcovplots(fitted.model, burnin = 0, thin = 1)
}


\arguments{

	\item{object}{object of class 'HLSM' obtained as an output from \code{LSM}, \code{HLSMrandomEF()} or \code{HLSMfixedEF()}
}
	\item{fitted.model}{model fit from LSM(), HLSMrandomEF() or HLSMfixedEF()}

\item{burnin}{numeric value to burn the chain for plotting the results from the 'HLSM'object }
	
	\item{thin}{a numeric thinning value}
}


\value{
 returns plot objects.
}


\author{Sam Adhikari & Tracy Sweet}

\examples{
#using advice seeking network of teachers in 15 schools
#to fit the data

#Random effect model#
priors = NULL
tune = NULL
initialVals = NULL
niter = 10

random.fit = HLSMrandomEF(Y = ps.advice.mat,FullX = ps.edge.vars.mat,
	initialVals = initialVals,priors = priors,
	tune = tune,tuneIn = FALSE,dd = 2,niter = niter)


HLSMcovplots(random.fit)
}
