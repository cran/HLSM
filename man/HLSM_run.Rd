\name{HLSMrandomEF}
\alias{HLSMrandomEF}
\alias{HLSMfixedEF}
\alias{LSM}
\alias{print.HLSM}
\alias{print.summary.HLSM}
\alias{summary.HLSM}
\alias{getIntercept}
\alias{getLS}
\alias{getLikelihood}
\alias{getBeta}

\title{Function to run the MCMC sampler in random effects latent space model, HLSMfixedEF for fixed effects model, or LSM for single network latent space model)
}

\description{
    Function to run the MCMC sampler to draw from the posterior distribution of intercept, slopes, and latent positions. HLSMrandomEF( ) fits random effects model; HLSMfixedEF( ) fits fixed effects model; LSM( ) fits single network model. 
}

\usage{

HLSMrandomEF(Y,edgeCov=NULL, receiverCov = NULL, senderCov = NULL, 
        FullX = NULL,initialVals = NULL, priors = NULL, tune = NULL,
        tuneIn = TRUE,dd=2, niter, verbose=TRUE)

HLSMfixedEF(Y,edgeCov=NULL, receiverCov = NULL, senderCov = NULL,
        FullX = NULL, initialVals = NULL, priors = NULL, tune = NULL,
        tuneIn = TRUE,dd=2, estimate.intercept=FALSE, niter, verbose=TRUE)
        
LSM(Y,edgeCov=NULL, receiverCov = NULL, senderCov = NULL, 
        FullX = NULL,initialVals = NULL, priors = NULL, tune = NULL,
        tuneIn = TRUE,dd=2, estimate.intercept=FALSE, niter, verbose=TRUE)

getBeta(object, burnin = 0, thin = 1)
getIntercept(object, burnin = 0, thin = 1)
getLS(object, burnin  = 0, thin = 1)
getLikelihood(object, burnin = 0, thin = 1)
}


\arguments{

\item{Y}{
	input outcome for different networks. Y can either be

        (i). list of sociomatrices for \code{K} different networks (Y[[i]] must be a matrix with named rows and columns)

        (ii). list of data frame with columns \code{Sender}, \code{Receiver} and \code{Outcome} for \code{K} different networks

        (iii). a dataframe with columns named as follows: \code{id} to identify network, \code{Receiver} for receiver nodes, \code{Sender} for sender nodes and finally, \code{Outcome} for the edge outcome.     
        Note that for LSM, Y must be an adjacency matrix.
}

\item{edgeCov}{
	data frame to specify edge level covariates with

        (i). a column for network id named \code{id}, 

        (ii). a column for sender node named \code{Sender},

        (iii). a column for receiver nodes named \code{Receiver}, and 

        (iv). columns for values of each edge level covariates.
}

\item{receiverCov}{
	a data frame to specify nodal covariates as edge receivers with

        (i.) a column for network id named \code{id}, 

        (ii.) a column \code{Node} for node names, and 

        (iii). the rest for respective node level covariates.
}

\item{senderCov}{
	a data frame to specify nodal covariates as edge senders with 

        (i). a column for network id named \code{id}, 

        (ii). a column \code{Node} for node names, and 

        (iii). the rest for respective node level covariates.
}

\item{FullX}{
	list of numeric arrays of dimension \code{n} by \code{n} by \code{p} of covariates for K different networks. When FullX is provided to the function, edgeCov, receiverCov and senderCov must be specified as NULL. 
}
		
\item{initialVals}{
	an optional list of values to initialize the chain. If \code{NULL} default initialization is used, else 
	\code{initialVals = list(ZZ, beta, intercept, alpha)}.
 
	For fixed effect model \code{beta} is a vector of length \code{p} and \code{intercept} is a vector of length 1.

	For random effect model \code{beta} is an array of dimension  \code{K} by \code{p}, and \code{intercept} is a vector of length \code{K}, where \code{p} is the number of covariates and \code{K} is the number of network.

	\code{ZZ} is an array of dimension \code{NN} by \code{dd}, where \code{NN} is the sum of nodes in all \code{K} networks.

}


\item{priors}{
      an optional list to specify the hyper-parameters for the prior distribution of the paramters.
        If priors = \code{NULL}, default value is used. Else,

    \code{priors=}

	\code{list(MuBeta,VarBeta,MuZ,VarZ,PriorA,PriorB)}
       
	\code{MuBeta} is a numeric vector of length PP + 1 specifying the mean of prior distribution for coefficients and intercept

	\code{VarBeta} is a numeric vector for the variance of the prior distribution of coefficients and intercept. Its length is same as that of MuBeta.

    \code{MuZ} is a numeric vector of length same as the dimension of the latent space, specifying the prior mean of the latent positions.
 
    \code{VarZ} is a numeric vector of length same as the dimension of the latent space, specifying diagonal of the variance covariance matrix of the prior of latent positions.

    \code{PriorA, PriorB} is a numeric variable to indicate the rate and scale parameters for the inverse gamma prior distribution of the hyper parameter of variance of slope and intercept
       }


\item{tune}{
    an optional list of tuning parameters for tuning the chain. If tune = \code{NULL}, default tuning is done. Else, 

	\code{tune = list(tuneBeta, tuneInt,tuneZ)}.

	\code{tuneBeta} and \code{tuneInt} have the same structure as \code{beta} and \code{intercept} in \code{initialVals}.

	 \code{ZZ} is a vector of length \code{NN}.
}
    \item{tuneIn}{
    a logical to indicate whether tuning is needed in the MCMC sampling. Default is \code{FALSE}.
}


\item{dd}{
    dimension of latent space.
}

\item{estimate.intercept}{

  When TRUE, the intercept will be estimated. If the variance of the latent positions are of interest, intercept=FALSE will allow users to obtain a unique variance. The intercept can also be inputed by the user. 
  }
  
\item{niter}{
    number of iterations for the MCMC chain.
}


\item{object}{
	object of class 'HLSM' returned by \code{HLSM()} or \code{HLSMfixedEF()}
}

\item{burnin}{
	numeric value to burn the chain while extracting results from the 'HLSM'object. Default is \code{burnin = 0}. 
}

\item{thin}{
	numeric value by which the chain is to be thinned while extracting results from the 'HLSM' object. Default is \code{thin = 1}. 
}
\item{verbose}{
logical value; TRUE results in messages during MCMC tuning
}
}


\value{
    Returns an object of class "HLSM". It is a list with following components:
    \item{draws}{
       list of posterior draws for each parameters.
    }
    \item{acc}{
      list of acceptance rates of the parameters.
    }
    \item{call}{
    the matched call.
 }

\item{tune}{
	final tuning values
}	

}

\details{
The \code{HLSMfixedEF} and \code{HLSMrandomEF} functions will not automatically assess thinning and burn-in. To ensure appropriate inference, see \code{HLSMdiag}. 
See also \code{LSM} for fitting network data from a single network.
}

\author{
    Sam Adhikari & Tracy Sweet
}

\references{Tracy M. Sweet, Andrew C. Thomas and Brian W. Junker (2013), "Hierarchical Network Models for Education Research: Hierarchical Latent Space Models", Journal of Educational and Behavorial Statistics.
}


\examples{

library(HLSM)
data(schoolsAdviceData)

#Set values for the inputs of the function
priors = NULL
tune = NULL
initialVals = NULL
niter = 10
\donttest{
lsm.fit = LSM(Y=School9Network,edgeCov=School9EdgeCov, 
senderCov=School9NodeCov, receiverCov=School9NodeCov, estimate.intercept=0, niter = niter)
}
}



