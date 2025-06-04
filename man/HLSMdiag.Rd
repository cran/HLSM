\name{HLSMdiag}
\alias{HLSMdiag}

\title{Function to conduct diagnostics the MCMC chain from a random effect HLSM (and HLSMfixedEF for fixed effects model)
}

\description{
    Function to compute and report diagnostic plots and statistics for a single or multiple HLSM objects.
}

\usage{

HLSMdiag(object, burnin = 0,
         diags = c('psrf', 'raftery', 'traceplot', 'autocorr'),
         col = 1:6, lty = 1)
}

\arguments{


\item{object}{
	object or list of objects of class 'HLSM' returned by \code{HLSMrandomEf()} or \code{HLSMfixedEF()}
}

\item{burnin}{
	numeric value to burn the chain while extracting results from the 'HLSM'object. Default is \code{burnin = 0}. 
}

\item{diags}{
  a character vector that is a subset of \code{c('psrf', 'raftery', 'traceplot', 'autocorr')}. Default returns all diagnostics. If only a single chain is supplied in object, \code{'psrf'} throws a warning if explicitly requested by user.
}

\item{col}{
  a character or integer vector specifying the colors for the traceplot and autocorr plot
}
\item{lty}{
  a character or integer vector specifying the linetype for the traceplot and autocorr plot
}
}



\value{
    Returns an object of class "HLSMdiag". It is a list that contains variable-level diagnostic tables from either or both of the raftery diagnostic and psrf diagnostic.
    \item{call}{
      the matched call.    
    }
    \item{raftery}{
       list of matrices of suggested niters, burnin, and thinning for each chain.
    }
    \item{psrf}{
      list containing \code{psrf}, a matrix of psrf estimates and upper limits for variable, and \code{mpsrf} the multivariate psrf estimate.
    }
}



\author{
    Christian Meyer
}


