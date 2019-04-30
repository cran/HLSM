\name{HLSMDiag}
\alias{HLSMDiag}
\alias{HLSMDiag}

\title{Diagnostic functions for  HLSM object

}
\description{
These functions can either take a single HLSM object as input to give you the optimal thinning, burn-in and chain length values along with traceplot, autocorrelation plots of the betas based on the selected options or can provide you with convergence diagnostic if you pass more one than one and less than five HLSM objects.
}
\usage{
HLSMDiag(...,burnin=0,thin=1,lag=500,type="0",varnum=1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{...}{object of class 'HLSM' obtained as an output from \code{HLSMrandomEF()} or   \code{HLSMfixedEF()}, you can either pass just one object or a list of objects if you wish to check for convergence as well.
    }

    \item{burnin}{numeric value to burn the chain for plotting the results from the 'HLSM'object. 
    }

    \item{thin}{a numeric thinning value.
    }
    
    \item{lag}{A numeric lag value to plot autocorrelation plot.
    }

    \item{type}{indicates the type of output 
      type="0" will give you just the optimal thinning, burn-in and chain length values,
      type="traceplot" will give you the optimal thinning, burn-in and chain length values along with traceplots of the betas,
      type="autocorr" will give you the optimal thinning, burn-in and chain length values along with autocorrelation plots of the betas.
      }
      \item{varnum}{
       The variable/beta number for which diagnostic needs to be calculated( only applicable to random effects model)
      }
}

\value{
 Returns Diagnostic plots and values for the HLSM objects.
}

\author{
Aditya Bhat
}



