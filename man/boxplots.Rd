\name{HLSMcovplots}
\alias{HLSMcovplots}
\title{
Diagnostic boxplots for HLSM objects
}
\description{
Plots the diagnostic boxplots for random or fixed effects model by automatically calculating the thinning and burnin value before plotting the boxplot of the betas.
}
\usage{
HLSMcovplots(model,burnin=0,thin=1)
}

\arguments{
  \item{model}{object of class 'HLSM' obtained as an output from \code{HLSMrandomEF()} or\code{HLSMfixedE()}.
  }
  
  \item{burnin}{numeric value to burn the chain for plotting the results from the 'HLSM'object. 
  }

  \item{thin}{a numeric thinning value.
  }
}

\value{
Returns a boxplot of the Betas for fixed effects and random effects objects.
}

\author{
Aditya Bhat 
}

