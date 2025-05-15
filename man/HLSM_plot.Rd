\name{HLSMcovplots}
\alias{HLSMcovplots}


\title{Function plot posterior summaries (boxplots) for regression coefficients
}

\description{
    Function plot posterior summaries (boxplots) for regression coefficients
    }

\usage{

HLSMcovplots(fitted.model, burnin=0, thin=1, verbose=TRUE)

}


\arguments{

\item{fitted.model}{
	Model fit using HLSM fitting function; should be a HLSM or LSM object

}

\item{burnin}{
	Amount of burnin to remove
}

\item{thin}{
 Amount to thin each chain
}
\item{verbose}{
logical to indicate whether message about order of covariates is included with plots
}

}

\details{
The plots show posterior means and 50 and 95 credible intervals.
}

\value{No return value, makes a plot in plotting window}

\author{
    Sam Adhikari & Tracy Sweet
}

\references{Tracy M. Sweet, Andrew C. Thomas and Brian W. Junker (2013), "Hierarchical Network Models for Education Research: Hierarchical Latent Space Models", Journal of Educational and Behavorial Statistics.
}






