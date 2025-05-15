\name{schoolsAdviceData}
\alias{ps.advice.mat}
\alias{ps.all.vars.mat}
\alias{ps.edge.vars.mat}
\alias{ps.school.vars.mat}
\alias{ps.teacher.vars.mat}
\alias{School9Network}
\alias{School9NodeCov}
\alias{School9EdgeCov}
\alias{ps.advice.df}
\alias{ps.edge.df}
\alias{ps.node.df}


\title{HLSM: Included Data Sets}
\description{
	Data set included with the HLSM package: network variables from Pitts
	and Spillane (2009).
}


\usage{
	ps.advice.mat
	ps.advice.df
	ps.all.vars.mat
	ps.edge.vars.mat
	ps.edge.df
	ps.school.vars.mat
	ps.teacher.vars.mat
	ps.node.df
	School9Network
	School9NodeCov
	School9EdgeCov
}

\format{
	ps.advice.mat: a list of 15 sociomatrices of advice seeking network, one for each school.
	
	ps.advice.df: a data frame of all ties.

	ps.all.vars.mat: a list of 15 arrays of all the covariates, one for each school.
	edge.vars.mat: a list of edge level covariates for 15 different school.
	
	ps.edge.df: a dataframe of all edge covariates.
	
	ps.school.vars.mat: a list of school level covariates for all 15 schools.

	ps.teacher.vars.mat: a list of node level covariates for all 15 schools.
	
	ps.node.df: a dataframe of all node covariates.

	ps.all.vars.mat: a single list of length 15 containing the covariates mentioned above.

  School9Network: a single adjacency matrix from School 9.
  
  School9NodeCov: a dataframe with node covariates
  
  School9EdgeCov: a dataframe with dyad-level covariates.
 }



\references{
 Pitts, V., & Spillane, J. (2009). "Using social network methods to study school leadership".International Journal of Research & Method in Education, 32, 185-207

 Sweet, T.M., Thomas, A.C., and Junker, B.W. (2013).
 "Hierarchical Network Models for Education Research: Hierarchical Latent Space Models". Journal of Educational and Behavorial Statistics.
}


\author{Sam Adhikari and Tracy Sweet}
\keyword{datasets}


