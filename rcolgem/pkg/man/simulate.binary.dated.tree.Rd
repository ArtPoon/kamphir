\name{simulate.binary.dated.tree}
\alias{simulate.binary.dated.tree}
\title{Simulate a binary dated tree using a coalescent model and given a demographic history}
\usage{
  simulate.binary.dated.tree(births, deaths, nonDemeDynamics,  t0, x0, sampleTimes, sampleStates,  migrations=NA,  parms=NA, fgyResolution = 2000, integrationMethod = 'rk4') 
}
\arguments{
\item{births}{ A vector or matrix of strings. These are evaluated as equations for the number of births within and between demes. Must have rownames and colnames corresponding to the names of demes.}
\item{deaths}{ A vector of strings. These are evaluated as equations for the rate that lineages in each deme are terminated.  Must have rownames corresponding to the names of demes.}
\item{nonDemeDynamics}{ A vector of strings. These are evaluated as equations for the rate of change of state variables that do not correspond to demes. Vector must have names of state variables. }
\item{t0}{ The time of origin of the process. Should be before root of genealogy. } 
\item{x0}{ A vector of initial conditions for the demographic process (the state of the system at time t0). Should include the name and value of all variables mentioned in births and nonDemeDynamics. }
\item{sampleTimes}{A vector of times for each sample. The names attribute will match tip.label of the generated tree. }
\item{sampleStates}{A matrix of states for each sample. There is one row for each sample, and row names should match the names of sampleTimes. Each column gives the probability that the sample is in the deme 1..m at the time of sampling. }
\item{migrations}{ A vector or matrix of strings. These are evaluated as equations for the number of migrations between demes. Must have rownames and colnames corresponding to the names of demes. Should be omitted if there is only one deme.}
\item{parms}{A list of parameters that will be accessible to differential equations specified in births, migrations and nonDemeDynamics.}
\item{fgyResolution}{Determines timestep of ODEs (larger values gives smaller time step)}
\item{integrationMethod}{Passed to lsoda. Higher precision may be obtained with methods such as \emph{adams} at some computational cost. }
}
\value{
	A binaryDatedTree, extends ape::phylo and includes heights of all nodes, time before most recent sample when node occurs. 
}
\references{
E. M. Volz, Complex population dynamics and the coalescent under neutrality, Genetics, January, 2012
}
\author{
Erik M Volz
}

