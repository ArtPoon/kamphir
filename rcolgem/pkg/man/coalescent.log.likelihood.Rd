\name{coalescent.log.likelihood}
\alias{coalescent.log.likelihood}
\title{
Compute the likelihood of a gene genealogy conditonal on a demographic history.
}
\description{
 The likelihood of a dated tree is calculated conditional on a demographic history, which consists of the timeseries of births, migrations, and population sizes in the population. 
}

\usage{
coalescent.log.likelihood( bdt, births, deaths, nonDemeDynamics
   ,  t0, x0,  migrations=NA,  parms=NA, fgyResolution = 2000
   , integrationMethod = 'rk4',  censorAtHeight=FALSE
   , forgiveAgtY=.2, returnTree=FALSE)
}

\arguments{
  \item{bdt}{A binaryDatedTree object which extends ape:phylo and includes dates of tips and internal nodes.}
\item{births}{ A vector or matrix of strings. These are evaluated as equations for the number of births within and between demes. Must have rownames and colnames corresponding to the names of demes.}
\item{deaths}{ A vector of strings. These are evaluated as equations for the rate that lineages in each deme are terminated.  Must have rownames corresponding to the names of demes.}
\item{nonDemeDynamics}{ A vector of strings. These are evaluated as equations for the rate of change of state variables that do not correspond to demes. Vector must have names of state variables. }
\item{t0}{ The time of origin of the process. Should be before root of genealogy. } 
\item{x0}{ A vector of initial conditions for the demographic process (the state of the system at time t0). Should include the name and value of all variables mentioned in births and nonDemeDynamics. }
\item{migrations}{ A vector or matrix of strings. These are evaluated as equations for the number of migrations between demes. Must have rownames and colnames corresponding to the names of demes. Should be omitted if there is only one deme.}
\item{parms}{A list of parameters that will be accessible to differential equations specified in births, migrations and nonDemeDynamics.}
\item{fgyResolution}{Determines timestep of ODEs (larger values gives smaller time step)}
\item{integrationMethod}{Passed to lsoda. Higher precision may be obtained with methods such as \emph{adams} at some computational cost. }
\item{censorAtHeight}{Boolean or numeric. If numeric, will ignore all parts of the tree older than this value when calculating the likelihood.}
\item{forgiveAgtY}{If this is zero, function will return -Inf if at any point the number of lineages exceeds population size. If this value is between zero and one, function will return -Inf if the number of lineages exceeds population size AND number of lineages is less than this value times sample size. This can be useful when fitting models where the number of lineages is close to population size, which often happens during early period of exponential growth. }
\item{returnTree}{ If TRUE, the return value is a list with two elements. The first element is the log likelihood, and the second element is a modified version of the binaryDatedTree that includes the estimated ancestral states along each lineage. Ancestral state information is contained in the matrices lstates and ustates, which has the same number of rows as the tree edge attribute. These respectively describe the state at the beginning and end of each edge (present to past). }
}
\details{
   See vignettes for  detailed usage. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
 The log likelihood (numeric). 
}
\references{
E. M. Volz, Complex population dynamics and the coalescent under neutrality, Genetics, January, 2012
}
\author{
Erik M Volz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{binaryDatedTree}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ coalescent }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
