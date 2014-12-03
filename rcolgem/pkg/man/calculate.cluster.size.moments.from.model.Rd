\name{calculate.cluster.size.moments.from.model}
\alias{calculate.cluster.size.moments.from.model}
\title{Calculate cluster size moments from model}
\usage{
  calculate.cluster.size.moments.from.model(sampleTime,
    sampleStates, maxTime = NA, minTime = NA,
    timeResolution = 50, discretizeRates = FALSE,
    fgyResolution = 100, integrationMethod = "adams")
}
\arguments{
  \item{sampleTime}{Time that sample is collected}

  \item{sampleStates}{(n X M) matrix of sample states for each taxon}

  \item{maxTime}{}

  \item{minTime}{}

  \item{timeResolution}{Determines time step for discrete approximations. Larger resolution gives smaller time step. }

  \item{discretizeRates}{if TRUE, a discrete approximation is used to calculate transition rates}

  \item{fgyResolution}{the number of time steps to use if discretizeRates=TRUE}

  \item{integrationMethod}{method passed to deSolve::ode}
}
\description{
  Calculate cluster moments of the distribution of cluster sizes (clade sizes), which change as a function of the threhold TMRCA used for defining  clades. Currently, only sampling at a single time point is supported.
}

