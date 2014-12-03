\name{calculate.cluster.size.moments.from.tree}
\alias{calculate.cluster.size.moments.from.tree}
\title{Calculate cluster size moments from tree}
\usage{
  calculate.cluster.size.moments.from.tree(bdt, heights)
}
\arguments{
  \item{bdt}{ a binary dated tree object}

  \item{heights}{vector numeric, heights at which to
  calculate cluster sizes}
}
\description{
  Calculate cluster moments of the distribution of cluster sizes (clade sizes), which change as a function of the threhold TMRCA used for defining  clades. Currently, only sampling at a single time point is supported.
}

