\name{binaryDatedTree}
\alias{binaryDatedTree}
\title{Create binary dated tree. The binaryDatedTree class includes heights for each node and other helper variables}
\usage{
  binaryDatedTree(phylo, sampleTimes, sampleStates=NULL
  , sampleStatesAnnotations=NULL)
}
\arguments{
	\item{phylo}{ape::phylo tree}
	\item{sampleTimes}{Named vector of sample times for each taxon. Names must correspond to phylo$tip.label}
	\item{sampleStates}{Matrix of discrete character information for each taxon. If NULL, this may be inferred from taxon label. Rownames must correspond to phylo$tip.label.  Colnames must correspond to names of demes used in subsequent analysis. Each row is a probability vector summing to 1 (the probability that taxon occupies each deme). }
	\item{sampleStatesAnnotations}{Vector of possible discrete character states for taxa. If inferring taxon state from label, this provides the possible matches for taxon annotations. The annotations are extracted from all text following '_' at the end of the taxon label. }
}
\description{
  Create binary dated tree binaryDatedTree class, includes
  heights for each node and other helper variables like maxSampleTime.
}
\value{
	A binaryDatedTree object. Includes numeric vector heights for each tip and internal node, the time before most recent sample when each node occurs. 
}
