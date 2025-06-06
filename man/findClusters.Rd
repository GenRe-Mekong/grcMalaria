% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public.R
\name{findClusters}
\alias{findClusters}
\title{Clustering}
\usage{
findClusters(
  ctx,
  sampleSet,
  clusterSet,
  minIdentity = 1,
  impute = TRUE,
  clusteringMethod = "louvain",
  minClusterSize = 10,
  resolution = 1,
  objectiveFunction = "CPM",
  beta = 0.01
)
}
\arguments{
\item{ctx}{The analysis context, created by intializeContext().}

\item{sampleSet}{The name of the sample set being used, which must have been previously created by selectSampleSet().}

\item{clusterSet}{Give the cluster set a name. In the example the clustering set is called "GMS".}

\item{minIdentity}{The minimal similarity level set for a pair of samples to be in a cluster. For example, "0.95" corresponds to at least 95 percent genetic barcode similarity.
Multiple similarity levels can be set at once, as in the example above, by putting number inside c() separated by , the default is 1.}

\item{impute}{To use imputed or filtered data. The default is TRUE.}

\item{clusteringMethod}{The clustering method. Methods available: "allNeighbours", "leiden" and "louvain". The default is "louvain".
The "allNeighbours" method clusters samples together that are above the set "minIdentity" threshold.
This method is less informative at low similarity levels, because each sample will be assigned to a single cluster.
The "louvain" and "leiden" methods are preferred. They use an algorithm to identify clusters within a network that 
are strongly connected to each other, and more weakly connected to other clusters.}

\item{minClusterSize}{To avoid creating very small clusters, one can set a minimum cluster size. The default is 10 samples.}

\item{resolution}{A tuning parameter for the "leiden" and "louvain" methods. 
Values are between 0 and 1; lowering the value produces larger, "looser" clusters. Default is 1.}

\item{objectiveFunction}{A parameter for the "leiden" method only. Possible values: "CPM" (Constant Potts Model) or "modularity". Default is "CPM".}

\item{beta}{A tuning parameter for the "leiden" method only. Affects the randomness in the algorithm. Default is 0.01.}
}
\description{
In order to plot maps of cluster prevalence and cluster sharing, first clusters of similar genetic background need to be assigned.
This function partitions the Plasmodium samples into clusters with a similar genetic background, based on the genetic barcode similarity.
This function was built on the igraph R package (Csárdi G, 2023), which computes a graph connecting sample pairs with S greater than a minimum threshold S(min) , subsequently partitioning the graph into clusters.
Functions providing cluster analyses include: mapClusterSharing(), and mapClusterPrevalence().
}
\examples{
\dontrun{
## Find clusters of similar genetic background ##
   findClusters(ctx, sampleSet="EBKK", clusterSet = "GMS",
                minIdentity = c(0.95, 0.80), impute = TRUE,
                clusteringMethod = "louvain")
}
}
