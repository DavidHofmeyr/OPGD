\name{OPGDcluster}
\alias{OPGDcluster}
\title{
  Optimal Projections for Clustering with Gaussian Mixture Models
}
\description{
  Obtains an optimal linear projection of a data matrix X, for clustering with a Gaussian Mixture Model
}
\usage{
  OPGDcluster(X, k, V0 = NULL, ndim = NULL, omega = NULL, standardise = FALSE)
}
\arguments{
  \item{X}{n x d matrix with observations row-wise.}
  \item{k}{integer number of clusters/mixture components.}
  \item{V0}{(optional) d x ndim matrix to be used as the initial projection to be optimised. If omitted then a heuristic
    similar to the LDA discriminant features based on the initial solution will be used for initialisation.}
  \item{ndim}{(optional) integer number of dimensions of the projection matrix. If omitted then this is set to the number of
    clusters minus one.}
  \item{omega}{(optional) non-negative numeric. Penalty coefficient to discourage projections which are not orthonormal. If
    omitted then this is set to n, the number of data.}
  \item{standardise}{(optional) logical. Should the data be standardised before applying clustering?}
}
\value{
a named list with fields
  \item{$V}{The optimal projection matrix.}
  \item{$X}{The data matrix passed to the function as argument. If standardise is set to TRUE then
            the output X will be the standardised version.}
  \item{$y}{The vector of cluster labels.}
  \item{$MU}{The matrix of component means (in the full dimensional space).}
  \item{$S}{The list of component covariance matrices (in the full dimensional space).}
  \item{$PI}{The vector of component mixing proportions.}
  \item{$k}{The number of clusters/components.}
  \item{$y0}{The clusters from the initial solution.}
  \item{$posterior}{The matrix of estimated probabilities of membership of each datum to each component/cluster.}
  \item{$V0}{The initial projection matrix.}
}
