\name{OPGD}
\alias{OPGD}
\title{
  Optimal Projections for Classification with Gaussian Class Densities
}
\description{
  Obtains an optimal linear projection of a data matrix X, for classification based on a Quadratic Discriminant type model.
  Predictions on test data can be obtained using the predict function.
}
\usage{
  OPGD(X, y, V0 = NULL, ndim = NULL, lambda = 0)
}
\arguments{
  \item{X}{n x d matrix with observations row-wise.}
  \item{y}{vector of class labels. The current implementation necessitates that these are the values 1, 2, ..., K, where K
  is the total number of classes.}
  \item{V0}{(optional) d x ndim matrix to be used as the initial projection to be optimised. If omitted then a heuristic
  similar to the LDA discriminant features will be used for initialisation.}
  \item{ndim}{(optional) integer number of dimensions of the projection matrix. If omitted then this is set to the number of
  classes minus one.}
  \item{lambda}{(optional) numeric in the interval [0, 1]. Allows the sharing of pooled information across classes, as in regularised
  discriminant analysis. Each class covariance matrix is given by (1-lambda)*Sigma_k + lambda*Sigma,
  where Sigma_k is an unbiased estimate of the class covariance and Sigma is a pooled estimate of the overall within class covariance.}
}
\value{
a named list with class OPGD and fields
  \item{$V}{The optimal projection matrix.}
  \item{$X}{The (mean-centered) data matrix passed to the function as argument.}
  \item{$y}{The class labels passed as argument.}
  \item{$MU}{The matrix of class means.}
  \item{$S}{The list of class covariance matrices.}
  \item{$PI}{The vector of class proportions.}
  \item{$nc}{Integer number of classes.}
}
