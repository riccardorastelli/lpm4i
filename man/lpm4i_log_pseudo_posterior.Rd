\name{lpm4i_log_pseudo_posterior}
\alias{lpm4i_log_pseudo_posterior}
\title{lpm4i_log_pseudo_posterior}
\description{
Calculates the log pseudo posterior associated to a configuration of the latent positions, for a given sequence of \code{K} binary undirected networks of \code{N} nodes.
}
\usage{
lpm4i_log_pseudo_posterior(adj, data, positions, eps = 0.001)
}
\arguments{
  \item{adj}{An array of size \code{NxNxK} where each of the \code{K} slices corresponds to the undirected binary adjacency matrix of a network.}
  \item{data}{A \code{NxK} matrix of real values representing nodal information. This corresponds to the response variable that is being modelled.}
  \item{positions}{An array of size \code{NxDxK} indicating the latent positions of each of the \code{N} nodes, at each of the \code{K} time frames. \code{D} is the number of dimensions for the latent space.}
  \item{eps}{A small non-negative value should be chosen. This corresponds to a small noise that is added to a denominator in order to avoid numerical problems. A small value will not affect the results.}
}
\value{
The value of the log pseudo posterior density.
}
\examples{
set.seed(12345)
N <- 10
K <- 3
D <- 2
adj <- array(0, c(N, N, K))
for (k in 1:K) for (i in 1:N) for (j in 1:N) if (i < j)
{
  adj[i,j,k] = sample(c(0,1), 1, TRUE)
  adj[j,i,k] = adj[i,j,k]
}
data <- matrix(rnorm(N*K), N, K)
positions <- array(rnorm(N*D*K), c(N, D, K))
temp <- seq(2, 0.001, -0.1)
lpm4i_log_pseudo_posterior(adj, data, positions)
}


