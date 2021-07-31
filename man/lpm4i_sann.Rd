\name{lpm4i_sann}
\alias{lpm4i_sann}
\title{lpm4i_sann}
\description{
Launches the simulated annealing algorithm to find the latent positions for the model, given a sequence of \code{K} binary undirected networks of \code{N} nodes.
}
\usage{
lpm4i_sann(adj, data, positions, temp, eps = 0.001, verbose = T)
}
\arguments{
  \item{adj}{An array of size \code{NxNxK} where each of the \code{K} slices corresponds to the undirected binary adjacency matrix of a network.}
  \item{data}{A \code{NxK} matrix of real values representing nodal information. This corresponds to the response variable that is being modelled.}
  \item{positions}{An array of size \code{NxDxK} indicating the initial values for the latent positions of each of the \code{N} nodes, at each of the \code{K} time frames. \code{D} is the number of dimensions for the latent space.}
  \item{temp}{A vector indicating the temperature value that should be used at each step of the simulated annealing procedure. The length of this vector implies the number of iterations for the procedure.}
  \item{eps}{A small non-negative value should be chosen. This corresponds to a small noise that is added to a denominator in order to avoid numerical problems. A small value will not affect the results.}
  \item{verbose}{True or false indicating whether extra comments should be printed out during the procedure.}
}
\value{
A list with components:
  \item{objective}{The values of the objective function at each of the steps of the optimisation procedure.}
  \item{positions}{The estimated latent positions of the nodes.}
  \item{elapsed_seconds}{Number of seconds required by the algorithm to complete the task.}
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
lpm4i_sann(adj, data, positions, temp)
}


