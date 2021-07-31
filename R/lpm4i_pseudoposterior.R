lpm4i_log_prior <- function(positions)
{
  res <- 0
  n_nodes <- dim(positions)[1]
  n_dimensions <- dim(positions)[2]
  n_time_frames <- dim(positions)[3]
  for (i in 1:n_nodes) for (d in 1:n_dimensions)
  {
    res = res + dnorm(positions[i,d,1], 0, 1, T)
    for (t in 2:n_time_frames) res = res + dnorm(positions[i,d,t], positions[i,d,t-1], 1, T)
  }
  res
}

lpm4i_log_pseudo_likelihood_per_node <- function(adj, data, positions, i, t, eps = 0.001)
{
  res <- 0
  A <- eps
  B <- 0
  for (j in 1:nrow(data)) if (i != j)
  {
    d_ijt <- 1 / sqrt(  sum( (positions[i,,t]-positions[j,,t])^2 )  )
    A = A + d_ijt * adj[i,j,t]
    B = B + d_ijt * adj[i,j,t] * data[j,t]
  }
  dnorm(data[i,t], B/A, sqrt(1/A), log = T)
}

lpm4i_log_pseudo_likelihood <- function(adj, data, positions, eps = 0.001)
{
  res <- 0
  for ( t in 1:ncol(data) ) for (i in 1:nrow(data)) 
  {
    res = res + lpm4i_log_pseudo_likelihood_per_node(adj, data, positions, i, t, eps)
  }
  res
}

lpm4i_log_pseudo_posterior <- function(adj, data, positions, eps = 0.001)
{
  res <- lpm4i_log_prior(positions)
  for ( t in 1:ncol(data) ) for (i in 1:nrow(data)) 
  {
    res = res + lpm4i_log_pseudo_likelihood_per_node(adj, data, positions, i, t, eps)
  }
  res
}
