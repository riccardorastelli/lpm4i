lpm4i_sann <- function(adj, data, positions, temp, eps = 0.001, verbose = T)
{
  start_time <- Sys.time()
  n_iterations <- length(temp)
  n_nodes <- dim(positions)[1]
  n_dimensions <- dim(positions)[2]
  n_time_frames <- dim(positions)[3]
  obj_current <- lpm4i_log_pseudo_posterior(adj, data, positions, eps)
  obj <- rep(0, n_iterations + 1)
  obj[1] = obj_current
  
  distances <- array(0, c(n_nodes,n_nodes,n_time_frames))
  A <- matrix(eps, n_nodes, n_time_frames)
  B <- matrix(0, n_nodes, n_time_frames)
  for (t in 1:n_time_frames) for (i in 1:n_nodes) for (j in 1:n_nodes) if (i != j)
  {
    distances[i,j,t] <- sqrt(  sum( (positions[i,,t]-positions[j,,t])^2)  )
    A[i,t] = A[i,t] + adj[i,j,t] / distances[i,j,t]
    B[i,t] = B[i,t] + data[j,t] * adj[i,j,t] / distances[i,j,t]
  }
  
  if (verbose) cat("\nSimulated annealing optimisation started\n")
  for (iter in 1:n_iterations) 
  {
    for (t in 1:n_time_frames) for (i in 1:n_nodes)
    {
      prop_position <- rep(0, n_dimensions)
      for (d in 1:n_dimensions) prop_position[d] = rnorm(1, positions[i,d,t], 1)
      
      prop_distances <- distances[i,,t]
      prop_A <- A[,t]
      prop_B <- B[,t]
      for (j in 1:n_nodes) if (i != j) prop_distances[j] <- sqrt(  sum( (prop_position-positions[j,,t])^2)  )
      for (j in 1:n_nodes) if (i != j) 
      {
        prop_A[i] = prop_A[i] + adj[i,j,t] / prop_distances[j] - adj[i,j,t] / distances[i,j,t]
        prop_B[i] = prop_B[i] + data[j,t] * adj[i,j,t] / prop_distances[j] - data[j,t] * adj[i,j,t] / distances[i,j,t]
      }
      for (j in 1:n_nodes) if (i != j) 
      {
        prop_A[j] = prop_A[j] + adj[i,j,t] / prop_distances[j] - adj[i,j,t] / distances[i,j,t]
        prop_B[j] = prop_B[j] + data[i,t] * adj[i,j,t] / prop_distances[j] - data[i,t] * adj[i,j,t] / distances[i,j,t]
      }
      
      delta_like <- 0
      for (j in 1:n_nodes)
      {
        delta_like = delta_like + dnorm(data[j,t], prop_B[j]/prop_A[j], sqrt(1/prop_A[j]), log = T) - 
          dnorm(data[j,t], B[j,t]/A[j,t], sqrt(1/A[j,t]), log = T)
      }
      
      delta_prior <- 0
      for (d in 1:n_dimensions)
      {
        if (t == 1) delta_prior = delta_prior + dnorm(prop_position[d], 0, 1, T) - dnorm(positions[i,d,t], 0, 1, T)
        if (t > 1)  delta_prior = delta_prior + dnorm(prop_position[d], positions[i,d,t-1], 1, T) - dnorm(positions[i,d,t], positions[i,d,t-1], 1, T)
        if (t < n_time_frames) delta_prior = delta_prior + dnorm(positions[i,d,t+1], prop_position[d], 1, T) - dnorm(positions[i,d,t+1], positions[i,d,t], 1, T)
      }
      
      delta <- delta_like + delta_prior
      
      if (  runif(1) < exp( delta / temp[iter] )  ) 
      {
        for (d in 1:n_dimensions) positions[i,d,t] = prop_position[d]
        obj_current = obj_current + delta
        for (j in 1:n_nodes)
        {
          distances[i,j,t] <- distances[j,i,t] <- prop_distances[j]
          A[j,t] = prop_A[j]
          B[j,t] = prop_B[j]
        }
        
      }
      
    }
    obj[iter+1] = obj_current
    if (verbose) if (iter %% 100 == 0) 
    {
      elapsed_seconds <- format( round(Sys.time() - start_time, 2), units = "secs")
      cat("Elapsed: ", elapsed_seconds, "\t\tEnd of iteration: ", iter, "\n")
    }
  }
  if (verbose) 
  {
    elapsed_seconds <- format( round(Sys.time() - start_time, 2), units = "secs")
    cat("Simulated annealing optimisation ended in ", elapsed_seconds, " seconds\n")
  }
  list(objective = obj, positions = positions, elapsed_seconds = elapsed_seconds)
}
