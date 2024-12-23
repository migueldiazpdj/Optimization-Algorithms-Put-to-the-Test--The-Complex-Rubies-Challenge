# ACO + Grid Search for a multidimensional problem.
setwd("C:/Users/Miguel/Desktop/OBA")
source("hill.R")    
source("blind.R")      
source("montecarlo.R") 

if (!require("rminer")) install.packages("rminer")
library(rminer)
if (!require("plotly")) install.packages("plotly")
library(plotly)

rubies <- function(x) {
  k <- 1  # Value proportionality constant
  c <- 5  # Penalty constant
  
  # Total value
  V <- k * (x[1]^2 + x[2]^2 + (3 - x[1] - x[2])^2) +
    k * sin(pi * x[1]) * cos(pi * x[2]) 
  
  # Penalty
  P <- c * abs(3 - x[1] - x[2])
  
  # Objective function
  T <- V - P
  
  return(T)
}


# ACO Implementation with Grid Search Refinement for 3D Rubies Functionn
aco_rubies <- function(max_ants = 1000, max_iter = 20, alpha = 1.8, beta = 2, rho = 0.1, q = 5) {
  n_bins <- 50  
  lower <- c(0, 0)  
  upper <- c(2, 2)  
  bins_x1 <- seq(lower[1], upper[1], length.out = n_bins)
  bins_x2 <- seq(lower[2], upper[2], length.out = n_bins)
  pheromones <- matrix(1, nrow = n_bins, ncol = n_bins) 
  
  evaluate <- function(x1_index, x2_index) {
    x1 <- bins_x1[x1_index]
    x2 <- bins_x2[x2_index]
    return(rubies(c(x1, x2)))
  }
  
  best_value <- -Inf
  best_solution <- c(NA, NA)
  values_over_time <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    solutions <- matrix(NA, nrow = max_ants, ncol = 2)
    values <- numeric(max_ants)
    
    for (ant in 1:max_ants) {
      probs_x1 <- (apply(pheromones, 1, sum)^alpha) * (1 / (1:n_bins))^beta
      probs_x1 <- probs_x1 / sum(probs_x1, na.rm = TRUE)
      if (any(is.na(probs_x1)) || sum(probs_x1) == 0) {
        probs_x1 <- rep(1 / n_bins, n_bins)
      }
      x1_index <- sample(1:n_bins, size = 1, prob = probs_x1)
      
      probs_x2 <- (pheromones[x1_index, ]^alpha) * (1 / (1:n_bins))^beta
      probs_x2 <- probs_x2 / sum(probs_x2, na.rm = TRUE)
      if (any(is.na(probs_x2)) || sum(probs_x2) == 0) {
        probs_x2 <- rep(1 / n_bins, n_bins)
      }
      x2_index <- sample(1:n_bins, size = 1, prob = probs_x2)
      
      solutions[ant, ] <- c(bins_x1[x1_index], bins_x2[x2_index])
      values[ant] <- evaluate(x1_index, x2_index)
      
      if (values[ant] > best_value) {
        best_value <- values[ant]
        best_solution <- solutions[ant, ]
      }
    }
    
    pheromones <- (1 - rho) * pheromones
    for (ant in 1:max_ants) {
      x1_index <- which.min(abs(bins_x1 - solutions[ant, 1]))
      x2_index <- which.min(abs(bins_x2 - solutions[ant, 2]))
      pheromones[x1_index, x2_index] <- pheromones[x1_index, x2_index] + q * values[ant]
    }
    pheromones <- pheromones / max(pheromones)
    values_over_time[iter] <- best_value
  }
  
  return(list(best_value = best_value, best_solution = best_solution))
}

# Grid Search Refinement
recursive_grid_search <- function(initial_solution, refinement_range, tolerance = 1e-6, max_iter = 50) {
  refined_solution <- initial_solution
  refined_value <- rubies(refined_solution)
  if (is.na(refined_value)) stop("Error: Initial solution returns NA.")
  
  iter <- 0
  while (iter < max_iter) {
    iter <- iter + 1
    search_x1 <- seq(refined_solution[1] - refinement_range, refined_solution[1] + refinement_range, length.out = 10)
    search_x2 <- seq(refined_solution[2] - refinement_range, refined_solution[2] + refinement_range, length.out = 10)
    search_values <- expand.grid(x1 = search_x1, x2 = search_x2)
    search_values$T <- apply(search_values, 1, function(row) {
      val <- rubies(c(row["x1"], row["x2"]))
      if (is.na(val)) return(-Inf)
      return(val)
    })
    best_index <- which.max(search_values$T)
    new_solution <- as.numeric(search_values[best_index, c("x1", "x2")])
    new_value <- search_values$T[best_index]
    if (abs(new_value - refined_value) < tolerance) break
    refined_solution <- new_solution
    refined_value <- new_value
  }
  return(list(best_split = refined_solution, best_value = refined_value))
}

# Run ACO
set.seed(2023)
aco_result <- aco_rubies(max_ants = 1000, max_iter = 20, alpha = 1.8, beta = 2.5, rho = 0.1, q = 5)
cat("ACO Best Split Point:", aco_result$best_solution, "\n")
cat("ACO Maximum Value of T(x):", aco_result$best_value, "\n")

# Refinement using Grid Search
refinement_range <- 0.1
refined_result <- recursive_grid_search(aco_result$best_solution, refinement_range)
cat("Refined Best Split Point:", refined_result$best_split, "\n")
cat("Refined Maximum Value of T(x):", refined_result$best_value, "\n")


# Plotting
x1 <- seq(0, 2, length.out = 100) 
x2 <- seq(0, 2, length.out = 100)
z <- outer(x1, x2, Vectorize(function(x1, x2) rubies(c(x1, x2))))

z <- t(z)  

# Optimum by ACO
aco_best_x <- 1.836735
aco_best_y <- 1.183673
aco_best_z <- rubies(c(aco_best_x, aco_best_y))  

# 3D graph with ACO optimum
fig <- plot_ly(
  x = x1, y = x2, z = z,
  type = "surface",
  colorscale = "Viridis"
) %>%
  add_trace(
    x = c(aco_best_x), 
    y = c(aco_best_y), 
    z = c(aco_best_z),
    type = "scatter3d", 
    mode = "markers",
    marker = list(color = "red", size = 5),
    name = "ACO Optimum"
  ) %>%
  layout(
    title = "3D Objective Function: Rubies with Corrected ACO Optimum",
    scene = list(
      xaxis = list(title = "x1"),
      yaxis = list(title = "x2"),
      zaxis = list(title = "T(x)")
    ),
    annotations = list(
      list(
        x = aco_best_x,
        y = aco_best_y,
        z = aco_best_z,
        text = paste0("ACO Optimum<br>x1 = ", round(aco_best_x, 4), 
                      "<br>x2 = ", round(aco_best_y, 4), 
                      "<br>T(x) = ", round(aco_best_z, 4)),
        showarrow = TRUE,
        arrowhead = 4,
        ax = 40,
        ay = -40
      )
    )
  )

# We sview the plot
fig
