# Here it is all the code for the basic optimization problem, to see the last algorithm explained in the video, and the more complex function view the AdvancedOptimizationProblem script.

# Performance Comparison of Monte Carlo, Hill Climbing, and Simulated Annealing for Maximizing T(x): Average Best Value Over Iterations. Also, tabu method.
setwd("C:/Users/Miguel/Desktop/OBA")
source("hill.R")       # Hill Climbing
source("blind.R")      # Blind Search
source("montecarlo.R") # Monte Carlo
install.packages("rminer")
library(rminer)        # Para manejo de promedios

# Objective function
rubies <- function(x) {
  k <- 1  # Constant proportionality
  c <- 10 # Penalty
  V <- k * x[1]^2 + k * (2 - x[1])^2  # Total value
  P <- c * abs(2 * x[1] - 2)          # Penalty
  T <- V - P                          # Objective function
  EV <<- EV + 1                       # Increased evaluations
  if (T > BEST) BEST <<- T            # Upgrade to best value
  if (EV <= MAXIT) F[EV] <<- BEST 
  return(-T)  # Negative for maximization
}

rchange1 <- function(par, lower, upper) {
  hchange(par, lower = lower, upper = upper, rnorm, mean = 0, sd = 0.1, round = FALSE)
}

rchange2 <- function(par) {
  hchange(par, lower = lower, upper = upper, rnorm, mean = 0, sd = 0.5, round = FALSE)
}

# Optimization parameters
Runs <- 100
D <- 1
MAXIT <- 1000
lower <- rep(0, D)
upper <- rep(2, D)

CHILL <- list(maxit = MAXIT, REPORT = 0)
CSANN <- list(maxit = MAXIT, temp = 1, trace = FALSE)
Methods <- c("monte carlo", "hill climbing", "simulated annealing")


RES <- vector("list", length(Methods))
for (m in 1:length(Methods)) {
  RES[[m]] <- matrix(nrow = MAXIT, ncol = Runs)
}


set.seed(2024)
for (R in 1:Runs) {
  s <- runif(D, lower, upper)  
  print(paste("Run:", R))      
  
  # Monte Carlo
  EV <- 0; BEST <- -Inf; F <- rep(NA, MAXIT)
  mcsearch(MAXIT, lower = lower, upper = upper, FUN = rubies)
  RES[[1]][, R] <- F
  print(paste("Monte Carlo BEST:", BEST))
  
  # Hill Climbing
  EV <- 0; BEST <- -Inf; F <- rep(NA, MAXIT)
  hclimbing(s, rubies, change = rchange1, lower = lower, upper = upper, control = CHILL, type = "min")
  RES[[2]][, R] <- F
  print(paste("Hill Climbing BEST:", BEST))
  
  # Simulated Annealing
  EV <- 0; BEST <- -Inf; F <- rep(NA, MAXIT)
  optim(s, rubies, method = "SANN", gr = NULL, control = CSANN)
  RES[[3]][, R] <- F
  print(paste("Simulated Annealing BEST:", BEST))
}


AV <- matrix(nrow = MAXIT, ncol = length(Methods))
for (m in 1:length(Methods)) {
  for (i in 1:MAXIT) {
    AV[i, m] <- mean(RES[[m]][i, ], na.rm = TRUE)
  }
}

# Plot
MIN <- min(AV)
MAX <- max(AV)
g1 <- seq(1, MAXIT, length.out = 1000)
plot(g1, AV[g1, 3], ylim = c(MIN, MAX), type = "l", lwd = 2,
     ylab = "Average Best", xlab = "Number of Evaluations")
lines(g1, AV[g1, 2], lwd = 2, lty = 2)
lines(g1, AV[g1, 1], lwd = 2, lty = 3)
legend("topright", legend = rev(Methods), lwd = 2, lty = 1:3)

# Another plot
plot(g1, AV[g1, 3], ylim = c(MIN, MAX), type = "l", lwd = 2, col = "blue",
     ylab = "Average Best Value", xlab = "Number of Evaluations",
     main = "Comparison of Optimization Methods", cex.main = 1.5, cex.lab = 1.2)
lines(g1, AV[g1, 2], lwd = 2, col = "green")
lines(g1, AV[g1, 1], lwd = 2, col = "red")
grid(nx = NULL, ny = NULL, col = "gray", lwd = 0.5)
legend("bottomright", legend = rev(Methods), col = c("red", "green", "blue"),
       lwd = 2, cex = 1.2, bg = "white", box.lty = 0)


# Results
max_values <- apply(AV, 2, max)
print(max_values)

iterations_to_max <- sapply(1:length(Methods), function(m) {
  apply(RES[[m]], 2, function(x) which.max(x))  
})
avg_iterations <- colMeans(iterations_to_max, na.rm = TRUE)
print(avg_iterations)

fastest_method <- Methods[which.min(avg_iterations)]
print(fastest_method)



# -------------------------------
# TABU SEARCH
library(tabuSearch)

rubies <- function(x) {
  k <- 1  
  c <- 10 
  V <- k * x[1]^2 + k * (2 - x[1])^2  
  P <- c * abs(2 * x[1] - 2)          
  T <- V - P                          
  return(T)  
}

intbin <- function(x, lower = 0, upper = 2, bits = 16) {
  bin_value <- sum(2^(which(rev(x == 1)) - 1))  
  real_value <- lower + (upper - lower) * bin_value / (2^bits - 1)  
  return(real_value)
}


rubies_binary <- function(x) {
  split_value <- intbin(x, lower = 0, upper = 2, bits = 16)  
  
  
  if (split_value < 0 || split_value > 2) {
    return(-Inf)  
  }
  
  f <- rubies(c(split_value)) 
  return(f)
}

# Parameters for Tabu Search
D <- 1  
bits <- 16 
size <- D * bits 
s <- sample(0:1, size, replace = TRUE) 


cat("Initial configuration:", intbin(s, lower = 0, upper = 2, bits = bits),
    "Objective value:", rubies_binary(s), "\n")


set.seed(2023)
result <- tabuSearch(size = size, iters = 100, objFunc = rubies_binary, config = s,
                     neigh = 2, listSize = 10, nRestarts = 1)

# Extract the best result
best_index <- which.max(result$eUtilityKeep)  
best_value <- result$eUtilityKeep[best_index]  
best_config <- intbin(result$configKeep[best_index, ], lower = 0, upper = 2, bits = bits) 

cat("Best configuration:", best_config, "\n")
cat("Maximum value:", best_value, "\n")

# Visualization, Plot
library(ggplot2)
iterations <- seq_along(result$eUtilityKeep)  
objective_values <- result$eUtilityKeep     

tabu_data <- data.frame(Iteration = iterations, ObjectiveValue = objective_values)

ggplot(tabu_data, aes(x = Iteration, y = ObjectiveValue)) +
  geom_line(color = "blue") +
  geom_point(color = "red", size = 1.5) +
  labs(title = "Tabu Search Performance for Rubies Optimization",
       x = "Iteration",
       y = "Objective Value (T(x))") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))


# -------------------------------
# Custom ACO Implementation for Rubies Optimization with Refinement
aco_rubies <- function(max_ants = 500, max_iter = 15, alpha = 1.8, beta = 2.5, rho = 0.05, q = 5) {
  n_bins <- 1000
  lower <- 0
  upper <- 2
  bins <- seq(lower, upper, length.out = n_bins)
  pheromones <- rep(1, n_bins)
  
  evaluate <- function(bin_index) {
    x <- bins[bin_index]
    return(rubies(c(x)))
  }
  
  best_value <- -Inf
  best_solution <- NA
  values_over_time <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    solutions <- numeric(max_ants)
    values <- numeric(max_ants)
    
    for (ant in 1:max_ants) {
      probs <- (pheromones^alpha) * ((1 / seq_len(n_bins))^beta)
      probs <- probs / sum(probs, na.rm = TRUE)
      
      chosen_bin <- sample(seq_len(n_bins), size = 1, prob = probs)
      
      solutions[ant] <- bins[chosen_bin]
      values[ant] <- evaluate(chosen_bin)
      
      if (values[ant] > best_value) {
        best_value <- values[ant]
        best_solution <- solutions[ant]
      }
    }
    
    pheromones <- (1 - rho) * pheromones
    for (ant in 1:max_ants) {
      bin_index <- which.min(abs(bins - solutions[ant]))
      pheromones[bin_index] <- pheromones[bin_index] + q * values[ant]
    }
    
    pheromones <- pmax(pheromones, 1e-10)
    pheromones <- pheromones / max(pheromones)
    
    # Track progress
    values_over_time[iter] <- best_value
    cat("Iteration:", iter, "| Best Value:", best_value, "| Best Solution:", best_solution, "\n")
  }
  
  # Plot
  plot(1:max_iter, values_over_time, type = "o", col = "blue", lwd = 2,
       main = "Convergence of ACO for Rubies", xlab = "Iteration", ylab = "Best Value Found",
       ylim = c(min(values_over_time) - 0.1, 2))
  abline(h = 2, col = "red", lty = 2)
  text(x = max_iter / 2, y = 2 - 0.1, labels = paste("Max Value =", round(best_value, 3)), col = "red")
  
  return(list(best_value = best_value, best_solution = best_solution))
}

# Run ACO
set.seed(2023)
aco_result <- aco_rubies(max_ants = 300, max_iter = 15, alpha = 1.8, beta = 2.5, rho = 0.05, q = 5)
cat("Best Split Point (ACO):", aco_result$best_solution, "\n")
cat("Maximum Value of T(x) (ACO):", aco_result$best_value, "\n")

# Recursive Grid Search
recursive_grid_search <- function(best_split, refinement_range, initial_bins = 500, max_recursions = 5, tolerance = 1e-6) {
  refined_best_split <- best_split
  refined_best_value <- -Inf
  
  for (i in 1:max_recursions) {
    
    fine_bins <- seq(refined_best_split - refinement_range, 
                     refined_best_split + refinement_range, 
                     length.out = initial_bins)
    
    
    fine_values <- sapply(fine_bins, function(x) rubies(c(x)))
    
    
    new_best_value <- max(fine_values)
    new_best_split <- fine_bins[which.max(fine_values)]
    
    
    if (abs(new_best_value - refined_best_value) < tolerance) {
      break
    }
    
    
    refined_best_value <- new_best_value
    refined_best_split <- new_best_split
    
    
    refinement_range <- refinement_range / 2
  }
  
  return(list(best_split = refined_best_split, best_value = refined_best_value))
}

# Initial Grid Search
best_split <- aco_result$best_solution
refinement_range <- 0.05  
refined_result <- recursive_grid_search(best_split, refinement_range)

cat("Refined Best Split Point:", refined_result$best_split, "\n")
cat("Refined Maximum Value of T(x):", refined_result$best_value, "\n")

# Plot the refined results
fine_bins <- seq(refined_result$best_split - refinement_range, 
                 refined_result$best_split + refinement_range, 
                 length.out = 1000)
fine_values <- sapply(fine_bins, function(x) rubies(c(x)))

plot(fine_bins, fine_values, type = "l", col = "blue", lwd = 2,
     main = "Recursive Grid Search Refinement of ACO Solution",
     xlab = "Split Point (grams)", ylab = "T(x) Value",
     ylim = c(min(fine_values) - 0.1, 2), xaxt = "n")
axis(1, at = seq(min(fine_bins), max(fine_bins), length.out = 5)) 
abline(h = 2, col = "red", lty = 2)  
abline(v = refined_result$best_split, col = "green", lty = 2)  
points(refined_result$best_split, refined_result$best_value, col = "black", pch = 19, cex = 1.5)  
text(refined_result$best_split, refined_result$best_value - 0.015,  
     labels = paste0("Best = ", round(refined_result$best_value, 6), 
                     "\nSplit = ", round(refined_result$best_split, 6)),
     pos = 1, col = "black", cex = 0.8)  
legend("topright", legend = c("T(x) Curve", "Max T(x) = 2", "Best Split Point"),
       col = c("blue", "red", "green"), lty = c(1, 2, 2), pch = c(NA, NA, 19), cex = 0.8, bty = "n")


#--------------------------------------------------#
#Graphical representation of best values for each method.
data <- data.frame(
  Algorithm = c("Blind Search", "Grid Search", "Montecarlo", "Hill Climbing",
                "SANN", "Tabu Search", "Evolutionary Algorithm", 
                "Differential Evolution", "Particle Swarm Optimization"),
  Category = c("Blind Search", "Blind Search", "Blind Search", "Local Search", 
               "Local Search", "Local Search", "Population Based", 
               "Population Based", "Population Based"),
  BestValue = c(1.9615, 1.5927, 1.9795, 1.9971, 1.9955, 1.9968, 
                1.9935, 1.9907, 1.9807)
)
library(ggplot2)

ggplot(data, aes(x = reorder(Algorithm, -BestValue), y = BestValue, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_text(aes(label = round(BestValue, 4)), vjust = -0.5, size = 3) +
  labs(title = "Best Value Comparison of Algorithms",
       x = "Algorithm",
       y = "Best Value",
       fill = "Category") +
  coord_cartesian(ylim = c(1.5, 2)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"))


# Deleting Grid Search for this plot, as we sais in the video
data_filtered <- subset(data, BestValue > 1.9)

# Crear el gráfico de puntos con líneas
ggplot(data_filtered, aes(x = reorder(Algorithm, -BestValue), y = BestValue, group = 1, color = Category)) +
  geom_line(size = 1, color = "grey") +             
  geom_point(size = 4) +                            
  geom_text(aes(label = round(BestValue, 4)),        
            vjust = -1, size = 3.5) +
  labs(title = "Best Value Comparison of Algorithms",
       x = "Algorithm",
       y = "Best Value",
       color = "Category") +
  coord_cartesian(ylim = c(1.95, 2)) +                
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"))


#-------------------------------------------
# Performance Comparison of Blind Search and Grid Search

rubies <- function(x) {
  k <- 1  
  c <- 10 
  V <- k * x[1]^2 + k * (2 - x[1])^2  
  P <- c * abs(2 * x[1] - 2)         
  T <- V - P                          
  EV <<- EV + 1                      
  if (T > BEST) BEST <<- T            
  if (EV <= MAXFN) F[EV] <<- BEST     
  return(T) 
}

# Blind Search Implementation
blind_search <- function(f, lower, upper, maxit) {
  best <- -Inf
  for (i in 1:maxit) {
    candidate <- runif(length(lower), lower, upper)
    val <- f(candidate)
    if (val > best) best <- val
  }
  return(best)
}

# Grid Search Implementation
grid_search <- function(f, lower, upper, steps) {
  grid_points <- seq(lower[1], upper[1], length.out = steps)
  best <- -Inf
  for (x in grid_points) {
    val <- f(c(x))
    if (val > best) best <- val
  }
  return(best)
}

# Unified Runner for Algorithms
crun_extended <- function(method, f, lower, upper, maxit, steps = 50) {
  if (method == "Blind") {
    return(blind_search(f, lower, upper, maxit))
  } else if (method == "Grid") {
    return(grid_search(f, lower, upper, steps))
  }
}

# Comparison Function
compare_algorithms <- function(Methods, f, lower, upper, Runs, MAXFN, maxit, steps, LIM) {
  VAL <- matrix(nrow = Runs, ncol = length(Methods))
  for (R in 1:Runs) {
    for (m in 1:length(Methods)) {
      EV <<- 0; BEST <<- -Inf; F <<- rep(NA, MAXFN)
      VAL[R, m] <- crun_extended(Methods[m], f, lower, upper, maxit, steps)
    }
  }
  
  summary_df <- data.frame(
    Method = Methods,
    Avg_Best = round(apply(VAL, 2, mean), digits = 4),
    Std_Dev = round(apply(VAL, 2, sd), digits = 4),
    Success_Rate = round(100 * apply(VAL, 2, function(x) sum(x > LIM)) / Runs, digits = 4)
  )
  print(summary_df)
  
  # Generate Plot
  par(mar = c(4.0, 4.0, 1.8, 0.6)) # Set plot margins
  barplot(
    summary_df$Avg_Best,
    names.arg = summary_df$Method,
    col = c("#1f78b4", "#33a02c"),
    main = "Performance Comparison: Blind vs Grid Search",
    ylab = "Average Best Value",
    xlab = "Method"
  )
  grid()
}

# Parameters
Runs <- 50
MAXFN <- 5000
maxit <- 500
steps <- 50
lower <- c(0)
upper <- c(2)
LIM <- 1.8 # Threshold for successes (we want our value to be as close to 2 as possible)
Methods <- c("Blind", "Grid")

# Run Comparison
compare_algorithms(
  Methods = Methods,
  f = rubies,
  lower = lower,
  upper = upper,
  Runs = Runs,
  MAXFN = MAXFN,
  maxit = maxit,
  steps = steps,
  LIM = LIM
)



#-------------------------------------------
# Modified compare2.R script for rubies function maximization 
library(genalg) 
library(DEoptim) 
library(pso) 

rubies <- function(x) {
  k <- 1 
  c <- 10 
  V <- k * x[1]^2 + k * (2 - x[1])^2 
  P <- c * abs(2 * x[1] - 2)          
  T <- V - P                          
  EV <<- EV + 1                       
  if (T > BEST) BEST <<- T            
  if (EV <= MAXFN) F[EV] <<- BEST     
  return(T)
}

# Auxiliary functions for optimization from Compare 2
crun <- function(method, f, lower, upper, LP, maxit) {
  if (method == "EA") {
    rbga(evalFunc = f, stringMin = lower, stringMax = upper, popSize = LP, iters = maxit * 1.5)
  } else if (method == "DE") {
    C = DEoptim.control(itermax = maxit, trace = FALSE, NP = LP)
    DEoptim(f, lower = lower, upper = upper, control = C)
  } else if (method == "PSO") {
    C = list(maxit = maxit, s = LP)
    psoptim(rep(NA, length(lower)), fn = f, lower = lower, upper = upper, control = C)
  }
}

successes <- function(x, LIM, type = "min") {
  if (type == "min") return(sum(x < LIM)) else return(sum(x > LIM))
}

ctest <- function(Methods, f, lower, upper, type = "min", Runs, D, MAXFN, maxit, LP, pdf, main, LIM) {
  RES <- vector("list", length(Methods))
  VAL <- matrix(nrow = Runs, ncol = length(Methods))
  ITER <- numeric(length(Methods)) 
  TIME <- numeric(length(Methods)) 
  ITER_TO_BEST <- numeric(length(Methods))
  
  for (m in 1:length(Methods)) { 
    RES[[m]] <- matrix(nrow = MAXFN, ncol = Runs)
  }
  set.seed(2023)
  for (R in 1:Runs) {
    for (m in 1:length(Methods)) {
      EV <<- 0; F <<- rep(NA, MAXFN) 
      if (type == "min") BEST <<- Inf else BEST <<- -Inf 
      start_time <- Sys.time() 
      suppressWarnings(crun(Methods[m], f, lower, upper, LP, maxit))
      end_time <- Sys.time() 
      RES[[m]][, R] <- F 
      VAL[R, m] <- F[MAXFN] 
      TIME[m] <- TIME[m] + as.numeric(difftime(end_time, start_time, units = "secs"))
      ITER[m] <- ITER[m] + EV 
      ITER_TO_BEST[m] <- ITER_TO_BEST[m] + which.max(F) 
    }
  }
  
  # Compute average F result per method
  AV <- matrix(nrow = MAXFN, ncol = length(Methods))
  SD <- matrix(nrow = MAXFN, ncol = length(Methods)) 
  for (m in 1:length(Methods)) {
    for (i in 1:MAXFN) {
      AV[i, m] <- mean(RES[[m]][i, ])
      SD[i, m] <- sd(RES[[m]][i, ])
    }
  }
  
  # Prepare summary table
  summary_df <- data.frame(
    Method = c("Evolutionary Algorithm", "Differential Evolution", "Particle Swarm Optimization"),
    Avg_Best = round(apply(VAL, 2, mean), digits = 4),
    Std_Dev = round(apply(VAL, 2, sd), digits = 4),
    Success_Rate = round(100 * apply(VAL, 2, successes, LIM, type) / Runs, digits = 4),
    Avg_Runtime = round(TIME / Runs, 4),
    Avg_Iterations = round(ITER / Runs, 0),
    Avg_Iterations_To_Best = round(ITER_TO_BEST / Runs, 0)
  )
  
  
  print(summary_df)
  
  # Plot results
  par(mar = c(4.0, 4.0, 1.8, 0.6)) 
  MIN <- min(AV); MAX <- max(AV)
  g1 <- seq(1, MAXFN, length.out = 500) # Grid for lines
  plot(g1, AV[g1, 1], ylim = c(MIN, MAX), type = "l", lwd = 2, col = "#1f78b4", main = main, ylab = "average best", xlab = "number of evaluations")
  for (i in 2:length(Methods)) lines(g1, AV[g1, i], lwd = 2, col = c("#1f78b4", "#33a02c", "#e31a1c")[i])
  grid()
  if (type == "min") position <- "topright" else position <- "bottomright"
  legend(position, legend = c("Evolutionary Algorithm", "Differential Evolution", "Particle Swarm Optimization"), lwd = 2, col = c("#1f78b4", "#33a02c", "#e31a1c"))
}

# Define optimization parameters
MAXFN <- 5000
EV <- 0; BEST <- -Inf; F <- rep(NA, MAXFN)

# Define methods
Methods <- c("EA", "DE", "PSO")

# Run optimization for rubies function
Runs <- 50
D <- 1
LP <- 50
maxit <- 100
lower <- rep(0, D)
upper <- rep(2, D)

ctest(
  Methods = Methods,
  f = rubies,
  lower = lower,
  upper = upper,
  type = "max",        
  Runs = Runs,
  D = D,
  MAXFN = MAXFN,
  maxit = maxit,
  LP = LP,
  pdf = "comp-rubies",  
  main = "Rubies Optimization Comparison",
  LIM = 1.8           
)
