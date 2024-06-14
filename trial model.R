library(deSolve)
# Function to calculate beta
load("C:\\Users\\tshir\\Downloads\\dataLondon.Rdata")
str(London)
tri<-London
measles_plot2 <- ggplot(London, aes(x = year, y = meas.inc)) +
  geom_point(color = "red") +
  geom_line(color = "red") +
  labs(title = "Measles in London: Annual incidence before vaccination",
       x = "Year",
       y = "Annual incidence (cases per 100,000 population)") +
  theme_minimal()
print(measles_plot2)
measles_plot <- ggplot(London[0:20,], aes(x = year, y = meas.inc)) +
  geom_point(color = "red") +
  geom_line(color = "red") +
  labs(title = "Measles in London: Annual incidence before vaccination",
       x = "Year",
       y = "Annual incidence (cases per 100,000 population)") +
  theme_minimal()
print(measles_plot)
measles_plot1 <- ggplot(London[21:46,], aes(x = year, y = meas.inc)) +
  geom_point(color = "red") +
  geom_line(color = "red") +
  labs(title = "Measles in London: Annual incidence after vaccination",
       x = "Year",
       y = "Annual incidence (cases per 100,000 population)") +
  theme_minimal()
print(measles_plot1)

beta_calc <- function(R_0, mu = 0.02 / 365.25, sigma = 1 / 8, gamma = 1 / 5) {
  R_0 / ((sigma / (mu + sigma)) * (1 / (mu + gamma)))
}

# Function to calculate S*
S.star <- function(beta, N, mu = 0.02 / 365.25, sigma = 1 / 8, gamma = 1 / 5) {
  N * ((mu + sigma) / beta) * ((mu + gamma) / sigma)
}

# Function to calculate E*
E.star <- function(beta, N, mu = 0.02 / 365.25, sigma = 1 / 8, gamma = 1 / 5) {
  mu * (N - S.star(beta, N, mu, sigma, gamma)) / (mu + sigma)
}

# Function to calculate I*
I.star <- function(beta, N, mu = 0.02 / 365.25, sigma = 1 / 8, gamma = 1 / 5) {
  (sigma / (mu + gamma)) * E.star(beta, N, mu, sigma, gamma)
}

# Function to calculate R*
R.star <- function(beta, N, mu = 0.02 / 365.25, sigma = 1 / 8, gamma = 1 / 5) {
  (gamma / mu) * I.star(beta, N, mu, sigma, gamma)
}

# Set initial conditions
initial_conditions <- c(S = 2000, E = 0, I = 1, R = 0)

N0 = initial_conditions[1] + initial_conditions[2] + initial_conditions[3] + initial_conditions[4]

#N0 <- 7696580
S.star(beta_calc(13),N0)+E.star(beta_calc(13),N0)+I.star(beta_calc(13),N0)+R.star(beta_calc(13),N0)

# Define the ODE system
seir_model <- function(t, y, parameters) {
  with(as.list(c(y, parameters)), {
    
    S <- y[1]
    E <- y[2]
    I <- y[3]
    R <- y[4]
    
    N = sum(y)
    dS <- nu * (1 - p) - (beta_calc(13) * S * I) / N - mu * S
    dE <- (beta_calc(13) * S * I) / N - (sigma + mu) * E
    dI <- sigma * E - (gamma + mu) * I
    dR <- nu * p + gamma * I - mu * R
    return(list(c(dS, dE, dI, dR)))
  })
}

seir <- function(t, y, params) {
  S <- y[1]
  E <- y[2]
  I <- y[3]
  R <- y[4]
  
  beta <- params["beta"]
  N <- params["N"]
  mu <- params["mu"]
  sigma <- params["sigma"]
  gamma <- params["gamma"]
  vacc_rate <- params["vacc_rate"]
  nu <- mu * N
  
  dSdt <- nu - (beta * S * I / N) - (mu * S) - (vacc_rate * S)
  dEdt <- (beta * S * I / N) - (mu * E) - (sigma * E)
  dIdt <- (sigma * E) - (mu * I) - (gamma * I)
  dRdt <- (gamma * I) - (mu * R) + (vacc_rate * S)
  
  return(list(c(dSdt, dEdt, dIdt, dRdt)))
}


# Set parameter values
parameters <- c(beta_calc(13), sigma = 1/8, gamma = 1/5, mu = 0.02, nu = 0.02*N0[[1]], p = 1)

# Set time range
time_range <- seq(0, 70, by = 1)

# Solve the ODE system
solution <- lsoda(initial_conditions, time_range, seir_model, parameters)


# Extract the solution
solution_matrix <- solution[, 2:5]
colnames(solution_matrix) <- c("S", "E", "I", "R")

# Plot the solution

matplot(time_range, solution_matrix, type = "l", xlab = "Time (years)", ylab = "Population",
        main = paste("SEIR Model : p =",parameters[length(parameters)]), lty = 1, lwd = 2, col = c("black", "green", "red", "blue"))
legend("top", legend = c("Susceptible", "Exposed", "Infected", "Recovered"),
       col = c("black", "green", "red", "blue"), lty = 1, bty = "n")

 
