#########################################
# GROUP 8, AIMS 2024, MMED
#########################################

#*****************************************
## 0. LOAD THE NECESSARY LIBRARIES
#*****************************************

# Remove all objects in memory
rm(list = ls())

# Libraries
library(deSolve)
library(tidyverse)


#*****************************************
# 1. LOADING AND EXPLORING THE DATA
#*****************************************

# Get the current directory
getwd()

# Define the path towards a new directory
path <- "/home/emmanuelnk/Desktop/MMED/TUTORIALS/"

# Set the working directory to the specified path
setwd(path)

# Loading the data
load("dataLondon.Rdata")

# Confirm the dataset is loaded by listing the objects in the environment
ls()

# Inspect the first 5 rows of the London dataset
head(London, 5)

# Shows last 5 rows
tail(London, 5) 

# Shows variable names
names(London) 

# Shape of our dataset:46X7
dim(London)

# Rows of our dataset
nrow(London)

# Columns of our dataset
ncol(London)

# Types (classes) of variables in our dataset
class(London$measles) # Nice: integer

# Let's inspect the biths dynamics
births_london <- ggplot(London, aes(x = year, y = births)) +
  geom_point(color = "violet") +
  geom_line(color = "violet") +
  labs(title = "Births in London",
       x = "Year",
       y = "Births") +
  theme_minimal()
print(births_london)
# Save the plot to a PDF file
ggsave(filename = "births_london_plot.pdf", plot = births_london, device = "pdf", path = path)

# Mean of the column births
mean(London$births)

# Let's display the measles dynamics
measles_london <- ggplot(London, aes(x = year, y = measles)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  labs(title = "Measles in London",
       x = "Year",
       y = "Measles cases") +
  theme_minimal()
print(measles_london)
# Save the plot to a PDF file
ggsave(filename = "measles_london_plot.pdf", plot = measles_london, device = "pdf", path = path)

# Mean of the column measles
mean(London$measles)


# Let's display the measles incidence plot
measles_plot <- ggplot(London, aes(x = year, y = meas.inc)) +
  geom_point(color = "magenta") +
  geom_line(color = "magenta") +
  labs(title = "Measles in London: Annual incidence",
       x = "Year",
       y = "Annual incidence (cases per 100,000 population)") +
  theme_minimal() +
  geom_vline(xintercept = 1967, color = "blue", linetype = "dashed", size = 1) # Add vertical line

print(measles_plot)
# Save the plot to a PDF file
ggsave(filename = "measles_incidence_vac_plot.pdf", plot = measles_plot, device = "pdf", path = path)


#**********************************************************
## 2. SEIR MODEL IMPLEMENTATION
#**********************************************************


# Function to calculate beta
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

# SEIR model function
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

# Function to simulate the SEIR model with adjusted tolerances
simulate_seir <- function(beta, N, vacc_rate, init, times) {
  # Define the parameters for the SEIR model, including beta, population size (N),
  # birth/death rate (mu), incubation rate (sigma), recovery rate (gamma), and vaccination rate (vacc_rate)
  params <- c(beta = beta, N = N, mu = 0.02 / 365.25, sigma = 1 / 8, gamma = 1 / 5, vacc_rate = vacc_rate)
  
  # Use the lsoda function to solve the SEIR differential equations
  # init(y): initial conditions for S, E, I, R
  # times: sequence of time points for which the solution is sought
  # seir: function that defines the SEIR model
  # params: parameters for the SEIR model
  # rtol and atol: relative and absolute tolerance for the ODE solver, set by default to 1e-6 for precision
  output <- data.frame(lsoda(init, times, seir, params, rtol = 1e-6, atol = 1e-6))
  
  # Convert the output to a data frame and return it
  return(output)
}


# Initial conditions and parameters
N0 <- 7696580 # Population in 1946
S0 <- 0.9 * N0 # 90% of the population
E0 <- 200 # Initial exposed population (arbitrary estimate)
I0 <- 22880 # Number of measles cases in 1946
R0 <- N0 - S0 - E0 - I0 # Remaining population

init <- c(S = S0, E = E0, I = I0, R = R0)
# Total number of days over 70 years.
times <- seq(0, 70 * 365, 7)


# Run the simulation for different vaccination rates
vacc_rates <- c(0, 0.1, 0.2, 0.5, 0.8) # 0%, 10%, 20%, 50%, 80% of the population is vaccinated
results <- list()

for (vacc_rate in vacc_rates) {
  results[[as.character(vacc_rate)]] <- simulate_seir(beta_calc(18), N0, vacc_rate, init, times)
}

# Plot the results
plot_results <- function(results, N0) {
  colors <- c("red", "green", "blue", "purple", "orange")
  i <- 1
  
  for (vacc_rate in names(results)) {
    result <- results[[vacc_rate]]
    if (i == 1) {
      plot(result$time / 365, result$I / N0, type = "l", col = colors[i], lwd = 2,
           xlab = "Time (years)", ylab = "Proportion of infected individuals",
           main = "SEIR Model with Varying Vaccination Rates")
    } else {
      lines(result$time / 365, result$I / N0, col = colors[i], lwd = 2)
    }
    lines(result$time / 365, result$S / N0, col = colors[i], lty = 2)
    i <- i + 1
  }
  
  legend("topright", legend = paste("Vaccination rate =", names(results)),
         col = colors, lty = 1, lwd = 2)
}

# Display the plot
plot_results(results, N0)



##########################################
## DISCUSSION 
#########################################
#This figure we obtenaid above illustrates how the proportion of infected individuals changes
#over time for different vaccination rates, using the SEIR model.

#The red line represents the scenario with no vaccination (0%), which shows frequent spikes in the proportion 
#of infected individuals. Without vaccination, the disease spreads quickly and extensively, 
#leading to high infection peaks. This indicates repeated outbreaks and high transmission rates, 
#which can overwhelm healthcare systems and cause significant morbidity and mortality.

#In the case of a 10% vaccination rate (green line), the infection peaks are lower compared
#to the no-vaccination scenario, but still significant outbreaks occur. A 10% vaccination rate
#helps reduce the spread of the disease, but it is not sufficient to prevent large-scale outbreaks.
#This suggests that low vaccination coverage has a limited impact on disease control.

#At a 20% vaccination rate (blue line), the infection peaks are further reduced, indicating a 
#noticeable impact of vaccination. Moderate vaccination coverage (30%) better controls the spread 
#of the disease, resulting in less severe infection peaks. 
#This implies that a moderate level of vaccination can significantly reduce disease 
#transmission and mitigate outbreaks.

#With a 50% vaccination rate (purple line), there are minimal infection peaks, and the proportion
#of infected individuals remains low throughout the period.
#A 60% vaccination rate effectively controls the disease, preventing large outbreaks 
#and maintaining low infection levels. This emphasizes the importance of high vaccination
#coverage in controlling disease spread and protecting public health.

#Finally, at a 80% vaccination rate (orange line), the infection levels are extremely low, 
#approaching eradication. With 90% vaccination, the disease is almost eliminated from the population, 
#resulting in very low infection rates. This shows that very high vaccination coverage 
#is crucial for eradicating the disease and preventing its re-emergence.
