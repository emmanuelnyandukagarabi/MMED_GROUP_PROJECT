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


