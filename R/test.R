rm(list = ls())

# Load functions
source("functions.R")

# Load Required Libraries
library(sandwich)
library(gtools)
library(foreach)
library(parallel)
library(doParallel)
library(doRNG)
library(plm)
library(tictoc)
library(HDcluster)

# Set Seed for Reproducibility
set.seed(1)

# Generate data

N <- 50
T <- 50
E <- matrix( rnorm(N*T), N, T)
U <- matrix( rnorm(N*T), N, T)
E[, 1] <- rnorm(N)
U[, 1] <- rnorm(N)


A <- rgamma(N, 1)
B <- rgamma(N, 1)

A_rep <- matrix(rep(A, T), N, T)
B_rep <- t(matrix(rep(B, N), T, N))


F <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 10)
H <- (0.5 * A_rep^10 + 0.5 * B_rep^10)^(1 / 5)


X <- H + U
Y <- X + F + E

id_code <- numeric(N * T)
for (i in 1:N) {
  for (t in 1:T) {
    id_code[((t - 1) * N + i):(t * N)] <- 10+i
  }
}

time_code <- numeric(N * T)
for (i in 1:N) {
  for (t in 1:T) {
    time_code[((t - 1) * N + 1):(t * N)] <- 10 + t
  }
}

vY <- c(Y)
vX <- c(X)

data <- data.frame(id_code, time_code, vY, vX)

formula <- vY ~ vX - 1
index = c("id_code", "time_code")
init <- 300
# Baseline estimate
est <- estimator_dc(formula, data, index, init = init)
ols <- est[["res"]]
G <- est[["G"]]
C <- est[["C"]]
coeff <- ols$coefficients
std_errors <- sqrt(vcovHC(ols, type = "HC0", method = "arellano"))

# Cross-fitted estimate
est_CF <- estimator_dc(formula, data, index, CF = TRUE, init = init)
ols_CF <- est_CF[["res"]]
coeff_CF <- ols_CF$coefficients
std_errors_CF <- sqrt(vcovHC(ols_CF, type = "HC0", method = "arellano")) * sqrt((N * T) / est_CF[["df"]])
