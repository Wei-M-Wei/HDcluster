# Clean R Script
rm(list = ls())

library(dplyr)
library(kableExtra)

DGP <- 2
rho <- 0.7
kappa <- 0.7

filename <- paste("results", "_DGP" , DGP, "_rho",rho,"_kappa", kappa ,".RData",sep="")

load(filename)

df <- results[, c(2,13:16,25:28,37:40,29:36)] %>% mutate(across(everything(), ~ round(., 3)))

# Output Results as a LaTeX Table
df %>% 
  kable("latex", align = "c", linesep = "", escape = FALSE)
