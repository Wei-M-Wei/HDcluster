# Clean R Script
rm(list = ls())

library(dplyr)
library(kableExtra)

DGP <- 2
rho <- 0.7
kappa <- 0.7

filename <- paste("results", "_DGP" , DGP, "_rho",rho,"_kappa", kappa ,".RData",sep="")

load(filename)

df <- results[, c(-1,-13,-14,-15,-16,-25:-48)] %>% mutate(across(everything(), ~ round(., 3)))

# Output Results as a LaTeX Table
df %>% 
  kable("latex", align = "c", linesep = "", escape = FALSE)
