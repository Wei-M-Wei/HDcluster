# Clean R Script
rm(list = ls())

library(dplyr)
library(kableExtra)

N = 50
T = 50

filename <- paste("results_simi", "_N" , N, "_T",T,".RData",sep="")

load(filename)

df <- results %>% mutate(across(everything(), ~ round(., 3)))

# Output Results as a LaTeX Table
df %>% 
  kable("latex", align = "c", linesep = "", escape = FALSE)
