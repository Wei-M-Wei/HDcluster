## HDcluster

This is the first version of the R package 'HDcluster'. 

To install and use the package, we recommend installing the package by
```{r }
install.packages("devtools")  # If not already installed
library(devtools)
install_github("Wei-M-Wei/HDcluster")
```
Once installed, load the package with
```{r }
library(HDcluster)
```

## Features
- **Main functionality**: Perform inference after discretizing unobserved heterogeneity in the panel data model, see the paper [^1].
- **Validation example**: An example 'test.R' is included. 'estimator_dc(formula, data, index)' is the main function.

## Additional resources
- **Replication code**: The repository includes replication code for all simulations and empirical applications.
- **Suggestions welcome**: Further improvements are planned, and we encourage feedback and suggestions to enhance the package.


## Example
```{r }
# data should contain the 'id', 'time', outcome variable Y and regressors X
data <- data.frame(id_code, time_code, vY, vX)

# Specify the regression formula, same as lm().
formula <- vY ~ vX - 1

# Specify the name of the 'id' and 'time' in your data
index = c("id_code", "time_code")

# Specify the initial iterations of kmeans cluster
init <- 300

# Baseline estimate, allows for cross-fitting
est <- estimator_dc(formula, data, index, init = init)

# We recommond having a look at the 'text.R' in the folder 'R'
```
A CRAN release is coming soon.

## Reference
[^1]: Beyhum, J., Mugnier, M. Inference after discretizing unobserved heterogeneity. [[arXiv:2502.09740
Search](https://arxiv.org/abs/2502.09740).](https://arxiv.org/abs/2412.07352)
