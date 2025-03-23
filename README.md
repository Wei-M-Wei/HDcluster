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
A CRAN release is coming soon.


## Features
- **Main functionality**: Perform inference after discretizing unobserved heterogeneity in the panel data model, see the paper [^1].
- **Validation example**: An example 'test.R' is included.

## Additional resources
- **Replication code**: The repository includes replication code for all simulations and empirical applications.
- **Suggestions welcome**: Further improvements are planned, and we encourage feedback and suggestions to enhance the package.
    
    
```

## Reference
[^1]: Beyhum, J., Mugnier, M. Inference after discretizing unobserved heterogeneity. [[arXiv:2502.09740
Search](https://arxiv.org/abs/2502.09740).](https://arxiv.org/abs/2412.07352)
