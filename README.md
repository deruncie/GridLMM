**GridLMM** is a package for fitting linear mixed models (LMMs) with multiple random effects.

The fitting process is optimized for repeated evaluation of the random effect
    model with different sets of fixed effects, ex. for GWAS analyses. The approximation
    is due to the use of a discrete grid of possible values for the random effect variance
    component proportions. We include functions for both frequentist and Bayesian GWAS, 
    (Restricted) Maximum Likelihood evaluation, Bayesian Posterior inference of variance components,
    and Lasso/Elastic Net fitting of high-dimensional models with random effects.
    
Please treat this as a Beta version and let me know of issues running the functions.    
    
### Installation:
```{r}
devtools::install_github('deruncie/GridLMM')
```

### Main functions:

- `GridLMM_ML`: estimates parameters of a LMM by (restricted) Maximum Likelihood
- `GridLMM_posterior`: Approximates the posterior distribution of the variance component
proportions of a LMM
- `GridLMM_GWAS`: Runs a GWAS with the error structure a LMM. 
By default, uses heuristics to efficiently sample the grid. 
Can run Wald tests (`method = 'REML'`), Likelihood ratio tests (`method = 'ML'`), or calculate Bayes Factors (`method = 'BF'`)
- `GridLMMnet`: Fit a multiple regression with LASSO / elastic net penalty and mixed model error term.

### Examples:

There is a vignette walking through the data format necessary for GridLMM and a few analyses using `GridLMM_GWAS_fast`
```{r}
vignette(topic = 'Running_GridLMM_GWAS',package='GridLMM')
```

