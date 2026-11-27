# Generate Simulated Panel Survey Data

A longitudinal dataset with explanatory variable age and dependent
variables income (numeric), ntrips (count) and satisfaction (ordinal) is
created. This function also allows to re-create the dataset `psurvey`
included in the package.

## Usage

``` r
generate_psurvey_data(n)
```

## Arguments

- n:

  number of sample units.

## Details

The example below shows how the dataset `psurvey` included in the
package has been generated.

## Value

A list containing:

- data :

  a data frame with covariates, random effects, predictors and sampled
  outcomes,

- sd_num :

  a vector of length three containing the true values of the standard
  deviations for the numeric covariate,

- beta_num :

  a list of true values (named vector) of \\\beta\\ parameters for
  numeric outcome for each cluster,

- beta_poi :

  a list of true values (named vector) of \\\beta\\ parameters for count
  outcome for each cluster,

- beta_ord :

  a list of true values (named vector) of \\\beta\\ parameters for
  ordinal outcome for each cluster,

- c_ord :

  a list of true values of intercepts specific to each level of the
  ordinal outcome for each cluster.

## Author

Implemented by Jan Vávra to demonstrate the use of the package on
artificial data in a vignette.

## Examples

``` r
### Generate psurvey data
set.seed(31415)
n <- 500
ldata <- generate_psurvey_data(n)
psurvey1 <- ldata$data[,c("i", "j", "age", "income", "ntrips", "satisfaction")]
psurvey_latent1 <- ldata$data[,c("i", "j", "g", "b_num", "b_poi", "b_ord")]

data("psurvey", package = "clustGLMM")
data("psurvey_latent", package = "clustGLMM")
all.equal(psurvey, psurvey1)
#> [1] TRUE
all.equal(psurvey_latent, psurvey_latent1)
#> [1] TRUE
```
