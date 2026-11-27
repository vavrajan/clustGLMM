# Generate an Artificial Longitudinal Dataset with Mixed-type Outcomes

Five different longitudinal outcomes are sampled according to the model
described by Vávra et al. (2024). Column `i` is the unit id, `j` is the
number of observation for `i`th unit. Numeric `ynum`, count `ypoi`,
binary `ybin`, ordinal `yord` and categorical `ycat` outcomes are
sampled according to predictors `pred_...` which are constructed from
the intercept, time covariate `x`, categorical covariate `f` and the
random effects `b_...`. The definition of the predictor can be adjusted.

## Usage

``` r
generate_longitudinal_mixed_type_data(n, n_i = 4, G = 2, timepar = "parallel", 
  catpar = "no", catsubjfix = TRUE, xin1 = 1, Sigma, 
  Kord = 5, Kcat = 4, kspec_bi_cat = FALSE)
```

## Arguments

- n:

  number of sample units.

- n_i:

  number of observations per unit.

- G:

  number of clusters, only `G=2` and `G=3` are predefined.

- timepar:

  a string declaring what is the parametrization of time covariate `x`.
  One of `no`, `parallel`, `cross`, `spline`, `mix`.

- catpar:

  a string declaring what is the parametrization of binary covariate
  `f`. One of `no`, `fixed`, `group`.

- catsubjfix:

  a logical value; should the values of categorical covariate be fixed
  for the same unit or alternate among observations?

- xin1:

  `1/xin1` is the width of observational window for one unit, so `xin1`
  declares how many times does the observational period of one unit fit
  into the overall data gathering period.

- Sigma:

  the covariance matrix for the random intercepts (of dimension 5, one
  for each outcome). Can be a list of covariance matrices, in such a
  case the model assumes group-specific covariance matrix.

- Kord:

  number of categories for ordinal outcome.

- Kcat:

  number of categories for categorical outcome (only up to 4
  categories).

- kspec_bi_cat:

  a logical value; should the random effects be specific to each level
  of categorical outcome? If so, the dimension of Sigma should be
  larger.

## Details

The function `generate_longitudinal_mixed_type_data` is limited in the
number of clusters `G` since it is quite laborious to define all the
true values of model parameters (e.g. effects \\\beta\\).

Time covariate `x` is sampled so that observations for one unit are
close to each other. First, half-length of the observational window
`xmid = 1/(2xin1)` is defined. Then, values for a unit are sampled by
uniform distribution on `[-xmid, xmid]` and randomly shifted by
multiples of `xmid`. Time covariate is limited to interval \[0, 1\].

The spline parametrization
([`library("splines")`](https://rdrr.io/r/base/library.html)) of the
time covariate is of `degree = 2` and has
`knots = c(0, 0.25, 0.5, 0.75, 1)`

Binary covariate `f` is sampled with probability 0.5 either once per
unit (`catsubjfix = TRUE`) or repeatedly for each observation
(`catsubjfix = FALSE`).

The group allocation indicator `g` is sampled uniformly on {1, ..., G}
for each unit.

Then the linear predictors for each observation and outcome are created
according to selected parametrization of time covariate `x` (`timepar`)
and parametrization of binary covariate `f` (`catpar`). The intercepts
are always group-specific which ensures at least some difference among
clusters.

Choices for `timepar`:

- no:

  covariate `x` does not modify the predictor,

- parallel:

  clusters have the same linear trend in time `x` but different
  intercepts,

- cross:

  clusters have different linear trends in time `x` (one is increasing,
  another decreasing),

- spline:

  each cluster has a different randomly generated spline
  parametrization,

- mix:

  outcomes have different selection of all of the above (currently,
  numeric is `spline`, count is `cross`, binary is `parallel`, ordinal
  is `spline`, categorical is `cross`).

Choices for `catpar`:

- no:

  covariate `f` does not modify the predictor,

- fixed:

  clusters have the same effect of `f`,

- group:

  clusters have different effects of `f`.

Numerical outcome is sampled from normal distribution with
group-specific standard deviations (from 0.5 to 1.0).

Count outcome is sampled from Poisson distribution with rate made by
[`exp()`](https://rdrr.io/r/base/Log.html) of the predictor.

Binary outcome (0 or 1) is sampled with probability computed by
[`plogis`](https://rdrr.io/r/stats/Logistic.html).

Ordinal outcome (ordered factor) is sampled by probabilities computed by
sequential differences of
[`plogis`](https://rdrr.io/r/stats/Logistic.html) applied to shifted
predictor.

Categorical outcome (unordered factor) is sampled by probabilities
computed by `softmax` function from several predictors specific for each
category; `softmax(x) = exp(x) / sum(exp(x))` where `x` is
`Kcat`-dimensional. By default, predictor for the first category is set
to 0 for identifiability purposes.

## Value

A list containing:

- data :

  a data frame with covariates, random effects, predictors and sampled
  outcomes,

- sd_num :

  a vector of true values of standard deviations for numeric outcome
  depending on the total number of clusters `G`,

- beta_num :

  a list of true values (named vector) of \\\beta\\ parameters for
  numeric outcome for each cluster,

- beta_poi :

  a list of true values (named vector) of \\\beta\\ parameters for count
  outcome for each cluster,

- beta_bin :

  a list of true values (named vector) of \\\beta\\ parameters for
  binary outcome for each cluster,

- beta_ord :

  a list of true values (named vector) of \\\beta\\ parameters for
  ordinal outcome for each cluster,

- c_ord :

  a list of true values of intercepts specific to each level of the
  ordinal outcome for each cluster,

- beta_cat :

  a list of true values (matrix with named columns) of \\\beta\\
  parameters for categorical outcome for each cluster.

## References

Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise
multivariate regression of mixed-type panel data. Stat Comput 34, 46
(2024). https://doi.org/10.1007/s11222-023-10304-5

## Author

Implemented by Jan Vávra for performing a simulation study in Vávra et
al. (2024).

## Examples

``` r
n <- 100      
n_i <- 4       
G <- 2         
timepar <- "mix" 
catpar <- "no"  
catsubjfix <- TRUE 
xin1 <- 1        
Kord <- 5        
Kcat <- 4        
kspec_bi_cat <- FALSE
corr <- matrix(c(1, -0.5, -0.5, -0.4, -0.4,
                 -0.5, 1, 0.3, 0.4, 0.5,
                 -0.5, 0.3, 1, 0.2, 0.3,
                 -0.4, 0.4, 0.2, 1, 0.1,
                 -0.4, 0.5, 0.3, 0.1, 1), byrow = TRUE, ncol = 5)
sdSigmas <- c(0.5, 0.5, 0.5, 0.5, 0.5)
Sigma <- diag(sdSigmas) %*% corr %*% diag(sdSigmas)

set.seed(31415)
ldata <- generate_longitudinal_mixed_type_data(n = n, n_i = n_i,
                     G = G, timepar = timepar,
                     catpar = catpar, catsubjfix = catsubjfix,
                     xin1 = xin1, Sigma = Sigma,
                     Kord = Kord, Kcat = Kcat, kspec_bi_cat = kspec_bi_cat)
data <- ldata$data
head(data)
#>   i         x j          bs1       bs2        bs3        bs4       bs5 f g
#> 1 1 0.1330718 1 0.6395799651 0.1416649 0.00000000 0.00000000 0.0000000 0 2
#> 2 1 0.3357378 2 0.2158566570 0.7253356 0.05880771 0.00000000 0.0000000 0 2
#> 3 1 0.4901114 3 0.0007822799 0.5379900 0.46122776 0.00000000 0.0000000 0 2
#> 4 1 0.9502223 4 0.0000000000 0.0000000 0.01982255 0.33875390 0.6414235 0 2
#> 5 2 0.2356964 1 0.5523042649 0.4444222 0.00000000 0.00000000 0.0000000 1 2
#> 6 2 0.5821669 2 0.0000000000 0.2253436 0.72064519 0.05401118 0.0000000 1 2
#>        b_num     b_poi      b_bin      b_ord      b_cat   pred_num   pred_poi
#> 1 -0.1229307 0.3223326  0.2887508 -0.5344187 -0.0766008  1.1114427  1.0561889
#> 2 -0.1229307 0.3223326  0.2887508 -0.5344187 -0.0766008  1.1770693  0.6508571
#> 3 -0.1229307 0.3223326  0.2887508 -0.5344187 -0.0766008  1.1770693  0.3421099
#> 4 -0.1229307 0.3223326  0.2887508 -0.5344187 -0.0766008 -0.2895885 -0.5781120
#> 5  0.7138559 0.4919251 -0.4880856 -0.5612637 -0.7358979  2.0128738  1.0205323
#> 6  0.7138559 0.4919251 -0.4880856 -0.5612637 -0.7358979  1.9436414  0.3275913
#>    pred_bin    pred_ord c_ord1 c_ord2 c_ord3 c_ord4 pred_cat1   pred_cat2
#> 1 1.5548945 -0.43483571   -2.5   -1.5   -0.5    0.5         0  0.65725553
#> 2 1.9602263 -0.58338757   -2.5   -1.5   -0.5    0.5         0  0.25192367
#> 3 2.2689735 -0.22675527   -2.5   -1.5   -0.5    0.5         0 -0.05682354
#> 4 3.1891954  0.02682738   -2.5   -1.5   -0.5    0.5         0 -0.97704541
#> 5 0.9833071 -0.53968731   -2.5   -1.5   -0.5    0.5         0 -0.20729059
#> 6 1.6762482  0.03684712   -2.5   -1.5   -0.5    0.5         0 -0.90023163
#>     pred_cat3  pred_cat4       ynum ypoi ybin yord ycat
#> 1 -0.81045712 -0.0766008  0.5326489    2    1    0    0
#> 2 -0.40512527 -0.0766008  0.6595827    1    0    2    0
#> 3 -0.09637806 -0.0766008  1.0163111    4    1    3    3
#> 4  0.82384381 -0.0766008 -0.7877999    0    1    4    3
#> 5 -1.26450511 -0.7358979  3.4585835    2    1    1    3
#> 6 -0.57156407 -0.7358979  2.3370678    3    1    0    0

# The same as the "longitudinal_mixed_type_data"
# which will be used in other examples
data("longitudinal_mixed_type_data")
```
