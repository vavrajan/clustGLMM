# Evaluate Clustering Probabilities and Deviance of Output of `clustGLMM`

To evaluate the clustering probabilities, all latent variables need to
be integrated out. No analytic closed form solution is available for the
integral in case of joint random effects for mixed-type outcomes.
Numerical approximations are obtained using Laplace's approximation or
Adaptive Gaussian Quadrature (AGQ). The samples of the model parameters
from
[`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md)
are sent to C for faster evaluation.

## Usage

``` r
clustering_probabilities_and_deviance(mcmc, id, data, dynamic_prob = FALSE, 
  start = 1, end = mcmc$iter, thin = 1, chains = mcmc$chains,
  howsave = mcmc$howsave,                                      
  tolerance = mcmc$tuning$double$tolerance, 
  maxiter = mcmc$tuning$integer$maxiter, 
  maxnrep = mcmc$tuning$integer$maxnrep, 
  NGQ = 1)
```

## Arguments

- mcmc:

  an object of class `"clustglmm"` (as returned by function
  [`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md)).

- id:

  a string, name of the column containing the unit ID. New units not
  used in MCMC sampling may be used here.

- data:

  a data.frame containing outcomes and regressors including the unit ID
  variable `id`.

- dynamic_prob:

  a Boolean; indicates if the clustering probabilities are computed
  dynamically. That is, if unit `i` is observed `ni` times, then `ni`
  different probabilities are computed, each based on the first `j`
  observations for `j in 1:ni`.

- start:

  an iteration number where to start.

- end:

  an iteration number where to end.

- thin:

  a thinning that should be used between `start` and `end`; e.g., if
  `thin = 10` only every 10th iteration is approximated; default is 1.

- chains:

  chain number(s) for which the approximation should be computed.

- howsave:

  a string declaring how the sampled chains should be returned as
  `draws`. One of:

  `"list"`

  :   `draws` is a list indexed by chain numbers, parameters will be
      saved in a structured list in which it is easy to orient,

  `"matrix"`

  :   `draws` is a large matrix, where one row corresponds to one
      iteration of one specific chain and columns are all model
      parameters,

  `"cmcmc"`

  :   `draws` is a list indexed by chain numbers and contains the
      condensed output of the inner C function (not recommended!).

- tolerance:

  tolerance used when computing the norm of the shift within
  Newton-Raphson.

- maxiter:

  maximum number of iterations allowed during Newton-Raphson.

- maxnrep:

  maximum number of repetitions of Newton-Raphson from different
  starting points.

- NGQ:

  number of points for Adaptive Gaussian Quadrature; if `NGQ = 1`, then
  simple Laplace's approximation is used. Do not use values higher than
  7, since there are `NGQ^totnran` summands with `totnran` denoting the
  dimension of the random effects.

## Details

An object of class `"clustglmm"` can be used regardless of the format,
i.e., `howsave = "list"` or `howsave = "matrix"`.

The data.frame `data` containing the outcomes and regressors does not
have to be the same as used for MCMC sampling. Even completely new units
are allowed here.

The indexing of probabilities is more complicated when
`dynamic_prob = TRUE`. In this case, for one unit there will be as many
probabilities as observations belong to this unit. If unit `i` is
observed 3 times, there will be 3 probabilities: one where only the
first observation is used, a second where the first two observations are
used and a third where all observations are used. With more available
data the classification towards one of the clusters should be more
pronounced. In the end, there will be as many probabilities as rows of
`data`. However, it does not mean that the rows correspond to each
other, e.g., when `data` is not ordered by unit IDs. The reason is that
we use a `for` cycle that goes through all units and then for a specific
unit through all its observations.

The code can handle the case where `data` is not sorted by ID. However,
for one specific unit the data are assumed to be ordered in time. This
is internally used as ordering for the dynamic probability calculation.
However, we recommend to insert data also sorted by ID so that pairing
rows with the probabilities obtained with `dynamic_prob = TRUE` is
easier.

Computations are usually expensive and a division of the workload into
several jobs where different chains or even subsamples are processed
simultaneously might be useful. This can be achieved by setting `start`,
`end` and `thin` differently for different jobs. However, you need to
divide this on your own and also gather the results back together. To do
this properly keep in mind that the used iteration numbers are created
by `seq(start, end, by = thin)`.

To approximate the `totnran`-dimensional integral over the random
effects, we use Adaptive Gaussian Quadrature. First, the Newton-Raphson
method is used to find the Taylor expansion of the log-pdf as a function
of random effects. If `NGQ = 1` then the integral is approximated by
evaluation only in the central point found by Newton-Raphson, which
corresponds to simple Laplace's approximation. If `NGQ = 2, 3, ..., 7`
then `hermite.he.quadrature.rules(n = NGQ, normalized = FALSE)` from
package gaussquad is used to find the roots and weights of the Hermite
quadrature rules. These are combined for every dimension of random
effects, which yields `NGQ^totnran` points at which the log-pdf needs to
be evaluated. So, think carefully whether higher approximation precision
with high `NGQ` is worth the expensive evaluation.

## Value

An object of class `"clustglmm"` is returned so that already existing
methods are applicable. Most of the outputs are taken from the given
`mcmc` object.

Depending on the choice of `howsave = "list"` or `howsave = "matrix"`,
the `$draws` element contains a list of outputs for different chains
(from `chains`) organized in either `list` or `matrix` format. Each
chain then contains the following quantities:

- `dev_i`:

  a vector of individual unit contribution to the model deviance,

- `deviance`:

  a scalar containing the model deviance (sum of all individual
  contributions),

- `pUig_int`:

  a group-specific vector of length equal to number of unique `id`
  labels containing clustering probabilities towards one of `mcmc$G`
  clusters. The mapping between ids and numbers is given via the item
  `$numbered_unique_ids`.

Outputs `$clustering` and `$draws` are derived from `pUig_int` by
averaging and evaluating, which cluster has the highest posterior
probability.

## References

Pinheiro, J.C., Chao, E.C. Efficient Laplacian and adaptive Gaussian
quadrature algorithms for multilevel generalized linear mixed models. J.
Comput. Graph. Stat. 15(1), 58–81 (2006).
https://doi.org/10.1198/106186006x96962

Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise
multivariate regression of mixed-type panel data. Stat Comput 34, 46
(2024). https://doi.org/10.1007/s11222-023-10304-5

## Author

Jan Vávra

## See also

[`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md)

## Examples

``` r
# \donttest{
### Get some longitudinal mixed type data
### here generated by generate_longitudinal_mixed_type_data()
data("longitudinal_mixed_type_data")

# Examples here demonstrate the basic use of the available functions on a toy dataset.
# For more thorough analyses see the vignettes.


###-----------------------------------------------------------------------------
###     Prepare inputs for clustGLMM()
###-----------------------------------------------------------------------------
id <- "i"

family <- c("num", "poi", "bin", "ord", "cat")
names(family) <- Ys <- c("ynum", "ypoi", "ybin", "yord", "ycat")

## Create list of formulae for each outcome
formula <- list()
for(y in Ys){
  formula[[y]] <- list()
  formula[[y]]$fixed <- ~ f
  formula[[y]]$group <- ~ x
  formula[[y]]$random <- ~ 1
  formula[[y]]$offset <- ""
}
# intercept term will be group-specific 
# since it appears in $group which has a preference over $fixed

## All other inputs need to be prepared,
## unless the default values are desired (varying, save, howsave, param, tuning).
## See the vignette for how these can be adjusted.

###-----------------------------------------------------------------------------
###     Calling clustGLMM()
###-----------------------------------------------------------------------------

Gmax <- 5     # maximal number of mixture components to be considered
nchains <- 1  # number of chains to be sampled 

## No inits specified
set.seed(123456789)
mcmc <- clustGLMM(formula = formula, id = id, family = family, 
                  data = longitudinal_mixed_type_data, 
                  G = Gmax, iter = 100, nchains = nchains)
#> Inits are missing, calculating...
#> Inits are ready.
#> Settings ready, about to start for ch in 1:nchains.
#> Inits for chain 1 are ready.
#> Triggering C function.
#> 
#> Chain: 1, (1 %)
Chain: 1, (2 %)
Chain: 1, (3 %)
Chain: 1, (4 %)
Chain: 1, (5 %)
Chain: 1, (6 %)
Chain: 1, (7 %)
Chain: 1, (8 %)
Chain: 1, (9 %)
Chain: 1, (10 %)
Chain: 1, (11 %)
Chain: 1, (12 %)
Chain: 1, (13 %)
Chain: 1, (14 %)
Chain: 1, (15 %)
Chain: 1, (16 %)
Chain: 1, (17 %)
Chain: 1, (18 %)
Chain: 1, (19 %)
Chain: 1, (20 %)
Chain: 1, (21 %)
Chain: 1, (22 %)
Chain: 1, (23 %)
Chain: 1, (24 %)
Chain: 1, (25 %)
Chain: 1, (26 %)
Chain: 1, (27 %)
Chain: 1, (28 %)
Chain: 1, (29 %)
Chain: 1, (30 %)
Chain: 1, (31 %)
Chain: 1, (32 %)
Chain: 1, (33 %)
Chain: 1, (34 %)
Chain: 1, (35 %)
Chain: 1, (36 %)
Chain: 1, (37 %)
Chain: 1, (38 %)
Chain: 1, (39 %)
Chain: 1, (40 %)
Chain: 1, (41 %)
Chain: 1, (42 %)
Chain: 1, (43 %)
Chain: 1, (44 %)
Chain: 1, (45 %)
Chain: 1, (46 %)
Chain: 1, (47 %)
Chain: 1, (48 %)
Chain: 1, (49 %)
Chain: 1, (50 %)
Chain: 1, (51 %)
Chain: 1, (52 %)
Chain: 1, (53 %)
Chain: 1, (54 %)
Chain: 1, (55 %)
Chain: 1, (56 %)
Chain: 1, (57 %)
Chain: 1, (58 %)
Chain: 1, (59 %)
Chain: 1, (60 %)
Chain: 1, (61 %)
Chain: 1, (62 %)
Chain: 1, (63 %)
Chain: 1, (64 %)
Chain: 1, (65 %)
Chain: 1, (66 %)
Chain: 1, (67 %)
Chain: 1, (68 %)
Chain: 1, (69 %)
Chain: 1, (70 %)
Chain: 1, (71 %)
Chain: 1, (72 %)
Chain: 1, (73 %)
Chain: 1, (74 %)
Chain: 1, (75 %)
Chain: 1, (76 %)
Chain: 1, (77 %)
Chain: 1, (78 %)
Chain: 1, (79 %)
Chain: 1, (80 %)
Chain: 1, (81 %)
Chain: 1, (82 %)
Chain: 1, (83 %)
Chain: 1, (84 %)
Chain: 1, (85 %)
Chain: 1, (86 %)
Chain: 1, (87 %)
Chain: 1, (88 %)
Chain: 1, (89 %)
Chain: 1, (90 %)
Chain: 1, (91 %)
Chain: 1, (92 %)
Chain: 1, (93 %)
Chain: 1, (94 %)
Chain: 1, (95 %)
Chain: 1, (96 %)
Chain: 1, (97 %)
Chain: 1, (98 %)
Chain: 1, (99 %)
Chain: 1, (100 %)

#> Sampling of chain 1 is completed.
#> Saving the last state of chain 1 is completed.
#> Saving chain 1 is completed.
#> All chains have been sampled.
                  
## Extract the last state
inits <- mcmc$last

## Initialize with inits = Continue in sampling from the last state
mcmc <- clustGLMM(formula = formula, id = id, family = family, 
                  data = longitudinal_mixed_type_data,
                  inits = inits, 
                  G = Gmax, iter = 100, nchains = nchains)
#> Settings ready, about to start for ch in 1:nchains.
#> Inits for chain 1 are ready.
#> Triggering C function.
#> 
#> Chain: 1, (1 %)
Chain: 1, (2 %)
Chain: 1, (3 %)
Chain: 1, (4 %)
Chain: 1, (5 %)
Chain: 1, (6 %)
Chain: 1, (7 %)
Chain: 1, (8 %)
Chain: 1, (9 %)
Chain: 1, (10 %)
Chain: 1, (11 %)
Chain: 1, (12 %)
Chain: 1, (13 %)
Chain: 1, (14 %)
Chain: 1, (15 %)
Chain: 1, (16 %)
Chain: 1, (17 %)
Chain: 1, (18 %)
Chain: 1, (19 %)
Chain: 1, (20 %)
Chain: 1, (21 %)
Chain: 1, (22 %)
Chain: 1, (23 %)
Chain: 1, (24 %)
Chain: 1, (25 %)
Chain: 1, (26 %)
Chain: 1, (27 %)
Chain: 1, (28 %)
Chain: 1, (29 %)
Chain: 1, (30 %)
Chain: 1, (31 %)
Chain: 1, (32 %)
Chain: 1, (33 %)
Chain: 1, (34 %)
Chain: 1, (35 %)
Chain: 1, (36 %)
Chain: 1, (37 %)
Chain: 1, (38 %)
Chain: 1, (39 %)
Chain: 1, (40 %)
Chain: 1, (41 %)
Chain: 1, (42 %)
Chain: 1, (43 %)
Chain: 1, (44 %)
Chain: 1, (45 %)
Chain: 1, (46 %)
Chain: 1, (47 %)
Chain: 1, (48 %)
Chain: 1, (49 %)
Chain: 1, (50 %)
Chain: 1, (51 %)
Chain: 1, (52 %)
Chain: 1, (53 %)
Chain: 1, (54 %)
Chain: 1, (55 %)
Chain: 1, (56 %)
Chain: 1, (57 %)
Chain: 1, (58 %)
Chain: 1, (59 %)
Chain: 1, (60 %)
Chain: 1, (61 %)
Chain: 1, (62 %)
Chain: 1, (63 %)
Chain: 1, (64 %)
Chain: 1, (65 %)
Chain: 1, (66 %)
Chain: 1, (67 %)
Chain: 1, (68 %)
Chain: 1, (69 %)
Chain: 1, (70 %)
Chain: 1, (71 %)
Chain: 1, (72 %)
Chain: 1, (73 %)
Chain: 1, (74 %)
Chain: 1, (75 %)
Chain: 1, (76 %)
Chain: 1, (77 %)
Chain: 1, (78 %)
Chain: 1, (79 %)
Chain: 1, (80 %)
Chain: 1, (81 %)
Chain: 1, (82 %)
Chain: 1, (83 %)
Chain: 1, (84 %)
Chain: 1, (85 %)
Chain: 1, (86 %)
Chain: 1, (87 %)
Chain: 1, (88 %)
Chain: 1, (89 %)
Chain: 1, (90 %)
Chain: 1, (91 %)
Chain: 1, (92 %)
Chain: 1, (93 %)
Chain: 1, (94 %)
Chain: 1, (95 %)
Chain: 1, (96 %)
Chain: 1, (97 %)
Chain: 1, (98 %)
Chain: 1, (99 %)
Chain: 1, (100 %)

#> Sampling of chain 1 is completed.
#> Saving the last state of chain 1 is completed.
#> Saving chain 1 is completed.
#> All chains have been sampled.
                       
## Note that sampling takes time, here we generate too short chains.
## Sample longer chains for more faithful posterior approximation.
## Parallelization for different chains is recommended.
## See the vignette for more examples. 

###-----------------------------------------------------------------------------
###     Sample completely new data
###-----------------------------------------------------------------------------
nnew <- 100
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

set.seed(27818)
ldatanew <- generate_longitudinal_mixed_type_data(n = nnew, n_i = 4,
                     G = 2, timepar = "cross", catpar = "fixed", 
                     catsubjfix = TRUE,
                     xin1 = 1, Sigma = Sigma,
                     Kord = 5, Kcat = 4, kspec_bi_cat = FALSE)
newdata <- ldatanew$data

###-----------------------------------------------------------------------------
###     Evaluate the clustering probabilities for newdata
###-----------------------------------------------------------------------------

probs <- clustering_probabilities_and_deviance(mcmc, id, newdata,
                                               dynamic_prob = FALSE,
                                               start = 1,
                                               end = 100,
                                               thin = 2,
                                               chains = 1:nchains,
                                               tolerance = 1e-7, 
                                               maxiter = 100,
                                               maxnrep = 20, 
                                               NGQ = 1)
#> 
#> Chain: 1, (2 %)
Chain: 1, (4 %)
Chain: 1, (6 %)
Chain: 1, (8 %)
Chain: 1, (10 %)
Chain: 1, (12 %)
Chain: 1, (14 %)
Chain: 1, (16 %)
Chain: 1, (18 %)
Chain: 1, (20 %)
Chain: 1, (22 %)
Chain: 1, (24 %)
Chain: 1, (26 %)
Chain: 1, (28 %)
Chain: 1, (30 %)
Chain: 1, (32 %)
Chain: 1, (34 %)
Chain: 1, (36 %)
Chain: 1, (38 %)
Chain: 1, (40 %)
Chain: 1, (42 %)
Chain: 1, (44 %)
Chain: 1, (46 %)
Chain: 1, (48 %)
Chain: 1, (50 %)
Chain: 1, (52 %)
Chain: 1, (54 %)
Chain: 1, (56 %)
Chain: 1, (58 %)
Chain: 1, (60 %)
Chain: 1, (62 %)
Chain: 1, (64 %)
Chain: 1, (66 %)
Chain: 1, (68 %)
Chain: 1, (70 %)
Chain: 1, (72 %)
Chain: 1, (74 %)
Chain: 1, (76 %)
Chain: 1, (78 %)
Chain: 1, (80 %)
Chain: 1, (82 %)
Chain: 1, (84 %)
Chain: 1, (86 %)
Chain: 1, (88 %)
Chain: 1, (90 %)
Chain: 1, (92 %)
Chain: 1, (94 %)
Chain: 1, (96 %)
Chain: 1, (98 %)
Chain: 1, (100 %)


###-----------------------------------------------------------------------------
###     Final clustering
###-----------------------------------------------------------------------------
certainty <- probs$certainty
clustering <- probs$clustering
# To be classified in the cluster with maximal clustering probability,
# it has to overcome the given threshold, otherwise remains unclassified.
threshold <- 0.8
clustering[certainty < threshold] <- 0

ch <- 1
newdata$clustering <- clustering[,ch]

newdata1 <- newdata[newdata$j==1,]
table(newdata1$clustering, newdata1$g)
#>    
#>      1  2
#>   0 10 18
#>   1 10 10
#>   3 26 26
# row 0: unclassified
# }
```
