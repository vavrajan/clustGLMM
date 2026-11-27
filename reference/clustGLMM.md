# Sample MCMC for a Mixture of GLMMs of Mixed-type

First, initial values are set based on separate univariate GLMMs for
each individual outcome. Then, an inner `C` function is called to do the
hard work of MCMC sampling (Metropolis within Gibbs sampling). Outputs
are processed and returned in required format (either in nicely
organized lists or one huge data.frame).

## Usage

``` r
clustGLMM(formula, id, family, data, G,  
          varying, save, param, tuning, inits,
          iter = 1000, nchains = 1, standardize = TRUE,
          howsave = c("data.frame", "list", "cmcmc"))
```

## Arguments

- formula:

  a structured list containing formulae for all modeled outcomes. Let
  `y` be a modeled outcome then `formula[[y]]` should be a list
  containing

  `$fixed`

  :   `formula` for effects common to all clusters,

  `$group`

  :   `formula` for cluster-specific effects,

  `$random`

  :   `formula` for `id`-specific effects,

  `$offset`

  :   a column name from `data` to be used as an offset; by default it
      is set to an empty string, i.e., "", which means no offset.

- id:

  a string, the name of the column containing the unit ID.

- family:

  a vector of distributional families named by the outcome name in
  `data`; see below for implemented options.

- data:

  a data.frame containing the outcomes, regressors, and the unit ID
  variable indicated in `id`.

- G:

  an integer, the number of mixture components (maximal number of
  clusters). Can be a vector of integers (of length `nchains`) to
  provide a different value for each chain.

- varying:

  a named vector of Boolean values; declares which parameters should be
  group-specific (`TRUE`) with the following names: `"prec_num"`,
  `"c_ord"`, `"InvSigma"`, `"InvQ"`, `"naY"`.

- save:

  a named vector of Boolean values; declares which model parameters
  should be returned (`TRUE`) and which should be dropped (`FALSE`). See
  also 'Details'.

- param:

  a named list of hyperparameters for prior distributions. See also
  'Details'.

- tuning:

  a list of integer valued and double-valued tuning parameters. See also
  'Details'.

- inits:

  a structured list of initial values; ideally the last generated values
  from a previous run. If `missing`, units are randomly assigned to
  components and parameters are initiated using values obtained based on
  separately estimated univariate models. See also 'Details'.

- iter:

  length of each sampled chain.

- nchains:

  number of chains to be sampled.

- standardize:

  a Boolean whether the numeric outcome and regressors should be scaled.

- howsave:

  a string declaring how the sampled chains should be returned as
  `draws`. One of:

  `"list"`

  :   `draws` is a list indexed by chain numbers, parameters will be
      saved in a structured list in which it is easy to orient,

  `"data.frame"`

  :   `draws` is a large data.frame, where one row corresponds to one
      iteration of one specific chain and columns are all model
      parameters,

  `"cmcmc"`

  :   `draws` is a list indexed by chain numbers and contains the
      condensed output of the inner C function (not recommended!).

## Details

Column `id` in the `data` data.frame contains different IDs. For
simplicity, these IDs are tied with numbers from 1 to total number of
unique IDs. Then, unit-specific quantities are referred to by these
integers instead of full IDs.

Outcomes in `data` are expected to belong to one of the following
families:

- numeric:

  `"num"` or `"gaussian"`; double values without any restriction,

- count:

  `"poi"` or `"poisson"`; integer values 0, 1, 2, ...,

- binary:

  `"bin"` or `"bernoulli"`; factor of two levels, first will be baseline
  (0), other (1),

- ordinal:

  `"ord"` or `"cumulative"`; (ordered) factor of at least 3 levels, if
  not apriori ordered, ordering will correspond to the one by
  `factor(, ordered = TRUE)`,

- categorical:

  `"cat"` or `"categorical"`; factor of at least 3 levels, first
  category will be baseline.

Only these 5 types of outcome families are supported. The link functions
cannot be changed. `family` is internally relabeled to contain only
`"num"`, `"poi"`, `"bin"`, `"ord"`, `"cat"` and reordered to follow this
order of types. The same ordering of outcomes is applied to the
covariance matrix `Sigma` of random effects.

A named vector `varying` contains a Boolean for each model parameter
that may be considered cluster-specific. The defaults are set by
[`default_varying`](https://vavrajan.github.io/clustGLMM/reference/default_varying.md).

A named vector `save` contains a Boolean for each model parameter
including some transformed parameters. For an overview on all model
parameters including the default values used, see
[`default_save`](https://vavrajan.github.io/clustGLMM/reference/default_varying.md).

A named list `param` contains the values of the hyperparameters needed
to set up the priors. The default values are set up by
[`default_param`](https://vavrajan.github.io/clustGLMM/reference/default_varying.md).
By adjusting the values of `ae` and `be` different levels of sparsity
can be achieved. The lower `ae/be` is, the sparser the results expected.

The list of initial values `inits` is convenient to supply in case one
is interested in continuing the sampling using the last sampled values
as initialization, e.g., if the chains have not yet converged. If
`inits` is missing, then a completely new set of initial values for each
of the sampled chains is generated in the following way:

1\. Units are randomly (uniformly) distributed to `G` components and `w`
is initiated with `1/G` for all components.

2\. For each numeric outcome, a regression model with combined `fixed`
and `group` formula is fit to set up initial values for fixed parameters
`beta_num_fixed`. Then, the model is fitted again (if possible) for each
component separately to determine the initial values for group-specific
effects `beta_num`. Depending on the choice of group-specificity for
`prec_num` the initial values are chosen as the inverse of the error
variance estimate of the corresponding model.

3\. Effects for count outcomes are initialized in a similar way using
`glm(, family = poisson(link = "log"))` instead.

4\. Effects for binary outcomes are initialized in a similar way using
`glm(, family = binomial(link = "logit"))` instead.

5\. Effects and ordered intercepts for ordinal outcomes are initialized
in a similar way using `polr(, Hess = TRUE)` from MASS instead.

6\. Effects for general categorical outcomes are initialized in a
similar way using `multinom()` from nnet instead.

7\. Matrix `InvSigma` is initialized by diagonal matrices with gamma
distributed (`rate = 1/init_sd_b^2`) diagonal.

8\. Matrix `InvQ` is initialized as inverse of (averaged and scaled)
`InvSigma`.

9\. Random effects are sampled randomly from a centered normal
distribution with `sd = param$init_sd_b`.

The algorithm also requires a set of tuning parameters in `tuning` that
tune the sampling. For clarity, the parameters are split into two
sublists depending on the expected type (`integer` or `double`). If
`tuning` is missing, then the default values are set by
[`default_tuning`](https://vavrajan.github.io/clustGLMM/reference/default_varying.md).

## Value

An object of class `"clustglmm"` where the most important elements
contained are:

- draws:

  a list indexed by chain with generated samples saved either in a large
  data.frame (`howsave = "data.frame"`) or a carefully structured list
  (`howsave = "list"`).

- param_names:

  individual parameter (column) names in a structured list by general
  parameter name.

- accept:

  a list indexed by chain with a `data.frame` of proposal acceptances
  with values in interval \[0, 1\] (could be in-between if multiple
  proposals per step are considered by `tuning$integer$times_proposal`).

In addition the object contains the inputs given to the function:
`iter`, `nchains`, `family` (sorted by type), `formula`, `G` (a vector
of Gs for all chains), `varying`, `save`, `howsave`, `param`, `tuning`,
`standardize`.

Furthermore elements consisting of newly defined quantities are:

- iterations:

  a list per each chain containing iterations available,

- modeGplus:

  a vector of the most frequent number of non-empty components per
  chain,

- sameclusters:

  a Boolean declaring whether the cluster labels in all sampled chains
  are the same (useful for plotting),

- clusters:

  a list of ordered non-empty components for each sampled chain
  (`modeGplus` most frequently appearing cluster labels),

- clustering:

  a matrix containing estimated cluster labels based on all sampled
  allocation indicators for each unit (row) and sampled chain (column),

- certainty:

  a matrix containing the probability of assignment into the cluster
  specified in `clustering` based on all sampled allocation indicators
  for each unit (row) and sampled chain (column),

- inits:

  a structured list of parameter values used as initial values,

- last:

  a structured list in the same format as `inits` to be used as initial
  values for additional sampling,

- InitType:

  string describing whether the initial values were given or created,

- dim_accept:

  a list containing a vector of `accept` dimensions per parameter for
  each chain,

- call:

  summary of what has been fitted,

- iter:

  length of each sampled chain,

- n:

  total number of units (unique IDs),

- numbered_unique_ids:

  numbers `1:n` named by unique values from column `id`,

- Nums:

  names of numerical outcomes,

- Pois:

  names of count outcomes,

- Bins:

  names of binary outcomes,

- Ords:

  names of ordinal outcomes,

- Cats:

  names of categorical outcomes,

- nfix:

  a vector named by outcomes that declares the dimension size of fixed
  effects common to all clusters,

- ngrp:

  a vector named by outcomes that declares the dimension size of
  group-specific effects,

- nran:

  a vector named by outcomes that declares the dimension size of random
  effects specific to each unit,

- noff:

  a vector named by outcomes that declares whether an offset is included
  in formula,

- totnran:

  total dimension of random effects, `totnran = sum(nran)`,

- lfixnames:

  a list of names of fixed effects common to all clusters for each
  modeled outcome,

- lgrpnames:

  a list of names of group-specific effects for each modeled outcome,

- lrannames:

  a list of names of random effects specific to each unit for each
  modeled outcome,

- loffnames:

  a list of names of offset variables for each modeled outcomes,

- fixnames:

  unlisted `lfixnames`,

- grpnames:

  unlisted `lgrpnames`,

- rannames:

  unlisted `lrannames`,

- offnames:

  unlisted `loffnames`,

- nY:

  number of outcomes per type,

- Kord:

  number of categories for each ordinal outcome,

- Kcat:

  number of categories for each general categorical outcome,

- settings:

  a list (for each chain) of matrices where each row belongs to one
  model parameter and describes its properties such as: 1. is it
  saved?, 2. is it group-specific?, 3. is any dimension size dependent
  on outcome?, 4. on which type of outcome it depends?, 5. is it point,
  vector or matrix? 5. what are the dimension sizes?, 6. is the matrix
  symmetric?, 7. does the matrix contain diagonal?, 8. what is the total
  size of the parameter (including group-specificity)?; very useful
  matrix for organizing the samples and writing automated plotting
  functions,

- yspecd1:

  a list of model parameters, which have their first dimension size
  dependent on outcomes and how large it is, e,g, fixed-effects
  `beta_num_fix` are vectors that can be of different size among
  outcomes under different formulae,

- yspecd2:

  a list of matrix-valued model parameters, which have their second
  dimension size dependent on outcomes and how large it is, e.g.
  fixed-effects `beta_cat_fix` for categorical outcome `y in Cats` are
  matrices (`Kcat[y]` \* `nfix[y]`) with both dimensions dependent on
  the outcome `y`,

- isYna:

  a vector of Booleans that determine which outcome values are not
  observed, outcomes are ordered by
  `Ys = c(Nums, Pois, Bins, Ords, Cats)`.

- centers:

  a named vector of centers of numeric outcomes and regressors.

- scales:

  a named vector of scales of numeric outcomes and regressors.

## References

Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise
multivariate regression of mixed-type panel data. Stat Comput 34, 46
(2024). https://doi.org/10.1007/s11222-023-10304-5

## Author

Jan Vávra

## See also

[`generate_longitudinal_mixed_type_data`](https://vavrajan.github.io/clustGLMM/reference/generate_longitudinal_mixed_type_data.md),
[`plot.clustglmm`](https://vavrajan.github.io/clustGLMM/reference/plot.clustglmm.md),
[`default_varying`](https://vavrajan.github.io/clustGLMM/reference/default_varying.md),
[`default_save`](https://vavrajan.github.io/clustGLMM/reference/default_varying.md),
[`default_param`](https://vavrajan.github.io/clustGLMM/reference/default_varying.md),
[`default_tuning`](https://vavrajan.github.io/clustGLMM/reference/default_varying.md)

## Examples

``` r

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
formula
#> $ynum
#> $ynum$fixed
#> ~f
#> <environment: 0x64f208573cd0>
#> 
#> $ynum$group
#> ~x
#> <environment: 0x64f208573cd0>
#> 
#> $ynum$random
#> ~1
#> <environment: 0x64f208573cd0>
#> 
#> $ynum$offset
#> [1] ""
#> 
#> 
#> $ypoi
#> $ypoi$fixed
#> ~f
#> <environment: 0x64f208573cd0>
#> 
#> $ypoi$group
#> ~x
#> <environment: 0x64f208573cd0>
#> 
#> $ypoi$random
#> ~1
#> <environment: 0x64f208573cd0>
#> 
#> $ypoi$offset
#> [1] ""
#> 
#> 
#> $ybin
#> $ybin$fixed
#> ~f
#> <environment: 0x64f208573cd0>
#> 
#> $ybin$group
#> ~x
#> <environment: 0x64f208573cd0>
#> 
#> $ybin$random
#> ~1
#> <environment: 0x64f208573cd0>
#> 
#> $ybin$offset
#> [1] ""
#> 
#> 
#> $yord
#> $yord$fixed
#> ~f
#> <environment: 0x64f208573cd0>
#> 
#> $yord$group
#> ~x
#> <environment: 0x64f208573cd0>
#> 
#> $yord$random
#> ~1
#> <environment: 0x64f208573cd0>
#> 
#> $yord$offset
#> [1] ""
#> 
#> 
#> $ycat
#> $ycat$fixed
#> ~f
#> <environment: 0x64f208573cd0>
#> 
#> $ycat$group
#> ~x
#> <environment: 0x64f208573cd0>
#> 
#> $ycat$random
#> ~1
#> <environment: 0x64f208573cd0>
#> 
#> $ycat$offset
#> [1] ""
#> 
#> 
# intercept term will be group-specific 
# since it appears in $group which has a preference over $fixed

## All other inputs need to be prepared,
## unless the default values are desired (varying, save, howsave, param, tuning).
## See the vignette for how these can be adjusted.

###-----------------------------------------------------------------------------
###     Calling clustGLMM()
###-----------------------------------------------------------------------------

Gmax <- 5     # maximal number of mixture components to be considered
nchains <- 2  # number of chains to be sampled 

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
#> Inits for chain 2 are ready.
#> Triggering C function.
#> 
#> Chain: 2, (1 %)
Chain: 2, (2 %)
Chain: 2, (3 %)
Chain: 2, (4 %)
Chain: 2, (5 %)
Chain: 2, (6 %)
Chain: 2, (7 %)
Chain: 2, (8 %)
Chain: 2, (9 %)
Chain: 2, (10 %)
Chain: 2, (11 %)
Chain: 2, (12 %)
Chain: 2, (13 %)
Chain: 2, (14 %)
Chain: 2, (15 %)
Chain: 2, (16 %)
Chain: 2, (17 %)
Chain: 2, (18 %)
Chain: 2, (19 %)
Chain: 2, (20 %)
Chain: 2, (21 %)
Chain: 2, (22 %)
Chain: 2, (23 %)
Chain: 2, (24 %)
Chain: 2, (25 %)
Chain: 2, (26 %)
Chain: 2, (27 %)
Chain: 2, (28 %)
Chain: 2, (29 %)
Chain: 2, (30 %)
Chain: 2, (31 %)
Chain: 2, (32 %)
Chain: 2, (33 %)
Chain: 2, (34 %)
Chain: 2, (35 %)
Chain: 2, (36 %)
Chain: 2, (37 %)
Chain: 2, (38 %)
Chain: 2, (39 %)
Chain: 2, (40 %)
Chain: 2, (41 %)
Chain: 2, (42 %)
Chain: 2, (43 %)
Chain: 2, (44 %)
Chain: 2, (45 %)
Chain: 2, (46 %)
Chain: 2, (47 %)
Chain: 2, (48 %)
Chain: 2, (49 %)
Chain: 2, (50 %)
Chain: 2, (51 %)
Chain: 2, (52 %)
Chain: 2, (53 %)
Chain: 2, (54 %)
Chain: 2, (55 %)
Chain: 2, (56 %)
Chain: 2, (57 %)
Chain: 2, (58 %)
Chain: 2, (59 %)
Chain: 2, (60 %)
Chain: 2, (61 %)
Chain: 2, (62 %)
Chain: 2, (63 %)
Chain: 2, (64 %)
Chain: 2, (65 %)
Chain: 2, (66 %)
Chain: 2, (67 %)
Chain: 2, (68 %)
Chain: 2, (69 %)
Chain: 2, (70 %)
Chain: 2, (71 %)
Chain: 2, (72 %)
Chain: 2, (73 %)
Chain: 2, (74 %)
Chain: 2, (75 %)
Chain: 2, (76 %)
Chain: 2, (77 %)
Chain: 2, (78 %)
Chain: 2, (79 %)
Chain: 2, (80 %)
Chain: 2, (81 %)
Chain: 2, (82 %)
Chain: 2, (83 %)
Chain: 2, (84 %)
Chain: 2, (85 %)
Chain: 2, (86 %)
Chain: 2, (87 %)
Chain: 2, (88 %)
Chain: 2, (89 %)
Chain: 2, (90 %)
Chain: 2, (91 %)
Chain: 2, (92 %)
Chain: 2, (93 %)
Chain: 2, (94 %)
Chain: 2, (95 %)
Chain: 2, (96 %)
Chain: 2, (97 %)
Chain: 2, (98 %)
Chain: 2, (99 %)
Chain: 2, (100 %)

#> Sampling of chain 2 is completed.
#> Saving the last state of chain 2 is completed.
#> Saving chain 2 is completed.
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
#> Inits for chain 2 are ready.
#> Triggering C function.
#> 
#> Chain: 2, (1 %)
Chain: 2, (2 %)
Chain: 2, (3 %)
Chain: 2, (4 %)
Chain: 2, (5 %)
Chain: 2, (6 %)
Chain: 2, (7 %)
Chain: 2, (8 %)
Chain: 2, (9 %)
Chain: 2, (10 %)
Chain: 2, (11 %)
Chain: 2, (12 %)
Chain: 2, (13 %)
Chain: 2, (14 %)
Chain: 2, (15 %)
Chain: 2, (16 %)
Chain: 2, (17 %)
Chain: 2, (18 %)
Chain: 2, (19 %)
Chain: 2, (20 %)
Chain: 2, (21 %)
Chain: 2, (22 %)
Chain: 2, (23 %)
Chain: 2, (24 %)
Chain: 2, (25 %)
Chain: 2, (26 %)
Chain: 2, (27 %)
Chain: 2, (28 %)
Chain: 2, (29 %)
Chain: 2, (30 %)
Chain: 2, (31 %)
Chain: 2, (32 %)
Chain: 2, (33 %)
Chain: 2, (34 %)
Chain: 2, (35 %)
Chain: 2, (36 %)
Chain: 2, (37 %)
Chain: 2, (38 %)
Chain: 2, (39 %)
Chain: 2, (40 %)
Chain: 2, (41 %)
Chain: 2, (42 %)
Chain: 2, (43 %)
Chain: 2, (44 %)
Chain: 2, (45 %)
Chain: 2, (46 %)
Chain: 2, (47 %)
Chain: 2, (48 %)
Chain: 2, (49 %)
Chain: 2, (50 %)
Chain: 2, (51 %)
Chain: 2, (52 %)
Chain: 2, (53 %)
Chain: 2, (54 %)
Chain: 2, (55 %)
Chain: 2, (56 %)
Chain: 2, (57 %)
Chain: 2, (58 %)
Chain: 2, (59 %)
Chain: 2, (60 %)
Chain: 2, (61 %)
Chain: 2, (62 %)
Chain: 2, (63 %)
Chain: 2, (64 %)
Chain: 2, (65 %)
Chain: 2, (66 %)
Chain: 2, (67 %)
Chain: 2, (68 %)
Chain: 2, (69 %)
Chain: 2, (70 %)
Chain: 2, (71 %)
Chain: 2, (72 %)
Chain: 2, (73 %)
Chain: 2, (74 %)
Chain: 2, (75 %)
Chain: 2, (76 %)
Chain: 2, (77 %)
Chain: 2, (78 %)
Chain: 2, (79 %)
Chain: 2, (80 %)
Chain: 2, (81 %)
Chain: 2, (82 %)
Chain: 2, (83 %)
Chain: 2, (84 %)
Chain: 2, (85 %)
Chain: 2, (86 %)
Chain: 2, (87 %)
Chain: 2, (88 %)
Chain: 2, (89 %)
Chain: 2, (90 %)
Chain: 2, (91 %)
Chain: 2, (92 %)
Chain: 2, (93 %)
Chain: 2, (94 %)
Chain: 2, (95 %)
Chain: 2, (96 %)
Chain: 2, (97 %)
Chain: 2, (98 %)
Chain: 2, (99 %)
Chain: 2, (100 %)

#> Sampling of chain 2 is completed.
#> Saving the last state of chain 2 is completed.
#> Saving chain 2 is completed.
#> All chains have been sampled.
                       
## Note that sampling takes time, here we generate too short chains.
## Sample longer chains for more faithful posterior approximation.
## Parallelization for different chains is recommended.
## See the vignette for more examples. 
```
