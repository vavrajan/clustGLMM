# Change the Output Formats of `clustGLMM`

These functions are used to process or rearrange data between different
formats. The main function `as_coda` converts the `draws` item of the
`clustglmm` class into a
[`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) object
compatible with functions from package coda. Other functions (usually
internally used) allow us to change between different saving formats
(`mcnc$howsave = "list"` or `mcmc$howsave = "data.frame"`).

## Usage

``` r
as_coda(mcmc, burnin=0, thin=1)
from_C_to_list(values, p, settings, yspecd1=NULL, yspecd2=NULL, d2spec = NULL, family)
from_C_to_matrix(values, p, settings, yspecd1=NULL, yspecd2=NULL, d2spec = NULL, family)
from_list_to_C(draws, iterations, p, settings, yspecd1=NULL, yspecd2=NULL,
               d2spec = NULL, family)
from_list_to_coda(mcmc, burnin=0, thin=1)
from_list_to_matrix(mcmc)
from_matrix_to_C(draws, iterations, p, settings)
from_matrix_to_coda(mcmc, burnin=0, thin=1)
from_matrix_to_list(mcmc)
```

## Arguments

- mcmc:

  an object of class `clustglmm`.

- values:

  a vector of values to be stored.

- p:

  a string containing the name of the parameter.

- settings:

  a matrix containing the information about model parameters.

- yspecd1:

  a list of outcome-specific sizes along the first dimension.

- yspecd2:

  a list of outcome-specific sizes along the second dimension.

- d2spec:

  a numeric vector containing the dimension sizes if the parameter has
  different second dimension for every element of the first dimension
  (applicable only for dynamic clustering probabilities).

- family:

  a vector of distributional families named by the outcome name in
  `data`.

- iterations:

  what iterations from a given chain to include.

- draws:

  either a list or a matrix; e.g., an element of `clustGLMM$draws` for a
  certain chain.

- burnin:

  the length of burn-in period (default 0), i.e., the number of initial
  samples to omit.

- thin:

  thinning rate (default 1); only every `thin`th value is kept (must be
  a positive integer).

## Details

Functions with `C` in the name are used internally in
[`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md)
or
[`clustering_probabilities_and_deviance`](https://vavrajan.github.io/clustGLMM/reference/clustering_probabilities_and_deviance.md)
to process data coming from C functions. These functions only convert
only draws from specific chain for a specific parameter, which explains
the arguments `values`, `p` or `draws`. They do not convert the whole
`clustglmm` object.

The functions that change between `howsave = "list"`,
`howsave = "data.frame"` and `coda` output go through every chain and
every saved parameter, thus, convert the whole `clustglmm` object.

Function `as_coda` is a wrapper function that triggers the respective
function `from_..._to_coda` depending on the current value of
`mcmc$howsave`.

## Value

For functions `...to_list`, the output is a carefully structured list.
Each chain `ch in 1:mcmc$chains` is saved in `mcmc$draws[[ch]]`. Then
one can select a specific parameter which may have even other sublists.
For example, group-specific effects for a numeric outcome would be
accessed by `mcmc$draws[[ch]]$beta_num[[g]][[mcmc$Nums[1]]]` where `g`
is the group.

For functions `...to_matrix`, the output is a separate list for each
chain. Each element of this list is a matrix with carefully named
columns forming a `data.frame`. The first column is `m` and contains the
iteration number to which the row belongs to. Other columns are model
parameters. For example, the first effect specific for group 2 for
numeric outcome called "ynum" would be named `"beta_num_ynum(2)[1]"`.
Note that the name of the variable has been pasted after the parameter
name, then follows a group number in parentheses (only if parameter
group-specific) and it ends with the dimensional coordinates in square
brackets.

For functions `...to_C`, the output is a long vector where all the
parameter values are stacked together in the order assumed in the C
functions.

For functions `...to_coda`, the output is a `"mcmc.list"` object from
coda. Then you can use this class and related methods for monitoring,
investigating, plotting. However, a minor warning: the number of
parameters is quite huge and can get overwhelming. Therefore, a set of
plotting functions designed for the different parameter blocks is
provided with this package, see
[`plot.clustglmm`](https://vavrajan.github.io/clustGLMM/reference/plot.clustglmm.md).

## References

Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise
multivariate regression of mixed-type panel data. Stat Comput 34, 46
(2024). https://doi.org/10.1007/s11222-023-10304-5

## Author

Jan Vávra

## See also

[`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md),
[`plot.clustglmm`](https://vavrajan.github.io/clustGLMM/reference/plot.clustglmm.md),
[`clustering_probabilities_and_deviance`](https://vavrajan.github.io/clustGLMM/reference/clustering_probabilities_and_deviance.md)

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
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ynum$group
#> ~x
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ynum$random
#> ~1
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ynum$offset
#> [1] ""
#> 
#> 
#> $ypoi
#> $ypoi$fixed
#> ~f
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ypoi$group
#> ~x
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ypoi$random
#> ~1
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ypoi$offset
#> [1] ""
#> 
#> 
#> $ybin
#> $ybin$fixed
#> ~f
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ybin$group
#> ~x
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ybin$random
#> ~1
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ybin$offset
#> [1] ""
#> 
#> 
#> $yord
#> $yord$fixed
#> ~f
#> <environment: 0x000001cfda58d8f8>
#> 
#> $yord$group
#> ~x
#> <environment: 0x000001cfda58d8f8>
#> 
#> $yord$random
#> ~1
#> <environment: 0x000001cfda58d8f8>
#> 
#> $yord$offset
#> [1] ""
#> 
#> 
#> $ycat
#> $ycat$fixed
#> ~f
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ycat$group
#> ~x
#> <environment: 0x000001cfda58d8f8>
#> 
#> $ycat$random
#> ~1
#> <environment: 0x000001cfda58d8f8>
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
#> Chain: 1, (1 %)Chain: 1, (2 %)Chain: 1, (3 %)Chain: 1, (4 %)Chain: 1, (5 %)Chain: 1, (6 %)Chain: 1, (7 %)Chain: 1, (8 %)Chain: 1, (9 %)Chain: 1, (10 %)Chain: 1, (11 %)Chain: 1, (12 %)Chain: 1, (13 %)Chain: 1, (14 %)Chain: 1, (15 %)Chain: 1, (16 %)Chain: 1, (17 %)Chain: 1, (18 %)Chain: 1, (19 %)Chain: 1, (20 %)Chain: 1, (21 %)Chain: 1, (22 %)Chain: 1, (23 %)Chain: 1, (24 %)Chain: 1, (25 %)Chain: 1, (26 %)Chain: 1, (27 %)Chain: 1, (28 %)Chain: 1, (29 %)Chain: 1, (30 %)Chain: 1, (31 %)Chain: 1, (32 %)Chain: 1, (33 %)Chain: 1, (34 %)Chain: 1, (35 %)Chain: 1, (36 %)Chain: 1, (37 %)Chain: 1, (38 %)Chain: 1, (39 %)Chain: 1, (40 %)Chain: 1, (41 %)Chain: 1, (42 %)Chain: 1, (43 %)Chain: 1, (44 %)Chain: 1, (45 %)Chain: 1, (46 %)Chain: 1, (47 %)Chain: 1, (48 %)Chain: 1, (49 %)Chain: 1, (50 %)Chain: 1, (51 %)Chain: 1, (52 %)Chain: 1, (53 %)Chain: 1, (54 %)Chain: 1, (55 %)Chain: 1, (56 %)Chain: 1, (57 %)Chain: 1, (58 %)Chain: 1, (59 %)Chain: 1, (60 %)Chain: 1, (61 %)Chain: 1, (62 %)Chain: 1, (63 %)Chain: 1, (64 %)Chain: 1, (65 %)Chain: 1, (66 %)Chain: 1, (67 %)Chain: 1, (68 %)Chain: 1, (69 %)Chain: 1, (70 %)Chain: 1, (71 %)Chain: 1, (72 %)Chain: 1, (73 %)Chain: 1, (74 %)Chain: 1, (75 %)Chain: 1, (76 %)Chain: 1, (77 %)Chain: 1, (78 %)Chain: 1, (79 %)Chain: 1, (80 %)Chain: 1, (81 %)Chain: 1, (82 %)Chain: 1, (83 %)Chain: 1, (84 %)Chain: 1, (85 %)Chain: 1, (86 %)Chain: 1, (87 %)Chain: 1, (88 %)Chain: 1, (89 %)Chain: 1, (90 %)Chain: 1, (91 %)Chain: 1, (92 %)Chain: 1, (93 %)Chain: 1, (94 %)Chain: 1, (95 %)Chain: 1, (96 %)Chain: 1, (97 %)Chain: 1, (98 %)Chain: 1, (99 %)Chain: 1, (100 %)
#> Sampling of chain 1 is completed.
#> Saving the last state of chain 1 is completed.
#> Saving chain 1 is completed.
#> Inits for chain 2 are ready.
#> Triggering C function.
#> 
#> Chain: 2, (1 %)Chain: 2, (2 %)Chain: 2, (3 %)Chain: 2, (4 %)Chain: 2, (5 %)Chain: 2, (6 %)Chain: 2, (7 %)Chain: 2, (8 %)Chain: 2, (9 %)Chain: 2, (10 %)Chain: 2, (11 %)Chain: 2, (12 %)Chain: 2, (13 %)Chain: 2, (14 %)Chain: 2, (15 %)Chain: 2, (16 %)Chain: 2, (17 %)Chain: 2, (18 %)Chain: 2, (19 %)Chain: 2, (20 %)Chain: 2, (21 %)Chain: 2, (22 %)Chain: 2, (23 %)Chain: 2, (24 %)Chain: 2, (25 %)Chain: 2, (26 %)Chain: 2, (27 %)Chain: 2, (28 %)Chain: 2, (29 %)Chain: 2, (30 %)Chain: 2, (31 %)Chain: 2, (32 %)Chain: 2, (33 %)Chain: 2, (34 %)Chain: 2, (35 %)Chain: 2, (36 %)Chain: 2, (37 %)Chain: 2, (38 %)Chain: 2, (39 %)Chain: 2, (40 %)Chain: 2, (41 %)Chain: 2, (42 %)Chain: 2, (43 %)Chain: 2, (44 %)Chain: 2, (45 %)Chain: 2, (46 %)Chain: 2, (47 %)Chain: 2, (48 %)Chain: 2, (49 %)Chain: 2, (50 %)Chain: 2, (51 %)Chain: 2, (52 %)Chain: 2, (53 %)Chain: 2, (54 %)Chain: 2, (55 %)Chain: 2, (56 %)Chain: 2, (57 %)Chain: 2, (58 %)Chain: 2, (59 %)Chain: 2, (60 %)Chain: 2, (61 %)Chain: 2, (62 %)Chain: 2, (63 %)Chain: 2, (64 %)Chain: 2, (65 %)Chain: 2, (66 %)Chain: 2, (67 %)Chain: 2, (68 %)Chain: 2, (69 %)Chain: 2, (70 %)Chain: 2, (71 %)Chain: 2, (72 %)Chain: 2, (73 %)Chain: 2, (74 %)Chain: 2, (75 %)Chain: 2, (76 %)Chain: 2, (77 %)Chain: 2, (78 %)Chain: 2, (79 %)Chain: 2, (80 %)Chain: 2, (81 %)Chain: 2, (82 %)Chain: 2, (83 %)Chain: 2, (84 %)Chain: 2, (85 %)Chain: 2, (86 %)Chain: 2, (87 %)Chain: 2, (88 %)Chain: 2, (89 %)Chain: 2, (90 %)Chain: 2, (91 %)Chain: 2, (92 %)Chain: 2, (93 %)Chain: 2, (94 %)Chain: 2, (95 %)Chain: 2, (96 %)Chain: 2, (97 %)Chain: 2, (98 %)Chain: 2, (99 %)Chain: 2, (100 %)
#> Sampling of chain 2 is completed.
#> Saving the last state of chain 2 is completed.
#> Saving chain 2 is completed.
#> All chains have been sampled.
                       
###-----------------------------------------------------------------------------
###     Changing formats
###-----------------------------------------------------------------------------

## How to change from "data.frame"="matrix" to "list":
mcmc <- from_matrix_to_list(mcmc)

## How to change from "list" to "data.frame"="matrix":
mcmc <- from_list_to_matrix(mcmc)

## How to change from "data.frame"="matrix" to mcmc from 'coda' package
codamcmc <- from_matrix_to_coda(mcmc, burnin = 10, thin = 1)

## How to change from any howsave type to mcmc from 'coda' package
codamcmc <- as_coda(mcmc, burnin = 10, thin = 1)
```
