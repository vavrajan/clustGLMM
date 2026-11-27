# Setting Default Values in `clustGLMM`

Sets default values for the hyperparameters of the priors
`default_param`, indicators of group-specificity `default_varying`,
integer and double valued tuning parameters `default_tuning` and
indicators of what variables should be saved `default_save`.

## Usage

``` r
default_param(sparse = TRUE, ...)
default_save(...)
default_tuning(experience = FALSE, ...)
default_varying(...)
```

## Arguments

- sparse:

  a logical, if `TRUE` the rate and shape hyperparameters of `e0` are
  set to values encouraging sparsity in the weights

- experience:

  a logical, indicating if the default values are changed to values
  proven suitable from empirical experience

- ...:

  named scalars to override default values

## Details

If arguments `param`, `save`, `tuning` and `varying` of
[`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md)
are not specified, these functions are used to create the default
values. They can also be used to create prototypical `param`, `varying`,
`tuning` and `save` argument specifications and then adjust them.

## Value

The output of `default_varying` is a named vector that contains a
Boolean for each model parameter that may be considered group-specific.
The default values are:

- `varying["prec_num"] = TRUE`:

  group-specific precisions for numeric outcomes,

- `varying["c_ord"] = TRUE`:

  group-specific ordered intercepts for ordinal outcomes,

- `varying["InvSigma"] = FALSE`:

  group-invariant combined precision matrix for all random effects,

- `varying["InvQ"] = FALSE`:

  group-invariant precision matrix used in prior for `InvSigma`,

- `varying["naY"] = FALSE`:

  group-invariant imputed missing outcome values.

The output of `default_save` is a named vector that contains a Boolean
for each model parameter including some transformed parameters
indicating if they are saved and included in the return object. The
default values are:

- `save["beta_num_fix"] = TRUE`:

  fixed effects for numeric outcomes,

- `save["beta_num"] = TRUE`:

  group-specific effects for numeric outcomes,

- `save["prec_num"] = TRUE`:

  precisions for numeric outcomes,

- `save["sd_num"] = TRUE`:

  standard deviations for numeric outcomes (transformed from
  `prec_num`),

- `save["var_num"] = TRUE`:

  variances for numeric outcomes (transformed from `prec_num`),

- `save["beta_poi_fix"] = TRUE`:

  fixed effects for count outcomes,

- `save["beta_poi"] = TRUE`:

  group-specific effects for count outcomes,

- `save["beta_bin_fix"] = TRUE`:

  fixed effects for binary outcomes,

- `save["beta_bin"] = TRUE`:

  group-specific effects for binary outcomes,

- `save["beta_ord_fix"] = TRUE`:

  fixed effects for ordinal outcomes,

- `save["beta_ord"] = TRUE`:

  group-specific effects for ordinal outcomes,

- `save["c_ord"] = TRUE`:

  ordered intercepts for ordinal outcomes,

- `save["a_ord"] = FALSE`:

  transformed ordered intercepts for ordinal outcomes into unrestricted
  variables,

- `save["pi_ord"] = TRUE`:

  transformed `a_ord` into probabilities,

- `save["beta_cat_fix"] = TRUE`:

  fixed effects for general categorical outcomes,

- `save["beta_cat"] = TRUE`:

  group-specific effects for general categorical outcomes,

- `save["InvSigma"] = TRUE`:

  combined precision matrix for all random effects,

- `save["Sigma"] = TRUE`:

  combined covariance matrix for all random effects (transformed from
  `InvSigma`),

- `save["sdSigma"] = TRUE`:

  standard deviations for all random effects (transformed from
  `InvSigma`),

- `save["corSigma"] = TRUE`:

  correlations between all random effects (transformed from `InvSigma`),

- `save["detInvSigma"] = FALSE`:

  determinant of `InvSigma` matrix,

- `save["InvQ"] = FALSE`:

  precision matrix used in prior for `InvSigma`,

- `save["Q"] = FALSE`:

  covariance matrix used in prior for `InvSigma` (transformed from
  `InvQ`),

- `save["detInvQ"] = FALSE`:

  determinant of `InvQ` matrix,

- `save["b"] = FALSE`:

  a matrix of random effects for all outcomes and units,

- `save["w"] = TRUE`:

  mixture weights of length `G` that sum up to 1,

- `save["ng"] = TRUE`:

  cluster occupancies of length `G` that sum up to number of units,

- `save["loglik"] = TRUE`:

  NOT a log-likelihood but a very crude approximation to it created as
  an easy to evaluate byproduct of full-conditional cluster allocation
  probabilities,

- `save["pUig"] = FALSE`:

  a matrix of full-conditional cluster allocation probabilities,

- `save["U"] = TRUE`:

  a vector of current allocations of units into clusters (values in
  `1, ..., G`),

- `save["Gplus"] = TRUE`:

  a number of non-empty components,

- `save["e0"] = TRUE`:

  precision parameter for prior of `w`,

- `save["naY"] = FALSE`:

  imputed missing outcome values.

The output of `default_param` is a named list that contains values of
the hyperparameters needed to set up the priors. The default values are:

- `param[["InvSigma_df"]] = 1`:

  prior degrees of freedom for Wishart prior for `InvSigma` (dimension
  of this matrix is automatically added to this value),

- `param[["InvQ_df"]] = 1`:

  prior degrees of freedom for Wishart prior for `InvQ` (dimension of
  this matrix is automatically added to this value),

- `param[["prec_num_shp"]] = 1`:

  prior shape for precision parameter `prec_num` for numeric outcomes,

- `param[["prec_num_rte"]] = 10`:

  prior rate for precision parameter `prec_num` for numeric outcomes,

- `param[["beta_num_fix_sd"]] = 1`:

  standard deviation for fixed effects for numeric outcomes,

- `param[["beta_num_sd"]] = 1`:

  standard deviation for group-specific effects for numeric outcomes,

- `param[["beta_poi_fix_sd"]] = 1`:

  standard deviation for fixed effects for count outcomes,

- `param[["beta_poi_sd"]] = 1`:

  standard deviation for group-specific effects for count outcomes,

- `param[["beta_bin_fix_sd"]] = 1`:

  standard deviation for fixed effects for binary outcomes,

- `param[["beta_bin_sd"]] = 1`:

  standard deviation for group-specific effects for binary outcomes,

- `param[["beta_ord_fix_sd"]] = 1`:

  standard deviation for fixed effects for ordinal outcomes,

- `param[["beta_ord_sd"]] = 1`:

  standard deviation for group-specific effects for ordinal outcomes,

- `param[["api_prior"]] = 1`:

  prior parameter for Dirichlet distribution of `pi_ord` (same for all
  levels),

- `param[["beta_cat_fix_sd"]] = 1`:

  standard deviation for fixed effects for general categorical outcomes,

- `param[["beta_cat_sd"]] = 1`:

  standard deviation for group-specific effects for general categorical
  outcomes,

- `param[["init_b_sd"]] = 0.01`:

  standard deviation for sampling initial values of random effects `b`,

- `param[["e0_shp"]] = 1`:

  prior shape of `e0`,

- `param[["e0_rte"]] = 100`:

  prior rate of `e0`,

- `param[["InvV"]] = 0.01`:

  diagonal element of the prior inverse scale matrix for `InvQ`.

By adjusting the values of `e0_shp` and `e0_rte` different levels of
sparsity can be achieved. The lower `e0_shp/e0_rte` is, the more sparse
results are expected. With `sparse = TRUE` the prior distribution is a
gamma distribution with `e0_shp` = 1 and `e0_rte` = 100. Otherwise
`e0_shp` = 4 and `e0_rte` = 1 will be used, which prefers non-sparse
solutions.

The output of `default_tuning` is a named list of parameters that tune
the sampling algorithm. For clarity, the parameters are split into two
sublists depending on the expected type (`integer` or `double`). The
default values are:

- `$integer$freq_proposal_update = 10`:

  every `freq_proposal_update` iterations change the proposal
  distributions for Metropolis steps,

- `$integer$times_proposal = 1`:

  how many times a new state is proposed within one Metropolis step,

- `$integer$maxiter = 25`:

  maximum iterations allowed during Newton-Raphson update of proposal
  distribution,

- `$integer$maxnrep = 100`:

  maximum number of repetitions of Newton-Raphson with a different start
  location,

- `$integer$kspec_bi_cat = 0`:

  Boolean as integer, if `TRUE` then each categorical outcome has a
  different set of random effects for each of `Kcat-1` levels,

- `$double$const_proposal_beta_poi_fix = 1.0`:

  constant by which the proposal direction for `beta_poi_fix` is
  multiplied (const_proposal^2 to variance),

- `$double$const_proposal_beta_poi = 1.0`:

  constant by which the proposal direction for `beta_poi` is multiplied
  (const_proposal^2 to variance),

- `$double$const_proposal_beta_bin_fix = 1.0`:

  constant by which the proposal direction for `beta_bin_fix` is
  multiplied (const_proposal^2 to variance),

- `$double$const_proposal_beta_bin = 1.0`:

  constant by which the proposal direction for `beta_bin` is multiplied
  (const_proposal^2 to variance),

- `$double$const_proposal_beta_ord_fix = 1.0`:

  constant by which the proposal direction for `beta_ord_fix` is
  multiplied (const_proposal^2 to variance),

- `$double$const_proposal_beta_ord = 1.0`:

  constant by which the proposal direction for `beta_ord` is multiplied
  (const_proposal^2 to variance),

- `$double$const_proposal_beta_cat_fix = 1.0`:

  constant by which the proposal direction for `beta_cat_fix` is
  multiplied (const_proposal^2 to variance),

- `$double$const_proposal_beta_cat = 1.0`:

  constant by which the proposal direction for `beta_cat` is multiplied
  (const_proposal^2 to variance),

- `$double$const_proposal_a_ord = 1.0`:

  constant by which the proposal direction for `a_ord` is multiplied
  (const_proposal^2 to variance),

- `$double$const_proposal_b = 1.0`:

  constant by which the proposal direction for `b` is multiplied
  (const_proposal^2 to variance),

- `$double$const_proposal_e0 = 1.0`:

  constant by which the proposal direction for `e0` is multiplied
  (const_proposal^2 to variance),

- `$double$df_proposal_beta_poi_fix = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `beta_poi_fix`,

- `$double$df_proposal_beta_poi = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `beta_poi`,

- `$double$df_proposal_beta_bin_fix = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `beta_bin_fix`,

- `$double$df_proposal_beta_bin = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `beta_bin`,

- `$double$df_proposal_beta_ord_fix = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `beta_ord_fix`,

- `$double$df_proposal_beta_ord = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `beta_ord`,

- `$double$df_proposal_beta_cat_fix = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `beta_cat_fix`,

- `$double$df_proposal_beta_cat = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `beta_cat`,

- `$double$df_proposal_a_ord = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `a_ord`,

- `$double$df_proposal_b = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `b`,

- `$double$df_proposal_e0 = Inf`:

  the degrees of freedom for multivariate Student-t proposal (`Inf`
  results in standard multivariate normal proposal) of `e0`,

- `$double$tolerance = 1e-07`:

  tolerance when computing the norm of the shift within Newton-Raphson.

With `experience = TRUE` the chosen values will be set to those proven
suitable from empirical experience.

## References

Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise
multivariate regression of mixed-type panel data. Stat Comput 34, 46
(2024). https://doi.org/10.1007/s11222-023-10304-5

## Author

Jan Vávra

## See also

[`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md)

## Examples

``` r
(param <- default_param())
#> $InvSigma_df
#> [1] 1
#> 
#> $InvQ_df
#> [1] 1
#> 
#> $prec_num_shp
#> [1] 1
#> 
#> $prec_num_rte
#> [1] 10
#> 
#> $beta_num_fix_sd
#> [1] 1
#> 
#> $beta_num_sd
#> [1] 1
#> 
#> $beta_poi_fix_sd
#> [1] 1
#> 
#> $beta_poi_sd
#> [1] 1
#> 
#> $beta_bin_fix_sd
#> [1] 1
#> 
#> $beta_bin_sd
#> [1] 1
#> 
#> $beta_ord_fix_sd
#> [1] 1
#> 
#> $beta_ord_sd
#> [1] 1
#> 
#> $api_prior
#> [1] 1
#> 
#> $beta_cat_fix_sd
#> [1] 1
#> 
#> $beta_cat_sd
#> [1] 1
#> 
#> $init_b_sd
#> [1] 0.01
#> 
#> $e0_shp
#> [1] 1
#> 
#> $e0_rte
#> [1] 100
#> 
#> $InvV
#> [1] 0.01
#> 
(param <- default_param(sparse = FALSE))
#> $InvSigma_df
#> [1] 1
#> 
#> $InvQ_df
#> [1] 1
#> 
#> $prec_num_shp
#> [1] 1
#> 
#> $prec_num_rte
#> [1] 10
#> 
#> $beta_num_fix_sd
#> [1] 1
#> 
#> $beta_num_sd
#> [1] 1
#> 
#> $beta_poi_fix_sd
#> [1] 1
#> 
#> $beta_poi_sd
#> [1] 1
#> 
#> $beta_bin_fix_sd
#> [1] 1
#> 
#> $beta_bin_sd
#> [1] 1
#> 
#> $beta_ord_fix_sd
#> [1] 1
#> 
#> $beta_ord_sd
#> [1] 1
#> 
#> $api_prior
#> [1] 1
#> 
#> $beta_cat_fix_sd
#> [1] 1
#> 
#> $beta_cat_sd
#> [1] 1
#> 
#> $init_b_sd
#> [1] 0.01
#> 
#> $e0_shp
#> [1] 4
#> 
#> $e0_rte
#> [1] 1
#> 
#> $InvV
#> [1] 0.01
#> 
(save <- default_save())
#> beta_num_fix     beta_num     prec_num       sd_num      var_num beta_poi_fix 
#>         TRUE         TRUE         TRUE         TRUE         TRUE         TRUE 
#>     beta_poi beta_bin_fix     beta_bin beta_ord_fix     beta_ord        c_ord 
#>         TRUE         TRUE         TRUE         TRUE         TRUE         TRUE 
#>        a_ord       pi_ord beta_cat_fix     beta_cat     InvSigma        Sigma 
#>        FALSE         TRUE         TRUE         TRUE         TRUE         TRUE 
#>      sdSigma     corSigma  detInvSigma         InvQ            Q      detInvQ 
#>         TRUE         TRUE        FALSE        FALSE        FALSE        FALSE 
#>            b            w           ng       loglik         pUig            U 
#>        FALSE         TRUE         TRUE         TRUE        FALSE         TRUE 
#>        Gplus           e0          naY 
#>         TRUE         TRUE        FALSE 
(tuning <- default_tuning())
#> $integer
#> $integer$freq_proposal_update
#> [1] 10
#> 
#> $integer$times_proposal
#> [1] 1
#> 
#> $integer$maxiter
#> [1] 25
#> 
#> $integer$maxnrep
#> [1] 100
#> 
#> $integer$kspec_bi_cat
#> [1] 0
#> 
#> 
#> $double
#> $double$const_proposal_beta_poi_fix
#> [1] 1
#> 
#> $double$const_proposal_beta_poi
#> [1] 1
#> 
#> $double$const_proposal_beta_bin_fix
#> [1] 1
#> 
#> $double$const_proposal_beta_bin
#> [1] 1
#> 
#> $double$const_proposal_beta_ord_fix
#> [1] 1
#> 
#> $double$const_proposal_beta_ord
#> [1] 1
#> 
#> $double$const_proposal_beta_cat_fix
#> [1] 1
#> 
#> $double$const_proposal_beta_cat
#> [1] 1
#> 
#> $double$const_proposal_a_ord
#> [1] 1
#> 
#> $double$const_proposal_b
#> [1] 1
#> 
#> $double$const_proposal_e0
#> [1] 1
#> 
#> $double$df_proposal_beta_poi_fix
#> [1] Inf
#> 
#> $double$df_proposal_beta_poi
#> [1] Inf
#> 
#> $double$df_proposal_beta_bin_fix
#> [1] Inf
#> 
#> $double$df_proposal_beta_bin
#> [1] Inf
#> 
#> $double$df_proposal_beta_ord_fix
#> [1] Inf
#> 
#> $double$df_proposal_beta_ord
#> [1] Inf
#> 
#> $double$df_proposal_beta_cat_fix
#> [1] Inf
#> 
#> $double$df_proposal_beta_cat
#> [1] Inf
#> 
#> $double$df_proposal_a_ord
#> [1] Inf
#> 
#> $double$df_proposal_b
#> [1] Inf
#> 
#> $double$df_proposal_e0
#> [1] Inf
#> 
#> $double$tolerance
#> [1] 1e-07
#> 
#> 
(tuning <- default_tuning(experience = TRUE))
#> $integer
#> $integer$freq_proposal_update
#> [1] 10
#> 
#> $integer$times_proposal
#> [1] 10
#> 
#> $integer$maxiter
#> [1] 25
#> 
#> $integer$maxnrep
#> [1] 100
#> 
#> $integer$kspec_bi_cat
#> [1] 0
#> 
#> 
#> $double
#> $double$const_proposal_beta_poi_fix
#> [1] 1
#> 
#> $double$const_proposal_beta_poi
#> [1] 1
#> 
#> $double$const_proposal_beta_bin_fix
#> [1] 1
#> 
#> $double$const_proposal_beta_bin
#> [1] 1
#> 
#> $double$const_proposal_beta_ord_fix
#> [1] 0.5
#> 
#> $double$const_proposal_beta_ord
#> [1] 0.5
#> 
#> $double$const_proposal_beta_cat_fix
#> [1] 1
#> 
#> $double$const_proposal_beta_cat
#> [1] 1
#> 
#> $double$const_proposal_a_ord
#> [1] 0.5
#> 
#> $double$const_proposal_b
#> [1] 0.5
#> 
#> $double$const_proposal_e0
#> [1] 1
#> 
#> $double$df_proposal_beta_poi_fix
#> [1] Inf
#> 
#> $double$df_proposal_beta_poi
#> [1] Inf
#> 
#> $double$df_proposal_beta_bin_fix
#> [1] Inf
#> 
#> $double$df_proposal_beta_bin
#> [1] Inf
#> 
#> $double$df_proposal_beta_ord_fix
#> [1] Inf
#> 
#> $double$df_proposal_beta_ord
#> [1] Inf
#> 
#> $double$df_proposal_beta_cat_fix
#> [1] Inf
#> 
#> $double$df_proposal_beta_cat
#> [1] Inf
#> 
#> $double$df_proposal_a_ord
#> [1] Inf
#> 
#> $double$df_proposal_b
#> [1] Inf
#> 
#> $double$df_proposal_e0
#> [1] Inf
#> 
#> $double$tolerance
#> [1] 1e-07
#> 
#> 
(varying <- default_varying())
#> prec_num    c_ord InvSigma     InvQ      naY 
#>     TRUE     TRUE    FALSE    FALSE    FALSE 
```
