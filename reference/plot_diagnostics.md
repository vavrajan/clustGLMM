# Plot MCMC Samples from `clustGLMM`

There are 4 main plots to be created:

- `ACF`:

  auto-correlation function,

- `ECDF`:

  empirical cumulative distribution function,

- `kerneldensity`:

  one-dimensional kernel density estimate,

- `traceplots`:

  parameter samples against the iteration number.

All of these for one specific parameter can be plotted by
`plot_diagnostics`. These functions above work only for a univariate
parameter. Additionally, `acceptance` draws cumulative means of
acceptance of Metropolis-Hastings proposals against the iteration number
for a parameter block. Plots for all elements of a multivariate
parameter can be created using `..._param` functions. A special function
is `plot_clusters` which allows to compare densities or ECDFs of
univariate group-specific parameter across non-empty clusters.

## Usage

``` r
plot_acceptance(mcmc, what, 
                burnin = 0, thin = 1, iterations, chains,
                dolegend = FALSE, legend.pos = "bottomright",
                COL, labcex = 1, whatLAB,
                adjust.ylim = 0)
plot_ACF(mcmc, scalar_samples, what, gspec, dimspec, yspec,
         burnin = 0, thin = 1, lag.max = 30, iterations, chains, 
         dolegend = FALSE, legend.pos = "topright",
         COL, move.width = 0.3, labcex = 1, whatLAB)
plot_ECDF(mcmc, scalar_samples, what, gspec, dimspec, yspec,
          burnin = 0, thin = 1, iterations, chains, 
          dolegend = FALSE, legend.pos = "bottomright",
          COL, labcex = 1, whatLAB)
plot_kerneldensity(mcmc, scalar_samples, what, gspec, dimspec, yspec, 
                   burnin = 0, thin = 1, iterations, chains, 
                   dolegend = FALSE, legend.pos = "topright", 
                   COL, labcex = 1, whatLAB)
plot_traceplots(mcmc, scalar_samples, what, gspec, dimspec, yspec,
                burnin = 0, thin = 1, iterations, chains, 
                dolegend = FALSE, legend.pos = "top",
                COL, labcex = 1, whatLAB)
plot_ng_trace_chain_split(mcmc, 
                          burnin = 0, thin = 1, iterations, chains,
                          COL, labcex = 1, setparmfrow = TRUE)
plot_diagnostics(mcmc, scalar_samples, what, gspec, dimspec, yspec,
                 burnin = 0, thin = 1, lag.max = 30, iterations, chains,
                 trueval, setparmfrow = TRUE, COL, move.width = 0.3, whatLAB)
plot_acceptance_param(mcmc, what = "e0", burnin = 0, thin = 1, iterations, chains, 
                      setparmfrow = TRUE, allclusters = FALSE, 
                      adjust.ylim = 0, which.b = 1:min(mcmc$n, 20))
plot_ACF_param(mcmc, what = "ng", yspec, burnin = 0, thin = 1, iterations,
               lag.max = 30, chains, setparmfrow = TRUE, allclusters = FALSE)
plot_ECDF_param(mcmc, what, yspec, burnin = 0, thin = 1, iterations, chains, 
                setparmfrow = TRUE, allclusters = FALSE)
plot_kerneldensity_param(mcmc, what, yspec, burnin = 0, thin = 1, iterations,
                         chains, setparmfrow = TRUE, allclusters = FALSE)
plot_traceplots_param(mcmc, what, yspec, burnin = 0, thin = 1, iterations,
                      chains, setparmfrow = TRUE, allclusters = FALSE)
plot_clusters(mcmc, what, dimspec, yspec, whichg,
              burnin = 0, thin = 1, iterations, chains,
              setparmfrow = TRUE, doKern = TRUE, doECDF = TRUE,
              COL, labcex = 0.8, whatLAB, wherelegend = "topright")
plot_clusters_kerneldensity_param(mcmc, what, whichg, yspec, burnin = 0,
                                  thin = 1, iterations, chains,
                                  setparmfrow = TRUE)
plot_clusters_ECDF_param(mcmc, what, whichg, yspec, burnin = 0, thin = 1,
                         iterations, chains, setparmfrow = TRUE)
```

## Arguments

- mcmc:

  an output of
  [`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md)
  accepted in both formats (list or matrix).

- scalar_samples:

  an output of `get_scalar_samples`, overrides any given `burnin` or
  `thin`.

- what:

  a name of the parameter to be extracted; see `names(mcmc$save)` for
  possible options. For `plot_acceptance` the full name of the parameter
  block has to be given, e.g. `"beta_poi_...(1)"`.

- gspec:

  a group specifier (values in `1:mcmc$G`); if not group-specific
  parameter this *should be missing*.

- dimspec:

  a dimension specifier; not specifying this is possible only when the
  parameter is scalar.

- yspec:

  an outcome specifier; if the parameter is not outcome-dependent then
  this *should be missing*.

- burnin:

  the length of burn-in period (default 0), i.e. number of initial
  samples to ignore for inference.

- thin:

  thinning rate (default 1); only every `thin`th value is used for
  inference (must be a positive integer).

- iterations:

  a vector of iterations numbers to extract, overrides any given
  `burnin` or `thin`.

- chains:

  a vector of chain numbers to be plotted.

- adjust.ylim:

  a probability that will shorten the `ylim` to
  `quantile(y, probs = c(adjust.ylim, 1-adjust.ylim))` where `y` stands
  for the y-axis values of the objects to be plotted.

- lag.max:

  maximal lag value to be plotted in ACF.

- trueval:

  a true value of the parameter to be plotted by vertical/horizontal
  line.

- setparmfrow:

  whether `par(mfrow = ..)` should be set automatically before creating
  all plots.

- allclusters:

  whether for group-specific parameters plot all the clusters (`TRUE`)
  or just those in `mcmc$clusters` (`FALSE`).

- which.b:

  a subset of `1:mcmc$n` giving the subject internal indices to be
  plotted; by default setting limited to first 20 subjects, since
  plotting all of them could produce too many plots.

- COL:

  colors for different chains; default values are set by
  `diverge_hcl(mcmc$nchains, c = 80, l = 70)`.

- move.width:

  ACF for different chains are jittered around the actual lag by
  `seq(-move.width, move.width, length.out = length(chains)+1)`

- whatLAB:

  the label for the parameter; default corresponds to the label for that
  parameter that would be used in *matrix* format of the output.

- labcex:

  `cex` value for label of plots.

- dolegend:

  a logical declaring whether legend should be added to the plot.

- legend.pos:

  position for the legend given by a string (see details of `legend`).

- whichg:

  a set of cluster labels to be compared in `plot_clusters`. By default
  `mcmc$clusters` is taken under consistency among selected chains. Can
  be a list in `plot_clusters_kerneldensity_param` and
  `plot_clusters_ECDF_param`.

- doKern:

  a logical declaring whether kernel density estimates should be plotted
  for comparing different clusters.

- doECDF:

  a logical declaring whether ECDFs should be plotted for comparing
  different clusters.

- wherelegend:

  a string that declares where the legend should be placed.

## Details

`scalar_samples` is a list (for each element of `chains`) containing
`data.frame` of three columns: `chain` declaring the chain number, `m`
declaring the iterations numbers and `y` containing chosen samples. This
could be created using
[`get_scalar_samples`](https://vavrajan.github.io/clustGLMM/reference/get_scalar_samples.md),
see its help file for more details. If missing,
[`get_scalar_samples`](https://vavrajan.github.io/clustGLMM/reference/get_scalar_samples.md)
is called anyway with the given specifications. If `iterations` is
missing, first, only iterations with `m > burnin` are selected, then,
from the remaining values only `seq(1, total_remaining_count, by=thin)`
are selected. Parameter is selected by filling/missing `what`, `gspec`,
`yspec`, `dimspec`. More details here:
[`get_scalar_samples`](https://vavrajan.github.io/clustGLMM/reference/get_scalar_samples.md).

Description of the plots:

- `ACF`:

  Autocorrelation function plots have the lag on x-axis \[0, `lag.max`\]
  and auto-correlations in interval (-1, 1) on y-axis. By default,
  `lag.max = min(30, mcmc$iter-1)`. Different chains have different
  colors and are slightly jittered for better visibility, see
  `move.width`.

- `ECDF`:

  Cumulative distribution function plots have the range of sampled
  values on x-axis and the cumulative probabilities in \[0, 1\] on
  y-axis. Different colors depict different chains of non-decreasing
  ECDFs.

- `kerneldensity`:

  Kernel density estimator plots have the range of sampled values on
  x-axis and the density value in \[0, max(density)\] on y-axis.
  Different colors depict different chains.

- `traceplot`:

  Traceplots have the range of the considered iteration numbers on
  x-axis and the sampled values of the parameter on y-axis. Different
  colors depict different chains. Stationarity of these plots is desired
  for reliable posterior approximation, if the plots are not stationary
  then consider continuation in sampling.

You can plot all 4 of these by `plot_diagnostics`. By default they are
produced in 2 by 2 layout, but if you would prefer some different
arrangement you can do that before calling the function and suppressing
`setparmfrow = FALSE`. Moreover, the `ACF` plot contains a legend for
the chains including the acceptance percentage, if the parameter is
sampled via Metropolis-Hastings proposal algorithm.

Using the functions above you can only plot a scalar model parameter
which is specified by filling/missing `what`, `gspec`, `yspec`,
`dimspec`. If you would like to plot results for all parameters that are
associated with parameter `what` then `..._param` functions can do this
for you automatically. By default, only cluster-specific parameters for
the (non-empty) clusters saved in `mcmc$clusters` are arranged. If the
labels across chains do not coincide, then chains are plotted
separately. The exception is made for `w` and `ng` for which always all
clusters are plotted. Moreover, there is special function
`plot_ng_trace_chain_split` for plotting traces of `ng` for all clusters
into one plot per chain (default by
[`plot.clustglmm`](https://vavrajan.github.io/clustGLMM/reference/plot.clustglmm.md)).
Using function
[`nice_nrow_ncol`](https://vavrajan.github.io/clustGLMM/reference/nice_nrow_ncol.md)
plots are arranged into a most effective grid. Moreover, parameters
which are symmetric matrices (with or without diagonal) are arranged as
expected and only the upper triangle is filled with plots. See examples
below.

Parameters `e0`, `b`, `beta_poi_fix`, `beta_poi`, `beta_bin_fix`,
`beta_bin`, `beta_ord_fix`, `beta_ord`, `a_ord`, `beta_cat_fix`,
`beta_cat` are sampled in blocks. Therefore, their acceptance is
gathered in `mcmc$accept` named after such blocks. In many cases the
corresponding name is much shorter, e.g., `b[i]` instead of `b[i,3]`,
`beta_poi_y(1)` instead of `beta_poi_y(1)[2]` or even `beta_cat_y(2)`
instead of `beta_poi_y(2)[2,3]`. Function `plot_acceptance` requires
these shorter names in the argument `what` (the name that appears in
columns of `mcmc$accept[[ch]]`), `yspec`, `dimspec`, `gspec` are no
longer required. It plots the cumulative mean of acceptance level
against the iteration number on percentage scale. The y-axis is by
default fixed to 0-100 percent, but could be shortened via `adjust.ylim`
probability. Since it depends on the plotted y-values through
`quantile(yvalues, probs = c(adjust.ylim, 1-adjust.ylim))` the final
`ylim` will be different among parameter blocks. To plot all cumulative
acceptance levels tied to the same parameter use `plot_acceptance_param`
or `plot(, which = "acceptance")`. Since `c_ord`, `a_ord` and `pi_ord`
are different transformations of the same parameters, the same
acceptance (saved only for `a_ord`) is plotted in these cases.

Special function in this collection is `plot_clusters` which allows to
plot kernel density estimates or empirical cumulative distribution
function of the same parameters for different cluster labels (`whichg`)
into one plot. There are also `plot_clusters_kerneldensity_param` and
`plot_clusters_ECDF_param` functions that create the plots automatically
for a set of given parameter `what`, however, the plotted cluster labels
only include those in `mcmc$clusters` of corresponding chain.

## Value

The functions do not return anything but produce one or more plots.

## References

Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise
multivariate regression of mixed-type panel data. Stat Comput 34, 46
(2024). https://doi.org/10.1007/s11222-023-10304-5

## Author

Jan Vávra

## See also

[`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md),
[`get_scalar_samples`](https://vavrajan.github.io/clustGLMM/reference/get_scalar_samples.md),
[`nice_nrow_ncol`](https://vavrajan.github.io/clustGLMM/reference/nice_nrow_ncol.md),
[`post_processing`](https://vavrajan.github.io/clustGLMM/reference/post_processing.md).

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
#> <environment: 0x623c2b161860>
#> 
#> $ynum$group
#> ~x
#> <environment: 0x623c2b161860>
#> 
#> $ynum$random
#> ~1
#> <environment: 0x623c2b161860>
#> 
#> $ynum$offset
#> [1] ""
#> 
#> 
#> $ypoi
#> $ypoi$fixed
#> ~f
#> <environment: 0x623c2b161860>
#> 
#> $ypoi$group
#> ~x
#> <environment: 0x623c2b161860>
#> 
#> $ypoi$random
#> ~1
#> <environment: 0x623c2b161860>
#> 
#> $ypoi$offset
#> [1] ""
#> 
#> 
#> $ybin
#> $ybin$fixed
#> ~f
#> <environment: 0x623c2b161860>
#> 
#> $ybin$group
#> ~x
#> <environment: 0x623c2b161860>
#> 
#> $ybin$random
#> ~1
#> <environment: 0x623c2b161860>
#> 
#> $ybin$offset
#> [1] ""
#> 
#> 
#> $yord
#> $yord$fixed
#> ~f
#> <environment: 0x623c2b161860>
#> 
#> $yord$group
#> ~x
#> <environment: 0x623c2b161860>
#> 
#> $yord$random
#> ~1
#> <environment: 0x623c2b161860>
#> 
#> $yord$offset
#> [1] ""
#> 
#> 
#> $ycat
#> $ycat$fixed
#> ~f
#> <environment: 0x623c2b161860>
#> 
#> $ycat$group
#> ~x
#> <environment: 0x623c2b161860>
#> 
#> $ycat$random
#> ~1
#> <environment: 0x623c2b161860>
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
plot(mcmc)

plot(mcmc, ng_trace_chain_split=FALSE)

                       
###-----------------------------------------------------------------------------
###     Creating plots
###-----------------------------------------------------------------------------

## Traceplots of number of non-empty clusters
# Plotting with defaults values:
# thin = 1, burnin = 0, chains = 1:mcmc$nchains
plot_traceplots(mcmc = mcmc, what = "Gplus")


## Traceplots of number of observations in cluster 1
plot_traceplots(mcmc = mcmc, what = "ng", gspec = 1)


# Omitting the burn-in-period
plot_traceplots(mcmc = mcmc, what = "ng", gspec = 1, burnin = 20)


# Thinned plots
plot_traceplots(mcmc = mcmc, what = "ng", gspec = 1, burnin = 20, thin = 5)


## First fixed beta coefficient for the first numerical outcome
plot_ACF(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1], dimspec = 1)

plot_ECDF(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1], dimspec = 1)

plot_kerneldensity(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1],
                   dimspec = 1, dolegend = TRUE)

plot_traceplots(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1],
                dimspec = 1)

plot_diagnostics(mcmc = mcmc, what = "beta_num_fix", yspec = mcmc$Nums[1], dimspec = 1)


## Acceptance level for beta coefficients for ordinal outcome in specific cluster
plot_acceptance(mcmc = mcmc, what = "beta_ord_yord(1)")


## Acceptance level for beta coefficients for all ordinal outcomes and clusters
plot(mcmc, "beta_ord", which = "acceptance")



## Correlation matrix between the random effects
# True correlation matrix used for generating the data
corr <- matrix(c(1, -0.5, -0.5, -0.4, -0.4,
                 -0.5, 1, 0.3, 0.4, 0.5,
                 -0.5, 0.3, 1, 0.2, 0.3,
                 -0.4, 0.4, 0.2, 1, 0.1,
                 -0.4, 0.5, 0.3, 0.1, 1), byrow = TRUE, ncol = 5)
                 
par(mfrow = c(mcmc$totnran-1, mcmc$totnran-1), mar = c(4, 4, 1, 1))
for(i in 1:(mcmc$totnran-1)){
  for(j in 2:(mcmc$totnran)){
    if(i < j){
      plot_traceplots(mcmc = mcmc, what = "corSigma", 
                      dimspec = c(i, j),
                      labcex = 0.6)
      abline(h = corr[i,j], col = "black", lty = 2)
    }else{
      plot(x = c(0, 1), y = c(0, 1), type = "n", xlab = "", ylab = "",
           xaxt = "n", yaxt = "n", bty = "n")
    }
  }
}

# or simply use
plot_traceplots_param(mcmc, "corSigma")

# but no true value displayed...

## Possible what arguments:
# "_param" - plots everything tied to the given what parameter
plot_traceplots_param(mcmc, "Gplus", thin = 5)

plot_traceplots_param(mcmc, "e0", thin = 5)

plot_traceplots_param(mcmc, "ng", thin = 5)

plot_traceplots_param(mcmc, "w", thin = 5)

plot_traceplots_param(mcmc, "sdSigma", thin = 5)

plot_traceplots_param(mcmc, "corSigma", thin = 5)

plot_traceplots_param(mcmc, "Sigma", thin = 5)

plot_traceplots_param(mcmc, "InvSigma", thin = 5)

plot_traceplots_param(mcmc, "prec_num", thin = 5)


plot_traceplots_param(mcmc, "beta_num_fix", thin = 5)

plot_traceplots_param(mcmc, "beta_num", thin = 5)


plot_traceplots_param(mcmc, "beta_poi_fix", thin = 5)

plot_traceplots_param(mcmc, "beta_poi", thin = 5)


plot_traceplots_param(mcmc, "beta_bin_fix", thin = 5)

plot_traceplots_param(mcmc, "beta_bin", thin = 5)


plot_traceplots_param(mcmc, "beta_ord_fix", thin = 5)

plot_traceplots_param(mcmc, "beta_ord", thin = 5)


plot_traceplots_param(mcmc, "c_ord", thin = 5)


plot_traceplots_param(mcmc, "beta_cat_fix", thin = 5)

plot_traceplots_param(mcmc, "beta_cat", thin = 5)





plot_traceplots_param(mcmc, "loglik", thin = 5)


## Other "_param" versions for different types of plots
plot_kerneldensity_param(mcmc, "beta_num")


plot_ECDF_param(mcmc, "prec_num")


plot_ACF_param(mcmc, "w")

plot_acceptance_param(mcmc, "c_ord")



## Comparing different cluster-specific parameters
## We recommend to do it for chains separately
ch <- 1
TAB <- table(mcmc$last[[ch]]$U)
nonempty <- as.numeric(names(TAB)[which(TAB>0)]) 
# non-empty clusters decided based on last sampled values
# by default, mcmc$clusters (based on all iterations) is taken

neff <- mcmc$ngrp[[mcmc$Nums[1]]]
oldpar <- par(mfrow = c(length(mcmc$Nums), neff), mar = c(4, 4, 1, 1))

for(y in mcmc$Nums){
  for(i in 1:neff){
    plot_clusters(mcmc = mcmc, what = "beta_num", 
                  yspec = y, dimspec = i, whichg = nonempty,
                  setparmfrow = FALSE, chains = ch,
                  doECDF = FALSE)
    # Now you can add anything to the plot
    # e.g. true values, etc.
  }
}

par(oldpar)

# or simply use
plot_clusters_kerneldensity_param(mcmc = mcmc, what = "beta_num", chains = c(ch))

# where cluster labels are taken from mcmc$clusters
```
