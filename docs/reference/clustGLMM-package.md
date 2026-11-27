# Model-Based Clustering of Mixed-Type Longitudinal Data

Provides tools for Bayesian estimation and inference for modelling
clusterwise multivariate regression models for numeric, count, binary,
ordinal and count outcomes observed repeatedly on the same units and
where possible relations among outcomes are captured through a joint
distribution of random effects. The clusters are defined through
cluster-specific parameters, which the analyst can choose, e.g., with
respect to the regression coefficients. In particular, the model
specification for each regression model via the formula is specific to
the outcome and consists of four parts: (1) fixed - regression
coefficients common to all clusters, (2) group - group-specific
regression coefficients, (3) random - random effects specific for each
unit, (3) offset - name of an offset variable (if needed). Estimation is
performed using MCMC sampling combining Gibbs and Metropolis-Hastings
steps. Post-processing tools allow to assess convergence and address
label switching and provide visual diagnostics. Units may be classified
based on sampled allocation indicators or by exploiting the posterior
distribution of the classification probabilities. For more details see
Vavra et al. (2024) \<doi:10.1007/s11222-023-10304-5\>.

## Details

The key function is
[`clustGLMM`](https://vavrajan.github.io/clustGLMM/reference/clustGLMM.md)
which creates MCMC samples. The data, formulae and other auxiliary
parameters need to be set up carefully before calling this function. The
examples and tutorials may be used for guidance. Other functions
contained in the package consist of:

- a function for generating longitudinal data of mixed-type:
  [`generate_longitudinal_mixed_type_data`](https://vavrajan.github.io/clustGLMM/reference/generate_longitudinal_mixed_type_data.md),

- functions to generate plots:
  [`plot_traceplots`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_traceplots_param`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_kerneldensity`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_kerneldensity_param`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_ACF`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_ACF_param`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_ECDF`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_ECDF_param`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_clusters`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_clusters_ECDF_param`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot_clusters_kerneldensity_param`](https://vavrajan.github.io/clustGLMM/reference/plot_diagnostics.md),
  [`plot.clustglmm`](https://vavrajan.github.io/clustGLMM/reference/plot.clustglmm.md),
  [`plot_num_vs_x_grouped`](https://vavrajan.github.io/clustGLMM/reference/plot_cat_vs_x_grouped.md),
  [`plot_cat_vs_x_grouped`](https://vavrajan.github.io/clustGLMM/reference/plot_cat_vs_x_grouped.md),

- functions for post-processing the samples:
  [`post_processing`](https://vavrajan.github.io/clustGLMM/reference/post_processing.md),
  [`permute_cluster_labels`](https://vavrajan.github.io/clustGLMM/reference/post_processing.md),

- functions for internal computations.

## Author

Jan Vavra \[aut, cre\] (ORCID:
\<https://orcid.org/0000-0003-0068-5356\>), Bettina Gruen \[aut\]
(ORCID: \<https://orcid.org/0000-0001-7265-4773\>), Gertraud
Malsiner-Walli \[aut\] (ORCID:
\<https://orcid.org/0000-0002-1213-4749\>), Arnost Komarek \[aut\]
(ORCID: \<https://orcid.org/0000-0001-8778-3762\>)

Maintainer: Jan Vavra \<vavraj@karlin.mff.cuni.cz\>

## References

Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise
multivariate regression of mixed-type panel data. Stat Comput 34, 46
(2024). https://doi.org/10.1007/s11222-023-10304-5
