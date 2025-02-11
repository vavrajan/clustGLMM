\name{post_processing}
\alias{post_processing}
\alias{permute_cluster_labels}
\title{
Post-process the clustGLMM sampled chains
}
\description{
Sampling MCMC chains via \code{\link{clustGLMM}} assumes underlying \code{Gmax} clusters and allows for learning the true number of underlying clusters. 
Many clusters are eventually empty and only some cluster labels survive. 
Given random initialization, different chains may end up with different labels of non-empty clusters including their total count. 
The function \code{post_processing} performs the post-processing suggested by Frühwirth-Schnatter (2011) and Malsiner-Walli et al (2016).
}
\usage{
post_processing(mcmc, cols, iter.max = 200, nstart = 30)
permute_cluster_labels(mcmc, perm)
}
\arguments{
  \item{mcmc}{
an output of \code{\link{clustGLMM}} accepted in both formats (list or matrix).
}
  \item{cols}{
a vector of column names of cluster-specific variables to be used for k-means clustering. 
If missing, all cluster-specific parameters are used. 
}
  \item{iter.max}{
the maximum number of iterations of \code{\link{kmeans}} algorithm allowed. 
}
  \item{nstart}{
how many random sets for \code{\link{kmeans}} should be chosen? 
}
  \item{perm}{
a list of permutation of cluster labels, one permutation for each chain.  
}
}
\details{
\code{post_processing} follows the post-processing suggested by Frühwirth-Schnatter (2011) and Malsiner-Walli et al (2016) to identify the number of underlying cluster and to resolve possible issues with label-switching. 
All post-processing is done separately for each chain, keep in mind that the results among chains may differ due to different initial conditions. 
The procedure can be divided into several steps:

0) The procedure works with \code{mcmc$draws} only under \code{howsave = "matrix"}, so if \code{howsave = "list"} then \code{mcmc} is converted to the \code{"matrix"} format, post-processed and converted back to "list".  

1) Take the most frequent \code{Gplus} value - \code{modeGplus} - the number of non-empty clusters across the sampled chain and find the posterior draws with \code{Gplus = modeGplus}. 
The indices are stored in \code{inds} and their length is \code{iters}. 

2) Target all cluster-specific parameters. 
To clarify these are parameters with \code{(g)} in their column name, parameters \code{w}, \code{U}, \code{pUig} do not count. 
If you prefer to use just a subset of those parameters, you have to specify the column names by \code{cols} input variable.

3) Create submatrices of sampled cluster-specific parameters. 
Index by index in \code{inds} gather the cluster-specific parameters from only non-empty clusters. 
Due to label-switching it may happen that different cluster labels correspond to the non-empty clusters at each index. 
After these operations it may happen that values in one column correspond to different clusters.

4) Prepare the data for k-means clustering to resolve the label-switching. 
Gather values for different clusters below each other. 
That is, if there are \code{L} different cluster-specific parameters and \code{iters} samples, then the matrix is rearranged from \code{iters × (L * modeGplus)} to \code{(iters * modeGplus) × L}. 

5) K-means clustering \code{\link{kmeans(., centers=modeGplus, iter.max = iter.max, nstart = nstart)}} is applied to matrix from 5). 
Ideally, rows corresponding to iteration \code{m} would be each clustered to a different mode. 
In such case, we would obtain a permutation for each row that aligns the meaning of the cluster labels across all iterations. 
Unfortunately, it may happen that some iterations do not result into clear permutation. 
In such case, this iteration has to be eliminated (replaced with \code{NA}), so the length of the chain used for inference is reduced from \code{iters} to possibly even lower number \code{Mrhos}. 

6) Regroup the matrix of cluster-specific data back to the format in 3) while also correcting for the permutations found in 5). 
This should result in a matrix where all the values in one column  correspond to the same cluster. 
This output matrix is called \code{clusterspec}.

7) Create the matrix of posterior draws by relabelling and permuting clustering-related parameters.

8) Update all other items saved in \code{mcmc}, see \code{\link{clustGLMM}}.
Especially, 
\code{clustering}, 
\code{certainty}, 
\code{G} (now \code{modeGplus}, possibly different among chains!), 
\code{iterations} (omit iterations with NA values created by not leading to permutation in step 5)), 
\code{settings} (changed dimensions for cluster-specific parameters), 
\code{param_names}, 
\code{last} (list of last sampled values - now from the last iterations with \code{modeGplus} clusters, properly permuted) 
and \code{call} (summary of the post-processing). 
Many items remain the same as in \code{mcmc}, e.g., even the initial values.

\code{permute_cluster_labels} simply takes the list of permutations \code{perm} and reshuffles the cluster labels as requested. 
Useful when you want to align the meaning of cluster labels across different chains. 
For example, an element of the \code{perm} list should be a vector \code{c(2, 3, 4, 1)}, if

cluster 2 is to become cluster 1,

cluster 3 is to become cluster 2,

cluster 4 is to become cluster 3,

cluster 1 is to become cluster 4.
}
\value{
Both functions return back list of class \code{\link{clustGLMM}}, the same structure as \code{mcmc}, but with modified items.
}
\references{
Frühwirth-Schnatter, S. (2011) Dealing with Label Switching under Model Uncertainty, John Wiley & Sons, chap 10, pp 213–239

Malsiner-Walli, G., Frühwirth-Schnatter, S., Grün, B.(2016) Model-based clustering based on sparse finite Gaussian mixtures. Statistics and Computing 26:303–324. https://doi.org/10.1007/s11222-014-9500-2

Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\author{
Jan Vávra
}
\seealso{
\code{\link{clustGLMM}}, \code{\link{plot.clustglmm}}, \code{\link{nice_nrow_ncol}}, \code{\link{post_processing}}, \code{\link{kmeans}}.
}
\examples{
###-----------------------------------------------------------------------------
###     Generate some longitudinal mixed-type data
###-----------------------------------------------------------------------------
corr <- matrix(c(1, -0.5, -0.5, -0.4, -0.4,
                 -0.5, 1, 0.3, 0.4, 0.5,
                 -0.5, 0.3, 1, 0.2, 0.3,
                 -0.4, 0.4, 0.2, 1, 0.1,
                 -0.4, 0.5, 0.3, 0.1, 1), byrow = TRUE, ncol = 5)
sdSigmas <- c(0.5, 0.5, 0.5, 0.5, 0.5)
Sigma <- diag(sdSigmas) %*% corr %*% diag(sdSigmas)

set.seed(31415)
n <- 500
ldata <- generate_longitudinal_mixed_type_data(n = n, n_i = 4,
                     G = 2, timepar = "cross", catpar = "fixed", 
                     catsubjfix = TRUE,
                     xin1 = 1, Sigma = Sigma,
                     Kord = 5, Kcat = 4, kspec_bi_cat = FALSE)
data <- ldata$data

###-----------------------------------------------------------------------------
###     Prepare inputs for clustGLMM
###-----------------------------------------------------------------------------
id <- "i"

family <- c("num", "poi", "bin", "ord", "cat")
names(family) <- Ys <- c("ynum", "ypoi", "ybin", "yord", "ycat")

## Create list of formulae for each outcome
formula <- list()
for(y in Ys){
  formula[[y]] <- list()
  formula[[y]]$fixed <- as.formula("~ f")
  formula[[y]]$group <- as.formula("~ x")
  formula[[y]]$random <- as.formula("~ 1")
  formula[[y]]$offset <- ""
}
formula
# intercept term will be group-specific 
# since it appears in $group which has a preference over $fixed

## All other inputs need to be prepared,
## unless the default values are desired (varying, save, howsave, param, tuning).
## See the tutorial for how these can be adjusted.


###-----------------------------------------------------------------------------
###     Sampling chains by clustGLMM
###-----------------------------------------------------------------------------

Gmax <- 5     # maximal number of mixture components to be considered
nchains <- 2  # number of chains to be sampled 

## No inits specified
mcmc <- clustGLMM(formula = formula, id = id, family = family, data = data, 
                  G = Gmax, iter = 1000, nchains = nchains)

# for the sake of an example we have too short iter
# use higher iter for proper inference

###-----------------------------------------------------------------------------
###     Post-processing the sampled chains
###-----------------------------------------------------------------------------

### First run the function post_processing
# takes some time with iter large....
mcmc <- post_processing(mcmc)
mcmc
summary(mcmc)

# All plotting functions should work even with these post-processed draws 
plot(mcmc, "beta_num", "clusters_kerneldensity", chains=1:nchains)

### Additionally permute so that the meaning of the number of clusters
### is kept the same across all chains

## first we have to recognize which chain has which meaning
neff <- mcmc$ngrp[["ynum"]]
opar <- par(no.readonly = TRUE)
par(mfrow = c(nchains,neff), mar = c(4,4,1,1))
for(ch in 1:nchains){
  TAB <- table(mcmc$last[[ch]]$U)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)]) 
  
  for(i in 1:neff){
    plot_clusters(mcmc = mcmc, what = "beta_num", 
                  yspec = "ynum", dimspec = i, whichg = nonempty,
                  setparmfrow = FALSE, chains = ch,
                  doECDF = FALSE)
    # true values
    for(gg in 1:2){
      abline(v = ldata$beta_num[[gg]][i], col = "black", lty = 2)
    }
  }
}
par(opar)

### Example how to change the permutation so that results across chains align
perm <- list()
perm[[1]] <- 1:mcmc$G[[1]] # change if needed
perm[[2]] <- rev(1:mcmc$G[[2]]) # change if needed

# function to permute meaning of cluster labels in mcmc according to perm
mcmc <- permute_cluster_labels(mcmc, perm)



###-----------------------------------------------------------------------------
###     Plotting post-processed chains
###-----------------------------------------------------------------------------
# Notice that now there are only modeGplus clusters and not original G

plot(mcmc, "w", "traceplots", thin = 10)
plot(mcmc, "w", "kerneldensity")
plot(mcmc, "beta_num", "clusters_kerneldensity")
plot(mcmc, "beta_poi", "clusters_kerneldensity")
plot(mcmc, "beta_bin", "clusters_kerneldensity")
plot(mcmc, "beta_ord", "clusters_kerneldensity")
plot(mcmc, "beta_cat", "clusters_kerneldensity")

}
