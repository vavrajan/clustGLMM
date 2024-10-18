\name{post_processing}
\alias{post_processing}
\alias{permute_cluster_labels}
\title{
Post-process the clustGLMM sampled chains
}
\description{
Sampling MCMC chains via \code{\link{clustGLMM}} assumes underlying \code{Gmax} clusters and allows for learning the true number of underlying clusters. Many clusters are eventually empty and only some cluster labels survive. Given random initialization, different chains may end up with different labels of non-empty clusters including their total count. The function \code{post_processing} performs the post-processing suggested by Frühwirth-Schnatter (2011) and Malsiner-Walli et al (2016).
}
\usage{
post_processing(mcmc, cols)
permute_cluster_labels(mcmc, perm)
}
\arguments{
  \item{mcmc}{
an output of \code{\link{clustGLMM}} accepted in both formats (list or matrix).
}
  \item{cols}{
a vector of column names of cluster-specific variables to be used for k-means clustering. If missing, all cluster-specific parameters are used. 
}
  \item{perm}{
a list of permutation of cluster labels, one permutation for each chain.  
}
}
\details{
\code{post_processing} follows the post-processing suggested by Frühwirth-Schnatter (2011) and Malsiner-Walli et al (2016) to identify the number of underlying cluster and to resolve possible issues with label-switching. All post-processing is done separately for each chain, keep in mind that the results among chains may differ due to different initial conditions. 
The procedure can be divided into several steps:

0) The procedure works with \code{mcmc$draws} only under \code{howsave = "matrix"}, so if \code{howsave = "list"} then \code{mcmc} is converted to the \code{"matrix"} format.  

1) Find the most frequent \code{Gplus} value - the number of non-empty clusters across the sampled chain and call it \code{modeGplus}.

2) Find the posterior draws with \code{Gplus = modeGplus}. The indices are stored in \code{inds} and their length is \code{iters}. 

3) Target all cluster-specific parameters. To clarify these are parameters with \code{(g)} in their column name, parameters \code{w}, \code{U}, \code{pUig} do not count. If you prefer to use just a subset of those parameters, you have to specify the column names by \code{cols} input variable.

4) Create submatrices of sampled cluster-specific parameters. Index by index in \code{inds} gather the cluster-specific parameters from only non-empty clusters. Due to label-switching it may happen that different cluster labels correspond to the non-empty clusters at each index. After these operations it may happen that values in one column correspond to different clusters.

5) Prepare the data for k-means clustering to resolve the label-switching. Gather values for different clusters below each other. That is, if there are \code{L} different cluster-specific parameters and \code{iters} samples, then the matrix is rearranged from \code{iters × (L * modeGplus)} to \code{(iters * modeGplus) × L}. 

6) K-means clustering \code{\link{kmeans}} is applied to matrix from 5). Ideally, rows corresponding to iteration \code{m} would be each clustered to a different mode. In such case we would obtain a permutation for each row that aligns the meaning of the cluster labels across all iterations. Unfortunately, it may happen that some iterations do not result into clear permutation. In such case, this iteration has to be eliminated (replaced with \code{NA}), so the length of the chain used for inference is reduced from \code{iters} to possibly even lower number \code{Mrhos}. 

7) Regroup the matrix of cluster-specific data back to the format in 4) while also correcting for the permutations found in 6). This should result in a matrix where all the values in one column  correspond to the same cluster. This output matrix is called \code{clusterspec}.

8) Relabel and reshuffled according to the same permutations also the values of other cluster-related parameters: \code{w}, \code{ng}, \code{U}, \code{pUig} and save them as \code{wGplus}, \code{ngGplus}, \code{UGplus}, \code{pUigGplus}.

These values are added to \code{mcmc}. If you wish to create plots with these values using the predefined plotting functions \code{\link{plot_all}} use \code{restrictedtoGplus = TRUE}.

\code{permute_cluster_labels} simply takes the list of permutations \code{perm} and reshuffles the cluster labels as requested. Useful when you want to align the meaning of cluster labels across different chains. For example, an element of the \code{perm} list should be a vector \code{c(2, 3, 4, 1)}, if

cluster 2 is to become cluster 1,

cluster 3 is to become cluster 2,

cluster 4 is to become cluster 3,

cluster 1 is to become cluster 4.
}
\value{
Both functions return back list of the same structure as \code{mcmc}. 

\code{post_processing} return \code{mcmc} without any changes, but includes additional elements: 
\item{modeGplus}{a list of the most frequent \code{Gplus} value for each chain.}
\item{iters}{a list of the number of samples with \code{Gplus} equal to the \code{modeGplus} for each chain, basically new \code{iter} value for each of the chains.}
\item{inds}{a list of indices of samples with \code{Gplus} equal to the \code{modeGplus} for each chain. Its length is \code{iters}.}
\item{Mrhos}{a list of number of samples out of \code{iters} that has been succesfully permuted to resolve label-switching.}
\item{clusterspec}{a list of matrices of cluster-specific parameters  where label-switching has been resolved. Cluster labels are now \code{1:modeGplus}.}

The following are cluster-relevant model parameters that do not fall into the set of cluster-specific parameters, they are non-empty only if also present in \code{mcmc}:
\item{wGplus}{a list of subsetted (by \code{inds}) and permuted samples of the prior cluster probabilities \code{w} for each chain.}
\item{ngGplus}{a list of subsetted (by \code{inds}) and permuted samples of the cluster occupancy numbers \code{ng} for each chain.}
\item{UGplus}{a list of subsetted (by \code{inds}) and permuted samples of the cluster allocation indicators \code{U} for each chain.}
\item{pUigGplus}{a list of subsetted (by \code{inds}) and permuted samples of the full-conditional clustering probabilities \code{pUig} for each chain.}

\code{permute_cluster_labels} returns \code{mcmc} without any additional elements, but only alters already existing elements. In particular, permutes all cluster-specific parameters as required. 
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
\code{\link{clustGLMM}}, \code{\link{plot_all}}, \code{\link{nice_nrow_ncol}}, \code{\link{post_processing}}.
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

# the most frequent Gplus value in each chain
mcmc$modeGplus
# number of states with modeGplus value
mcmc$iters
# number of states (out of total iters) not resulting in a permutation
mcmc$Mrhos
# (potentially reshuffled) cluster-specific parameters
# also reduced notation to 1:modeGplus
# different chains may have still permuted meaning compared to other chains
dim(mcmc$clusterspec[[1]])
head(mcmc$clusterspec[[1]])
# NA rows are those where the k-means clustering failed to find a permuation

# indexes of states with current Gplus=modeGplus
mcmc$inds[[1]]

# cluster related parameters permuted to fit clusterspec matrix
head(mcmc$wGplus[[1]])
head(mcmc$ngGplus[[1]])
head(mcmc$UGplus[[1]][,1:20])

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
perm[[1]] <- 1:mcmc$modeGplus[[1]] # change if needed
perm[[2]] <- rev(1:mcmc$modeGplus[[2]]) # change if needed

# function to permute meaning of cluster labels in mcmc according to perm
mcmc <- permute_cluster_labels(mcmc, perm)



###-----------------------------------------------------------------------------
###     Plotting post-processed chains
###-----------------------------------------------------------------------------

# for cluster-specific parameters use restrictedtoGplus = T
# this will use mcmc$clusterspec or mcmc$...Gplus values instead of those
# saved in mcmc$draws

thin <- 10

## probabilites w
for(ch in 1:mcmc$nchains){
  par(mfrow = c(mcmc$modeGplus[[ch]],4))
  for(g in 1:mcmc$modeGplus[[ch]]){
    plot_traceplots(mcmc = mcmc, what = "w", 
                    dimspec = g, 
                    whatLAB = paste0("w(",g,")"),
                    thin = thin, chains = ch, labcex = 0.6, 
                    restrictedtoGplus = TRUE)
    plot_kerneldensity(mcmc = mcmc, what = "w",
                       dimspec = g,  
                       whatLAB = paste0("w(",g,")"),
                       thin = thin, chains = ch, labcex = 0.6, 
                       restrictedtoGplus = TRUE)
    plot_ECDF(mcmc = mcmc, what = "w", 
              dimspec = g, 
              whatLAB = paste0("w(",g,")"),
              thin = thin, chains = ch, labcex = 0.6, 
              restrictedtoGplus = TRUE)
    plot_ACF(mcmc = mcmc, what = "w", 
             dimspec = g, 
             whatLAB = paste0("w(",g,")"),
             thin = thin, lag.max = 30, chains = ch, labcex = 0.6, 
             restrictedtoGplus = TRUE)
  }
}

## beta_num
nfr <- length(mcmc$lgrpnames$ynum)
for(ch in 1:mcmc$nchains){
  par(mfrow = c(mcmc$modeGplus[[ch]],nfr))
  for(g in 1:mcmc$modeGplus[[ch]]){
    # traceplots
    for(i in 1:nfr){
      plot_traceplots(mcmc = mcmc, what = "beta_num", 
                      gspec = g, yspec = "ynum", dimspec = i,
                      thin = thin, chains = ch, labcex = 0.6, 
                      restrictedtoGplus = TRUE)
    }
  }
}

par(opar)
## other similarly...
## plot_..._param not (yet) implemented for post-processed chains.
}
