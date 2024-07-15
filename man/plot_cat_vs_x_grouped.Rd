\name{plot_cat_vs_x_grouped}
\alias{plot_cat_vs_x_grouped}
\alias{plot_num_vs_x_grouped}
\title{
Plotting clustered longitudinal (panel) data 
}
\description{
Functions that work with longitudinal data and plot the outcomes with respect to some grouping, e.g. clustering found by \code{\link{clusGLMM}}. 
}
\usage{
plot_num_vs_x_grouped(data, y = "y", x = "x", group = "g", id = "i", units, 
                      xbreaks, add_lowess = TRUE, legend_placement = "topright")
plot_cat_vs_x_grouped(data, y = "y", x = "x", group = "g", 
                      xbreaks, margin = 0.05, reverse = FALSE)
}
\arguments{
  \item{data}{
a \code{data.frame} containing \code{y}, \code{x}, \code{group}, \code{id}.
}
  \item{y}{
a name of the outcome to be plotted.
}
  \item{x}{
a name of numeric regressor. 
}
  \item{group}{
a name of the clustering variable that divides different ids into groups.
}
  \item{id}{
a name of the id variable.
}
  \item{units}{
a vector of unit ids to be plotted, if missing all ids are plotted.
}
  \item{xbreaks}{
a set of breaks along which to cut regressor \code{x}, if missing \code{x} is divided equidistantly into 10 categories.
}
  \item{add_lowess}{
a boolean whether to add lowess curves for each group label.
}
  \item{margin}{
how far from each other should different groups be displayed.
}
  \item{reverse}{
a boolean whether to reverse the order of ordinal categories; if \code{TRUE} 0 on top, if \code{FALSE} 0 at the bottom.
}
  \item{legend_placement}{
a string describing the position of legend, see \code{\link{legend}}.
}
}
\details{
\code{plot_num_vs_x_grouped} spaghetti-plots the outcome \code{y} wrt covariate \code{x} coloured by \code{group}, grey for level 0 (unclassified). Only units included in \code{units} are plotted. 

\code{plot_cat_vs_x_grouped} divides the data by covariate \code{x} into groups divided by \code{xbreaks}. Then, for each group label in \code{group} (including level 0 = unclassified) a spineplot capturing the evolution of proportions of a categorical outcome \code{y} wrt the covariate \code{x}. 
}
\value{
The functions do not return anything but produce a plot.
}
\references{
Vávra, J., Komárek, A., Grün, B., Malsiner-Walli, G. Clusterwise multivariate regression of mixed-type panel data. Stat Comput 34, 46 (2024). https://doi.org/10.1007/s11222-023-10304-5
}
\author{
Jan Vávra
}
\seealso{
\code{\link{clusGLMM}}
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
###     Prepare inputs for clusGLMM
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
## unless the default values are desired (spec, whatsave, howsave, param, tuning).
## See the tutorial for how these can be adjusted.

###-----------------------------------------------------------------------------
###     Final call of clusGLMM
###-----------------------------------------------------------------------------

Gmax <- 5     # maximal number of mixture components to be considered
chains <- 2   # number of chains to be sampled 

## No inits specified
mcmc <- clusGLMM(formula = formula, id = id, family = family, data = data, 
                 G = Gmax, iter = 1500, chains = chains)
                       
###-----------------------------------------------------------------------------
###     Clustering based on not post-processed data (using sampled indicators)
###-----------------------------------------------------------------------------
### Take the sampled U indicators
threshold <- 0.5 # if not exceeded, unit remains unclassified

# posterior mode of Gplus=number of nonemtpy components
modeGplus <- as.numeric(tapply(mcmc$all$Gplus, mcmc$all$chain, 
                               function(x){names(table(x))[which.max(table(x))]}))

clusters <- list()
for(ch in 1:chains){
  TAB <- table(unlist(mcmc$all[mcmc$all$chain == ch,grep("U\\\\[",colnames(mcmc$all))]))
  TAB <- TAB[order(TAB, decreasing = TRUE)]
  clusters[[ch]] <- as.numeric(names(TAB)[1:modeGplus[ch]])
}
TAB # clusters are labelled from 0 to mcmc$G-1
clusters # the labels of non-empty clusters in each chain 

Us <- mcmc$all[,c("chain", paste0("U[",1:n,"]"))]
clusteringU <- matrix(-1, nrow = as.numeric(n), ncol = mcmc$chains)
certaintyU <- matrix(-1, nrow = as.numeric(n), ncol = mcmc$chains)
for(i in 1:n){
  TAB <- table(Us[,i+1], Us$chain)
  # the cluster is the one with the highest frequency
  clusteringU[i,] <- as.numeric(rownames(TAB)[apply(TAB, 2, which.max)])
  # the highest frequency ratio
  certaintyU[i,] <- apply(TAB/mcmc$iter, 2, max)
}
# units with frequency ratio smaller than threshold remain unclassified (group 0)
clusteringU[certaintyU < threshold] <- 0

# proposed clustering by each chain separately
for(ch in 1:mcmc$chains){
  data[,paste0("clustering",ch)] <- clusteringU[data[,id],ch]
  data[,paste0("certainty",ch)] <- certaintyU[data[,id],ch]
}

data1 <- data[data$j==1,]
# Confusion matrices for each chain separately:
table(data1$clustering1, data1$g)
# We can also see the suitable permutation for labels 

###-----------------------------------------------------------------------------
###     Plotting the resulting clustering
###-----------------------------------------------------------------------------

# Plotting numeric outcomes
plot_num_vs_x_grouped(data, y = "ynum", x = "x", group = "g", id = "i", units = 1:n)
plot_num_vs_x_grouped(data, y = "ynum", x = "x", group = "clustering1", id = "i", units = 1:n)

# Plotting binary, ordinal, categorical outcomes
plot_cat_vs_x_grouped(data, y = "ybin", x = "x", group = "g")
plot_cat_vs_x_grouped(data, y = "ybin", x = "x", group = "g")
plot_cat_vs_x_grouped(data, y = "ybin", x = "x", group = "clustering1")
plot_cat_vs_x_grouped(data, y = "yord", x = "x", group = "g", reverse = TRUE)
plot_cat_vs_x_grouped(data, y = "yord", x = "x", group = "g", reverse = FALSE)
plot_cat_vs_x_grouped(data, y = "yord", x = "x", group = "clustering1")
plot_cat_vs_x_grouped(data, y = "ycat", x = "x", group = "g")
plot_cat_vs_x_grouped(data, y = "ycat", x = "x", group = "clustering1")


}
