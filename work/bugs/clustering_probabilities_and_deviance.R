library("clustGLMM")
data("longitudinal_mixed_type_data")
id <- "i"

family <- c("num", "poi", "bin", "ord", "cat")
names(family) <- Ys <- c("ynum", "ypoi", "ybin", "yord", "ycat")

formula <- list()
for(y in Ys){
  formula[[y]] <- list()
  formula[[y]]$fixed <- ~ f
  formula[[y]]$group <- ~ x
  formula[[y]]$random <- ~ 1
  formula[[y]]$offset <- ""
}

Gmax <- 5     
nchains <- 2  

set.seed(123456789)
mcmc <- clustGLMM(formula = formula, id = id, family = family, 
                  data = longitudinal_mixed_type_data, 
                  G = Gmax, iter = 100, nchains = nchains)
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

probs <- clustering_probabilities_and_deviance(mcmc, id, newdata,
                                               dynamic_prob = FALSE,
                                               start = 1,
                                               end = 100,
                                               thin = 1,
                                               chains = 1:nchains,
                                               tolerance = 1e-7, 
                                               maxiter = 100,
                                               maxnrep = 20, 
                                               NGQ = 1)

class(probs)
methods(class = "clustglmm")
print(probs)
summary(probs)
predict(probs, newdata = newdata[1:10,])
plot(probs, chains = 1)
plot(probs, what = "deviance")
plot(probs, what = "ng")

