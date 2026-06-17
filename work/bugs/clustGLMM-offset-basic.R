library("clustGLMM")
id <- "i"
family <- c("num", "poi", "bin")
names(family) <- Ys <- c("ynum", "ypoi", "ybin")

formula <- list()
for(y in Ys){
  formula[[y]] <- list()
  formula[[y]]$fixed <- ~ f
  formula[[y]]$group <- ~ 1
  formula[[y]]$random <- ~ 1
  formula[[y]]$offset <- ""
}

Gmax <- 10     
nchains <- 1   

data("longitudinal_mixed_type_data")
longitudinal_mixed_type_data$f <- as.integer(longitudinal_mixed_type_data$f == 1)
set.seed(123456789)
mcmc0 <- suppressMessages(clustGLMM(formula = formula, id = id, family = family, 
                                   data = longitudinal_mixed_type_data, 
                                   G = Gmax, iter = 100, nchains = nchains,
                                   standardize = FALSE))

longitudinal_mixed_type_data1 <- longitudinal_mixed_type_data
longitudinal_mixed_type_data1$offsety <- 0
formula[["ynum"]]$offset <- "offsety"
longitudinal_mixed_type_data1$offsetp <- 0
formula[["ypoi"]]$offset <- "offsetp"
longitudinal_mixed_type_data1$offsetb <- 0
formula[["ybin"]]$offset <- "offsetb"

set.seed(123456789)
mcmc <- suppressMessages(clustGLMM(formula = formula, id = id, family = family, 
                                   data = longitudinal_mixed_type_data1, 
                                   G = Gmax, iter = 100, nchains = nchains,
                                   standardize = FALSE))
identical(mcmc$inits, mcmc0$inits)
all.equal(mcmc, mcmc0)
elements <- setdiff(names(mcmc),
                    c("formula", "noff", "loffnames", "offnames", "centers", "scales"))
identical(mcmc[elements], mcmc0[elements])
