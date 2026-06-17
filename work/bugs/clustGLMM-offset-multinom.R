library("clustGLMM")
id <- "i"
family <- "cat"
names(family) <- Ys <- "ycat"

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
longitudinal_mixed_type_data1$offsetcat <- matrix(0, nrow = nrow(longitudinal_mixed_type_data),
                                                  ncol = nlevels(longitudinal_mixed_type_data$ycat))
longitudinal_mixed_type_data1$offsetcat <- rep(0, nrow = nrow(longitudinal_mixed_type_data))

formula[["ycat"]]$offset <- "offsetcat"


set.seed(123456789)
mcmc <- suppressMessages(clustGLMM(formula = formula, id = id, family = family, 
                                   data = longitudinal_mixed_type_data1, 
                                   G = Gmax, iter = 100, nchains = nchains,
                                   standardize = FALSE))
identical(mcmc, mcmc0)
