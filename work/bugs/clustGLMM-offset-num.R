settings <- expand.grid(offsety = c(FALSE, TRUE),
                        standardize = c(FALSE, TRUE))

library("clustGLMM")
id <- "i"
family <- "num"
names(family) <- Ys <- "ynum"

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
longitudinal_mixed_type_data1 <- longitudinal_mixed_type_data
longitudinal_mixed_type_data1$offsety <- c(2, 0)
longitudinal_mixed_type_data1$ynum1 <- with(longitudinal_mixed_type_data1,
                                           ynum + offsety)

longitudinal_mixed_type_data2 <- longitudinal_mixed_type_data
longitudinal_mixed_type_data2$offsety <- c(1, 0)
longitudinal_mixed_type_data2$ynum1 <- with(longitudinal_mixed_type_data2,
                                           ynum + offsety)

checks <- sapply(seq_len(nrow(settings)), function(i) {
    offsety <- settings[i, "offsety"]
    standardize <- settings[i, "standardize"]

    formula1 <- formula
    family1 <- family
    if (offsety) {
        formula1[["ynum"]]$offset <- "offsety"
        names(family1)[1] <- "ynum1"
        names(formula1)[1] <- "ynum1"
    }
    
    set.seed(123456789)
    mcmc <- suppressMessages(clustGLMM(formula = formula1, id = id, family = family1, 
                                       data = longitudinal_mixed_type_data1, 
                                       G = Gmax, iter = 100, nchains = nchains,
                                       standardize = standardize))
    set.seed(123456789)
    mcmc0 <- suppressMessages(clustGLMM(formula = formula1, id = id, family = family1, 
                                        data = longitudinal_mixed_type_data2, 
                                        G = Gmax, iter = 100, nchains = nchains,
                                        standardize = standardize))
    c(inits_tol = isTRUE(all.equal(mcmc0$inits, mcmc$inits, tolerance = 10^-6)),
      inits = isTRUE(all.equal(mcmc0$inits, mcmc$inits, tolerance = 10^-12)),
      all = isTRUE(all.equal(mcmc0, mcmc)))
})

settings$valid_init <- checks["inits_tol",]
settings$identical_init <- checks["inits",]
settings$identical_all <- checks["all",]
settings

