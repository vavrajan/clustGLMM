
test_that("clustGLMM starts with the same initial values regardless of offset", {
  
  expect_equal({
    set.seed(123456789)
    mcmc <- suppressMessages(clustGLMM(
      formula = list(y = list(fixed = ~ 0, group = ~ x, random = ~ 1)), 
      id = "id", family = c(y = "num"), 
      data = data.frame(id = 1:10,
                        y = 1:10,
                        x = rep(c(0,1), each = 5),
                        o = rep(c(0,-10), each = 5)),
      G = 2, nchains = 1, iter = 10, standardize = TRUE))
    mcmc$inits
  }, {
    set.seed(123456789)
    mcmc <- suppressMessages(clustGLMM(
      formula = list(y = list(fixed = ~ 0, group = ~ x, random = ~ 1, offset = "o")), 
      id = "id", family = c(y = "num"), 
      data = data.frame(id = 1:10,
                        y = 11:20,
                        x = rep(c(0,1), each = 5),
                        o = rep(-10, 10)),
      G = 2, nchains = 1, iter = 10, standardize = TRUE))
    mcmc$inits
  })
  
  expect_equal({
    set.seed(123456789)
    mcmc <- suppressMessages(clustGLMM(
      formula = list(y = list(fixed = ~ 0, group = ~ x, random = ~ 1)), 
      id = "id", family = c(y = "poi"), 
      data = data.frame(id = 1:10,
                        y = 1:10,
                        x = rep(c(0,1), each = 5),
                        o = rep(c(0,-10), each = 5)),
      G = 2, nchains = 1, iter = 10, standardize = TRUE))
    c(mcmc$inits[[1]]$beta_poi$y[[1]]["x"],
      mcmc$inits[[1]]$beta_poi$y[[2]]["x"])
  }, {
    set.seed(123456789)
    mcmc <- suppressMessages(clustGLMM(
      formula = list(y = list(fixed = ~ 0, group = ~ x, random = ~ 1, offset = "o")), 
      id = "id", family = c(y = "poi"), 
      data = data.frame(id = 1:10,
                        y = c(1:10) * 2,
                        x = rep(c(0,1), each = 5),
                        o = rep(-log(2), 10)),
      G = 2, nchains = 1, iter = 10, standardize = TRUE))
    c(mcmc$inits[[1]]$beta_poi$y[[1]]["x"],
      mcmc$inits[[1]]$beta_poi$y[[2]]["x"])
  })
  
  expect_warning(suppressMessages(clustGLMM(
    formula = list(y = list(fixed = ~ 0, group = ~ x, random = ~ 1, offset = "o")), 
    id = "id", family = c(y = "cat"), 
    data = data.frame(id = 1:10,
                      y = rep(1:5, 2),
                      x = rep(c(0,1), each = 5),
                      o = rep(-log(2), 10)),
    G = 2, nchains = 1, iter = 10, standardize = TRUE)))
})
