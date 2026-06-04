test_that("A zero offset is automatically added to formula if missing.", {

  expect_equal({
    form <- list(y = list(fixed = ~ 0, group = ~ x, random = ~ 1))
    mcmc <- suppressMessages(clustGLMM(
      formula = form, id = "id", family = c(y = "num"), 
      data = data.frame(id = 1:10,
                        y = 1:10,
                        x = rep(c(0,1), each = 5),
                        o = rep(c(0,-10), each = 5)),
      G = 2, nchains = 1, iter = 10, standardize = TRUE))
    mcmc$formula
  }, 
  {
    list(y = list(fixed = ~ 0, group = ~ x, random = ~ 1, offset = ""))
  })
})

test_that("Even when not standardized, if intercept in both fixed and group, it becomes g-specific", {
  expect_equal({
    form <- list(y = list(fixed = ~ 1, group = ~ x, random = ~ 1))
    mcmc <- suppressMessages(clustGLMM(
      formula = form, id = "id", family = c(y = "num"), 
      data = data.frame(id = 1:10,
                        y = 1:10,
                        x = rep(c(0,1), each = 5),
                        o = rep(c(0,-10), each = 5)),
      G = 2, nchains = 1, iter = 10, standardize = FALSE))
    list(mcmc$lfixnames,
         mcmc$lgrpnames)
  }, 
  {
    list(list(y = character()),
         list(y = c("(Intercept)", "x")))
  })
  
})

test_that("When not standardized, if intercept not group, it stays not group-specific.", {
  
  expect_equal({
    form <- list(y = list(fixed = ~ 1, group = ~ x+0, random = ~ 1))
    mcmc <- suppressMessages(clustGLMM(
      formula = form, id = "id", family = c(y = "num"), 
      data = data.frame(id = 1:10,
                        y = 1:10,
                        x = rep(c(0,1), each = 5),
                        o = rep(c(0,-10), each = 5)),
      G = 2, nchains = 1, iter = 10, standardize = FALSE))
    list(mcmc$lfixnames,
         mcmc$lgrpnames)
  }, 
  {
    list(list(y = "(Intercept)"),
         list(y = "x"))
  })
  
})

test_that("When standardized, if intercept not group, it becomes group-specific.", {
  
  expect_warning(
    expect_equal({
      form <- list(y = list(fixed = ~ 1, group = ~ x-1, random = ~ 1))
      mcmc <- suppressMessages(clustGLMM(
        formula = form, id = "id", family = c(y = "num"), 
        data = data.frame(id = 1:10,
                          y = 1:10,
                          x = rep(c(0,1), each = 5),
                          o = rep(c(0,-10), each = 5)),
        G = 2, nchains = 1, iter = 10, standardize = TRUE))
      list(mcmc$lfixnames,
           mcmc$lgrpnames)
    }, 
    {
      list(list(y = character()),
           list(y = c("(Intercept)", "x")))
    })
  )
  
})

test_that("When standardized, if all fixed, remains fixed, even intercept does not change to group-specific.", {

    expect_equal({
      form <- list(y = list(fixed = ~ x, group = ~ 0, random = ~ 1))
      mcmc <- suppressMessages(clustGLMM(
        formula = form, id = "id", family = c(y = "num"), 
        data = data.frame(id = 1:10,
                          y = 1:10,
                          x = rep(c(0,1), each = 5),
                          o = rep(c(0,-10), each = 5)),
        G = 2, nchains = 1, iter = 10, standardize = TRUE))
      list(mcmc$lfixnames,
           mcmc$lgrpnames)
    }, 
    {
      list(list(y = c("(Intercept)", "x")),
           list())
    })
  
})
