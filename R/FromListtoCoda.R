FromListtoCoda <- function(mcmc, burnin=0, thin=1){
  mcmc <- FromListtoMatrix(mcmc)
  codamcmc <- FromMatrixtoCoda(mcmc, burnin=burnin, thin=thin)
  return(codamcmc)
}