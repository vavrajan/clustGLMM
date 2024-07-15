as_coda <- function(mcmc, burnin=0, thin=1){
  if(mcmc$howsave == "matrix"){
    codamcmc <- FromMatrixtoCoda(mcmc, burnin=burnin, thin=thin)
  }
  
  if(mcmc$howsave == "list"){
    codamcmc <- FromListtoCoda(mcmc, burnin=burnin, thin=thin)
  }
  
  return(codamcmc)
}