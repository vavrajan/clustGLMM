as_coda <- function(mcmc, burnin=0, thin=1){
  if(mcmc$howsave == "matrix"){
    codamcmc <- from_matrix_to_coda(mcmc, burnin=burnin, thin=thin)
  }
  
  if(mcmc$howsave == "list"){
    codamcmc <- from_list_to_coda(mcmc, burnin=burnin, thin=thin)
  }
  
  return(codamcmc)
}