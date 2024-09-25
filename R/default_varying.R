default_varying <- function(){
  varying <- c(T, T, 
               F, F, # common covariance matrix of random effects
               F) 
  names(varying) <- c("prec_num", "c_ord", "InvSigma", "InvQ", "naY")
  return(varying)
}