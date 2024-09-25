default_save <- function(){
  save <- c(T, T,        # beta_num_fix, beta_num,
            T, T, T,     # prec_num, sd_num, var_num
            T, T,        # beta_poi_fix, beta_poi
            T, T,        # beta_bin_fix, beta_bin
            T, T,        # beta_ord_fix, beta_ord
            T, F, T,     # c_ord, a_ord, pi_ord, 
            T, T,        # beta_cat_fix, beta_cat
            T,T,T,T,F,   # InvSigma, Sigma, sdSigma, corSigma, detInvSigma
            F, F, F,     # InvQ, Q, detInvQ
            F,           # b
            T, T,        # w, ng
            T, F, T,     # loglik, pUig, U
            T, T,        # Gplus, e0
            F            # naY
  )
  names(save) <- c("beta_num_fix", "beta_num",
                   "prec_num", "sd_num", "var_num",
                   "beta_poi_fix", "beta_poi",
                   "beta_bin_fix", "beta_bin",
                   "beta_ord_fix", "beta_ord",
                   "c_ord", "a_ord", "pi_ord",
                   "beta_cat_fix", "beta_cat",
                   "InvSigma", "Sigma", "sdSigma", "corSigma", "detInvSigma",
                   "InvQ", "Q", "detInvQ",
                   "b",
                   "w", "ng",
                   "loglik", "pUig", "U",
                   "Gplus", "e0",
                   "naY")
  return(save)
}