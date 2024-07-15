default_tuning <- function(){
  tuning <- list()
  tuning$integer <- list(
    freq_proposal_update = 10,
    times_proposal = 10,
    maxiter = 25,
    maxnrep = 100,
    kspec_bi_cat = 0
  )
  tuning$double <- list(
    const_proposal_beta_poi_fix = 1,
    const_proposal_beta_poi = 1,
    const_proposal_beta_bin_fix = 1,
    const_proposal_beta_bin = 1,
    const_proposal_beta_ord_fix = 0.5,
    const_proposal_beta_ord = 0.5,
    const_proposal_beta_cat_fix = 1,
    const_proposal_beta_cat = 1,
    const_proposal_a_ord = 0.5,
    const_proposal_b = 0.5,
    const_proposal_e0 = 1,
    tolerance = 1e-7
  )
  return(tuning)
}