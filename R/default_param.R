default_param <- function(){
  param <- list(InvSigma_df = 1,
                InvQ_df = 1,
                prec_num_shp=1, prec_num_rte=10,
                beta_num_fix_sd = 1,
                beta_num_sd = 1,
                beta_poi_fix_sd = 1,
                beta_poi_sd = 1,
                beta_bin_fix_sd = 1,
                beta_bin_sd = 1,
                beta_ord_fix_sd = 1,
                beta_ord_sd = 1,
                api_prior = 1,
                beta_cat_fix_sd = 1,
                beta_cat_sd = 1,
                init_b_sd = 0.01,
                e0_shp = 1, e0_rte = 100,
                InvV = 0.01)
  return(param)
}