default_param <- function(){
  param <- list(nu_0 = 1,
                nu_1 = 1,
                gamma_a=1, gamma_b=10,
                sd_beta_num_fix = 1,
                sd_beta_num = 1,
                sd_beta_poi_fix = 1,
                sd_beta_poi = 1,
                sd_beta_bin_fix = 1,
                sd_beta_bin = 1,
                sd_beta_ord_fix = 1,
                sd_beta_ord = 1,
                api_prior = 1,
                sd_beta_cat_fix = 1,
                sd_beta_cat = 1,
                init_sd_b = 0.01,
                ae = 1, be = 100,
                InvV = 0.01)
  return(param)
}