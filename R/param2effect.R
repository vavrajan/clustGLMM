param2effect <- function(param, 
                         settings,
                         lfixnames,
                         lgrpnames,
                         lrannames,
                         loffnames){
  if(grepl("^beta", param)){
    sparam <- strsplit(param, "_")
    str_extract(param, "(?<=\\[)\\s*\\d+")
    if(length(sparam[[1]]) == 3){
      # fix
      
    }else{
      # grp
      
    }
  }
  
  if(grepl("Sigma", param)){
    if(grepl("^sd", param)){
      
    }else{
      if(!grepl("^det", param)){
        
      }
    }
  }
}

# or directly
param_names <- colnames(mcmc$draws[[ch]])
effect_names <- rep(NA, length(param_names))
names(effect_names) <- param_names

for(y in names(mcmc$family)){
  k <- ifelse(mcmc$family[y] == "cat", mcmc$Kcat[y], 1)
  # fixed
  pattern <- paste0("^beta_", mcmc$family[y], "_fix_", y)
  bnames <- param_names[grep(pattern, param_names)]
  if(length(bnames) > 0){
    effect_names[bnames] <- rep(mcmc$lfixnames[[y]], each = k)
  }
  
  # group
  pattern <- paste0("^beta_", mcmc$family[y], "_", y)
  bnames <- param_names[grep(pattern, param_names)]
  if(length(bnames) > 0){
    effect_names[bnames] <- rep(mcmc$lgrpnames[[y]], each = k)
  }
}
# random
sdnames <- param_names[grep("^sdSigma", param_names)]
refnames <- unlist(mcmc$lrannames)
refnames <- paste0(names(refnames), "_", refnames)
effect_names[bnames] <- refnames
#TODO other Sigma parameters
