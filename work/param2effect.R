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

## fixed and group effects parameters
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
## random effects parameters
#sdSigma
sdnames <- param_names[grep("^sdSigma", param_names)]
refnames <- unlist(mcmc$lrannames)
refnames <- paste0(rep(names(mcmc$family), times = lapply(mcmc$lrannames, length)), 
                   "_", refnames)
effect_names[sdnames] <- refnames

#corSigma, Sigma, InvSigma
mat_names <- param_names[grep("Sigma\\[(\\d+),(\\d+)\\]", param_names)]
positions <- regexec("Sigma\\[(\\d+),(\\d+)\\]", mat_names)
extracted_numbers <- matrix(unlist(regmatches(mat_names, positions)), nrow = 3)
coordinates <- cbind(as.numeric(extracted_numbers[2,]), as.numeric(extracted_numbers[3,]))
coord_names <- matrix(refnames[coordinates], ncol = 2)
coord_names <- paste0(coord_names[,1], "|", coord_names[,2])
effect_names[mat_names] <- coord_names

## individual original ids
bis <- grep("^b\\[", param_names)
if(length(bis) > 0){
  effect_names[bis] <- paste0(rep(names(mcmc$numbered_unique_ids), each = length(refnames)),
                              ", ",
                              rep(refnames,length(mcmc$numbered_unique_ids)))
}

Us <- grep("^U\\[", param_names)
if(length(Us) > 0){
  effect_names[Us] <- names(mcmc$numbered_unique_ids)
}

nas <- grep("^naY_", param_names)
matisYna <- matrix(mcmc$isYna, ncol = length(mcmc$family))
na_names <- c()
for(iy in 1:length(mcmc$family)){
  na_names <- c(na_names, 
                names(mcmc$numbered_unique_ids)[data[matisYna[,iy], id]])
}
effect_names[nas] <- na_names

pUigs <- grep("^pUig\\(", param_names)
if(length(pUigs) > 0){
  effect_names[pUigs] <- rep(names(mcmc$numbered_unique_ids), mcmc$G[ch])
}

pUigs <- grep("^pUig_int\\(", param_names)
if(length(pUigs) > 0){
  effect_names[pUigs] <- rep(names(mcmc$numbered_unique_ids), mcmc$G[ch])
}

dis <- grep("^dev_i\\[", param_names)
if(length(dis) > 0){
  effect_names[dis] <- names(mcmc$numbered_unique_ids)
}

