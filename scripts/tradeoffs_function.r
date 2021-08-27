#################################################
# R helper function from the trade-off analysis #
#################################################


## helper functions
### this function reverse variables that need to be minimized
revv <- function(vec){
  return(- vec + (max(vec) + min(vec)))
}

# derive desirability as weighted average
desirab <- function(linpred, direction, importance, quantiles = c(0.1, 0.5, 0.9), sum = FALSE){
  for(i in 1:dim(linpred)[3]){ # reverse the functions that need to
    if(direction[i] == "minimize"){
      linpred[,,i] <- t(apply(linpred[,,i], 1, revv))
    }
  }
  
  avg_fw <- apply(linpred[,,1:8], c(1, 2), function(x) sum(importance[1:8] * x) / sum(importance[1:8]))
  if(sum){
    avg_fs <- rowSums(avg_fw)
    avg_fqq <- quantile(avg_fs, probs = quantiles)
    avg_fout <- data.frame(LCI = avg_fqq[1], Median = avg_fqq[2], UCI = avg_fqq[3], type = "functioning")
  }
  else{
    avg_fqq <- apply(avg_fw, 2, quantile, probs = quantiles) 
    avg_fout <- data.frame(LCI = avg_fqq[1,], Median = avg_fqq[2,], UCI = avg_fqq[3,], type = "functioning")
  }
  
  
  avg_dw <- apply(linpred[,,9:16], c(1, 2), function(x) sum(importance[9:16] * x) / sum(importance[9:16]))
  if(sum){
    avg_ds <- rowSums(avg_dw)
    avg_dqq <- quantile(avg_ds, probs = quantiles)
    avg_dout <- data.frame(LCI = avg_dqq[1], Median = avg_dqq[2], UCI = avg_dqq[3], type = "diversity")
  }
  else{
    avg_dqq <- apply(avg_dw, 2, quantile, probs = quantiles) 
    avg_dout <- data.frame(LCI = avg_dqq[1,], Median = avg_dqq[2,], UCI = avg_dqq[3,], type = "diversity")
  }
  avg_out <- rbind(avg_fout, avg_dout)
  
  
  return(avg_out)
}

# function to compute the percent of variance explained by the different covariates in the model
get_sd <- function(post, name = "Cstock"){
  # grab fixed effects
  fixx <- post[,grep(paste0("b_",name), names(post))]
  # the sd for speccomb
  comp_sd <- apply(fixx[,1:7],1,function(b) sd(b[as.numeric(brm_dat$speccomb)]))
  # the sd for edge
  edge_sd <- abs(fixx[,9])
  # the sd for prox
  prox_sd <- abs(fixx[,8])
  # the sd for density
  # dens_sd <- abs(fixx[,10])
  # the sd of the resid
  res_sd <- resid_all[,name]
  sd.all <- cbind(comp_sd, edge_sd, prox_sd, res_sd)
  out <- broom.mixed::tidyMCMC(100 * sd.all / rowSums(sd.all), conf.int = TRUE, conf.level = 0.8)
  out$variable <- name
  out$fraction <- c("composition","edge","proximity", "residual")
  return(out)
}