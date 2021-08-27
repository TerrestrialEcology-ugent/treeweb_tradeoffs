# hypervolume stuff

## load libraries
library(plyr)
library(tidyverse)
library(brms)

## load data
func <- read.csv("~/Documents/PostDoc_Ghent/Synthesis_stuff/treeweb_synthesis/data/synthesis_responsedata_std.csv", sep = " ")
plots <- read.csv("~/Documents/PostDoc_Ghent/Synthesis_stuff/treeweb_synthesis/data/synthesis_expldata_raw.csv", sep = " ")

## subset the data to focus on a limited subset of indicators
resp <- func[,c("C_stock","CN","Biomass","Cover","Arth_div","Herbivory","Bird_smi","Predation")]

## fit a separate model to the different tree composition
comp <- unique(plots$speccomb)
m <- list()
for(k in comp){
  subs <- resp[which(plots$speccomb == k),]
  m[[k]] <- brm(mvbind(C_stock,CN,Biomass,Cover, Arth_div, Herbivory, Bird_smi, Predation) ~ 1, data = subs)
}

# a function to grab hypervolume
hypervolume <- function(covariance){
  dim <- ncol(covariance)
  eig <- eigen(nearPD(covariance)$mat)$values # use nearPD to ensure ne negative eigenvalues
  sf <- qchisq(0.95, df = dim)
  ax <- sqrt(sf * eig)
  volume <- (2 / dim) * (pi ** (dim / 2)) / factorial((dim / 2) - 1) * prod(ax)
  return(volume)
}

# a function to grab the an estimated variance-covariance matrix from a brms object
cov_brms <- function(model, fun = "median"){
  post <- posterior_samples(model)
  # grab estimate of sd
  sd_var <- apply(post[,grep("sigma",names(post))], 2, match.fun(fun))
  # grab estimate of correlation
  r_var <- apply(post[,grep("rescor",names(post))], 2, match.fun(fun))
  # put this in a matrix format
  R <- data.frame(X = names(r_var),Cor = r_var)
  R <- separate(R, X, c("drop","X1","X2"))
  R$X1 <- factor(R$X1,levels = c("Cstock","CN","Biomass","Cover","Arthdiv",
                                 "Herbivory","Birdsmi","Predation"))
  R$X2 <- factor(R$X2,levels = c("Cstock","CN","Biomass","Cover","Arthdiv",
                                 "Herbivory","Birdsmi","Predation"))
  
  R2 <- dcast(R,X1~X2,value.var="Cor")
  rownames(R2) <- R2[,1]
  R2 <- R2[,-1]
  R2 <- cbind(NA, R2)
  colnames(R2)[1] <- "Cstock"
  R2 <- rbind(R2,NA)
  rownames(R2)[nrow(R2)] <- "Predation"
  R2 <- as.matrix(R2)
  R2[lower.tri(R2)] <- t(R2)[lower.tri(R2)]
  diag(R2) <- sd_var
  
 return(R2)
}

# apply this to the different models
ldply(m, function(x) hypervolume(cov_brms(x)))

# compare this with empirical covariance matrices
cc <- list()
for(k in comp){
  subs <- resp[which(plots$speccomb == k),]
  cc[[k]] <- cov(subs)
}

ldply(cc,function(x) hypervolume(x))
