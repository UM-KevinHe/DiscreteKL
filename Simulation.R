library(mvtnorm)
library(discSurv)
library(survival)
library(tidyr)
library(Rcpp)
library(parallel)
library(devtools)
#devtools::install_github("UM-KevinHe/DiscreteKL")
library(DiscreteKL)
###### this simulation code can be applied with parallel computing
##### detect max number of cores on this machine
max.cores = parallel::detectCores()
##### by default, we will use floor(0.5*max.cores) to conduct the analysis 
####### internal sample size and censoring rate
####### censoring rate is set as 40%
####### sample size is 150 or 200
n_local = 150
cens_rate = 40
#### number of simulation replicates
#here we set it as 50, note that in the paper we replicated the simulations 500 times
rep = 50

loglik <- as.data.frame(matrix(rep(0, 4*rep), rep, 4))
names(loglik) <- c("KL_logit", "prior", "local", "stacked")
eta <- as.data.frame(matrix(rep(0, 1*rep), rep, 1))
names(eta) <- c("KL_prior")
beta_KL <- as.data.frame(matrix(rep(0, 20*rep), rep, 20))
beta_internal <- as.data.frame(matrix(rep(0, 20*rep), rep, 20))

n_local = n_local
n_prior = 10000
n_test = 1000
p_local = 10
p_prior = 10

if (cens_rate==40){
  cens_upper = 11
}
if (cens_rate==20){
  cens_upper = 1000
}

Z.char_prior <- paste0('Z', 1:p_prior)
Z.char_local <- paste0('Z', 1:p_local)

#prior and local data parameters
local_beta <- c(2,-1,2,3,-1,4,-1,3,4,-1)
external_beta <- local_beta
day_effect <- c(-6.0, -6.0, -6.0, -4.5, -4.5,
                -4.5, -3.0, -3.0, -1.5, -1.5)

#simulate prior data assuming local_beta is the beta for true model
prior_data <- sim.disc(local_beta, day_effect, n_prior, Z.char_prior, cens_upper)
1-mean(prior_data$status)
X_prior <- dplyr::select(prior_data, -c("time", "status", "Z9", "Z10"))
X_prior <- as.matrix(X_prior)

betap=discSurv_logit(prior_data$time, X_prior, prior_data$status)
prior_beta <- as.vector(c(betap$beta_v[1:8],0,0))
prior_day_effect <- as.vector(betap$beta_t)
betap$iter

local_beta <- as.matrix(local_beta)
external_beta <- as.matrix(external_beta)
prior_beta <- as.matrix(prior_beta)

data_internal_list <- vector(mode="list", length=rep)

KL_EST <- function(x){
  KL_prior <- kl_logit(prior_day_effect, prior_beta, 0, 10, 0.1, x)
  return(KL_prior)
}

for (i in 1:rep){
  set.seed(i)
  data_internal_list[[i]] <- sim.disc(local_beta, day_effect, n_local, Z.char_local, cens_upper)
}

models_KL <- mclapply(data_internal_list, KL_EST, mc.cores = floor(0.5*max.cores))

#################################
for (i in 1:rep){
set.seed(i)
  if(class(models_KL[[i]])!="try-error"){
#############Prior
local_data <- data_internal_list[[i]]
1-mean(local_data$status)
external_data <- sim.disc(external_beta, day_effect, n_test, Z.char_local, cens_upper)

Z_local <- dplyr::select(local_data, -c("time", "status"))
Z_external <- dplyr::select(external_data, -c("time", "status"))

#prior
loglik$prior[i] <- DiscLoglik_logit(external_data$time, Z_external, external_data$status,prior_day_effect,prior_beta)

#local
estLocal <- discSurv_logit(local_data$time, Z_local, local_data$status)
loglik$local[i] <- DiscLoglik_logit(external_data$time, Z_external, external_data$status,estLocal$beta_t,estLocal$beta_v)
beta_internal[i,] <- as.numeric(c(estLocal$beta_v,estLocal$beta_t))

#KL
#KL_prior <- kl_lik(prior_day_effect, prior_beta, 0, 10, 0.25, local_data)
KL_prior <- models_KL[[i]]
estKL <- KL_prior$model
eta$KL_prior[i] <- KL_prior$eta
loglik$KL_logit[i] <- DiscLoglik_logit(external_data$time, Z_external, external_data$status,estKL$beta_t,estKL$beta_v)
beta_KL[i,] <- as.numeric(c(estKL$beta_v,estKL$beta_t))

#stacked regression: only beta from the prior model: estimate day effect
Z_local_stacked <- data.frame(Z_stacked = as.matrix(Z_local[,1:p_prior])%*%as.matrix(prior_beta))
Z_external_stacked <- data.frame(Z_stacked = as.matrix(Z_external[,1:p_prior])%*%as.matrix(prior_beta))
estStacked <- discSurv_logit(local_data$time, Z_local_stacked, local_data$status)
loglik$stacked[i] <- DiscLoglik_logit(external_data$time, Z_external_stacked, external_data$status,estStacked$beta_t,estStacked$beta_v)
  }
}

deviance <- -loglik
deviance <- deviance/n_test

cols <- c("#fcd575", "#6C71C2", "#abcea6", "#b3c6e8")
library(ggplot2)
mr <- data.frame(Method = c(
  rep("KL", rep),
  rep("External", rep),
  rep("Internal", rep),
  rep("Stacked", rep)),
  loglik = c(
    deviance$KL_logit,
    deviance$prior,
    deviance$local,
    deviance$stacked))
mr$Method <- factor(mr$Method,
                    levels = c('Internal','KL', 'Stacked','External'),ordered = TRUE)
g <- ggplot(mr, aes(x=Method, y=loglik, fill=Method)) + 
  geom_boxplot() + labs(y = "Predictive Deviance", x = "Method")+  scale_colour_manual(
    values = cols,
    aesthetics = c("fill")
  )+ theme_bw() + theme(panel.border = element_blank(), 
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(size = 14))+theme(axis.text = element_text(size = 14))+theme(legend.position = "none")+theme(axis.title.x = element_blank())
print(g)





