################################################################################
# Expanded Simulation Study: Comparing Logit vs cauchit Link Functions
# 
# This code extends the original simulation study to examine what happens
# when using the cauchit link function (a non-canonical link) instead of
# the logit link function (the canonical link for binomial outcomes).
################################################################################

## loading libraries
library(pbapply)
library(mvtnorm)
library(boot)
library(ggplot2)

# Simulation parameters
n_sims <- 10000
n_sample <- 200
n_boot <- 10

## Ueeful functions:
ilogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

################################################################################
# Main estimation function with link parameter
################################################################################

ITC_estimates <- function(Y_all, X_all, S_all, IPD_for_Study0=TRUE, M=10000, boot_n=1000, outcome_link="logit"){
  
  n1 <- sum(S_all==1)
  n0 <- sum(S_all==0)
  n <- n1 + n0
  
  # calculating overlap of first covariate:
  d1dens <- density(X_all[S_all==1,], 
                    from = -10, 
                    to = 10)
  d2dens <- density(X_all[S_all==0,], 
                    from = -10, 
                    to = 10)
  joint <- pmin(d1dens$y, d2dens$y)
  
  df2 <- data.frame(x = rep(d1dens$x, 3), 
                    y = c(d1dens$y, d2dens$y, joint),
                    Data = rep(c("D1", "D2", "overlap"), each = length(d1dens$x)))
  
  overlap <- 2*sum(joint) / sum(d1dens$y, d2dens$y)
  
  # g-function is the log odds (this defines the estimand scale, kept as logit):
  g_function <- function(p){log(p/(1-p))}
  
  #############################################
  # Naive estimate
  mu1_naive <- mean(Y_all[S_all==1])
  mu0_naive <- mean(Y_all[S_all==0])
  ATC_naive <- g_function(mu1_naive) - g_function(mu0_naive)
  
  SE_g_mu0 <- summary(glm(Y_all[S_all==0]~1, family="binomial"))$coef[2]
  SE_g_mu1 <- summary(glm(Y_all[S_all==1]~1, family="binomial"))$coef[2]
  
  naive_function <- function(data, indices) {
    temp <- data[indices,]
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu1_naive <- mean(temp[temp[,"S_all"]==1,"Y_all"])
    ATC_naive <- g_function(mu1_naive) - g_function(mu0_naive)
    return(c(g_function(mu0_naive), g_function(mu1_naive), ATC_naive))
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all), 
                       statistic=naive_function, R=boot_n,  
                       strata=S_all,
                       parallel = "multicore") 
  
  SE_g_mu0_naive <- sd(boot_samples$t[,1], na.rm=TRUE)
  SE_g_mu1_naive <- sd(boot_samples$t[,2], na.rm=TRUE)
  SE_ATC <- sd(boot_samples$t[,3], na.rm=TRUE)
  
  ATC_naive_CI <- c(boot_samples$t0[3] - abs(qnorm(0.025))*sd(boot_samples$t[,3], na.rm=TRUE), 
                    boot_samples$t0[3] + abs(qnorm(0.025))*sd(boot_samples$t[,3], na.rm=TRUE))
  
  #############################################    
  if(!IPD_for_Study0){
    M <- 100000
    X_0_mean <- apply(X_all[S_all==0,], 2, mean)
    X_0_sd <- apply(X_all[S_all==0,], 2, sd)
    X_all <- rbind(X_all[S_all==1,], rmvnorm(M, X_0_mean, diag(X_0_sd^2)))
    S_all <- c(S_all[S_all==1], rep(0, M))
    n1 <- sum(S_all==1)
    n0 <- sum(S_all==0)
    n <- n1 + n0
  }
  
  #############################################
  # Inverse odds weighting estimate
  # Note: Propensity score model always uses logit (this is for study membership)
  ps_mod <- glm(S_all ~ X_all, family="binomial")
  ps <- predict(ps_mod, type="response")
  w_IOW <- (1-ps)*S_all/ps + (1-S_all)
  
  mu1_IOW1 <- (1/n0)*sum(w_IOW[1:n1]*Y_all[1:n1])
  mu0_IOW1 <- mu0_naive
  ATC_IOW1 <- g_function(mu1_IOW1) - g_function(mu0_IOW1)
  
  mu1_IOW2 <- (1/sum(w_IOW[1:n1]))*sum(w_IOW[1:n1]*Y_all[1:n1])
  mu0_IOW2 <- mu0_naive
  ATC_IOW2 <- g_function(mu1_IOW2) - g_function(mu0_IOW2)
  
  IOW1_function <- function(data, indices){
    temp <- data[indices,]
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[,grep("X_all",colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    mu1_IOW1 <- (1/sum(temp[,"S_all"]==0))*sum(w_IOW[temp[,"S_all"]==1]*temp[,"Y_all"][temp[,"S_all"]==1])
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    ATC_IOW1 <- g_function(mu1_IOW1) - g_function(mu0_naive)
    return(ATC_IOW1)
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=IOW1_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  ATC_IOW1_CI <- c(boot_samples$t0[1] - abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE), 
                   boot_samples$t0[1] + abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE))
  
  IOW2_function <- function(data, indices){
    temp <- data[indices,]
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[,grep("X_all",colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    mu1_IOW2 <- (1/sum(w_IOW[temp[,"S_all"]==1]))*sum(w_IOW[temp[,"S_all"]==1]*temp[,"Y_all"][temp[,"S_all"]==1])
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    ATC_IOW2 <- g_function(mu1_IOW2) - g_function(mu0_naive)
    return(ATC_IOW2)
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=IOW2_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  ATC_IOW2_CI <- c(boot_samples$t0[1] - abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE), 
                   boot_samples$t0[1] + abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE))
  
  #############################################
  # 3.3 Entropy balancing (matching-adjusted indirect comparison)
  
  objfn <- function(a1, X){ sum(exp(X %*% a1)) }
  gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
  cov_names <- paste0("X.", 1:dim(X_all)[2])
  AC.IPD <- data.frame(y=Y_all[S_all==1], X= X_all[S_all==1,])
  BC.ALD <- data.frame(matrix(apply(X_all[S_all==0,], 2, mean), 1,))
  colnames(BC.ALD) <- paste0("mean.", cov_names)
  X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                  as.matrix((BC.ALD[, c(paste("mean.", cov_names, sep=""))])), '-')
  
  gamma <- optim(par = rep(0, dim(X.EM.0)[2]), 
                 fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
  wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
  ESS <- ((sum(wt_EB))^2)/sum(wt_EB^2)
  
  mu1_EB <- sum(wt_EB*Y_all[S_all==1])
  mu0_EB <- mu0_naive
  ATC_EB <- g_function(mu1_EB) - g_function(mu0_EB)
  
  ATC_EB_function <- function(data, indices){
    temp <- data[indices,]
    AC.IPD <- data.frame(y=temp[,"Y_all"][temp[,"S_all"]==1], X= temp[temp[,"S_all"]==1, grep("X_all", colnames(temp))])
    BC.ALD <- data.frame(matrix(apply(temp[temp[,"S_all"]==0, grep("X_all", colnames(temp))], 2, mean), 1,))
    objfn <- function(a1, X){ sum(exp(X %*% a1)) }
    gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
    cov_names <- paste0("X.X_all.", 1:dim(X_all)[2])
    colnames(BC.ALD) <- paste0("mean.", cov_names)
    X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                    as.matrix((BC.ALD[, c(paste("mean.", cov_names, sep=""))])), '-')
    
    gamma <- optim(par = rep(0, dim(X.EM.0)[2]), 
                   fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
    wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
    
    mu1_EB <- sum(wt_EB*temp[,"Y_all"][temp[,"S_all"]==1])
    mu0_naive <- mean(temp[temp[,"S_all"]==0, "Y_all"])
    ATC_EB <- g_function(mu1_EB) - g_function(mu0_naive)
    return(ATC_EB)
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_EB_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  
  ATC_EB_CI <- c(boot_samples$t0[1] - abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE), 
                 boot_samples$t0[1] + abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE))
  
  #############################################
  # 3.4 Doubly robust augmented weighting estimators
  
  ################################
  # The G-computation estimator
  # NOTE: This uses the outcome_link parameter
  
  outcome_model <- glm(y ~ ., 
                       data = data.frame(y=c(Y_all[S_all==1]), 
                                         x=(X_all[S_all==1,])), 
                       family = binomial(link = outcome_link))
  
  Y1_hat <- predict(outcome_model, newdata = data.frame(x=X_all[S_all==0,]), 
                    type = "response")
  
  mu1_GCOMP <- (1/n0)*sum(Y1_hat)
  mu0_GCOMP1 <- mu0_naive
  ATC_GCOMP <- g_function(mu1_GCOMP) - g_function(mu0_GCOMP1)
  
  # Create a version of ATC_GCOMP_function that uses the specified link
  ATC_GCOMP_function <- function(data, indices){
    temp <- data[indices,]
    outcome_model <- glm(y ~ ., 
                         data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                           x=temp[temp[,"S_all"]==1, grep("X_all", colnames(temp))]), 
                         family = binomial(link = outcome_link))
    
    Y1_hat <- predict(outcome_model, newdata = data.frame(x=temp[temp[,"S_all"]==0, grep("X_all", colnames(temp))]), 
                      type = "response")
    
    mu1_GCOMP <- (1/sum(temp[,"S_all"]==0))*sum(Y1_hat)
    mu0_naive <- mean(temp[temp[,"S_all"]==0, "Y_all"])
    ATC_GCOMP <- g_function(mu1_GCOMP) - g_function(mu0_naive)
    return(ATC_GCOMP)
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_GCOMP_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  
  ATC_GCOMP_CI <- c(boot_samples$t0[1] - abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE), 
                    boot_samples$t0[1] + abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE))
  
  ################################
  # DR1 (augmented inverse odds weighting estimator):
  # Using Inverse odds weights
  # NOTE: This uses the outcome_link parameter
  
  Y1_hat_all <- predict(outcome_model, newdata = data.frame(x=X_all), 
                        type = "response")
  
  mu1_DR1 <- (1/n0)*sum(w_IOW[1:n1]*(Y_all[1:n1] - Y1_hat_all[1:n1])) + 
    (1/n0)*sum(Y1_hat_all[(n1+1):n])
  mu0_DR1 <- mu0_naive
  ATC_DR1 <- g_function(mu1_DR1) - g_function(mu0_DR1)
  
  ATC_DR1_function <- function(data, indices){
    temp <- data[indices,]
    outcome_model <- glm(y ~ ., 
                         data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                           x=temp[temp[,"S_all"]==1, grep("X_all", colnames(temp))]), 
                         family = binomial(link = outcome_link))
    
    Y1_hat_all <- predict(outcome_model, newdata = data.frame(x=temp[, grep("X_all", colnames(temp))]), 
                          type = "response")
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[, grep("X_all", colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    mu1_DR1 <- (1/sum(temp[,"S_all"]==0))*sum(w_IOW[temp[,"S_all"]==1]*(temp[,"Y_all"][temp[,"S_all"]==1] - Y1_hat_all[temp[,"S_all"]==1])) + 
      (1/sum(temp[,"S_all"]==0))*sum(Y1_hat_all[temp[,"S_all"]==0])
    
    mu0_naive <- mean(temp[temp[,"S_all"]==0, "Y_all"])
    ATC_DR1 <- g_function(mu1_DR1) - g_function(mu0_naive)
    return(ATC_DR1)
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR1_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  
  ATC_DR1_CI <- c(boot_samples$t0[1] - abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE), 
                  boot_samples$t0[1] + abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE))
  
  ################################
  # DR2 (augmented inverse odds weighting estimator with normalized weights)
  # NOTE: This uses the outcome_link parameter
  
  mu1_DR2 <- (1/sum(w_IOW[1:n1]))*sum(w_IOW[1:n1]*(Y_all[1:n1] - Y1_hat_all[1:n1])) + 
    (1/n0)*sum(Y1_hat_all[(n1+1):n])
  mu0_DR2 <- mu0_naive
  ATC_DR2 <- g_function(mu1_DR2) - g_function(mu0_DR2)
  
  ATC_DR2_function <- function(data, indices){
    temp <- data[indices,]
    
    outcome_model <- glm(y ~ ., 
                         data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                           x=temp[temp[,"S_all"]==1, grep("X_all", colnames(temp))]), 
                         family = binomial(link = outcome_link))
    
    Y1_hat_all <- predict(outcome_model, newdata = data.frame(x=temp[, grep("X_all", colnames(temp))]), 
                          type = "response")
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[, grep("X_all", colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    mu1_DR2 <- (1/sum(w_IOW[temp[,"S_all"]==1]))*sum(w_IOW[temp[,"S_all"]==1]*(temp[,"Y_all"][temp[,"S_all"]==1] - Y1_hat_all[temp[,"S_all"]==1])) + 
      (1/sum(temp[,"S_all"]==0))*sum(Y1_hat_all[temp[,"S_all"]==0])
    
    mu0_naive <- mean(temp[temp[,"S_all"]==0, "Y_all"])
    ATC_DR2 <- g_function(mu1_DR2) - g_function(mu0_naive)
    return(ATC_DR2)
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR2_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  
  ATC_DR2_CI <- c(boot_samples$t0[1] - abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE), 
                  boot_samples$t0[1] + abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE))
  
  ################################
  # DR3: Augmented MAIC (entropy balancing) estimator
  # NOTE: This uses the outcome_link parameter
  
  mu1_DR3 <- (1/sum(wt_EB[1:n1]))*sum(wt_EB[1:n1]*(Y_all[1:n1] - Y1_hat_all[1:n1])) + 
    (1/n0)*sum(Y1_hat_all[(n1+1):n])
  mu0_DR3 <- mu0_naive
  ATC_DR3 <- g_function(mu1_DR3) - g_function(mu0_DR3)
  
  ATC_DR3_function <- function(data, indices){
    temp <- data[indices,]
    
    outcome_model <- glm(y ~ ., 
                         data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                           x=temp[temp[,"S_all"]==1, grep("X_all", colnames(temp))]), 
                         family = binomial(link = outcome_link))
    
    Y1_hat_all <- predict(outcome_model, newdata = data.frame(x=temp[, grep("X_all", colnames(temp))]), 
                          type = "response")
    
    AC.IPD <- data.frame(y=temp[,"Y_all"][temp[,"S_all"]==1], X= temp[temp[,"S_all"]==1, grep("X_all", colnames(temp))])
    BC.ALD <- data.frame(matrix(apply(temp[temp[,"S_all"]==0, grep("X_all", colnames(temp))], 2, mean), 1,))
    objfn <- function(a1, X){ sum(exp(X %*% a1)) }
    gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
    cov_names <- paste0("X.X_all.", 1:dim(X_all)[2])
    colnames(BC.ALD) <- paste0("mean.", cov_names)
    X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                    as.matrix((BC.ALD[, c(paste("mean.", cov_names, sep=""))])), '-')
    
    gamma <- optim(par = rep(0, dim(X.EM.0)[2]), 
                   fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
    wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
    
    mu1_DR3 <- (1/sum(wt_EB[temp[,"S_all"]==1]))*sum(wt_EB[temp[,"S_all"]==1]*(temp[,"Y_all"][temp[,"S_all"]==1] - Y1_hat_all[temp[,"S_all"]==1])) + 
      (1/sum(temp[,"S_all"]==0))*sum(Y1_hat_all[temp[,"S_all"]==0])
    
    mu0_naive <- mean(temp[temp[,"S_all"]==0, "Y_all"])
    ATC_DR3 <- g_function(mu1_DR3) - g_function(mu0_naive)
    return(ATC_DR3)
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR3_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  
  ATC_DR3_CI <- c(boot_samples$t0[1] - abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE), 
                  boot_samples$t0[1] + abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE))
  
  ################################
  # 3.5 Other augmented weighting estimators
  # DR4: Weighted G-computation with IOW weights
  # NOTE: This uses the outcome_link parameter
  
  outcome_model_w_IOW <- glm(y ~ ., 
                             data = data.frame(y=c(Y_all[S_all==1]), 
                                               x=(X_all[S_all==1,])), 
                             family = binomial(link = outcome_link),
                             weights=w_IOW[S_all==1])
  
  Y1_hat_w_IOW <- predict(outcome_model_w_IOW, newdata = data.frame(x=X_all[S_all==0,]), 
                          type = "response")
  
  # Note: the n0 here refers to the n0 from the parent environment and this is intentional
  # The stratified resampling (strata=S_all) implies that each bootstrap sample will have 
  # the same number of S=0 and S=1 observations as the original data.
  mu1_DR4 <- (1/n0)*sum(Y1_hat_w_IOW)
  mu0_DR4 <- mu0_naive
  ATC_DR4 <- g_function(mu1_DR4) - g_function(mu0_DR4)
  
  ATC_DR4_function <- function(data, indices){
    temp <- data[indices,]
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[, grep("X_all", colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    outcome_model_w_IOW <- glm(y ~ ., 
                               data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                                 x=temp[temp[,"S_all"]==1, grep("X_all", colnames(temp))]), 
                               family = binomial(link = outcome_link),
                               weights=w_IOW[temp[,"S_all"]==1])
    
    Y1_hat_w_IOW <- predict(outcome_model_w_IOW, newdata = data.frame(x=temp[temp[,"S_all"]==0, grep("X_all", colnames(temp))]), 
                            type = "response")
    
    mu1_DR4 <- (1/sum(temp[,"S_all"]==0))*sum(Y1_hat_w_IOW)
    mu0_naive <- mean(temp[temp[,"S_all"]==0, "Y_all"])
    ATC_DR4 <- g_function(mu1_DR4) - g_function(mu0_naive)
    return(ATC_DR4)
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR4_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  
  ATC_DR4_CI <- c(boot_samples$t0[1] - abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE), 
                  boot_samples$t0[1] + abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE))
  
  ################################
  # DR5: Weighted G-computation with MAIC weights
  # NOTE: This uses the outcome_link parameter
  
  outcome_model_wt_EB <- glm(y ~ ., 
                             data = data.frame(y=c(Y_all[S_all==1]), 
                                               x=(X_all[S_all==1,])), 
                             family = binomial(link = outcome_link),
                             weights=wt_EB)
  
  Y1_hat_wt_EB <- predict(outcome_model_wt_EB, newdata = data.frame(x=X_all[S_all==0,]), 
                          type = "response")
  
  mu1_DR5 <- (1/n0)*sum(Y1_hat_wt_EB)
  mu0_DR5 <- mu0_naive
  ATC_DR5 <- g_function(mu1_DR5) - g_function(mu0_DR5)
  
  ATC_DR5_function <- function(data, indices){
    temp <- data[indices,]
    
    AC.IPD <- data.frame(y=temp[,"Y_all"][temp[,"S_all"]==1], X= temp[temp[,"S_all"]==1, grep("X_all", colnames(temp))])
    BC.ALD <- data.frame(matrix(apply(temp[temp[,"S_all"]==0, grep("X_all", colnames(temp))], 2, mean), 1,))
    objfn <- function(a1, X){ sum(exp(X %*% a1)) }
    gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
    cov_names <- paste0("X.X_all.", 1:dim(X_all)[2])
    colnames(BC.ALD) <- paste0("mean.", cov_names)
    X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                    as.matrix((BC.ALD[, c(paste("mean.", cov_names, sep=""))])), '-')
    
    gamma <- optim(par = rep(0, dim(X.EM.0)[2]), 
                   fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
    wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
    
    outcome_model_wt_EB <- glm(y ~ ., 
                               data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                                 x=temp[temp[,"S_all"]==1, grep("X_all", colnames(temp))]), 
                               family = binomial(link = outcome_link),
                               weights=wt_EB[temp[,"S_all"]==1])
    
    Y1_hat_wt_EB <- predict(outcome_model_wt_EB, 
                            newdata = data.frame(x=temp[temp[,"S_all"]==0, grep("X_all", colnames(temp))]), 
                            type = "response")
    
    mu1_DR5 <- (1/sum(temp[,"S_all"]==0))*sum(Y1_hat_wt_EB)
    mu0_naive <- mean(temp[temp[,"S_all"]==0, "Y_all"])
    ATC_DR5 <- g_function(mu1_DR5) - g_function(mu0_naive)
    return(ATC_DR5)
  }
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR5_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  
  ATC_DR5_CI <- c(boot_samples$t0[1] - abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE), 
                  boot_samples$t0[1] + abs(qnorm(0.025))*sd(boot_samples$t[,1], na.rm=TRUE))
  
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  # Return results
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  all_ATC_estimates <- list(
    overlap=c(overlap, ESS),
    ATC_naive=c(ATC_naive, n0), 
    ATC_IOW1=c(ATC_IOW1, n1), 
    ATC_IOW2=c(ATC_IOW2, 1), 
    ATC_EB=c(ATC_EB, 1), 
    ATC_GCOMP=c(ATC_GCOMP, 1), 
    ATC_DR1=c(ATC_DR1, 1), 
    ATC_DR2=c(ATC_DR2, 1), 
    ATC_DR3=c(ATC_DR3, 1), 
    ATC_DR4=c(ATC_DR4, 1), 
    ATC_DR5=c(ATC_DR5, 1))
  
  all_ATC_CI <- list(
    ATC_naive_CI=ATC_naive_CI, 
    ATC_IOW1_CI=ATC_IOW1_CI, 
    ATC_IOW2_CI=ATC_IOW2_CI, 
    ATC_EB_CI=ATC_EB_CI, 
    ATC_GCOMP_CI=ATC_GCOMP_CI, 
    ATC_DR1_CI=ATC_DR1_CI, 
    ATC_DR2_CI=ATC_DR2_CI, 
    ATC_DR3_CI=ATC_DR3_CI, 
    ATC_DR4_CI=ATC_DR4_CI, 
    ATC_DR5_CI=ATC_DR5_CI)
  
  return(list(all_ATC_estimates=all_ATC_estimates, all_ATC_CI=all_ATC_CI))
}


################################################################################
# Simulation study function with link parameter
################################################################################

simstudy_Kang_Schafer <- function(scenario=1, IPD=TRUE, n=1000, n_boot=500, 
                                  findtrueATC=FALSE, outcome_link="logit"){
  
  X <- cbind(rnorm(n,0,1), rnorm(n,0,1), rnorm(n,0,1), rnorm(n,0,1))
  X1 <- X[,1]; X2 <- X[,2]; X3 <- X[,3]; X4 <- X[,4]
  Z1 <- scale(exp(X1/2))
  Z2 <- scale(X2^2)
  Z3 <- scale((X1*X3+0.6)^3)
  Z4 <- scale((X2+X4+20)^2)
  Z <- cbind(Z1, Z2, Z3, Z4)
  
  if(scenario %in% c(1,2)){prob_of_STUDY0 <- ilogit(-1*X1 + 0.5*X2 - 0.25*X3 - 0.5*X4)}
  if(scenario %in% c(3,4)){prob_of_STUDY0 <- ilogit(-1*Z1 + 0.5*Z2 - 0.25*Z3 - 0.5*Z4)}
  
  odds_STUDY0 <- prob_of_STUDY0/(1-prob_of_STUDY0)
  
  # warning message:
  if(quantile(prob_of_STUDY0, c(0.1)) < 0.01 | quantile(prob_of_STUDY0, c(0.9)) > 0.99){ 
    print("very low overlap")
  }
  
  STUDY0_id <- unlist(lapply(1:n, function(q){
    sample(x=c(0,1), size=1, prob=c(1-prob_of_STUDY0[q], prob_of_STUDY0[q]))
  }))
  STUDY1_id <- 1 - STUDY0_id
  S <- STUDY1_id
  
  # generating outcomes:
  beta0 <- 0
  delta <- 1.5
  delta2 <- -0.5
  beta1 <- c(1, -1.5, 0.5, -0.5)
  
  if(scenario %in% c(1,3)){prob_y <- ilogit(beta0 + c(beta1 %*% t(X)) + delta*S + delta2*S*X1)}
  if(scenario %in% c(2,4)){prob_y <- ilogit(beta0 + c(beta1 %*% t(Z)) + delta*S + delta2*S*Z1)}
  
  if(findtrueATC){
    S_all <- S
    n1 <- sum(S_all==1)
    n0 <- sum(S_all==0)
    n <- n1 + n0
    
    # calculating overlap of first covariate:
    overlap_plot <- list()
    for(j in 1:4){
      d1dens <- density(X[S_all==1, j], from = -10, to = 10)
      d2dens <- density(X[S_all==0, j], from = -10, to = 10)
      joint <- pmin(d1dens$y, d2dens$y)
      
      df2 <- data.frame(x = rep(d1dens$x, 3), 
                        y = c(d1dens$y, d2dens$y, joint),
                        Data = rep(c("S=1", "S=0", "The overlap"), each = length(d1dens$x)))
      
      overlap <- 2*sum(joint) / sum(d1dens$y, d2dens$y)
      print(overlap)
      
      overlap_plot[[j]] <- ggplot(df2, aes(x, y, fill = Data)) + 
        geom_area(position = position_identity(), color = "black") +
        scale_fill_brewer(palette = "Pastel2") +
        theme_bw() + xlim(-4, 4) +
        annotate("text", x=0, y=0.1, label=paste("Overlap =", round(overlap, 2))) +
        annotate("text", x=3, y=0.35, label=paste("X[", j, "]", sep=""), cex=6) +
        theme(legend.position="none") + ylab("") + xlab("")
    }
    print(c(mean(S_all==0), mean(S_all==1)))
    
    if(scenario %in% c(1,3)){
      p0 <- mean(ilogit(beta0 + c(beta1 %*% t(X[S==0,])) + delta*0 + delta2*0*X1[S==0]))
      p1 <- mean(ilogit(beta0 + c(beta1 %*% t(X[S==0,])) + delta*1 + delta2*1*X1[S==0]))
      true_OR <- log((p1/(1-p1))/(p0/(1-p0)))
    }
    
    if(scenario %in% c(2,4)){
      p0 <- mean(ilogit(beta0 + c(beta1 %*% t(Z[S==0,])) + delta*0 + delta2*0*Z1[S==0]))
      p1 <- mean(ilogit(beta0 + c(beta1 %*% t(Z[S==0,])) + delta*1 + delta2*1*Z1[S==0]))
      true_OR <- log((p1/(1-p1))/(p0/(1-p0)))
    }
    
    return(true_OR)
  }
  
  if(!findtrueATC){
    y_all <- unlist(lapply(1:n, function(q){
      sample(x=c(0,1), size=1, prob=c(1-prob_y[q], prob_y[q]))
    }))
    
    Y_all <- c(y_all[STUDY1_id==1], y_all[STUDY1_id==0])
    X_all <- rbind(X[STUDY1_id==1,], X[STUDY1_id==0,])
    S_all <- c(S[STUDY1_id==1], S[STUDY1_id==0])
    
    estimates <- ITC_estimates(Y_all, X_all, S_all, IPD_for_Study0=IPD, 
                               boot_n=n_boot, outcome_link=outcome_link)
    
    return(estimates)
  }
}



################################################################################
# Function to process simulation results
################################################################################
process_simulation_results <- function(simstudy_results, trueATC, method_names){
  
  # bias
  bias <- apply(simstudy_results[seq(1, 22, 2),] - trueATC, 1, function(z){ mean(z[is.finite(z)]) })
  
  # MC SE of bias
  MCSE_bias <- apply(simstudy_results[seq(1, 22, 2),], 1, function(z){ 
    sd(z[is.finite(z)]) / sqrt(length(z[is.finite(z)])) 
  })
  
  # std error
  ESE <- apply(simstudy_results[seq(1, 22, 2),], 1, function(z) sd(z[is.finite(z)]))
  
  # MC SE of ESE
  MCSE_ESE <- apply(simstudy_results[seq(1, 22, 2),], 1, function(z){ 
    z_finite <- z[is.finite(z)]
    sd(z_finite) / sqrt(2 * (length(z_finite) - 1))
  })
  
  # coverage
  CIcov <- apply((simstudy_results[-c(1:22),][seq(1, 20, 2),] < trueATC) & 
                   (simstudy_results[-c(1:22),][seq(2, 20, 2),] > trueATC), 1, 
                 function(z) mean(z[is.finite(z)]))
  
  # MC SE of coverage
  MCSE_CIcov <- apply((simstudy_results[-c(1:22),][seq(1, 20, 2),] < trueATC) & 
                        (simstudy_results[-c(1:22),][seq(2, 20, 2),] > trueATC), 1, 
                      function(z){
                        z_finite <- z[is.finite(z)]
                        p <- mean(z_finite)
                        sqrt(p * (1 - p) / length(z_finite))
                      })
  
  # CI width
  CIwidth <- apply(abs((simstudy_results[-c(1:22),][seq(2, 20, 2),]) - 
                         (simstudy_results[-c(1:22),][seq(1, 20, 2),])), 1, 
                   function(z) mean(z[is.finite(z)]))
  
  results <- cbind(bias[-1], ESE[-1], CIcov, CIwidth)
  rownames(results) <- method_names
  
  return(list(results = results, 
              MCSE_bias = MCSE_bias, 
              MCSE_ESE = MCSE_ESE, 
              MCSE_CIcov = MCSE_CIcov))
}
#################################################################


################################################################################
# Run simulation study
################################################################################

# Method names for output tables
method_names <- c(
  "1. The naive estimator", 
  "2. IOW", 
  "3. Normalized IOW",
  "4. MAIC",
  "5. G-computation",
  "6. DR augmented IOW" ,
  "7. DR augmented normalized IOW",
  "8. DR augmented MAIC",
  "9. Weighted G-computation (normalized IOW weights)",
  "10. Weighted G-computation (MAIC weights)"
)



################################################################################
# LOGIT LINK SIMULATIONS (original)
################################################################################

cat("\n\n========================================\n")
cat("LOGIT LINK SIMULATIONS\n")
cat("========================================\n\n")

# True ATC values (same for both links since data generation uses logit)
set.seed(123)
trueATC1 <- simstudy_Kang_Schafer(1, n=10000000, findtrueATC=TRUE)
round(trueATC1,3)
set.seed(123)
trueATC2 <- simstudy_Kang_Schafer(2, n=10000000, findtrueATC=TRUE)
set.seed(123)
trueATC3 <- simstudy_Kang_Schafer(3, n=10000000, findtrueATC=TRUE)
set.seed(123)
trueATC4 <- simstudy_Kang_Schafer(4, n=10000000, findtrueATC=TRUE)

cat("True ATC values:\n")
cat("Scenario 1:", round(trueATC1,3), "\n")
cat("Scenario 2:", round(trueATC2,3), "\n")
cat("Scenario 3:", round(trueATC3,3), "\n")
cat("Scenario 4:", round(trueATC4,3), "\n\n")

# Scenario 1 - Logit
cat("Running Scenario 1 with LOGIT link...\n")
set.seed(123)
simstudy1_logit <- pbapply(cbind(1:n_sims), 1, function(z) {
  unlist(simstudy_Kang_Schafer(scenario=1, n=n_sample, n_boot=n_boot, outcome_link="logit"))
})
KS1_logit <- process_simulation_results(simstudy1_logit, trueATC1, method_names)

round(max(KS1_logit$MCSE_bias),4)
round(max(KS1_logit$MCSE_ESE),4)
round(max(KS1_logit$MCSE_CIcov),4)

# Scenario 2 - Logit
cat("Running Scenario 2 with LOGIT link...\n")
set.seed(123)
simstudy2_logit <- pbapply(cbind(1:n_sims), 1, function(z) {
  unlist(simstudy_Kang_Schafer(scenario=2, n=n_sample, n_boot=n_boot, outcome_link="logit"))
})
KS2_logit <- process_simulation_results(simstudy2_logit, trueATC2, method_names)

round(max(KS2_logit$MCSE_bias),4)
round(max(KS2_logit$MCSE_ESE),4)
round(max(KS2_logit$MCSE_CIcov),4)

# Scenario 3 - Logit
cat("Running Scenario 3 with LOGIT link...\n")
set.seed(123)
simstudy3_logit <- pbapply(cbind(1:n_sims), 1, function(z) {
  unlist(simstudy_Kang_Schafer(scenario=3, n=n_sample, n_boot=n_boot, outcome_link="logit"))
})
KS3_logit <- process_simulation_results(simstudy3_logit, trueATC3, method_names)

round(KS3_logit$results,3)

round(max(KS3_logit$MCSE_bias),4)
round(max(KS3_logit$MCSE_ESE),4)
round(max(KS3_logit$MCSE_CIcov),4)

# Scenario 4 - Logit
cat("Running Scenario 4 with LOGIT link...\n")
set.seed(123)
simstudy4_logit <- pbapply(cbind(1:n_sims), 1, function(z) {
  unlist(simstudy_Kang_Schafer(scenario=4, n=n_sample, n_boot=n_boot, outcome_link="logit"))
})
KS4_logit <- process_simulation_results(simstudy4_logit, trueATC4, method_names)
round(max(KS4_logit$MCSE_bias),4)
round(max(KS4_logit$MCSE_ESE),4)
round(max(KS4_logit$MCSE_CIcov),4)

################################################################################
# cauchit LINK SIMULATIONS
################################################################################

cat("\n\n========================================\n")
cat("cauchit LINK SIMULATIONS\n")
cat("========================================\n\n")

# Scenario 1 - cauchit
cat("Running Scenario 1 with cauchit link...\n")
set.seed(123)
simstudy1_cauchit <- pbapply(cbind(1:n_sims), 1, function(z) {
  unlist(simstudy_Kang_Schafer(scenario=1, n=n_sample, n_boot=n_boot, outcome_link="cauchit"))
})
KS1_cauchit <- process_simulation_results(simstudy1_cauchit, trueATC1, method_names)
round(max(KS1_cauchit$MCSE_bias),4)
round(max(KS1_cauchit$MCSE_ESE),4)
round(max(KS1_cauchit$MCSE_CIcov),4)

# Scenario 2 - cauchit
cat("Running Scenario 2 with cauchit link...\n")
set.seed(123)
simstudy2_cauchit <- pbapply(cbind(1:n_sims), 1, function(z) {
  unlist(simstudy_Kang_Schafer(scenario=2, n=n_sample, n_boot=n_boot, outcome_link="cauchit"))
})
KS2_cauchit <- process_simulation_results(simstudy2_cauchit, trueATC2, method_names)
round(max(KS2_cauchit$MCSE_bias),4)
round(max(KS2_cauchit$MCSE_ESE),4)
round(max(KS2_cauchit$MCSE_CIcov),4)

# Scenario 3 - cauchit
cat("Running Scenario 3 with cauchit link...\n")
set.seed(123)
simstudy3_cauchit <- pbapply(cbind(1:n_sims), 1, function(z) {
  unlist(simstudy_Kang_Schafer(scenario=3, n=n_sample, n_boot=n_boot, outcome_link="cauchit"))
})
KS3_cauchit <- process_simulation_results(simstudy3_cauchit, trueATC3, method_names)
round(max(KS3_cauchit$MCSE_bias),4)
round(max(KS3_cauchit$MCSE_ESE),4)
round(max(KS3_cauchit$MCSE_CIcov),4)

# Scenario 4 - cauchit
cat("Running Scenario 4 with cauchit link...\n")
set.seed(123)
simstudy4_cauchit <- pbapply(cbind(1:n_sims), 1, function(z) {
  unlist(simstudy_Kang_Schafer(scenario=4, n=n_sample, n_boot=n_boot, outcome_link="cauchit"))
})
KS4_cauchit <- process_simulation_results(simstudy4_cauchit, trueATC4, method_names)
round(max(KS4_cauchit$MCSE_bias),4)
round(max(KS4_cauchit$MCSE_ESE),4)
round(max(KS4_cauchit$MCSE_CIcov),4)

################################################################################
# Print Results
################################################################################

library(xtable)

cat("\n\n========================================\n")
cat("RESULTS COMPARISON\n")
cat("========================================\n\n")

# Scenario 1
cat("\n--- SCENARIO 1: Both models correctly specified ---\n")
cat("\nLogit Link:\n")
print(round(KS1_logit$results, 3))
xtable(round(KS1_logit$results, 3), digits=3)
cat("\ncauchit Link:\n")
print(round(KS1_cauchit$results, 3))
xtable(round(KS1_cauchit$results, 3), digits=3)

# Scenario 2
cat("\n--- SCENARIO 2: Outcome model misspecified, PS model correct ---\n")
cat("\nLogit Link:\n")
print(round(KS2_logit$results, 3))
xtable(round(KS2_logit$results, 3), digits=3)
cat("\ncauchit Link:\n")
print(round(KS2_cauchit$results, 3))
xtable(round(KS2_cauchit$results, 3), digits=3)

# Scenario 3
cat("\n--- SCENARIO 3: Outcome model correct, PS model misspecified ---\n")
cat("\nLogit Link:\n")
print(round(KS3_logit$results, 3))
xtable(round(KS3_logit$results, 3), digits=3)
cat("\ncauchit Link:\n")
print(round(KS3_cauchit$results, 3))
xtable(round(KS3_cauchit$results, 3), digits=3)

# Scenario 4
cat("\n--- SCENARIO 4: Both models misspecified ---\n")
cat("\nLogit Link:\n")
print(round(KS4_logit$results, 3))
xtable(round(KS4_logit$results, 3), digits=3)
cat("\ncauchit Link:\n")
print(round(KS4_cauchit$results, 3))
xtable(round(KS4_cauchit$results, 3), digits=3)

