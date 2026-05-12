################################################################################
ilogit<- function(x){exp(x)/(1+exp(x))}
logit<-function(x){log(x/(1-x))}

library(pbapply)
library(mvtnorm)
library(boot)
library(ggplot2)
library(xtable)


ITC_estimates <- function(Y_all, X_all, S_all, IPD_for_Study0=TRUE, M=10000, boot_n=1000){
  
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
  overlap
  n1/n0
  
  # g-function is the log odds:
  g_function <- function(p){log(p/(1-p))}
  
  #############################################
  # Naive estimate
  mu1_naive <- mean(Y_all[S_all==1])
  mu0_naive <- mean(Y_all[S_all==0])
  ATC_naive =  g_function(mu1_naive) - g_function(mu0_naive)
  ATC_naive
  SE_g_mu0 <- summary(glm(Y_all[S_all==0]~1, family="binomial"))$coef[2]
  SE_g_mu1 <- summary(glm(Y_all[S_all==1]~1, family="binomial"))$coef[2]
  
  naive_function <- function(data, indices) {
    temp <- data[indices,]
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu1_naive <- mean(temp[temp[,"S_all"]==1,"Y_all"])
    ATC_naive <- g_function(mu1_naive) - g_function(mu0_naive)
    return(c(g_function(mu0_naive),g_function(mu1_naive),ATC_naive))}
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all), 
                       statistic=naive_function, R=boot_n,  
                       strata=S_all,
                       parallel = "multicore") 
  
  boot_samples$t0
  c(g_function(mu0_naive), g_function(mu1_naive), ATC_naive)
  
  SE_g_mu0_naive <- sd(boot_samples$t[,1],na.rm=TRUE)
  SE_g_mu0_naive
  SE_g_mu0
  
  SE_g_mu1_naive <- sd(boot_samples$t[,2],na.rm=TRUE)
  SE_g_mu1_naive
  SE_g_mu1
  
  SE_ATC <- sd(boot_samples$t[,3],na.rm=TRUE)
  SE_ATC
  SE_ATC2 <- sqrt(SE_g_mu0_naive^2  + SE_g_mu1_naive^2)
  SE_ATC2
  
  SE_ATC3 <- sqrt(SE_g_mu0^2  + SE_g_mu1^2)
  SE_ATC3
  
  
  c(ATC_naive-abs(qnorm(0.025))*SE_ATC, ATC_naive+abs(qnorm(0.025))*SE_ATC)
  
  ATC_naive_CI<- c(boot_samples$t0[3]-abs(qnorm(0.025))*sd(boot_samples$t[,3],na.rm=TRUE), 
                   boot_samples$t0[3]+abs(qnorm(0.025))*sd(boot_samples$t[,3],na.rm=TRUE))
  
  ATC_naive_CI
  ##########
  # naive_function <- function(data, indices) {
  #   temp <- data[indices,]
  #   ATC_naive =  g_function(mean(temp[temp[,"S_all"]==1, "Y_all"])) - 
  #     g_function(mean(temp[temp[,"S_all"]==0, "Y_all"]))
  #   return(c(ATC_naive,g_function(mean(temp[temp[,"S_all"]==0, "Y_all"]))))}
  # 
  # boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all), 
  #                      statistic=naive_function, R=boot_n,  
  #                      strata=S_all, parallel = "multicore") 
  # 
  # SE_ATC_0<-sd(boot_samples$t[,1])
  # SE_ATC_0
  # 
  # SE_g_mu0_0 <-sd(boot_samples$t[,2])
  # 
  # ATC_naive_CI_0 <- quantile(c(boot_samples$t),c(0.025,0.975))
  # ATC_naive_CI_0
  # 
  #############################################    
  if(!IPD_for_Study0){
    M = 100000
    X_0_mean <- apply(X_all[S_all==0,],2,mean)
    X_0_sd <- apply(X_all[S_all==0,],2,sd)
    X_all <- rbind(X_all[S_all==1,], rmvnorm(M, X_0_mean, diag(X_0_sd^2)))
    S_all <-  c(S_all[S_all==1], rep(0,M))
    n1 <- sum(S_all==1)
    n0 <- sum(S_all==0)
    n <- n1 + n0
  }
  
  #############################################
  # Inverse odds weighting estimate
  ps_mod <- glm(S_all~X_all, family="binomial")
  ps <- predict(ps_mod, type="response")
  w_IOW <- (1-ps)*S_all/ps + (1-S_all)
  
  # if X is discrete, these should be relatively close:
  #coef(ps_mod)[1]
  #log(mean(S_all[X_all==0])/mean(1-S_all[X_all==0]))
  # But....
  coef(ps_mod)[1]
  # does not equal this:
  log(mean(S_all==0)/mean(S_all==1))
  
  
  mu1_IOW1 <- (1/n0)*sum(w_IOW[1:n1]*Y_all[1:n1])
  mu0_IOW1 <- mu0_naive
  ATC_IOW1 =  g_function(mu1_IOW1) - g_function(mu0_IOW1)
  ATC_IOW1
  
  mu1_IOW2 <- (1/sum(w_IOW[1:n1]))*sum(w_IOW[1:n1]*Y_all[1:n1])
  mu0_IOW2 <- mu0_naive
  ATC_IOW2 =  g_function(mu1_IOW2) - g_function(mu0_IOW2)
  ATC_IOW2
  
  
  IOW1_function <- function(data, indices){
    temp <- data[indices,]
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[,grep("X_all",colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    mu1_IOW1 <- (1/sum(temp[,"S_all"]==0))*sum(w_IOW[temp[,"S_all"]==1]*temp[,"Y_all"][temp[,"S_all"]==1])
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    ATC_IOW1 =  g_function(mu1_IOW1) - g_function(mu0_naive)
    return(ATC_IOW1)
  }
  
  
  
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=IOW1_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  ATC_IOW1_CI<- c(boot_samples$t0[1]-abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE), 
                   boot_samples$t0[1]+abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE))
  

  
  IOW2_function <- function(data, indices){
    temp <- data[indices,]
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[,grep("X_all",colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    
    mu1_IOW2 <- (1/sum(w_IOW[temp[,"S_all"]==1]))*sum(w_IOW[temp[,"S_all"]==1]*temp[,"Y_all"][temp[,"S_all"]==1])
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu0_IOW2 <- mu0_naive
    ATC_IOW2 =  g_function(mu1_IOW2) - g_function(mu0_IOW2)
    return(ATC_IOW2)
  }
  
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=IOW2_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  ATC_IOW2_CI <- c(boot_samples$t0[1]-abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE), 
                   boot_samples$t0[1]+abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE))
  
  
  #print("ATC_IOW2_CI")
  #print(ATC_IOW2_CI)
  
  #############################################
  # 3.3 Entropy balancing (matching-adjusted indirect comparison)
  
  # ...equivalent to minimizing the objective function:
  objfn <- function(a1, X){ sum(exp(X %*% a1)) }
  gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
  cov_names <- paste0("X.", 1:dim(X_all)[2])
  AC.IPD <- data.frame(y=Y_all[S_all==1], X= X_all[S_all==1,])
  BC.ALD <- data.frame(matrix(apply(X_all[S_all==0,],2,mean),1,))
  colnames(BC.ALD)<- paste0("mean.",cov_names)
  X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                  as.matrix((BC.ALD[,c(paste("mean.",cov_names, sep=""))])), '-')
  
  gamma <- optim(par = rep(0,dim(X.EM.0)[2]), 
                 fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
  wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
  ESS <- ((sum(wt_EB))^2)/sum(wt_EB^2)
  #print(ESS)
  mu1_EB <- sum(wt_EB*Y_all[S_all==1])
  mu0_EB <- mu0_naive
  ATC_EB <- g_function(mu1_EB) - g_function(mu0_EB)
  ATC_EB
  
  
  ##
  ATC_EB_function <- function(data, indices){
    temp <- data[indices,]
    AC.IPD <- data.frame(y=temp[,"Y_all"][temp[,"S_all"]==1], X= temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))])
    BC.ALD <- data.frame(matrix(apply(temp[temp[,"S_all"]==0,grep("X_all",colnames(temp))],2,mean),1,))
    objfn <- function(a1, X){ sum(exp(X %*% a1)) }
    gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
    cov_names <- paste0("X.X_all.", 1:dim(X_all)[2])
    colnames(BC.ALD)<- paste0("mean.",cov_names)
    X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                    as.matrix((BC.ALD[,c(paste("mean.",cov_names, sep=""))])), '-')
    
    gamma <- optim(par = rep(0,dim(X.EM.0)[2]), 
                   fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
    wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
    
    mu1_EB <- sum(wt_EB*temp[,"Y_all"][temp[,"S_all"]==1])
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu0_EB <- mu0_naive
    ATC_EB <- g_function(mu1_EB) - g_function(mu0_EB)
    return(ATC_EB)
  }
  
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_EB_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 

  ATC_EB_CI <- c(boot_samples$t0[1]-abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE), 
                 boot_samples$t0[1]+abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE))
  
  #plot(wt_EB[S_all==1]/sum(wt_EB[S_all==1])~c(w_IOW[S_all==1]/sum(w_IOW[S_all==1]))); abline(0,1)
  
  #############################################
  # 3.4 Doubly robust augmented weighting estimators
  
  ################################
  # The G-computation estimator
  
  # The G-computation estimator for the ATC contrasts the average of potential counterfactual outcomes under the active intervention with the average of observed outcomes for the external control. 
  
  outcome_model <- glm(y ~ ., 
                       data = data.frame(y=c(Y_all[S_all==1]), 
                                         x=(X_all[S_all==1,])), 
                       family = "binomial")	
  
  Y1_hat <- (predict(outcome_model, newdata = data.frame(x=X_all[S_all==0,]), 
                     type = "response"))
  
  mu1_GCOMP <- (1/n0)*sum(Y1_hat)
  mu0_GCOMP1 <- mu0_naive
  ATC_GCOMP <-  g_function(mu1_GCOMP) - g_function(mu0_GCOMP1)
  ATC_GCOMP
  
  
  
  ###
  ATC_GCOMP_function <- function(data, indices){
    temp <- data[indices,]
    outcome_model <- glm(y ~ ., 
                         data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                           x=temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))]), 
                         family = "binomial")	
    
    
    Y1_hat <- (predict(outcome_model, newdata = data.frame(x=temp[temp[,"S_all"]==0,grep("X_all",colnames(temp))]), 
                       type = "response"))
    
    mu1_GCOMP <- (1/sum(temp[,"S_all"]==0))*sum(Y1_hat)
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu0_GCOMP1 <- mu0_naive
    ATC_GCOMP <-  g_function(mu1_GCOMP) - g_function(mu0_GCOMP1)
    return(ATC_GCOMP)}
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_GCOMP_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 

  ATC_GCOMP_CI <- c(boot_samples$t0[1]-abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE), 
                    boot_samples$t0[1]+abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE))
  
  
  ################################
  # DR1 (augmented inverse odds weighting estimator):
  # Using Inverse odds weights
  
  Y1_hat_all <- (predict(outcome_model, newdata = data.frame(x=X_all), 
                         type = "response"))
  
  mu1_DR1 <- (1/n0)*sum(w_IOW[1:n1]*(Y_all[1:n1]  - Y1_hat_all[1:n1])) + 
    (1/n0)*sum(Y1_hat_all[(n1+1):n])
  mu0_DR1 <- mu0_naive
  ATC_DR1 <-  g_function(mu1_DR1) - g_function(mu0_DR1)
  ATC_DR1
  
  
  ###
  ATC_DR1_function <- function(data, indices){
    temp <- data[indices,]
    outcome_model <- glm(y ~ ., 
                         data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                           x=temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))]), 
                         family = "binomial")	
    
    
    Y1_hat_all <- (predict(outcome_model, newdata = data.frame(x=temp[,grep("X_all",colnames(temp))]), 
                           type = "response"))
    
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[,grep("X_all",colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    
    mu1_DR1 <- (1/sum(temp[,"S_all"]==0))*sum(w_IOW[temp[,"S_all"]==1]*(temp[,"Y_all"][temp[,"S_all"]==1]  - Y1_hat_all[temp[,"S_all"]==1])) + 
      (1/sum(temp[,"S_all"]==0))*sum(Y1_hat_all[temp[,"S_all"]==0])
    
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu0_DR1 <- mu0_naive
    ATC_DR1 <-  g_function(mu1_DR1) - g_function(mu0_DR1)
    return(ATC_DR1)}
  
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR1_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
  
  ATC_DR1_CI <- c(boot_samples$t0[1]-abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE), 
                  boot_samples$t0[1]+abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE))
  
  
  ################################
  # DR2 (augmented inverse odds weighting estimator with Normalizing the weights so that they sum to one and to ensure bounded):
  
  mu1_DR2 <- (1/sum(w_IOW[1:n1]))*sum(w_IOW[1:n1]*(Y_all[1:n1]  - Y1_hat_all[1:n1])) + 
    (1/n0)*sum(Y1_hat_all[(n1+1):n])
  mu0_DR2 <- mu0_naive
  ATC_DR2 <-  g_function(mu1_DR2) - g_function(mu0_DR2)
  ATC_DR2
  
  ###
  ATC_DR2_function <- function(data, indices){
    temp <- data[indices,]
    
    outcome_model <- glm(y ~ ., 
                         data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                           x=temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))]), 
                         family = "binomial")	
    
    
    Y1_hat_all <- (predict(outcome_model, newdata = data.frame(x=temp[,grep("X_all",colnames(temp))]), 
                           type = "response"))
    
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[,grep("X_all",colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    mu1_DR2 <- (1/sum(w_IOW[temp[,"S_all"]==1]))*sum(w_IOW[temp[,"S_all"]==1]*(temp[,"Y_all"][temp[,"S_all"]==1]  - Y1_hat_all[temp[,"S_all"]==1])) + 
      (1/sum(temp[,"S_all"]==0))*sum(Y1_hat_all[temp[,"S_all"]==0])
    
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu0_DR2 <- mu0_naive
    ATC_DR2 <-  g_function(mu1_DR2) - g_function(mu0_DR2)
    return(ATC_DR2)}
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR2_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
 
  ATC_DR2_CI <- c(boot_samples$t0[1]-abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE), 
                  boot_samples$t0[1]+abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE))
  
   ################################
  # DR3
  # Our novel contribution is combining the entropy balancing-based MAIC approach
  
  mu1_DR3 <- (1/sum(wt_EB[1:n1]))*sum(wt_EB[1:n1]*(Y_all[1:n1]  - Y1_hat_all[1:n1])) + 
    (1/n0)*sum(Y1_hat_all[(n1+1):n])
  mu0_DR3 <- mu0_naive
  ATC_DR3 <-  g_function(mu1_DR3) - g_function(mu0_DR3)
  ATC_DR3
  
  
  ###
  ATC_DR3_function <- function(data, indices){
    temp <- data[indices,]
    
    outcome_model <- glm(y ~ ., 
                         data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                           x=temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))]), 
                         family = "binomial")	
    
    
    Y1_hat_all <- (predict(outcome_model, newdata = data.frame(x=temp[,grep("X_all",colnames(temp))]), 
                           type = "response"))
    
    AC.IPD <- data.frame(y=temp[,"Y_all"][temp[,"S_all"]==1], X= temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))])
    BC.ALD <- data.frame(matrix(apply(temp[temp[,"S_all"]==0,grep("X_all",colnames(temp))],2,mean),1,))
    objfn <- function(a1, X){ sum(exp(X %*% a1)) }
    gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
    cov_names <- paste0("X.X_all.", 1:dim(X_all)[2])
    colnames(BC.ALD)<- paste0("mean.",cov_names)
    X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                    as.matrix((BC.ALD[,c(paste("mean.",cov_names, sep=""))])), '-')
    
    gamma <- optim(par = rep(0,dim(X.EM.0)[2]), 
                   fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
    wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
    
    
    mu1_DR3 <- (1/sum(wt_EB[temp[,"S_all"]==1]))*sum(wt_EB[temp[,"S_all"]==1]*(temp[,"Y_all"][temp[,"S_all"]==1]  - Y1_hat_all[temp[,"S_all"]==1])) + 
      (1/sum(temp[,"S_all"]==0))*sum(Y1_hat_all[temp[,"S_all"]==0])
    
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu0_DR3 <- mu0_naive
    ATC_DR3 <-  g_function(mu1_DR3) - g_function(mu0_DR3)
    return(ATC_DR3)}
  
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR3_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 

  ATC_DR3_CI <- c(boot_samples$t0[1]-abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE), 
                  boot_samples$t0[1]+abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE))

  ################################
  # 3.5 Other augmented weighting estimators
  
  # One approach claimed to be doubly robust consists of 
  # G-computation based on the predictions of a weighted outcome model:
  
  outcome_model_w_IOW <- glm(y ~ ., 
                             data = data.frame(y=c(Y_all[S_all==1]), 
                                               x=(X_all[S_all==1,])), 
                             family = "binomial",
                             weights=w_IOW[S_all==1])	
  
  Y1_hat_w_IOW <- (predict(outcome_model_w_IOW, newdata = data.frame(x=X_all[S_all==0,]), 
                           type = "response"))
  
  mu1_DR4 <- (1/n0)*sum(Y1_hat_w_IOW)
  mu0_DR4 <- mu0_naive
  ATC_DR4 <-  g_function(mu1_DR4) - g_function(mu0_DR4)
  ATC_DR4
  
  ###
  ATC_DR4_function <- function(data, indices){
    temp <- data[indices,]
    
    ps_mod <- glm(temp[,"S_all"]~., data=temp[,grep("X_all",colnames(temp))], family="binomial")
    ps <- predict(ps_mod, type="response")
    w_IOW <- (1-ps)*temp[,"S_all"]/ps + (1-temp[,"S_all"])
    
    outcome_model_w_IOW <- glm(y ~ ., 
                               data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                                 x=temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))]), 
                               family = "binomial",
                               weights=w_IOW[temp[,"S_all"]==1])
    
    Y1_hat_w_IOW <- (predict(outcome_model_w_IOW, newdata = data.frame(x=temp[temp[,"S_all"]==0,grep("X_all",colnames(temp))]), 
                             type = "response"))
    
    mu1_DR4 <- (1/n0)*sum(Y1_hat_w_IOW)
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu0_DR4 <- mu0_naive
    ATC_DR4 <-  g_function(mu1_DR4) - g_function(mu0_DR4)
    return(ATC_DR4)}
  
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR4_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
 
  ATC_DR4_CI <- c(boot_samples$t0[1]-abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE), 
                  boot_samples$t0[1]+abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE))
  
  
  ################################
  # finally:
  
  outcome_model_wt_EB <- glm(y ~ ., 
                             data = data.frame(y=c(Y_all[S_all==1]), 
                                               x=(X_all[S_all==1,])), 
                             family = "binomial",
                             weights=wt_EB)	
  
  Y1_hat_wt_EB <- (predict(outcome_model_wt_EB, newdata = data.frame(x=X_all[S_all==0,]), 
                           type = "response"))
  
  mu1_DR5 <- (1/n0)*sum(Y1_hat_wt_EB)
  mu0_DR5 <- mu0_naive
  ATC_DR5 <-  g_function(mu1_DR5) - g_function(mu0_DR5)
  ATC_DR5
  
  ###
  ATC_DR5_function <- function(data, indices){
    temp <- data[indices,]
    
    AC.IPD <- data.frame(y=temp[,"Y_all"][temp[,"S_all"]==1], X= temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))])
    BC.ALD <- data.frame(matrix(apply(temp[temp[,"S_all"]==0,grep("X_all",colnames(temp))],2,mean),1,))
    objfn <- function(a1, X){ sum(exp(X %*% a1)) }
    gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
    cov_names <- paste0("X.X_all.", 1:dim(X_all)[2])
    colnames(BC.ALD)<- paste0("mean.",cov_names)
    X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                    as.matrix((BC.ALD[,c(paste("mean.",cov_names, sep=""))])), '-')
    
    gamma <- optim(par = rep(0,dim(X.EM.0)[2]), 
                   fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
    wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
    
    outcome_model_wt_EB <- glm(y ~ ., 
                               data = data.frame(y=c(temp[,"Y_all"][temp[,"S_all"]==1]), 
                                                 x=temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))]), 
                               family = "binomial",
                               weights=wt_EB[temp[,"S_all"]==1])
    
    Y1_hat_wt_EB <- (predict(outcome_model_wt_EB, 
                             newdata = data.frame(x=temp[temp[,"S_all"]==0,grep("X_all",colnames(temp))]), 
                             type = "response"))
    
    
    mu1_DR5 <- (1/n0)*sum(Y1_hat_wt_EB)
    mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
    mu0_DR5 <- mu0_naive
    ATC_DR5 <-  g_function(mu1_DR5) - g_function(mu0_DR5)
    return(ATC_DR5)}
  
  
  boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                       statistic=ATC_DR5_function, R=boot_n,  
                       strata=S_all, parallel = "multicore") 
 
  ATC_DR5_CI <- c(boot_samples$t0[1]-abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE), 
                  boot_samples$t0[1]+abs(qnorm(0.025))*sd(boot_samples$t[,1],na.rm=TRUE))
 
   #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
  
  all_ATC_estimates <- list(
    overlap=c(overlap,ESS),
    ATC_naive=c(ATC_naive,n0), 
    ATC_IOW1=c(ATC_IOW1,n1), 
    ATC_IOW2=c(ATC_IOW2,1), 
    ATC_EB=c(ATC_EB,1), 
    ATC_GCOMP=c(ATC_GCOMP,1), 
    ATC_DR1=c(ATC_DR1,1), 
    ATC_DR2=c(ATC_DR2,1), 
    ATC_DR3=c(ATC_DR3,1), 
    ATC_DR4=c(ATC_DR4,1), 
    ATC_DR5=c(ATC_DR5,1))
  
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
    
  
  
  
  return(list(all_ATC_estimates=all_ATC_estimates, all_ATC_CI=all_ATC_CI))}



################################################################################
################################################################################
simstudy_Kang_Schafer <- function(scenario=1, IPD=TRUE, n=1000, n_boot=500, findtrueATC=FALSE){
  
  X <- cbind(rnorm(n,0,1),rnorm(n,0,1),rnorm(n,0,1),rnorm(n,0,1))
  X1<-X[,1];X2<-X[,2];X3<-X[,3];X4<-X[,4];
  Z1 <- scale(exp(X1/2))
  Z2 <- scale(X2^2)

  Z3 <- scale((X1*X3+0.6)^3)
  Z4 <- scale((X2+X4+20)^2)
  Z <- cbind(Z1,Z2,Z3,Z4)
  
  if(scenario%in%c(1,2)){prob_of_STUDY0 <- (ilogit(-1*X1+0.5*X2-0.25*X3-0.5*X4))}
  if(scenario%in%c(3,4)){prob_of_STUDY0 <- (ilogit(-1*Z1+0.5*Z2-0.25*Z3-0.5*Z4))}
  
  odds_STUDY0 <- (prob_of_STUDY0/(1-prob_of_STUDY0))
  #hist(prob_of_STUDY0)
  # warning message:
  if( quantile(prob_of_STUDY0, c(0.1))<0.01 | quantile(prob_of_STUDY0, c(0.9))>0.99){ 
    print("very low overlap")
  }
  
  STUDY0_id <- unlist(lapply(1:n, 
                             function(q){
                               sample(x=c(0,1), size=1, 
                                      prob=c(1-prob_of_STUDY0[q], 
                                             prob_of_STUDY0[q]))}))
  STUDY1_id <- 1-STUDY0_id
  # S indicates study (S=1 for SAT and S=0 ECA)
  S <- STUDY1_id
  
  # generating outcomes:
  beta0 <- 0;
  delta <- 1.5
  delta2 <- -0.5
  beta1 <- c(1, -1.5, 0.5, -0.5)
  if(scenario%in%c(1,3)){prob_y <- ilogit(beta0 + c(beta1%*%t(X))  + delta*S + delta2*S*X1)}
  if(scenario%in%c(2,4)){prob_y <- ilogit(beta0 + c(beta1%*%t(Z))  + delta*S + delta2*S*Z1)}
  
  if(findtrueATC){
    
 
    S_all<- S
    n1 <- sum(S_all==1)
    n0 <- sum(S_all==0)
    n <- n1 + n0
    # calculating overlap of first covariate:
    overlap_plot<-list()
   for(j in 1:4){
    d1dens <- density(X[S_all==1,j], 
                      from = -10, 
                      to = 10)
    d2dens <- density(X[S_all==0,j], 
                      from = -10, 
                      to = 10)
    joint <- pmin(d1dens$y, d2dens$y)
    
    df2 <- data.frame(x = rep(d1dens$x, 3), 
                      y = c(d1dens$y, d2dens$y, joint),
                      Data = rep(c("S=1", "S=0", "The overlap"), each = length(d1dens$x)))
    
    

    
    overlap <- 2*sum(joint) / sum(d1dens$y, d2dens$y)
    print(overlap)
    
    overlap_plot[[j]] <- ggplot(df2, aes(x, y, fill = Data)) + 
      geom_area(position = position_identity(), color = "black") +
      scale_fill_brewer(palette = "Pastel2") +
      theme_bw()+xlim(-4,4)+
      annotate("text",x=0,y=0.1,label=paste("Overlap =",round(overlap,2)))+
      annotate("text",x=3,y=0.35,label=(paste("X[",j,"]",sep="")), cex=6)+
      theme(legend.position="none")+ylab("")+xlab("")
    
    }
    print(c(mean(S_all==0), mean(S_all==1)))
    n1/n0
    #plot_grid(plotlist = overlap_plot, nrow = 2)
    
    if(scenario%in%c(1,3)){
      p0 <- mean(ilogit(beta0 + c(beta1%*%t(X[S==0,]))  + delta*0 + delta2*0*X1[S==0]))
      p1 <- mean(ilogit(beta0 + c(beta1%*%t(X[S==0,]))  + delta*1 + delta2*1*X1[S==0]))
      true_OR <- log((p1/(1-p1))/(p0/(1-p0)))
    }
    
    if(scenario%in%c(2,4)){
      p0 <- mean(ilogit(beta0 + c(beta1%*%t(Z[S==0,]))  + delta*0 + delta2*0*Z1[S==0]))
      p1 <- mean(ilogit(beta0 + c(beta1%*%t(Z[S==0,]))  + delta*1 + delta2*1*Z1[S==0]))
      true_OR <- log((p1/(1-p1))/(p0/(1-p0)))
    }
    
    return(true_OR)
    
    }
  if(!findtrueATC){
    #hist(prob_y)
    y_all <- unlist(lapply(1:n, 
                           function(q){
                             sample(x=c(0,1), size=1, 
                                    prob=c(1-prob_y[q], prob_y[q]))}))
    
    Y_all <- c(y_all[STUDY1_id==1],y_all[STUDY1_id==0])
    X_all <- rbind(X[STUDY1_id==1,],X[STUDY1_id==0,])
    S_all <- c(S[STUDY1_id==1],S[STUDY1_id==0])
    
    estimates <- ITC_estimates(Y_all, X_all, S_all, IPD_for_Study0=IPD, boot_n=n_boot)
   # round(unlist(lapply(estimates$all_ATC_estimates,function(z)z[1])),2)
    
    return(estimates)}
  }




###########

# n = 1000

###########



#####
# Kang-Schafer Example Simulation study 1
set.seed(123)
trueATC1 <- simstudy_Kang_Schafer(1,n=1000000, findtrueATC=TRUE)
#
set.seed(123)
simstudy1 <- pbapply(cbind(1:10000), 1, 
                     function(z) unlist(simstudy_Kang_Schafer(scenario=1, n=1000, n_boot=100)))


# emperical SE of g(mu1) compared to mean of boot SE of g(mu1)
# cbind(apply(simstudy1[43:nrow(simstudy1),][seq(1,20,2),],1, function(z) sd(z,na.rm=TRUE)), rowMeans(simstudy1[63:nrow(simstudy1),][seq(1,20,2),],na.rm=TRUE))

# bias
bias1 <- apply(simstudy1[seq(1,22,2),]-trueATC1,1,function(z){ mean(z,na.rm=TRUE)} )

# MC SE of estimate
MCSE_bias1 <- apply(simstudy1[seq(1,22,2),]-trueATC1,1,function(z){ (1/(length(z)*(length(z)-1) ))*sum(z^2) })
max(MCSE_bias1, na.rm=TRUE)


# std error
ESE1 <- apply(simstudy1[seq(1,22,2),],1,function(z) sd(z, na.rm=TRUE))
# coverage
CIcov1 <- apply((simstudy1[-c(1:22),][seq(1,20,2),]<trueATC1)&(simstudy1[-c(1:22),][seq(2,20,2),]>trueATC1), 1, function(z) mean(z,na.rm=TRUE))
# CI width
CIwidth1 <- apply(abs((simstudy1[-c(1:22),][seq(2,20,2),])-(simstudy1[-c(1:22),][seq(1,20,2),])), 1, function(z) mean(z,na.rm=TRUE))

library(xtable)
KS1results <- cbind(bias1[-1],ESE1[-1],CIcov1,CIwidth1)
rownames(KS1results) <- c("1. The naive estimator", 
                          "2. IOW with weights from modeling", 
                          "3. IOW with normalized weights from modeling",
                          "4. MAIC",
                          "5. G-computation",
                          "6. DR with ``modeling'' IOW weights",
                          "7. DR with normalized ``modeling'' IOW weights",
                          "8. DR with MAIC weights",
                          "9. Augmented ``weighted G-computation'' with normalized ``modeling'' IOW weights",
                          "10. Augmented ``weighted G-computation'' with MAIC weights")
xtable(KS1results, digits=3)


# Kang-Schafer Example Simulation study 2
set.seed(123)
trueATC2 <- simstudy_Kang_Schafer(2,n=1000000, findtrueATC=TRUE)
set.seed(123)
simstudy2 <- pbapply(cbind(1:10000),1, function(z) unlist(simstudy_Kang_Schafer(scenario=2, n=1000, n_boot=100)))
# bias
bias2 <- apply(simstudy2[seq(1,22,2),]-trueATC2,1,function(z){ mean(z,na.rm=TRUE)} )
round(bias2,2)
# std error
ESE2 <- apply(simstudy2[seq(1,22,2),],1,function(z) sd(z, na.rm=TRUE))
# coverage
CIcov2 <- apply((simstudy2[-c(1:22),][seq(1,20,2),]<trueATC2)&(simstudy2[-c(1:22),][seq(2,20,2),]>trueATC2), 1, function(z) mean(z,na.rm=TRUE))
# CI width
CIwidth2 <- apply(abs((simstudy2[-c(1:22),][seq(2,20,2),])-(simstudy2[-c(1:22),][seq(1,20,2),])), 1, function(z) mean(z,na.rm=TRUE))


library(xtable)
KS2results <- cbind(bias2[-1],ESE2[-1],CIcov2,CIwidth2)
rownames(KS2results) <- c("1. The naive estimator", 
                          "2. IOW with weights from modeling", 
                          "3. IOW with normalized weights from modeling",
                          "4. MAIC",
                          "5. G-computation",
                          "6. DR with ``modeling'' IOW weights",
                          "7. DR with normalized ``modeling'' IOW weights",
                          "8. DR with MAIC weights",
                          "9. Augmented ``weighted G-computation'' with normalized ``modeling'' IOW weights",
                          "10. Augmented ``weighted G-computation'' with MAIC weights")
xtable(KS2results, digits=3)



# Kang-Schafer Example Simulation study 3
set.seed(123)
trueATC3 <- simstudy_Kang_Schafer(3,n=1000000, findtrueATC=TRUE)
set.seed(123)
simstudy3 <- pbapply(cbind(1:10000),1, function(z) unlist(simstudy_Kang_Schafer(scenario=3,n=1000, n_boot=100)))
# bias
bias3<-apply(simstudy3[seq(1,22,2),]-trueATC3,1,function(z){ mean(z,na.rm=TRUE)} )
round(bias3,2)
# std error
ESE3<-apply(simstudy3[seq(1,22,2),],1,function(z) sd(z, na.rm=TRUE))
# coverage
CIcov3<-apply((simstudy3[-c(1:22),][seq(1,20,2),]<trueATC3)&(simstudy3[-c(1:22),][seq(2,20,2),]>trueATC3), 1, function(z) mean(z,na.rm=TRUE))
# CI width
CIwidth3<-apply(abs((simstudy3[-c(1:22),][seq(2,20,2),])-(simstudy3[-c(1:22),][seq(1,20,2),])), 1, function(z) mean(z,na.rm=TRUE))

library(xtable)
KS3results <- cbind(bias3[-1],ESE3[-1],CIcov3,CIwidth3)
rownames(KS3results) <- c("1. The naive estimator", 
                          "2. IOW with weights from modeling", 
                          "3. IOW with normalized weights from modeling",
                          "4. MAIC",
                          "5. G-computation",
                          "6. DR with ``modeling'' IOW weights",
                          "7. DR with normalized ``modeling'' IOW weights",
                          "8. DR with MAIC weights",
                          "9. Augmented ``weighted G-computation'' with normalized ``modeling'' IOW weights",
                          "10. Augmented ``weighted G-computation'' with MAIC weights")
xtable(KS3results,digits=3)


# Kang-Schafer Example Simulation study 4
set.seed(123)
trueATC4 <- simstudy_Kang_Schafer(4,n=1000000, findtrueATC=TRUE)
set.seed(123)
simstudy4 <- pbapply(cbind(1:10000),1, function(z) unlist(simstudy_Kang_Schafer(scenario=4, n=1000, n_boot=100)))
# bias
bias4<-apply(simstudy4[seq(1,22,2),]-trueATC4,1,function(z){ mean(z,na.rm=TRUE)} )
# std error
ESE4<-apply(simstudy4[seq(1,22,2),],1,function(z) sd(z, na.rm=TRUE))
# coverage
CIcov4<-apply((simstudy4[-c(1:22),][seq(1,20,2),]<trueATC4)&(simstudy4[-c(1:22),][seq(2,20,2),]>trueATC4), 1, function(z) mean(z,na.rm=TRUE))
# CI width
CIwidth4<-apply(abs((simstudy4[-c(1:22),][seq(2,20,2),])-(simstudy4[-c(1:22),][seq(1,20,2),])), 1, function(z) mean(z,na.rm=TRUE))


library(xtable)
KS4results <- cbind(bias4[-1],ESE4[-1],CIcov4,CIwidth4)
rownames(KS4results) <- c("1. The naive estimator", 
                          "2. IOW with weights from modeling", 
                          "3. IOW with normalized weights from modeling",
                          "4. MAIC",
                          "5. G-computation",
                          "6. DR with ``modeling'' IOW weights",
                          "7. DR with normalized ``modeling'' IOW weights",
                          "8. DR with MAIC weights",
                          "9. Augmented ``weighted G-computation'' with normalized ``modeling'' IOW weights",
                          "10. Augmented ``weighted G-computation'' with MAIC weights")
xtable(KS4results,digits=3)





###########

# n = 200

###########



#####
# Kang-Schafer Example Simulation study 1
set.seed(123)
trueATC1 <- simstudy_Kang_Schafer(1,n=1000000, findtrueATC=TRUE)
#
set.seed(123)
simstudy1_200 <- pbapply(cbind(1:10000),1, function(z) unlist(simstudy_Kang_Schafer(scenario=1, n=200, n_boot=100)))
# bias
bias1_200 <- apply(simstudy1_200[seq(1,22,2),]-trueATC1,1,function(z){ mean(z,na.rm=TRUE)} )

# MC SE of estimate
MCSE_bias1_200 <- apply(simstudy1_200[seq(1,22,2),]-trueATC1,1,function(z){ (1/(length(z)*(length(z)-1) ))*sum(z^2) })
max(MCSE_bias1, na.rm=TRUE)



# std error
ESE1_200 <- apply(simstudy1_200[seq(1,22,2),],1,function(z){ sd(z,na.rm=TRUE)} )
MCSE_ESE1_200 <- apply(simstudy1_200[seq(1,22,2),],1,function(z){ (1/(length(z)*(length(z)-1) ))*(sum((z-trueATC1)^2-sd(z))) })
max(MCSE_ESE1_200, na.rm=TRUE)


# coverage
CIcov1_200 <- apply((simstudy1_200[-c(1:22),][seq(1,20,2),]<trueATC1)&(simstudy1_200[-c(1:22),][seq(2,20,2),]>trueATC1), 1, function(z) mean(z,na.rm=TRUE))
# CI width
CIwidth1_200 <- apply(abs((simstudy1_200[-c(1:22),][seq(2,20,2),])-(simstudy1_200[-c(1:22),][seq(1,20,2),])), 1, function(z) mean(z,na.rm=TRUE))

library(xtable)
KS1results_200 <- cbind(bias1_200[-1],ESE1_200[-1],CIcov1_200,CIwidth1_200)
rownames(KS1results_200) <- c("1. The naive estimator", 
                          "2. IOW with weights from modeling", 
                          "3. IOW with normalized weights from modeling",
                          "4. MAIC",
                          "5. G-computation",
                          "6. DR with ``modeling'' IOW weights",
                          "7. DR with normalized ``modeling'' IOW weights",
                          "8. DR with MAIC weights",
                          "9. Augmented ``weighted G-computation'' with normalized ``modeling'' IOW weights",
                          "10. Augmented ``weighted G-computation'' with MAIC weights")
xtable(KS1results_200, digits=3)


# Kang-Schafer Example Simulation study 2
set.seed(123)
trueATC2 <- simstudy_Kang_Schafer(2,n=1000000, findtrueATC=TRUE)
set.seed(123)
simstudy2_200 <- pbapply(cbind(1:10000),1, function(z){ print(z); unlist(simstudy_Kang_Schafer(scenario=2, n=200, n_boot=100))})
# bias
bias2_200 <- apply(simstudy2_200[seq(1,22,2),]-trueATC2,1,function(z){ mean(z,na.rm=TRUE)} )
round(bias2,2)

# MC SE of estimate
MCSE_bias2_200 <- apply(simstudy2_200[seq(1,22,2),]-trueATC2,1,function(z){ (1/(length(z)*(length(z)-1) ))*sum(z^2) })
max(MCSE_bias2, na.rm=TRUE)

# std error
ESE2_200 <- apply(simstudy2_200[seq(1,22,2),],1,function(z){ sd(z,na.rm=TRUE)} )
MCSE_ESE2_200 <- apply(simstudy2_200[seq(1,22,2),],1,function(z){ (1/(length(z)*(length(z)-1) ))*(sum((z-trueATC2)^2-sd(z))) })
max(MCSE_ESE2_200, na.rm=TRUE)

# coverage
CIcov2_200 <- apply((simstudy2_200[-c(1:22),][seq(1,20,2),]<trueATC2)&(simstudy2_200[-c(1:22),][seq(2,20,2),]>trueATC2), 1, function(z) mean(z,na.rm=TRUE))
# CI width
CIwidth2_200 <- apply(abs((simstudy2_200[-c(1:22),][seq(2,20,2),])-(simstudy2_200[-c(1:22),][seq(1,20,2),])), 1, function(z) mean(z,na.rm=TRUE))


library(xtable)
KS2results_200 <- cbind(bias2_200[-1],ESE2_200[-1],CIcov2_200,CIwidth2_200)
rownames(KS2results_200) <- c("1. The naive estimator", 
                          "2. IOW with weights from modeling", 
                          "3. IOW with normalized weights from modeling",
                          "4. MAIC",
                          "5. G-computation",
                          "6. DR with ``modeling'' IOW weights",
                          "7. DR with normalized ``modeling'' IOW weights",
                          "8. DR with MAIC weights",
                          "9. Augmented ``weighted G-computation'' with normalized ``modeling'' IOW weights",
                          "10. Augmented ``weighted G-computation'' with MAIC weights")
xtable(KS2results_200, digits=3)






# Kang-Schafer Example Simulation study 3
set.seed(123)
trueATC3 <- simstudy_Kang_Schafer(3,n=1000000, findtrueATC=TRUE)
set.seed(123)
simstudy3_200 <- pbapply(cbind(1:10000),1, function(z) unlist(simstudy_Kang_Schafer(scenario=3,n=200, n_boot=100)))
# bias
bias3_200<-apply(simstudy3_200[seq(1,22,2),]-trueATC3,1,function(z){ mean(z,na.rm=TRUE)} )
round(bias3_200,2)

# MC SE of estimate
MCSE_bias3_200 <- apply(simstudy3_200[seq(1,22,2),]-trueATC3,1,function(z){ (1/(length(z)*(length(z)-1) ))*sum(z^2) })
max(MCSE_bias3_200, na.rm=TRUE)

# std error
ESE3_200<-apply(simstudy3_200[seq(1,22,2),],1,function(z){ sd(z,na.rm=TRUE)} )
# coverage
CIcov3_200<-apply((simstudy3_200[-c(1:22),][seq(1,20,2),]<trueATC3)&(simstudy3_200[-c(1:22),][seq(2,20,2),]>trueATC3), 1, function(z) mean(z,na.rm=TRUE))
# CI width
CIwidth3_200<-apply(abs((simstudy3_200[-c(1:22),][seq(2,20,2),])-(simstudy3_200[-c(1:22),][seq(1,20,2),])), 1, function(z) mean(z,na.rm=TRUE))

library(xtable)
KS3results_200 <- cbind(bias3_200[-1],ESE3_200[-1],CIcov3_200,CIwidth3_200)
rownames(KS3results_200) <- c("1. The naive estimator", 
                          "2. IOW with weights from modeling", 
                          "3. IOW with normalized weights from modeling",
                          "4. MAIC",
                          "5. G-computation",
                          "6. DR with ``modeling'' IOW weights",
                          "7. DR with normalized ``modeling'' IOW weights",
                          "8. DR with MAIC weights",
                          "9. Augmented ``weighted G-computation'' with normalized ``modeling'' IOW weights",
                          "10. Augmented ``weighted G-computation'' with MAIC weights")
xtable(KS3results_200, digits=3)


# Kang-Schafer Example Simulation study 4
set.seed(123)
trueATC4 <- simstudy_Kang_Schafer(4,n=1000000, findtrueATC=TRUE)
set.seed(123)
simstudy4_200 <- pbapply(cbind(1:10000),1, function(z) unlist(simstudy_Kang_Schafer(scenario=4, n=200, n_boot=100)))
# bias
bias4_200<-apply(simstudy4_200[seq(1,22,2),]-trueATC4,1,function(z){ mean(z,na.rm=TRUE)} )

# MC SE of estimate
MCSE_bias4_200 <- apply(simstudy4_200[seq(1,22,2),]-trueATC4,1,function(z){ (1/(length(z)*(length(z)-1) ))*sum(z^2) })
max(MCSE_bias4_200, na.rm=TRUE)

# std error
ESE4_200<-apply(simstudy4_200[seq(1,22,2),],1,function(z){ sd(z,na.rm=TRUE)} )
# coverage
CIcov4_200<-apply((simstudy4_200[-c(1:22),][seq(1,20,2),]<trueATC4)&(simstudy4_200[-c(1:22),][seq(2,20,2),]>trueATC4), 1, function(z) mean(z,na.rm=TRUE))
# CI width
CIwidth4_200<-apply(abs((simstudy4_200[-c(1:22),][seq(2,20,2),])-(simstudy4_200[-c(1:22),][seq(1,20,2),])), 1, function(z) mean(z,na.rm=TRUE))


library(xtable)
KS4results_200 <- cbind(bias4_200[-1],ESE4_200[-1],CIcov4_200,CIwidth4_200)
rownames(KS4results_200) <- c("1. The naive estimator", 
                          "2. IOW with weights from modeling", 
                          "3. IOW with normalized weights from modeling",
                          "4. MAIC",
                          "5. G-computation",
                          "6. DR with ``modeling'' IOW weights",
                          "7. DR with normalized ``modeling'' IOW weights",
                          "8. DR with MAIC weights",
                          "9. Augmented ``weighted G-computation'' with normalized ``modeling'' IOW weights",
                          "10. Augmented ``weighted G-computation'' with MAIC weights")
xtable(KS4results_200,digits=3)




