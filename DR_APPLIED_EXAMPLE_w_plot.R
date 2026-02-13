library(dplyr)
library(boot)
library(survival)
library(MAIC)
library(ggplot2)
library(survminer)
library(flextable)
library(officer)

# set seed to give reproducible example
set.seed(1894)

boot_n<-10000

# g-function is the log-odds:
g_function <- function(p){log(p/(1-p))}

#### Intervention data

# Read in ADaM data and rename variables of interest

adsl <- read.csv(system.file("extdata", "adsl.csv", package = "MAIC", mustWork = TRUE))
adrs <- read.csv(system.file("extdata", "adrs.csv", package = "MAIC", mustWork = TRUE))
adtte <- read.csv(system.file("extdata", "adtte.csv", package = "MAIC", mustWork = TRUE))

adsl <- adsl %>% # Data containing the matching variables
  mutate(SEX=ifelse(SEX=="Male", 1, 0)) # Coded 1 for males and 0 for females

adrs <- adrs %>% # Response data
  filter(PARAM=="Response") %>%
  transmute(USUBJID, ARM, response=AVAL)

adtte <- adtte %>% # Time to event data (overall survival)
  filter(PARAMCD=="OS") %>%
  mutate(Event=1-CNSR) %>% #Set up coding as Event = 1, Censor = 0
  transmute(USUBJID, ARM, Time=AVAL, Event)

# Combine all intervention data
intervention_input <- adsl %>%
  full_join(adrs, by=c("USUBJID", "ARM")) %>%
  full_join(adtte, by=c("USUBJID", "ARM"))
head(intervention_input)

# List out matching covariates
match_cov <- c("AGE",
               "SEX",
               "SMOKE",
               "ECOG0")

# Baseline aggregate data for the comparator population
target_pop <- read.csv(system.file("extdata", "aggregate_data.csv",
                                   package = "MAIC", mustWork = TRUE))

# Renames target population cols to be consistent with match_cov
match_cov
names(target_pop)
target_pop_standard <- target_pop %>%
  #EDIT
  dplyr::rename(N=N,
                Treatment=ARM,
                AGE=age.mean,
                SEX=prop.male,
                SMOKE=prop.smoke,
                ECOG0=prop.ecog0
  ) %>%
  transmute(N, Treatment, AGE, SEX, SMOKE, ECOG0)

target_pop_standard

# Simulate response data based on the known proportion of responders
comparator_n <- target_pop$N # total number of patients in the comparator data
comparator_prop_events <- 0.4 # proportion of responders
# Calculate number with event
# Use round() to ensure we end up with a whole number of people
# number without an event = Total N - number with event to ensure we keep the same number of patients
n_with_event <- round(comparator_n*comparator_prop_events, digits = 0)
comparator_binary <- data.frame("response"= c(rep(1, n_with_event), rep(0, comparator_n - n_with_event)))

n0 <- dim(comparator_binary)[1]
n1 <- length(intervention_input$response)
Y_all <- unlist(c(intervention_input$response,comparator_binary))
S_all <- c(rep(1, n1), rep(0,n0))
X_all <- rbind(cbind(intervention_input%>%select(AGE, SEX, SMOKE, ECOG0)),
               cbind(AGE=rep(NA, n0), 
                     SEX=rep(NA, n0), 
                     SMOKE=rep(NA, n0), 
                     ECOG0=rep(NA, n0)))

X_all$AGE_SQ<-(X_all$AGE)^2

n1 <- sum(S_all==1)
n <- n1 + n0
n1/n0

#############################################
# naive estimate
mu1_naive <- mean(Y_all[S_all==1])
mu0_naive <- mean(Y_all[S_all==0])
ATC_naive =  g_function(mu1_naive) - g_function(mu0_naive)
ATC_naive
SE_g_mu1 <- sqrt(1/(length((Y_all[S_all==1]))*mu1_naive*(1-mu1_naive)))
SE_g_mu0 <- sqrt(1/(length((Y_all[S_all==0]))*mu0_naive*(1-mu0_naive)))

(sqrt((mu0_naive*(1-mu0_naive)/length((Y_all[S_all==0])))))

naive_function <- function(data, indices) {
  temp <- data[indices,]
  return( g_function(mean(temp[temp[,"S_all"]==1, "Y_all"])))}

boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all), 
                     statistic=naive_function, R=boot_n,  
                     strata=S_all,
                     parallel = "multicore") 

SE_ATC <- sqrt(sd(boot_samples$t,na.rm=TRUE)^2  + SE_g_mu0^2)
SE_ATC
ATC_naive_CI_boot<- c(ATC_naive-abs(qnorm(0.025))*SE_ATC, ATC_naive+abs(qnorm(0.025))*SE_ATC)
round(ATC_naive_CI_boot,3)

SE_ATC <- sqrt(SE_g_mu1^2  + SE_g_mu0^2)
ATC_naive_CI<- c(ATC_naive-abs(qnorm(0.025))*SE_ATC, ATC_naive+abs(qnorm(0.025))*SE_ATC)
round(ATC_naive_CI,3)

round((c(ATC_naive, ATC_naive_CI)),3)

round(exp(c(ATC_naive, ATC_naive_CI)),3)
# 5.318 3.880 7.289
# compare with published result:
# 5.318 (3.888 to 7.275)  https://roche.github.io/MAIC/articles/MAIC.html

#############################################
# 3.3 Entropy balancing (matching-adjusted indirect comparison)

# ...equivalent to minimizing the objective function:
objfn <- function(a1, X){ sum(exp(X %*% a1)) }
gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
cov_names <- paste0("X.", colnames(X_all))
AC.IPD <- data.frame(y=Y_all[S_all%in%c(1)], X= X_all[S_all==1,])

BC.ALD <- data.frame(cbind(target_pop$age.mean, 
                           target_pop$prop.male,
                           target_pop$prop.smoke,
                           target_pop$prop.ecog0,
                           target_pop$age.mean^2 + target_pop$age.sd^2))

colnames(BC.ALD)<- paste0("mean.",cov_names)

library(maicChecks)
maicLP(AC.IPD[,-1], BC.ALD) # Checks if AD is within the convex hull of IPD using lp-solve

X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                as.matrix((BC.ALD[,c(paste("mean.",cov_names, sep=""))])), '-')

gamma <- optim(par = rep(0,dim(X.EM.0)[2]), 
               fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))

mu1_EB <- sum(wt_EB*Y_all[S_all==1])
mu0_EB <- mu0_naive
ATC_EB <- g_function(mu1_EB) - g_function(mu0_EB)
ATC_EB

ATC_EB_function <- function(data, indices){
  temp <- data[indices,]
  AC.IPD <- data.frame(y=temp[,"Y_all"][temp[,"S_all"]==1], 
                       X= temp[temp[,"S_all"]==1,grep("X_all",colnames(temp))])
  BC.ALD <- data.frame(cbind(target_pop$age.mean, 
                             target_pop$prop.male,
                             target_pop$prop.smoke,
                             target_pop$prop.ecog0,
                             target_pop$age.mean^2 + target_pop$age.sd^2))
  objfn <- function(a1, X){ sum(exp(X %*% a1)) }
  gradfn <- function(a1, X){ colSums(sweep(X, 1, exp(X %*% a1), "*")) }
  cov_names <- paste0("X.X_all.", colnames(X_all))
  colnames(BC.ALD)<- paste0("mean.",cov_names)
  X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                  as.matrix((BC.ALD[,c(paste("mean.",cov_names, sep=""))])), '-')
  
  gamma <- optim(par = rep(0,dim(X.EM.0)[2]), 
                 fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
  wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
  
  mu1_EB <- sum(wt_EB*temp[,"Y_all"][temp[,"S_all"]==1])
  mu0_naive <- mean(temp[temp[,"S_all"]==0,"Y_all"])
  return(g_function(mu1_EB))
}

set.seed(123)
boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                     statistic=ATC_EB_function, R=boot_n,  
                     strata=S_all, parallel = "multicore") 
sd(boot_samples$t,na.rm=TRUE)
SE_ATC <- sqrt(sd(boot_samples$t,na.rm=TRUE)^2  + SE_g_mu0^2)
SE_ATC
ATC_EB_CI <- c(ATC_EB-abs(qnorm(0.025))*SE_ATC,ATC_EB+abs(qnorm(0.025))*SE_ATC)

round((c(ATC_EB, ATC_EB_CI)),3)

round(exp(c(ATC_EB, ATC_EB_CI)),3)
# 3.787 2.497 5.742
# compare with published result:
# 3.787 (2.558 to 5.605)  https://roche.github.io/MAIC/articles/MAIC.html

# Store MAIC SE for forest plot
SE_ATC_EB <- SE_ATC

#############################################    
# Simulate M individual values from target population
#############################################    
library(multinma)
M = 10000
set.seed(123)
out2<- add_integration(data.frame(Y_all=NA), 
                       AGE = distr(qnorm, mean=target_pop$age.mean, sd=target_pop$age.sd), 
                       SEX = distr(qbern, prob=target_pop$prop.male),
                       SMOKE = distr(qbern, prob=target_pop$prop.smoke),
                       ECOG0 = distr(qbern, prob=target_pop$prop.ecog0),
                       cor =  cor(X_all[S_all==1,c("AGE", "SEX", "SMOKE", "ECOG0")]),
                       n_int = M)

x_star <-  cbind(unlist(out2$.int_AGE), 
                 unlist(out2$.int_SEX),
                 unlist(out2$.int_SMOKE),
                 unlist(out2$.int_ECOG0))
dim((X_all[S_all==1,]))
dim((X_all[S_all==0,]))
dim(x_star)
# add squared age values:
x_star <- cbind(x_star,(x_star[,1])^2)
dim(x_star)
#add names
colnames(x_star)<-colnames(X_all[S_all==1,])

n_with_event <- round(M*comparator_prop_events, digits = 0)
Y_all <-  c(Y_all[S_all==1], c(rep(1, n_with_event), rep(0, M - n_with_event))) 

X_all <- rbind(X_all[S_all==1,], x_star)
S_all <-  c(S_all[S_all==1], rep(0,M))

dim((X_all[S_all==1,]))
dim((X_all[S_all==0,]))
n1 <- sum(S_all==1)
n0 <- sum(S_all==0)
n <- n1 + n0

################################
# The G-computation estimator

# The G-computation estimator for the ATC contrasts the average of potential counterfactual outcomes
# under the active intervention with the average of observed outcomes for the external control. 

outcome_model <- glm(y ~ ., 
                     data = data.frame(y=c(Y_all[S_all==1]), 
                                       x=(X_all[S_all==1,])), 
                     family = binomial(link="logit"))	

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
  return(g_function(mu1_GCOMP))}
set.seed(123)
boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                     statistic=ATC_GCOMP_function, R=boot_n,  
                     strata=S_all, parallel = "multicore") 
sd(boot_samples$t,na.rm=TRUE)
SE_ATC <- sqrt(sd(boot_samples$t,na.rm=TRUE)^2  + SE_g_mu0^2)
ATC_GCOMP_CI <- c(ATC_GCOMP-abs(qnorm(0.025))*SE_ATC,ATC_GCOMP+abs(qnorm(0.025))*SE_ATC)
round((c(ATC_GCOMP, ATC_GCOMP_CI)),3)
round(exp(c(ATC_GCOMP, ATC_GCOMP_CI)),3)

# Store G-comp SE for forest plot
SE_ATC_GCOMP <- SE_ATC

################################
# DR3
# Our contribution is augmenting the entropy balancing-based MAIC approach with an outcome model
data_for_outcome_model <- data.frame(y=c(Y_all[S_all==1]), 
                                     X_all[S_all==1,])
colnames(data_for_outcome_model)<-c("y", colnames(X_all))
outcome_model <- glm(y ~ ., 
                     data = data_for_outcome_model, 
                     family = "binomial")	
Y1_hat_all <- (predict(outcome_model, newdata = data.frame(X_all), 
                       type = "response"))

mu1_DR3 <- (1/sum(wt_EB))*sum(wt_EB*(Y_all[S_all==1]  - Y1_hat_all[S_all==1])) + 
  (1/n0)*sum(Y1_hat_all[S_all==0])
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
  cov_names <- paste0("X.X_all.", colnames(X_all))
  colnames(BC.ALD)<- paste0("mean.",cov_names)
  X.EM.0 <- sweep(cbind(as.matrix((AC.IPD[, cov_names]))), 2, 
                  as.matrix((BC.ALD[,c(paste("mean.",cov_names, sep=""))])), '-')
  
  gamma <- optim(par = rep(0,dim(X.EM.0)[2]), 
                 fn = objfn, gr = gradfn, X = X.EM.0, method = "BFGS")$par
  wt_EB <- exp(X.EM.0 %*% gamma)/sum(exp(X.EM.0 %*% gamma))
  
  mu1_DR3 <- (1/sum(wt_EB))*sum(wt_EB*(temp[,"Y_all"][temp[,"S_all"]==1]  -
                                         Y1_hat_all[temp[,"S_all"]==1])) + 
    (1/sum(temp[,"S_all"]==0))*sum(Y1_hat_all[temp[,"S_all"]==0])    
  return( g_function(mu1_DR3))}

set.seed(123)
boot_samples <- boot(data=data.frame(Y_all=Y_all, S_all=S_all, X_all=X_all), 
                     statistic=ATC_DR3_function, R=boot_n,  
                     strata=S_all, parallel = "multicore") 
sd(boot_samples$t,na.rm=TRUE)
SE_ATC <- sqrt(sd(boot_samples$t,na.rm=TRUE)^2  + SE_g_mu0^2)
SE_ATC
ATC_DR3_CI <- c(ATC_DR3-abs(qnorm(0.025))*SE_ATC,ATC_DR3+abs(qnorm(0.025))*SE_ATC)
ATC_DR3_CI

round((c(ATC_DR3, ATC_DR3_CI)),3)

round(exp(c(ATC_DR3, ATC_DR3_CI)),3)

# Store DR SE for forest plot
SE_ATC_DR3 <- SE_ATC

#############################################
# Normalized Inverse Odds Weighting (Hajek type)
#############################################

# Fit propensity score model: Pr(S=1 | X) on concatenated data
# (intervention IPD + M simulated external control profiles)
ps_data <- data.frame(
  S      = S_all,
  AGE    = X_all[, "AGE"],
  SEX    = X_all[, "SEX"],
  SMOKE  = X_all[, "SMOKE"],
  ECOG0  = X_all[, "ECOG0"],
  AGE_SQ = X_all[, "AGE_SQ"]
)

ps_model <- glm(S ~ AGE + SEX + SMOKE + ECOG0 + AGE_SQ,
                data = ps_data, family = binomial(link = "logit"))
summary(ps_model)

# Predicted propensity scores for intervention subjects
e_hat_all <- predict(ps_model, type = "response")
e_hat_intervention <- e_hat_all[S_all == 1]

# Inverse odds weights: w_i = (1 - e_i) / e_i
iow_raw <- (1 - e_hat_intervention) / e_hat_intervention

# Normalized weights (sum to 1)
iow_normalized <- iow_raw / sum(iow_raw)

cat("\n--- IOW Weight Summary ---\n")
cat("Min:    ", round(min(iow_raw), 4), "\n")
cat("Median: ", round(median(iow_raw), 4), "\n")
cat("Mean:   ", round(mean(iow_raw), 4), "\n")
cat("Max:    ", round(max(iow_raw), 4), "\n")
cat("Sum:    ", round(sum(iow_raw), 2), "\n")

# Effective sample size
ESS_IOW <- (sum(iow_raw))^2 / sum(iow_raw^2)
cat("Effective sample size (IOW): ", round(ESS_IOW, 2), "\n")

# Histogram of normalized IOW weights is included in the combined plot above

# --- Covariate balance diagnostics ---

X_intervention <- intervention_input %>% select(AGE, SEX, SMOKE, ECOG0)

# Unweighted means in intervention SAT
means_intervention_unweighted <- colMeans(X_intervention)

# IOW-weighted means in intervention SAT
means_intervention_IOW <- colSums(iow_normalized * as.matrix(X_intervention))

# IOW-weighted SD of age
weighted_var_age_IOW <- sum(iow_normalized * (X_intervention$AGE - means_intervention_IOW["AGE"])^2)
weighted_sd_age_IOW <- sqrt(weighted_var_age_IOW)

# Target population means
means_target <- c(
  AGE   = target_pop$age.mean,
  SEX   = target_pop$prop.male,
  SMOKE = target_pop$prop.smoke,
  ECOG0 = target_pop$prop.ecog0
)



## Combined histogram of MAIC (EB) and normalized IOW weights
weights_df <- data.frame(
  Weight = c(as.numeric(wt_EB), as.numeric(iow_normalized)),
  Method = factor(
    c(rep("MAIC (entropy balancing)", length(wt_EB)),
      rep("Normalized IOW", length(iow_normalized))),
    levels = c("Normalized IOW", "MAIC (entropy balancing)")
  )
)

p_hist <- ggplot(weights_df, aes(x = Weight)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "black", linewidth = 0.3) +
  facet_wrap(~ Method, scales = "free", ncol = 2) +
  labs(x = "Weight", y = "Frequency") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

print(p_hist)

cat("\n====================================================================\n")
cat("            COVARIATE BALANCE DIAGNOSTICS (IOW)\n")
cat("====================================================================\n\n")

balance_table <- data.frame(
  Covariate = c("Age (mean)", "Age (SD)", "Sex (prop male)",
                "ECOG0 (prop)", "Smoking (prop)"),
  Intervention_Unweighted = c(
    round(means_intervention_unweighted["AGE"], 2),
    round(sd(X_intervention$AGE), 2),
    round(means_intervention_unweighted["SEX"], 2),
    round(means_intervention_unweighted["ECOG0"], 2),
    round(means_intervention_unweighted["SMOKE"], 2)
  ),
  External_Control_Target = c(
    round(target_pop$age.mean, 2),
    round(target_pop$age.sd, 2),
    round(target_pop$prop.male, 2),
    round(target_pop$prop.ecog0, 2),
    round(target_pop$prop.smoke, 2)
  ),
  Intervention_IOW_Weighted = c(
    round(means_intervention_IOW["AGE"], 2),
    round(weighted_sd_age_IOW, 2),
    round(means_intervention_IOW["SEX"], 2),
    round(means_intervention_IOW["ECOG0"], 2),
    round(means_intervention_IOW["SMOKE"], 2)
  )
)
print(balance_table, row.names = FALSE)

# Standardized mean differences before and after IOW weighting
smd_before <- c(
  AGE   = (means_intervention_unweighted["AGE"] - means_target["AGE"]) /
    sqrt((sd(X_intervention$AGE)^2 + target_pop$age.sd^2) / 2),
  SEX   = (means_intervention_unweighted["SEX"] - means_target["SEX"]) /
    sqrt(means_target["SEX"] * (1 - means_target["SEX"])),
  SMOKE = (means_intervention_unweighted["SMOKE"] - means_target["SMOKE"]) /
    sqrt(means_target["SMOKE"] * (1 - means_target["SMOKE"])),
  ECOG0 = (means_intervention_unweighted["ECOG0"] - means_target["ECOG0"]) /
    sqrt(means_target["ECOG0"] * (1 - means_target["ECOG0"]))
)

smd_after <- c(
  AGE   = (means_intervention_IOW["AGE"] - means_target["AGE"]) /
    sqrt((sd(X_intervention$AGE)^2 + target_pop$age.sd^2) / 2),
  SEX   = (means_intervention_IOW["SEX"] - means_target["SEX"]) /
    sqrt(means_target["SEX"] * (1 - means_target["SEX"])),
  SMOKE = (means_intervention_IOW["SMOKE"] - means_target["SMOKE"]) /
    sqrt(means_target["SMOKE"] * (1 - means_target["SMOKE"])),
  ECOG0 = (means_intervention_IOW["ECOG0"] - means_target["ECOG0"]) /
    sqrt(means_target["ECOG0"] * (1 - means_target["ECOG0"]))
)

cat("\n--- Standardized Mean Differences ---\n")
smd_table <- data.frame(
  Covariate = c("Age", "Sex", "Smoking", "ECOG0"),
  SMD_Before = round(unname(smd_before), 4),
  SMD_After_IOW = round(unname(smd_after), 4)
)
print(smd_table, row.names = FALSE)

# --- ATC estimation: Normalized IOW (Hajek type) ---

mu1_IOW_Hajek <- sum(iow_normalized * intervention_input$response)
ATC_IOW_Hajek <- g_function(mu1_IOW_Hajek) - g_function(mu0_naive)
cat("\nATC_IOW_Hajek (log-OR):", round(ATC_IOW_Hajek, 3), "\n")

# Bootstrap for SE
ATC_IOW_Hajek_function <- function(data, indices) {
  temp <- data[indices, ]
  
  ps_model_b <- glm(S ~ AGE + SEX + SMOKE + ECOG0 + AGE_SQ,
                    data = temp, family = binomial(link = "logit"))
  
  e_hat_b <- predict(ps_model_b, type = "response")
  e_hat_int_b <- e_hat_b[temp$S == 1]
  iow_raw_b <- (1 - e_hat_int_b) / e_hat_int_b
  iow_norm_b <- iow_raw_b / sum(iow_raw_b)
  
  Y_int_b <- temp$Y[temp$S == 1]
  mu1_b <- sum(iow_norm_b * Y_int_b)
  
  # Bound to avoid log-odds issues at boundaries
  mu1_b <- max(min(mu1_b, 0.9999), 0.0001)
  
  return(g_function(mu1_b))
}

set.seed(123)
boot_data_IOW <- data.frame(
  Y      = Y_all,
  S      = S_all,
  AGE    = X_all[, "AGE"],
  SEX    = X_all[, "SEX"],
  SMOKE  = X_all[, "SMOKE"],
  ECOG0  = X_all[, "ECOG0"],
  AGE_SQ = X_all[, "AGE_SQ"]
)

boot_samples_IOW <- boot(data = boot_data_IOW,
                         statistic = ATC_IOW_Hajek_function,
                         R = boot_n,
                         strata = S_all,
                         parallel = "multicore")

SE_mu1_IOW <- sd(boot_samples_IOW$t, na.rm = TRUE)
SE_ATC_IOW <- sqrt(SE_mu1_IOW^2 + SE_g_mu0^2)
ATC_IOW_Hajek_CI <- c(ATC_IOW_Hajek - qnorm(0.975) * SE_ATC_IOW,
                      ATC_IOW_Hajek + qnorm(0.975) * SE_ATC_IOW)

cat("SE(g(mu_0^1)):      ", round(SE_mu1_IOW, 3), "\n")
cat("SE(ATC_IOW_Hajek):  ", round(SE_ATC_IOW, 3), "\n")

round(c(ATC_IOW_Hajek, ATC_IOW_Hajek_CI), 3)
round(exp(c(ATC_IOW_Hajek, ATC_IOW_Hajek_CI)), 3)

#############################################
# Forest plot comparing all estimators
#############################################

# Assemble results into a data frame
# Use the analytical SE for naive CI (consistent with the LaTeX writeup)
SE_ATC_naive <- sqrt(SE_g_mu1^2 + SE_g_mu0^2)

forest_data <- data.frame(
  Estimator = factor(
    c("Naive", "MAIC (EB)", "IOW (Hajek)", "G-computation", "DR augmented MAIC"),
    levels = c("DR augmented MAIC", "G-computation", "IOW (Hajek)", "MAIC (EB)", "Naive")
  ),
  Estimate  = c(ATC_naive, ATC_EB, ATC_IOW_Hajek, ATC_GCOMP, ATC_DR3),
  SE        = c(SE_ATC_naive, SE_ATC_EB, SE_ATC_IOW, SE_ATC_GCOMP, SE_ATC_DR3),
  CI_lower  = c(ATC_naive_CI[1], ATC_EB_CI[1], ATC_IOW_Hajek_CI[1], ATC_GCOMP_CI[1], ATC_DR3_CI[1]),
  CI_upper  = c(ATC_naive_CI[2], ATC_EB_CI[2], ATC_IOW_Hajek_CI[2], ATC_GCOMP_CI[2], ATC_DR3_CI[2])
)

# Labels for the right-hand side annotation
forest_data$label <- paste0(
  sprintf("%.3f", forest_data$Estimate),
  " (",
  sprintf("%.3f", forest_data$CI_lower),
  ", ",
  sprintf("%.3f", forest_data$CI_upper),
  ")"
)

# --- Log-odds ratio scale ---
p_logOR <- ggplot(forest_data, aes(x = Estimate, y = Estimator)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.2, linewidth = 0.6) +
  geom_point(size = 4, shape = 20) +
  geom_text(aes(x = max(forest_data$CI_upper) + 0.15, label = label),
            hjust = 0, size = 3.2) +
  labs(x = "ATC (marginal log-odds ratio)",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 11)
  ) +
  coord_cartesian(xlim = c(min(forest_data$CI_lower) - 0.1,
                           max(forest_data$CI_upper) + 1.2))

print(p_logOR)



# --- Print summary table ---
cat("\n====================================================================\n")
cat("               SUMMARY OF ALL ATC ESTIMATES\n")
cat("====================================================================\n\n")

summary_table <- data.frame(
  Estimator = c("Naive", "MAIC (EB)", "IOW (Hajek)", "G-computation", "DR augmented MAIC"),
  logOR     = round(c(ATC_naive, ATC_EB, ATC_IOW_Hajek, ATC_GCOMP, ATC_DR3), 3),
  SE        = round(c(SE_ATC_naive, SE_ATC_EB, SE_ATC_IOW, SE_ATC_GCOMP, SE_ATC_DR3), 3),
  CI_lower_logOR = round(c(ATC_naive_CI[1], ATC_EB_CI[1], ATC_IOW_Hajek_CI[1], ATC_GCOMP_CI[1], ATC_DR3_CI[1]), 3),
  CI_upper_logOR = round(c(ATC_naive_CI[2], ATC_EB_CI[2], ATC_IOW_Hajek_CI[2], ATC_GCOMP_CI[2], ATC_DR3_CI[2]), 3),
  OR        = round(exp(c(ATC_naive, ATC_EB, ATC_IOW_Hajek, ATC_GCOMP, ATC_DR3)), 3),
  CI_lower_OR = round(exp(c(ATC_naive_CI[1], ATC_EB_CI[1], ATC_IOW_Hajek_CI[1], ATC_GCOMP_CI[1], ATC_DR3_CI[1])), 3),
  CI_upper_OR = round(exp(c(ATC_naive_CI[2], ATC_EB_CI[2], ATC_IOW_Hajek_CI[2], ATC_GCOMP_CI[2], ATC_DR3_CI[2])), 3)
)
print(summary_table, row.names = FALSE)







#############################################
# Covariate balance table (xtable / LaTeX output)
# Replicates the LaTeX Table 1 structure with
# an additional column for IOW-weighted intervention SAT
#############################################

library(xtable)

# --- Compute MAIC-weighted statistics ---
# wt_EB (normalized, sum to 1) was computed on the original intervention data
# before X_all/S_all were redefined for the simulated external control

# MAIC ESS
ESS_MAIC <- 1 / sum(wt_EB^2)

# MAIC-weighted means (wt_EB already sums to 1)
maic_mean_age   <- sum(wt_EB * intervention_input$AGE)
maic_mean_sex   <- sum(wt_EB * intervention_input$SEX)
maic_mean_ecog  <- sum(wt_EB * intervention_input$ECOG0)
maic_mean_smoke <- sum(wt_EB * intervention_input$SMOKE)

# MAIC-weighted SD of age
maic_sd_age <- sqrt(sum(wt_EB * (intervention_input$AGE - maic_mean_age)^2))

# --- Compute IOW-weighted statistics ---
# iow_normalized (sum to 1) was computed earlier

iow_mean_age   <- sum(iow_normalized * intervention_input$AGE)
iow_mean_sex   <- sum(iow_normalized * intervention_input$SEX)
iow_mean_ecog  <- sum(iow_normalized * intervention_input$ECOG0)
iow_mean_smoke <- sum(iow_normalized * intervention_input$SMOKE)

# IOW-weighted SD of age
iow_sd_age <- sqrt(sum(iow_normalized * (intervention_input$AGE - iow_mean_age)^2))

# --- Format helper ---
fmt2 <- function(x) sprintf("%.2f", x)

# --- Build the table data frame ---
balance_df <- data.frame(
  Covariate = c(
    "Age in years (mean; SD)",
    "Sex (proportion male)",
    "ECOG (proportion status 1)",
    "Smoking (proportion smokers)"
  ),
  Intervention = c(
    paste0(fmt2(mean(intervention_input$AGE)), "; ", fmt2(sd(intervention_input$AGE))),
    fmt2(mean(intervention_input$SEX)),
    fmt2(mean(intervention_input$ECOG0)),
    fmt2(mean(intervention_input$SMOKE))
  ),
  External = c(
    paste0(fmt2(target_pop$age.mean), "; ", fmt2(target_pop$age.sd)),
    fmt2(target_pop$prop.male),
    fmt2(target_pop$prop.ecog0),
    fmt2(target_pop$prop.smoke)
  ),
  IOW = c(
    paste0(fmt2(iow_mean_age), "; ", fmt2(iow_sd_age)),
    fmt2(iow_mean_sex),
    fmt2(iow_mean_ecog),
    fmt2(iow_mean_smoke)
  ),
  MAIC = c(
    paste0(fmt2(maic_mean_age), "; ", fmt2(maic_sd_age)),
    fmt2(maic_mean_sex),
    fmt2(maic_mean_ecog),
    fmt2(maic_mean_smoke)
  ),
  stringsAsFactors = FALSE
)

# --- Generate LaTeX code via xtable ---

# Custom column headers with multiline information
# We use sanitize.text.function to allow raw LaTeX in headers

col_headers <- c(
  "\\textbf{Covariate}",
  paste0("\\textbf{Intervention SAT} \\\\ ($n_1=", n1, "$)"),
  paste0("\\textbf{External control} \\\\ ($n_0=", comparator_n, "$)"),
  paste0("\\textbf{IOW-weighted} \\\\ \\textbf{intervention SAT} \\\\ ($\\textrm{ESS}=", fmt2(ESS_IOW), "$)"),
  paste0("\\textbf{MAIC-weighted} \\\\ \\textbf{intervention SAT} \\\\ ($\\textrm{ESS}=", fmt2(ESS_MAIC), "$)")
)

xt <- xtable(balance_df,
             caption = paste0(
               "Summary statistics of the four baseline covariates identified as ",
               "imbalanced prognostic factors, before and after weighting using MAIC ",
               "(entropy balancing) and normalized inverse odds weighting (IOW). ",
               "The standard deviation of age in the weighted columns is ",
               "$\\sqrt{\\sum_{i=1}^{n_{1}} v_{i} (X_{1,i} - \\sum_{i=1}^{n_{1}} v_{i} X_{1,i})^2}$, ",
               "where $X_{1,i}$ and $v_i$ are the age and the weight, respectively, ",
               "for subject $i=1,\\dots, n_1$ in the intervention SAT."
             ),
             label = "tab:balance",
             align = c("l", "l", "c", "c", "c", "c"))

# Print LaTeX to console
cat("\n% ---- LaTeX table code ----\n")
print(xt,
      include.rownames = FALSE,
      sanitize.text.function = identity,  # allow LaTeX commands in cells
      sanitize.colnames.function = function(x) {
        # Replace the auto-generated column names with our custom headers
        col_headers
      },
      hline.after = c(-1, 0, nrow(balance_df)),
      table.placement = "!htb",
      caption.placement = "bottom",
      booktabs = FALSE)
