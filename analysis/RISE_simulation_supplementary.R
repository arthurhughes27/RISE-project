# Load libraries
library(ggplot2)
library(dplyr)
library(quantreg)
library(Rsurrogate)
library(SurrogateRank)
library(MASS)
library(reshape2)
library(gridExtra)
library(cowplot)
library(grid)
library(tidyr)

# Set seed
set.seed(23102024)

# Set parameters
alpha <- 0.05 # nominal significance level
beta <- 0.2 # desired power for univariate surrogate test = (1-beta)*100% 
s0_mean <- 0 # invalid surrogate mean value in untreated group
s0_sd <- 1 # invalid surrogate  sd value in untreated group
s1_mean <- 0 # invalid surrogate mean value in treated group
s1_sd <- 1 # invalid surrogate sd value in treated group
y0_mean <- 0 # primary response mean value in untreated group
y0_sd <- 1 # primary response sd value in untreated group
y1_mean <- 3 # primary response mean value in treated group
y1_sd <- 1 # primary response sd value in treated group
valid_sigma <- 10 # level of perturbation for valid surrogates (lower = stronger surrogates)
n <- 50 # total sample size
n0 <- n / 2 # treated sample size 
n1 <- n / 2 # untreated sample size 
p <- 500 # total number of predictors
prop_valid <- 0.1 # proportion of predictors which are valid surrogates
p_valid <- prop_valid * p # number of valid surrogates
p_invalid = (1 - prop_valid) * p # number of invalid surrogates
n_sim <- 500 # number of simulations
corr = 0 # correlation - off-diagonal element in surrogate covariance matrix

# Data generation function
gen.data <- function(n1, n0, p, prop_valid, valid_sigma, corr, mode = "simple") {
  if (mode == "simple"){ # DG1
    
    # calculate the number of valid and invalid surrogates
    p_valid = prop_valid * p 
    p_invalid = (1 - prop_valid) * p
    
    # generate primary responses from multivariate normal distributions
    y1 <- rnorm(n1, y1_mean, y1_sd) # treated response 
    y0 <- rnorm(n0, y0_mean, y0_sd) # untreated repsonse
    
    # generate candidate surrogates 
    mm <- runif(p_invalid, min = 0.5, max = 2.5) # generate random invalid surrogate means from uniform distribution
    ss <- runif(p_invalid, min = 0.5, max = 2)  # generate random invalid surrogate sds from uniform distribution
    Sigma_invalid = matrix(corr, nrow = p_invalid, ncol= p_invalid) # Invalid surrogate covariance matrix
    diag(Sigma_invalid) = ss 
    
    if (prop_valid != 0){ # if valid surrogates required
      
      # valid surrogate covariance matrix
      Sigma_valid = matrix(corr*valid_sigma, nrow = p_valid, ncol= p_valid) 
      diag(Sigma_valid) = rep(valid_sigma, p_valid)
      
      # valid surrogates in treated group by perturbing primary response
      s1.valid <- matrix(y1, nrow = n1, ncol = p_valid, byrow = TRUE) + 
        mvrnorm(n = n1, mu = rep(0, p_valid), Sigma = Sigma_valid)
      
      # valid surrogates in untreated group by perturbing primary response
      s0.valid <- matrix(y0, nrow = n0, ncol = p_valid, byrow = TRUE) + 
        mvrnorm(n = n0, mu = rep(0, p_valid), Sigma = Sigma_valid)
      
      # invalid surrogates
      if (prop_valid != 1){ # if invalid surrogates required
        s1.invalid <- mvrnorm(n = n1, mu = mm, Sigma = Sigma_invalid)
        s0.invalid <- mvrnorm(n = n0, mu = mm, Sigma = Sigma_invalid)
        # bind candidates together 
        s1 <- cbind(s1.valid, s1.invalid)
        s0 <- cbind(s0.valid, s0.invalid)
      } else {
        s1 = s1.valid 
        s0 = s0.valid 
      }
      
    } else { # if no valid surrogates required 
      # invalid surrogates
      s1 <- mvrnorm(n = n1, mu = mm, Sigma = Sigma_invalid)
      s0 <- mvrnorm(n = n0, mu = mm, Sigma = Sigma_invalid)
    }
    # store hypothesis truths
    hyp <- c(rep("null false", p_valid), rep("null true", p_invalid))
    
    # return generated data as a list
    return(list(y1 = y1, y0 = y0, s1 = s1, s0 = s0, hyp = hyp))
  } else if (mode == "complex") { # DG2
    # calculate the number of valid and invalid surrogates
    p_valid = prop_valid * p 
    p_invalid = (1 - prop_valid) * p
    
    # generate primary responses from multivariate normal distributions
    y1 <- rnorm(n1, y1_mean, y1_sd) # treated response 
    y0 <- rnorm(n0, y0_mean, y0_sd) # untreated repsonse
    
    # generate candidate surrogates 
    lambda = runif(p_invalid, min = 0.5, max = 2.5)
    if (prop_valid != 0){ #if valid surrogates required
      # valid surrogate covariance matrix
      Sigma_valid = matrix(corr*valid_sigma, nrow = p_valid, ncol= p_valid) 
      diag(Sigma_valid) = rep(valid_sigma, p_valid)
      s1.valid = matrix(y1^3, nrow = n1, ncol = p_valid, byrow = TRUE) + 
        mvrnorm(n = n1, mu = rep(0, p_valid), Sigma = Sigma_valid)
      
      s0.valid = matrix(y0^3, nrow = n0, ncol = p_valid, byrow = TRUE) + 
        mvrnorm(n = n0, mu = rep(0, p_valid), Sigma = Sigma_valid)
      
      if (prop_valid != 1){ #if invalid surrogates required
        s0.invalid = sapply(lambda, function(rate) rexp(n0, rate))
        s1.invalid = sapply(lambda, function(rate) rexp(n1, rate))
        s1 <- cbind(s1.valid, s1.invalid)
        s0 <- cbind(s0.valid, s0.invalid)
      } else {
        s1 = s1.valid 
        s0 = s0.valid  
      }
    } else {
      s0 = sapply(lambda, function(rate) rexp(n0, rate))
      s1 = sapply(lambda, function(rate) rexp(n1, rate))
    }
    hyp <- c(rep("null false", p_valid), rep("null true", p_invalid))
    return(list(y1 = y1, y0 = y0, s1 = s1, s0 = s0, hyp = hyp))
  }
}


# Function to estimate the ground truth using asymptotic properties 
calc.truth <- function(p, prop_valid, valid_sigma, corr, mode = "simple") {
  # generate a dataset with a large sample size 
  dd <- gen.data(n1 = 10000, n0 = 10000, p, prop_valid, valid_sigma, corr = corr, mode = mode)
  
  # find the asymptotic U statistic for Y
  uy <- (10000 * 10000)^(-1) * wilcox.test(dd$y1, dd$y0)$statistic
  
  # find the asymptotic U statistic for each surrogate
  us <- numeric(length(which(dd[["hyp"]] == "null false")))
  for (j in which(dd[["hyp"]] == "null false")) {
    us[j] <- (10000 * 10000)^(-1) * wilcox.test(dd$s1[, j], dd$s0[, j])$statistic
  }
  
  # take the truth as the mean difference
  delta_mean <- mean(uy - us)
  delta_sd <- sd(uy - us)
  return(list(uy_true = uy, 
              us_true = us, 
              delta_true = delta_mean, 
              delta_true_sd = delta_sd))
}

# truth <- calc.truth(p, prop_valid, valid_sigma, corr, mode = "simple")

# Function to simulate results (performance metrics)
simulate_results <- function(n1, n0, p, prop_valid, n_sim, valid_sigma, corr, mode = "simple") {
  p_unadjusted <- matrix(0, nrow = n_sim, ncol = p) # store unadjusted p_values for each candidate and for each simulation
  fpr_unadjusted <- numeric(n_sim) # store unadjusted false positive rates per simulation
  fdr_unadjusted <- numeric(n_sim) # store unadjusted proportion of false positives per simulation
  tpr_unadjusted <- numeric(n_sim) # store unadjusted true positive rates (empirical power) per simulation
  ppv_unadjusted <- numeric(n_sim) # store unadjusted positive predictive values per simulation
  
  p_bonf <- matrix(0, nrow = n_sim, ncol = p) # store bonferroni corrected p-values
  fpr_bonf <- numeric(n_sim)
  fdr_bonf <- numeric(n_sim)
  tpr_bonf <- numeric(n_sim)
  ppv_bonf <- numeric(n_sim)
  
  p_bh <- matrix(0, nrow = n_sim, ncol = p) # store benjamini-hochberg corrected p-values
  fpr_bh <- numeric(n_sim)
  fdr_bh <- numeric(n_sim)
  tpr_bh <- numeric(n_sim)
  ppv_bh <- numeric(n_sim)
  
  p_by <- matrix(0, nrow = n_sim, ncol = p) # store benjamini-yekutieli corrected p-values
  fpr_by <- numeric(n_sim)
  fdr_by <- numeric(n_sim)
  tpr_by <- numeric(n_sim)
  ppv_by <- numeric(n_sim)
  
  for (k in 1:n_sim) { # for loop to do many simulations
    # generate data
    data <- gen.data(n1, n0, p, prop_valid, valid_sigma, corr, mode = mode)
    
    # Estimate u_y
    u_y_estimated <- test.surrogate(
      yone = data$y1, yzero = data$y0, sone = data$y1, szero = data$y0, epsilon = 0.1
    )$u.y
    # u_y_true = truth$uy_true
    # choose epsilon as maximum of estimated value of UY - 0.5 and 0
    eps <- max(0, u_y_estimated - 0.5)
    
    for (j in 1:p) {  # for each candidate surrogate 
      ss.test <- test.surrogate(
        yone = data$y1, 
        yzero = data$y0, 
        sone = data$s1[, j], # jth candidate surrogate 
        szero = data$s0[, j], 
        epsilon = eps)
      # compute unadjusted p-value for jth candidate surrogate
      p_unadjusted[k, j] <- pnorm(ss.test$delta.estimate, ss.test$epsilon.used, ss.test$sd.delta)
    }
    # correct p-values per-simulation with bonferroni, BH, and BY. 
    p_bonf[k, ] = p.adjust(p_unadjusted[k, ], method = "bonferroni")
    p_bh[k, ] = p.adjust(p_unadjusted[k, ], method = "BH")
    p_by[k, ] = p.adjust(p_unadjusted[k, ], method = "BY")
    
    
    p_values_unadjusted = p_unadjusted[k, ]
    # Classify test outcomes for uncorrected p-values
    TP_unadjusted <- sum(data[["hyp"]] == "null false" & p_values_unadjusted < alpha) # true positives
    FP_unadjusted <- sum(data[["hyp"]] == "null true" & p_values_unadjusted < alpha) # false positives
    TN_unadjusted <- sum(data[["hyp"]] == "null true" & p_values_unadjusted >= alpha) # true negatives
    FN_unadjusted <- sum(data[["hyp"]] == "null false" & p_values_unadjusted >= alpha) # false negatives
    
    # repeat for each correction method
    p_values_bonf = p_bonf[k, ]
    TP_bonf <- sum(data[["hyp"]] == "null false" & p_values_bonf < alpha)
    FP_bonf <- sum(data[["hyp"]] == "null true" & p_values_bonf < alpha)
    TN_bonf <- sum(data[["hyp"]] == "null true" & p_values_bonf >= alpha)
    FN_bonf <- sum(data[["hyp"]] == "null false" & p_values_bonf >= alpha)
    
    p_values_bh = p_bh[k, ]
    TP_bh <- sum(data[["hyp"]] == "null false" & p_values_bh < alpha)
    FP_bh <- sum(data[["hyp"]] == "null true" & p_values_bh < alpha)
    TN_bh <- sum(data[["hyp"]] == "null true" & p_values_bh >= alpha)
    FN_bh <- sum(data[["hyp"]] == "null false" & p_values_bh >= alpha)
    
    p_values_by = p_by[k, ]
    TP_by <- sum(data[["hyp"]] == "null false" & p_values_by < alpha)
    FP_by <- sum(data[["hyp"]] == "null true" & p_values_by < alpha)
    TN_by <- sum(data[["hyp"]] == "null true" & p_values_by >= alpha)
    FN_by <- sum(data[["hyp"]] == "null false" & p_values_by >= alpha)
    
    
    if (prop_valid == 0){ # if no true surrogates we only calculate the false positive rate
      fpr_unadjusted[k] <- FP_unadjusted / (FP_unadjusted + TN_unadjusted) # FPR = FP/(FP + TN)                 # Proportion of false positives
      fdr_unadjusted[k] <- 0  # set the other metrics as 0
      tpr_unadjusted[k] <- 0                 
      ppv_unadjusted[k] <- 0
    } else if (prop_valid != 0){  # if any true positives, calculate FPR, PFP, TPR, PPV
      fpr_unadjusted[k] <- FP_unadjusted / (FP_unadjusted + TN_unadjusted)  # False positive rate = FP/(FP + TN)   
      fdr_unadjusted[k] <- FP_unadjusted / (TP_unadjusted + FP_unadjusted)  # Proportion of False Positives = FP/(FP + TP) 
      tpr_unadjusted[k] <- TP_unadjusted / (TP_unadjusted + FN_unadjusted)  # True Positive Rate (Sensitivity) = TP/(TP + FN)
      ppv_unadjusted[k] <- TP_unadjusted / (TP_unadjusted + FP_unadjusted)  # Positive Predictive Value = TP/(TP + FP)
      
      # repeat for all the correction methods
      fpr_bonf[k] <- FP_bonf / (FP_bonf + TN_bonf)                 
      fdr_bonf[k] <- FP_bonf / (TP_bonf + FP_bonf)                  
      tpr_bonf[k] <- TP_bonf / (TP_bonf + FN_bonf)                
      ppv_bonf[k] <- TP_bonf / (TP_bonf + FP_bonf)                 
      
      fpr_bh[k] <- FP_bh / (FP_bh + TN_bh)                
      fdr_bh[k] <- FP_bh / (TP_bh + FP_bh)                 
      tpr_bh[k] <- TP_bh / (TP_bh + FN_bh)                  
      ppv_bh[k] <- TP_bh / (TP_bh + FP_bh)                
      
      fpr_by[k] <- FP_by / (FP_by + TN_by)                 
      fdr_by[k] <- FP_by / (TP_by + FP_by)                
      tpr_by[k] <- TP_by / (TP_by + FN_by)                 
      ppv_by[k] <- TP_by / (TP_by + FP_by)                 
    }
    # print progress within simulation runs
    # print(paste0("Sigma = ", valid_sigma, ", p = ", p, ", n = ", n, " : simulation ", k, " of ", n_sim))
  }
  
  # Where we have NA in our metrics due to division by 0 
  # (i.e. in particular, where we do not identify any positives)
  # replace with a 0 
  fpr_unadjusted[is.nan(fpr_unadjusted)] <- 0
  fdr_unadjusted[is.nan(fdr_unadjusted)] <- 0
  tpr_unadjusted[is.nan(tpr_unadjusted)] <- 0
  ppv_unadjusted[is.nan(ppv_unadjusted)] <- 0
  
  fpr_bonf[is.nan(fpr_bonf)] <- 0
  fdr_bonf[is.nan(fdr_bonf)] <- 0
  tpr_bonf[is.nan(tpr_bonf)] <- 0
  ppv_bonf[is.nan(ppv_bonf)] <- 0
  
  fpr_bh[is.nan(fpr_bh)] <- 0
  fdr_bh[is.nan(fdr_bh)] <- 0
  tpr_bh[is.nan(tpr_bh)] <- 0
  ppv_bh[is.nan(ppv_bh)] <- 0
  
  fpr_by[is.nan(fpr_by)] <- 0
  fdr_by[is.nan(fdr_by)] <- 0
  tpr_by[is.nan(tpr_by)] <- 0
  ppv_by[is.nan(ppv_by)] <- 0
  
  # calculate the mean metrics across simulations
  avg_fpr_unadjusted <- mean(fpr_unadjusted)
  avg_fdr_unadjusted <- mean(fdr_unadjusted)
  avg_tpr_unadjusted <- mean(tpr_unadjusted)
  avg_ppv_unadjusted <- mean(ppv_unadjusted)
  
  avg_fpr_bonf <- mean(fpr_bonf)
  avg_fdr_bonf <- mean(fdr_bonf)
  avg_tpr_bonf <- mean(tpr_bonf)
  avg_ppv_bonf <- mean(ppv_bonf)
  
  avg_fpr_bh <- mean(fpr_bh)
  avg_fdr_bh <- mean(fdr_bh)
  avg_tpr_bh <- mean(tpr_bh)
  avg_ppv_bh <- mean(ppv_bh)
  
  avg_fpr_by <- mean(fpr_by)
  avg_fdr_by <- mean(fdr_by)
  avg_tpr_by <- mean(tpr_by)
  avg_ppv_by <- mean(ppv_by)
  
  # store the metrics in lists per correction method
  metrics_unadjusted = list("avg_fpr" = avg_fpr_unadjusted,
                            "avg_fdr" = avg_fdr_unadjusted,
                            "avg_tpr" = avg_tpr_unadjusted,
                            "avg_ppv" = avg_ppv_unadjusted)
  
  metrics_bonf = list("avg_fpr" = avg_fpr_bonf,
                      "avg_fdr" = avg_fdr_bonf,
                      "avg_tpr" = avg_tpr_bonf,
                      "avg_ppv" = avg_ppv_bonf)
  
  metrics_bh = list("avg_fpr" = avg_fpr_bh,
                    "avg_fdr" = avg_fdr_bh,
                    "avg_tpr" = avg_tpr_bh,
                    "avg_ppv" = avg_ppv_bh)
  
  metrics_by = list("avg_fpr" = avg_fpr_by,
                    "avg_fdr" = avg_fdr_by,
                    "avg_tpr" = avg_tpr_by,
                    "avg_ppv" = avg_ppv_by)
  
  metrics = list("metrics_unadjusted" = metrics_unadjusted, 
                 "metrics_bonf" = metrics_bonf,
                 "metrics_bh" = metrics_bh,
                 "metrics_by" = metrics_by)
  
  # print progress of simulation number
  # print(paste0("p = ", p, ", n = ", n, " : simulation complete!"))
  
  p_values = list("p_unadjusted" = p_unadjusted, 
                  "p_bonf" = p_bonf,
                  "p_bh" = p_bh,
                  "p_by" = p_by)
  # output of function are the p-value lists and the performance metrics
  return(list("p_values" = p_values, "metrics"= metrics))
}



##--Figure S1 :  distribution of p-values at the boundary value under the null--##
alpha <- 0.05 # nominal significance level
beta <- 0.2 # desired power for univariate surrogate test = (1-beta)*100% 
s0_mean <- 0 # invalid surrogate mean value in untreated group
s0_sd <- 1 # invalid surrogate  sd value in untreated group
s1_mean <- 0 # invalid surrogate mean value in treated group
s1_sd <- 1 # invalid surrogate sd value in treated group
y0_mean <- 0 # primary response mean value in untreated group
y0_sd <- 1 # primary response sd value in untreated group
y1_mean <- 3 # primary response mean value in treated group
y1_sd <- 1 # primary response sd value in treated group
valid_sigma <- 10 # level of perturbation for valid surrogates (lower = stronger surrogates)
n <- 50 # total sample size
n0 <- n / 2 # treated sample size 
n1 <- n / 2 # untreated sample size 
p <- 500 # total number of predictors
prop_valid <- 0 # proportion of predictors which are valid surrogates
p_valid <- prop_valid * p # number of valid surrogates
p_invalid = (1 - prop_valid) * p # number of invalid surrogates
n_sim <- 1000 # number of simulations
corr = 0 # correlation - off-diagonal element in surrogate covariance matrix

p_figs1 = matrix(0)

for (j in 1:n_sim){
    res <- simulate_results(n1 = n / 2, 
                            n0 = n / 2, 
                            p = p, 
                            prop_valid = 0, 
                            n_sim = 1, 
                            valid_sigma,
                            corr = corr, 
                            mode = "simple")
    # Append p-values to results 
    p_figs1 = append(p_figs1,res[["p_values"]][["p_unadjusted"]][1,] %>% as.matrix())
    
    print(paste0("simulation ", j, " of ", n_sim, " complete." ))
  }

p_figs1 =  data.frame("p" = p_figs1[-1])

save(p_figs1, file = "./output/simulation_results_final/metrics_figures1.Rdata")
metrics_figures1 = get(load("./output/simulation_results_final/metrics_figures1.Rdata"))

fig_s1 = metrics_figures1 %>% ggplot(aes(x = p)) +
  geom_histogram(binwidth = 0.03, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of P-values under the null hypothesis", x = "P-value", y = "Frequency") +
  xlim(0,1) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = 'white', color = 'white'))

fig_s1

ggsave(plot = fig_s1,
       "./output/figures/main/figures1.pdf", 
       units = "cm", 
       width = 35, 
       height = 7)

##--Figure S2-- Empirical power and proportion of false positives 
## as a function of surrogate strength
## under different correction methods for a fixed sample size--##

## Step 1 - find sigmas corresponding to desired surrogate strengths
## We want average surrogate strengths from U_S = 0.55 to U_Y (to 2 d.p.)
## We do this by trial and error

# round(mean(calc.truth(p, prop_valid, valid_sigma = 244, corr, mode = "simple")$us_true),2)
# round(mean(calc.truth(p, prop_valid, valid_sigma = 68, corr, mode = "simple")$us_true),2)
# round(mean(calc.truth(p, prop_valid, valid_sigma = 30, corr, mode = "simple")$us_true),2)
# round(mean(calc.truth(p, prop_valid, valid_sigma = 15, corr, mode = "simple")$us_true),2)
# round(mean(calc.truth(p, prop_valid, valid_sigma = 9, corr, mode = "simple")$us_true),2)
# round(mean(calc.truth(p, prop_valid, valid_sigma = 5.5, corr, mode = "simple")$us_true),2)
# round(mean(calc.truth(p, prop_valid, valid_sigma = 3, corr, mode = "simple")$us_true),2)
# round(mean(calc.truth(p, prop_valid, valid_sigma = 1.8, corr, mode = "simple")$us_true),2)
# round(mean(calc.truth(p, prop_valid, valid_sigma = 0.65, corr, mode = "simple")$us_true),2)
# round(mean(calc.truth(p, prop_valid, valid_sigma = 0.01, corr, mode = "simple")$us_true),2)

# Set up the values to simulate
# Set parameters
alpha <- 0.05 # nominal significance level
beta <- 0.2 # desired power for univariate surrogate test = (1-beta)*100% 
s0_mean <- 0 # invalid surrogate mean value in untreated group
s0_sd <- 1 # invalid surrogate  sd value in untreated group
s1_mean <- 0 # invalid surrogate mean value in treated group
s1_sd <- 1 # invalid surrogate sd value in treated group
y0_mean <- 0 # primary response mean value in untreated group
y0_sd <- 1 # primary response sd value in untreated group
y1_mean <- 3 # primary response mean value in treated group
y1_sd <- 1 # primary response sd value in treated group
n <- 50 # total sample size
n0 <- n / 2 # treated sample size 
n1 <- n / 2 # untreated sample size 
p <- 500 # total number of predictors
prop_valid <- 0.1 # proportion of predictors which are valid surrogates
p_valid <- prop_valid * p # number of valid surrogates
p_invalid = (1 - prop_valid) * p # number of invalid surrogates
n_sim <- 500 # number of simulations
corr = 0 # correlation - off-diagonal element in surrogate covariance matrix
sigma_grid = c(0.01,0.65,1.8,3,5.5,9,15,30,68,244) # values of sigma

# Initialise an empty list for p-values and dataframe for metrics
p_values = list()
metrics_figs2 <- data.frame(
  n = integer(), 
  p = integer(), 
  sigma_S = numeric(),
  avg_us = numeric(),
  avg_fpr_unadjusted = numeric(),
  avg_fdr_unadjusted = numeric(),
  avg_tpr_unadjusted = numeric(),
  avg_ppv_unadjusted = numeric(),
  avg_fpr_bonf = numeric(),
  avg_fdr_bonf = numeric(),
  avg_tpr_bonf = numeric(),
  avg_ppv_bonf = numeric(),
  avg_fpr_bh = numeric(),
  avg_fdr_bh = numeric(),
  avg_tpr_bh = numeric(),
  avg_ppv_bh = numeric(),
  avg_fpr_by = numeric(),
  avg_fdr_by = numeric(),
  avg_tpr_by = numeric(),
  avg_ppv_by = numeric()
)

# Perform simulations over grid
for (sigma in sigma_grid){
  avg_us = round(mean(calc.truth(p, prop_valid, valid_sigma = sigma, corr, mode = "simple")$us_true),2)
  res <- simulate_results(n1 = n / 2, 
                          n0 = n / 2, 
                          p = p, 
                          prop_valid, 
                          n_sim = n_sim, 
                          valid_sigma = sigma, 
                          corr = corr, 
                          mode = "simple")
      
      # p_values = append(p_values,list(res[["p_values"]][["p_unadjusted"]]))
      
      # Append results to data frame
      metrics_figs2 <- rbind(metrics_figs2, data.frame(
        n = n, p = p,
        sigma_S = sigma,
        avg_us = avg_us,
        avg_fpr_unadjusted = res[["metrics"]][["metrics_unadjusted"]][["avg_fpr"]], 
        avg_fdr_unadjusted = res[["metrics"]][["metrics_unadjusted"]][["avg_fdr"]], 
        avg_tpr_unadjusted = res[["metrics"]][["metrics_unadjusted"]][["avg_tpr"]],
        avg_ppv_unadjusted = res[["metrics"]][["metrics_unadjusted"]][["avg_ppv"]],
        avg_fpr_bonf = res[["metrics"]][["metrics_bonf"]][["avg_fpr"]], 
        avg_fdr_bonf = res[["metrics"]][["metrics_bonf"]][["avg_fdr"]], 
        avg_tpr_bonf = res[["metrics"]][["metrics_bonf"]][["avg_tpr"]],
        avg_ppv_bonf = res[["metrics"]][["metrics_bonf"]][["avg_ppv"]],
        avg_fpr_bh = res[["metrics"]][["metrics_bh"]][["avg_fpr"]], 
        avg_fdr_bh = res[["metrics"]][["metrics_bh"]][["avg_fdr"]], 
        avg_tpr_bh = res[["metrics"]][["metrics_bh"]][["avg_tpr"]],
        avg_ppv_bh = res[["metrics"]][["metrics_bh"]][["avg_ppv"]],
        avg_fpr_by = res[["metrics"]][["metrics_by"]][["avg_fpr"]], 
        avg_fdr_by = res[["metrics"]][["metrics_by"]][["avg_fdr"]], 
        avg_tpr_by = res[["metrics"]][["metrics_by"]][["avg_tpr"]],
        avg_ppv_by = res[["metrics"]][["metrics_by"]][["avg_ppv"]]
      ))
      paste0("Simulation ", which(sigma_grid == sigma), " of ", length(sigma_grid), " complete.")
}

save(metrics_figs2, file = "./output/simulation_results_final/metrics_figures2.Rdata")

metrics_figs2 = get(load("./output/simulation_results_final/metrics_figures2.Rdata"))

# Pivot metrics to long format 
metrics_figs2_fdr_long <- metrics_figs2  %>%
  pivot_longer(cols = c(avg_fdr_unadjusted, 
                        avg_fdr_bonf, 
                        avg_fdr_bh, 
                        avg_fdr_by), 
               names_to = "Correction", 
               values_to = "Value") %>% 
  dplyr::select(avg_us, Correction, Value)

metrics_figs2_tpr_long <- metrics_figs2  %>%
  pivot_longer(cols = c(avg_tpr_unadjusted, 
                        avg_tpr_bonf, 
                        avg_tpr_bh, 
                        avg_tpr_by), 
               names_to = "Correction", 
               values_to = "Value") %>% 
  dplyr::select(avg_us, Correction, Value)

# Plot
FDR_figs2 = ggplot(metrics_figs2_fdr_long, aes(x = avg_us, y = Value, linetype = Correction)) +
  geom_line(color = "red") +
  labs(x = "U_S", y = "Proportion of False Positives", linetype = "Correction") +
  scale_linetype_manual(values = c("avg_fdr_unadjusted" = "solid", 
                                   "avg_fdr_bonf" = "dashed",
                                   "avg_fdr_bh" = "dotted",
                                   "avg_fdr_by" = "dotdash"),
                        labels = c("avg_fdr_unadjusted" = "Unadjusted", 
                                   "avg_fdr_bonf" = "Bonferroni",
                                   "avg_fdr_bh" = "B-H",
                                   "avg_fdr_by" = "B-Y")) + 
  theme_minimal(base_size = 18) +
  theme(plot.background = element_rect(fill = 'white', color = 'white'))

TPR_figs2 = ggplot(metrics_figs2_tpr_long, aes(x = avg_us, y = Value, linetype = Correction)) +
  geom_line(color = "red") +
  labs(x = "U_S", y = "Proportion of False Positives", linetype = "Correction") +
  scale_linetype_manual(values = c("avg_tpr_unadjusted" = "solid", 
                                   "avg_tpr_bonf" = "dashed",
                                   "avg_tpr_bh" = "dotted",
                                   "avg_tpr_by" = "dotdash"),
                        labels = c("avg_tpr_unadjusted" = "Unadjusted", 
                                   "avg_tpr_bonf" = "Bonferroni",
                                   "avg_tpr_bh" = "B-H",
                                   "avg_tpr_by" = "B-Y")) + 
  theme_minimal(base_size = 18) +
  theme(plot.background = element_rect(fill = 'white', color = 'white'))

## FPR and TPR in the same grid
TPR_patchwork_figs2 = TPR_figs2  + 
  ylim(0,1) + 
  theme(legend.position = "none",
        title = element_blank())

FDR_patchwork_figs2 = FDR_figs2  + 
  ylim(0,1) + 
  theme(legend.position = "none",
        title = element_blank())

legend = get_legend(FDR_figs2)

figures2 = grid.arrange(
  arrangeGrob(TPR_patchwork_figs2, FDR_patchwork_figs2, legend, 
              ncol = 3, widths = c(1, 1, 0.3)),  # Adjust widths as needed
  top = textGrob("Empirical Power and Proportion of False Positives as a function of surrogate strength - multiple testing corrections",
                 gp = gpar(fontsize = 17, fontface = "bold"))  # Adjust fontsize and fontface
)

ggsave(plot = figures2,
       "./output/figures/main/figures2.pdf", 
       units = "cm", 
       width = 35, 
       height = 18)


##--Figure s3 - PFP and empirical power as a function of correlation--##
##-- for a fixed sample size and surrogate strength--##
# Set parameters
alpha <- 0.05 # nominal significance level
beta <- 0.2 # desired power for univariate surrogate test = (1-beta)*100% 
s0_mean <- 0 # invalid surrogate mean value in untreated group
s0_sd <- 1 # invalid surrogate sd value in untreated group
s1_mean <- 0 # invalid surrogate mean value in treated group
s1_sd <- 2 # invalid surrogate sd value in treated group
y0_mean <- 0 # primary response mean value in untreated group
y0_sd <- 1 # primary response sd value in untreated group
y1_mean <- 3 # primary response mean value in treated group
y1_sd <- 1 # primary response sd value in treated group
valid_sigma <- 1.8 # set sigma such that average surrogate strength = 0.9
n <- 50 # fixed total sample size
n0 <- n / 2 # treated sample size 
n1 <- n / 2 # untreated sample size 
p <- 500 # total number of predictors
prop_valid <- 0.1 # proportion of predictors which are valid surrogates
p_valid <- prop_valid * p # number of valid surrogates
p_invalid = (1 - prop_valid) * p # number of invalid surrogates
n_sim <- 500 # number of simulations

corr_grid = seq(0,0.5,0.1) # grid of correlation values

metrics_figs3 = matrix(0,nrow = length(corr_grid)*n_sim, ncol = 3)

for (i in 1:length(corr_grid)){
  for (j in 1:n_sim){
    res <- simulate_results(n1 = n / 2, 
                            n0 = n / 2, 
                            p = p, 
                            prop_valid, 
                            n_sim = 1, 
                            valid_sigma,
                            corr = corr_grid[i],
                            mode = "simple")
    
    # Append results to data frame
    
    metrics_figs3[(i-1)*(n_sim) + j,1] = corr_grid[i]
    metrics_figs3[(i-1)*(n_sim) + j,2] = res[["metrics"]][["metrics_unadjusted"]][["avg_fdr"]]
    metrics_figs3[(i-1)*(n_sim) + j,3] = res[["metrics"]][["metrics_unadjusted"]][["avg_tpr"]]
    
    
    print(paste0("Correlation : ", corr_grid[i], ", simulation ", j, " of ", n_sim, " complete." ))
  }
}

metrics_figs3 = data.frame("correlation" = metrics_figs3[,1], "fdr" = metrics_figs3[,2], "tpr" = metrics_figs3[,3])

metrics_figs3 %>% ggplot(aes(x = as.factor(correlation), y = fdr))


fdr_plot_s3 = metrics_figs3 %>% ggplot(aes(x = as.factor(correlation), y = fdr, fill = as.factor(correlation))) + 
                           geom_violin() +
                           ylim(0,1) + 
                           ylab("Proportion of False Positives") +
                           xlab("Inter-predictor correlation") + 
                           ggtitle("Proportion of False Positives as a function of correlation strength") +
                           theme_minimal(base_size = 18) +
                           theme(plot.background = element_rect(fill = 'white', color = 'white'),
                                 plot.title = element_blank(),
                                 legend.position = "none")

tpr_plot_s3 = metrics_figs3 %>% ggplot(aes(x = as.factor(correlation), y = tpr, fill = as.factor(correlation))) + 
  geom_violin() +
  ylim(0,1) + 
  ylab("Empirical Power") +
  xlab("Inter-predictor correlation") + 
  ggtitle("Proportion of False Positives as a function of correlation strength") +
  theme_minimal(base_size = 18) +
  theme(plot.background = element_rect(fill = 'white', color = 'white'),
        plot.title = element_blank(),
        legend.position = "none")

figures3 = grid.arrange(
  arrangeGrob(tpr_plot_s3, fdr_plot_s3,  
              ncol = 2, widths = c(1, 1)),  # Adjust widths as needed
  top = textGrob("Empirical Power and Proportion of False Positives as a function of correlation",
                 gp = gpar(fontsize = 20, fontface = "bold"))  # Adjust fontsize and fontface
)

save(metrics_figs3, file = "./output/simulation_results_final/metrics_figures3.Rdata")

ggsave(plot = figures3,
       "./output/figures/main/figures3.pdf", 
       units = "cm", 
       width = 35, 
       height = 18)

##--Figure s4 - boxplots of FPR as a function of sample size in the uncorrelated setting, complex model--##

# Set parameters
alpha <- 0.05 # nominal significance level
beta <- 0.2 # desired power for univariate surrogate test = (1-beta)*100% 
s0_mean <- 0 # invalid surrogate mean value in untreated group
s0_sd <- 1 # invalid surrogate sd value in untreated group
s1_mean <- 0 # invalid surrogate mean value in treated group
s1_sd <- 2 # invalid surrogate sd value in treated group
y0_mean <- 0 # primary response mean value in untreated group
y0_sd <- 1 # primary response sd value in untreated group
y1_mean <- 3 # primary response mean value in treated group
y1_sd <- 1 # primary response sd value in treated group
p <- 500 # total number of predictors
prop_valid <- 0 # proportion of predictors which are valid surrogates
p_valid <- prop_valid * p # number of valid surrogates
p_invalid = (1 - prop_valid) * p # number of invalid surrogates
n_sim <- 500 # number of simulations

n_grid = seq(10,100,10)

valid_sigma = 10
truth <- calc.truth(p, prop_valid, valid_sigma, corr, mode = "complex")
mean(truth$us_true)

metrics_fig4 = matrix(0,nrow = length(n_grid)*n_sim, ncol = 2)

for (i in 1:length(n_grid)){
  for (j in 1:n_sim){
    res <- simulate_results(n1 = n_grid[i] / 2, 
                            n0 = n_grid[i] / 2, 
                            p = p, 
                            prop_valid = 0, 
                            n_sim = 1, 
                            valid_sigma,
                            corr = 0,
                            mode = "complex")
    
    # p_values = append(p_values,list(res[["p_values"]][["p_unadjusted"]]))
    
    # Append results to data frame
    
    metrics_fig4[(i-1)*(n_sim) + j,1] = n_grid[i]
    metrics_fig4[(i-1)*(n_sim) + j,2] = res[["metrics"]][["metrics_unadjusted"]][["avg_fpr"]]
    print(paste0("Sample size : ", n_grid[i], ", simulation ", j, " of ", n_sim, " complete." ))
  }
}

metrics_fig4 = data.frame("n" = metrics_fig4[,1],
                           "FPR" = metrics_fig4[,2])

save(metrics_fig4, file = "./output/simulation_results_final/metrics_figure4.Rdata")

fig4 = metrics_fig4 %>% ggplot(aes(x = as.factor(n), 
                                     y = FPR, 
                                     fill = as.factor(n))) +
  geom_boxplot() +
  labs(title = "Observed FPR across different sample sizes - complex data generation",
       x = "Sample size",
       y = "FPR",
       fill = "Sample Size") + 
  ylim(0,0.15) +
  geom_hline(yintercept = 0.05, colour = "red", linetype = "dashed", linewidth = 1.1) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = 'white', color = 'white'),
        legend.position = "none")

fig4

ggsave(plot = fig4,
       "./output/figures/main/figure4.pdf", 
       units = "cm", 
       width = 40, 
       height = 25)

##--Figure s5 - violin plots of FPR as a function of correlation--##

# Set parameters
alpha <- 0.05 # nominal significance level
beta <- 0.2 # desired power for univariate surrogate test = (1-beta)*100% 
s0_mean <- 0 # invalid surrogate mean value in untreated group
s0_sd <- 1 # invalid surrogate sd value in untreated group
s1_mean <- 0 # invalid surrogate mean value in treated group
s1_sd <- 2 # invalid surrogate sd value in treated group
y0_mean <- 0 # primary response mean value in untreated group
y0_sd <- 1 # primary response sd value in untreated group
y1_mean <- 3 # primary response mean value in treated group
y1_sd <- 1 # primary response sd value in treated group
valid_sigma <- 10 # level of perturbation for valid surrogates (lower = stronger surrogates)
n <- 50 # total sample size
n0 <- n / 2 # treated sample size 
n1 <- n / 2 # untreated sample size 
p <- 500 # total number of predictors
prop_valid <- 0 # proportion of predictors which are valid surrogates
p_valid <- prop_valid * p # number of valid surrogates
p_invalid = (1 - prop_valid) * p # number of invalid surrogates
n_sim <- 500 # number of simulations

corr_grid = seq(0,0.5,0.1)

metrics_fig2 = matrix(0,nrow = length(corr_grid)*n_sim, ncol = 2)

for (i in 1:length(corr_grid)){
  for (j in 1:n_sim){
    res <- simulate_results(n1 = n / 2, 
                            n0 = n / 2, 
                            p = p, 
                            prop_valid = 0, 
                            n_sim = 1, 
                            valid_sigma,
                            corr = corr_grid[i], 
                            mode = "simple")
    
    # p_values = append(p_values,list(res[["p_values"]][["p_unadjusted"]]))
    
    # Append results to data frame
    
    metrics_fig2[(i-1)*(n_sim) + j,1] = corr_grid[i]
    metrics_fig2[(i-1)*(n_sim) + j,2] = res[["metrics"]][["metrics_unadjusted"]][["avg_fpr"]]
    print(paste0("Correlation : ", corr_grid[i], ", simulation ", j, " of ", n_sim, " complete." ))
  }
}

metrics_fig2 = data.frame("correlation" = metrics_fig2[,1],
                          "FPR" = metrics_fig2[,2])

save(metrics_fig2, file = "./output/simulation_results_final/metrics_figure2.Rdata")

# metrics_fig2 = get(load("./output/simulation_results_final/simulation/metrics_figure2.Rdata"))

fig2 = metrics_fig2 %>% ggplot(aes(x = as.factor(correlation), 
                                   y = FPR, 
                                   fill = as.factor(correlation))) +
  geom_violin() +
  labs(title = "Observed FPR across different correlation levels",
       x = "Inter-predictor correlation",
       y = "FPR",
       fill = "Correlation") + 
  geom_hline(yintercept = 0.05, colour = "red", linetype = "dashed", linewidth = 1.1) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = 'white', color = 'white'),
        legend.position = "none")

fig2

ggsave(plot = fig2,
       "./output/figures/main/figure2.pdf", 
       units = "cm", 
       width = 40, 
       height = 25)
                         