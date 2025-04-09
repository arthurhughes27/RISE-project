#' Calculate Delta: Difference in Rank-based Statistics for Two Outcomes
#'
#' This function estimates the difference (`delta`) between two rank-based statistics
#' (e.g., Wilcoxon statistics or paired ranks) for a primary outcome and a surrogate,
#' under either an independent or paired design.
#'
#' @param full.data Optional. A data matrix or data frame with three columns:
#'   - Column 1: Primary outcome (y)
#'   - Column 2: Surrogate variable (s)
#'   - Column 3: Treatment indicator (0 or 1)
#'   If `full.data` is provided, individual vectors (`yone`, `yzero`, `sone`, `szero`) are ignored.
#' 
#' @param yone Numeric vector of outcome values in the treatment group.
#' @param yzero Numeric vector of outcome values in the control group.
#' @param sone Numeric vector of surrogate values in the treatment group.
#' @param szero Numeric vector of surrogate values in the control group.
#' @param mode Character. Either "independent" (default) or "paired" to specify
#'        the study design.
#'
#' @return A list with the following elements:
#'   - `u.y`: Rank-based test statistic for the primary outcome
#'   - `u.s`: Rank-based test statistic for the surrogate
#'   - `delta.estimate`: Estimated difference between outcome and surrogate statistics
#'   - `sd.u.y`: Standard deviation of the outcome statistic
#'   - `sd.u.s`: Standard deviation of the surrogate statistic
#'   - `sd.delta`: Standard error of the delta estimate

delta.calculate.paired = function(full.data = NULL, yone = NULL, yzero = NULL, sone = NULL, 
                           szero = NULL, mode = "independent") 
{
  # If full.data is provided, extract outcome and surrogate values by treatment group
  if (!is.null(full.data)) {
    yone  = full.data[full.data[, 3] == 1, 1]  # outcome, treatment group
    yzero = full.data[full.data[, 3] == 0, 1]  # outcome, control group
    sone  = full.data[full.data[, 3] == 1, 2]  # surrogate, treatment group
    szero = full.data[full.data[, 3] == 0, 2]  # surrogate, control group
  }
  
  if (mode == "independent") {
    # Use Wilcoxon rank-sum test for independent samples
    test.y = wilcox.test(yone, yzero, exact = FALSE)
    test.s = wilcox.test(sone, szero, exact = FALSE)
    
    # Normalize test statistics to estimate U statistics
    n1.f = length(yone)
    n0.f = length(yzero)
    u.y = (n1.f * n0.f)^(-1) * test.y$statistic
    u.s = (n1.f * n0.f)^(-1) * test.s$statistic
    delta.estimate = u.y - u.s
    
    # Estimate the variance using Hajek projection variance estimators
    m.count = n1.f
    n.count = n0.f
    
    # Variance components from treated group comparisons
    V10.Xi.Y = sapply(yone, var.wil, b = yzero)
    V10.Xi.S = sapply(sone, var.wil, b = szero)
    
    # Variance components from control group comparisons (flip ranks)
    V01.Yj.Y = sapply(yzero, var.wil, b = yone, flip = TRUE)
    V01.Yj.S = sapply(szero, var.wil, b = sone, flip = TRUE)
    
    # Covariance matrix components from treated group
    s10.11.YY = var(V10.Xi.Y)
    s10.12.YS = cov(V10.Xi.Y, V10.Xi.S)
    s10.22.SS = var(V10.Xi.S)
    s10.21.SY = cov(V10.Xi.Y, V10.Xi.S)  # same as s10.12
    
    # Covariance matrix components from control group
    s01.11.YY = var(V01.Yj.Y)
    s01.12.YS = cov(V01.Yj.Y, V01.Yj.S)
    s01.22.SS = var(V01.Yj.S)
    s01.21.SY = cov(V01.Yj.Y, V01.Yj.S)
    
    # Combine into covariance matrices
    S10 = matrix(c(s10.11.YY, s10.12.YS, s10.21.SY, s10.22.SS), nrow = 2)
    S01 = matrix(c(s01.11.YY, s01.12.YS, s01.21.SY, s01.22.SS), nrow = 2)
    
    # Total variance matrix
    S.mat = (1/m.count) * S10 + (1/n.count) * S01
    
    # Standard deviations for outcome and surrogate
    sd.y = sqrt(S.mat[1, 1])
    sd.s = sqrt(S.mat[2, 2])
    
    # Linear combination L = (1, -1) to compute variance of delta = u.y - u.s
    L = t(as.matrix(c(1, -1)))
    sd.est = sqrt(L %*% S.mat %*% t(L))
    
  } else if (mode == "paired") {
    # For paired samples, use within-pair comparisons
    n = length(yone)
    
    # Proportion of pairs where yone > yzero (similar to concordance)
    u.y = mean(ifelse(yone > yzero, 1, 0))
    u.s = mean(ifelse(sone > szero, 1, 0))
    delta.estimate = u.y - u.s
    
    # Estimate variances
    sd.y = sd(ifelse(yone > yzero, 1, 0))
    sd.s = sd(ifelse(sone > szero, 1, 0))
    
    # Difference in indicators for each pair
    di = ifelse(yone > yzero, 1, 0) - ifelse(sone > szero, 1, 0)
    var_delta = var(di) / n
    sd.est = sqrt(var_delta)
  }
  
  # Return a list of results
  return(list(
    u.y = as.numeric(u.y),
    u.s = as.numeric(u.s),
    delta.estimate = as.numeric(delta.estimate),
    sd.u.y = sd.y,
    sd.u.s = sd.s,
    sd.delta = as.numeric(sd.est)
  ))
}
