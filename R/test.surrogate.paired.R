#' Rank-Based test for Trial-level Surrogacy
#'
#' This function tests whether a candidate surrogate marker explains the treatment effect
#' using a rank-based non-inferiority framework. It computes the difference in treatment effects 
#' on the outcome and surrogate and checks whether the loss is within an acceptable margin.
#'
#' @param full.data Optional data matrix with 3 columns: outcome, surrogate, and treatment indicator (1 = treatment, 0 = control).
#' @param yone Outcome values for the treatment group (ignored if \code{full.data} is provided).
#' @param yzero Outcome values for the control group.
#' @param sone Surrogate values for the treatment group.
#' @param szero Surrogate values for the control group.
#' @param epsilon Optional non-inferiority margin. If \code{NULL}, it will be estimated using \code{power.want.s}.
#' @param power.want.s Desired power for the surrogate effect (default is 0.7). Used to estimate \code{epsilon} if not provided.
#' @param u.y.hyp Optional hypothesized value for treatment effect on the outcome (\code{u.y}).
#'        If provided, \code{epsilon} is computed as \code{u.y.hyp - u.s.power}.
#' @param alpha Significance level for the one-sided test (default is 0.05).
#' @param mode Either \code{"independent"} (default) or \code{"paired"} depending on study design.
#' @param test Type of test procedure to conduct. Default is \code{"non.inferiority"}, where a candidate is 
#'        declared a valid trial-level surrogate if the treatment effect on it is at least at great as the 
#'        treatment effect on the outcome. The alternative is \code{"two.one.sided"}, which performs two one sided 
#'        non-inferiority tests to determine if the treatment effect is within a symmetric interval of width 2 epsilon
#'        around the treatment effect on the outcome. 
#' 
#' 
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{u.y}}{Estimated rank-based treatment effect on the outcome.}
#'   \item{\code{u.s}}{Estimated rank-based treatment effect on the surrogate.}
#'   \item{\code{delta.estimate}}{Estimated difference in treatment effects: \code{u.y - u.s}.}
#'   \item{\code{sd.u.y}}{Standard deviation of \code{u.y}.}
#'   \item{\code{sd.u.s}}{Standard deviation of \code{u.s}.}
#'   \item{\code{sd.delta}}{Standard deviation of \code{delta.estimate}.}
#'   \item{\code{ci.delta}}{One-sided confidence interval upper bound for \code{delta.estimate}.}
#'   \item{\code{p.delta}}{p-value for validity of trial-level surrogacy.}
#'   \item{\code{epsilon.used}}{Non-inferiority threshold used in the test.}
#'   \item{\code{is.surrogate}}{\code{TRUE} if the surrogate passes the test, else \code{FALSE}.}
#' }
#'
#' @examples
#' set.seed(123)
#' yone <- rnorm(50, mean = 1)
#' yzero <- rnorm(50, mean = 0)
#' sone <- rnorm(50, mean = 0.8)
#' szero <- rnorm(50, mean = 0)
#' test.surrogate(yone = yone, yzero = yzero, sone = sone, szero = szero)
#'
#' @export


test.surrogate = function (full.data = NULL, yone = NULL, yzero = NULL, sone = NULL, 
                           szero = NULL, epsilon = NULL, power.want.s = 0.7, u.y.hyp = NULL,
                           alpha = 0.05, mode = "independent", test = "non.inferiority") 
{
  # Validity checks
  if (mode == "paired"){
    if (length(yone) != length(yzero)){
      stop("You have requested the paired setting mode
           but the lengths of yone and yzero do not match.
           Ensure that the kth element of yzero is matched
           to the kth element of yone.")
    } else if (length(sone) != length(szero)){
      stop ("You have requested the paired setting mode
           but the lengths of sone and szero do not match.
           Ensure that the kth element of szero is matched
           to the kth element of sone.")
    }
  }
  
  # Extract variables from full.data if provided
  if (!is.null(full.data)) {
    yone = full.data[full.data[, 3] == 1, 1]
    yzero = full.data[full.data[, 3] == 0, 1]
    sone = full.data[full.data[, 3] == 1, 2]
    szero = full.data[full.data[, 3] == 0, 2]
  }
  
  # Compute treatment effects and standard errors using delta.calculate()
  dd = delta.calculate(yone = yone, yzero = yzero, sone = sone, 
                       szero = szero, mode = mode)
  
  if (test == "non.inferiority"){
    # Compute upper bound of one-sided (1 - alpha) CI for delta = u.y - u.s
    z.alpha = qnorm(1 - alpha)
    ci.delta = c(-1, dd$delta.estimate + z.alpha * dd$sd.delta)  # CI = (-∞, upper bound]
    
    # If epsilon is not provided, estimate it based on power.want.s
    if (is.null(epsilon) & mode == "independent") {
      n1 = length(yone)
      n0 = length(yzero)
      
      # Null SD under Wilcoxon/Mann-Whitney assumptions
      sd.null = sqrt((n1 + n0 + 1) / (12 * n1 * n0))
      
      # Adjust quantiles for two-sided test and power calculation
      z.alpha.2 = qnorm(1 - (alpha / 2))
      u.s.power = 0.5 - (qnorm(1 - power.want.s) - z.alpha.2) * sd.null
      
      # Estimate epsilon as the difference between observed u.y (or hypothesized) and u.s under power
      if (is.null(u.y.hyp)) {
        epsilon = dd$u.y - u.s.power
      } else {
        epsilon = u.y.hyp - u.s.power
      }
    } else if (is.null(epsilon) & mode == "paired"){
      # If epsilon is not provided, estimate it based on power.want.s
      n = length(yone)
      
      # Null SD under paired-setting assumptions
      sd.null = 1/(2*sqrt(n))
      
      # Adjust quantiles for two-sided test and power calculation
      z.alpha.2 = qnorm(1 - (alpha / 2))
      u.s.power = 0.5 - (qnorm(1 - power.want.s) - z.alpha.2) * sd.null
      
      # Estimate epsilon as the difference between observed u.y (or hypothesized) and u.s under power
      if (is.null(u.y.hyp)) {
        epsilon = dd$u.y - u.s.power
      } else {
        epsilon = u.y.hyp - u.s.power
      }
    }
    
    # Decision rule for non-inferiority:
    # If upper bound of CI is less than epsilon, surrogate is acceptable
    is.surrogate = ci.delta[2] < epsilon
    
    # Compute a corresponding p-value 
    p = pnorm(dd$delta.estimate, 
              epsilon, 
              dd$sd.delta)
    
  } else if (test == "two.one.sided"){
    n = length(yone)
    # Compute (1 - 2*alpha) CI for delta = u.y - u.s
    z.alpha = qnorm(1 - alpha)
    ci.delta = c(dd$delta.estimate - z.alpha * dd$sd.delta, dd$delta.estimate + z.alpha * dd$sd.delta)
    
    # Calculate p-value corresponding to null: delta > epsilon
    p1 = pnorm(dd$delta.estimate, 
               epsilon, 
               dd$sd.delta)
    
    # Calculate p-value corresponding to null: delta < -epsilon
    p2 = 1-pnorm(dd$delta.estimate, 
                 -epsilon, 
                 dd$sd.delta)
    
    p = max(p1,p2)
  }
  
  # Return all relevant quantities
  return(list(
    u.y = as.numeric(dd$u.y),
    u.s = as.numeric(dd$u.s),
    delta.estimate = as.numeric(dd$delta.estimate),
    sd.u.y = dd$sd.u.y,
    sd.u.s = dd$sd.u.s,
    sd.delta = as.numeric(dd$sd.delta),
    ci.delta = ci.delta,
    p.delta = p,
    epsilon.used = epsilon,
    is.surrogate = is.surrogate
  ))
}
