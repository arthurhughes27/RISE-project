#' Applies the Rank-based Identification of high-dimensional SurrogatE markers
#' (RISE) screening process to identify a set of candidate surrogates for a given response.
#' 
#' @param Y Numeric vector containing the primary response values.
#' @param X Numeric matrix or dataframe where each column represents a surrogate candidate. 
#' Rows must correspond to those in \code{Y}.
#' @param A Binary vector indicating the treatment values. 
#' Rows must correspond to those in \code{Y}.
#' @param reference String specifying the reference level (i.e., untreated group) in \code{A}.
#' @param alpha Numeric value for the desired significance level for screening. Default is 0.05. 
#' @param power_desired Numeric value between 0 and 1 specifying the desired power for the individual markers 
#' to detect a treatment effect. Either this or \code{epsilon} must be supplied.
#' @param epsilon Numeric value between 0 and 1 specifying the margin for non-inferiority testing. 
#' Either this or \code{power_desired} must be supplied.
#' @param p_correction Character string describing the multiple testing correction to be applied to the raw p-values. 
#' Options are: \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, 
#' \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}. Default is \code{"none"}.
#' @param cores Integer specifying the number of cores to use for parallel computation. Default is 1. 
#' @param mode Either \code{"independent"} (default) or \code{"paired"} depending on study design.
#' @param test Type of test procedure to conduct. Default is \code{"non.inferiority"}, where a candidate is 
#'        declared a valid trial-level surrogate if the treatment effect on it is at least at great as the 
#'        treatment effect on the outcome. The alternative is \code{"two.one.sided"}, which performs two one sided 
#'        non-inferiority tests to determine if the treatment effect is within a symmetric interval of width 2 epsilon
#'        around the treatment effect on the outcome. 
#' @author Arthur Hughes
#' 
#' @import dplyr 
#' @import pbmcapply
#' @import SurrogateRank
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{results}: A dataframe summarizing the results of the surrogate screening process. 
#'   Includes columns:
#'   \itemize{
#'     \item \code{marker}: The name of the surrogate marker.
#'     \item \code{epsilon}: The non-inferiority margin used in the surrogacy test.
#'     \item \code{delta}: The estimate of the treatment effect mediated by the surrogate marker.
#'     \item \code{sd}: The standard deviation of the delta estimate.
#'     \item \code{ci_lower}: Lower bound of the confidence interval for \code{delta}.
#'     \item \code{ci_upper}: Upper bound of the confidence interval for \code{delta}.
#'     \item \code{p_unadjusted}: Raw p-value for the surrogacy test.
#'     \item \code{p_adjusted}: P-value after applying the specified multiple testing correction.
#'   }
#'   \item \code{significant_markers}: A character vector of markers with adjusted p-values below \code{alpha}.
#'   \item \code{weights}: A dataframe containing the significant markers and their corresponding weights 
#'   (i.e., \code{delta} estimates) for constructing a weighted surrogate.
#' }


rise_screen = function(Y, 
                       X,
                       A,
                       reference = NULL, 
                       alpha = 0.05,
                       power_desired = NULL,
                       epsilon = NULL,
                       p_correction = "none", 
                       cores = 1,
                       mode = "independent", 
                       test = "non.inferiority"){
  
  # validity checks
  n = length(Y)
  P = ncol(X)
  
  if(nrow(X) != n){
    stop("Number of rows in X does not match length of Y")
  }
  
  if(length(A) != n){
    stop("Length of A does not match length of Y")
  }
  
  if(is.null(power_desired) & is.null(epsilon)){
    stop("either power_desired or epsilon argument must be specified")
  }
  
  # Coercion of data
  X = as.matrix(X)
  Y = as.matrix(Y)
  A = as.character(A)
  A = factor(A, levels = c(setdiff(unique(A), reference), reference), labels = c(1,0))  
  
  # function to perform surrogate test in parallel
  surrogate_test = function(df, ind){
    yone = Y[which(A == 1)]
    yzero = Y[which(A == 0)]
    sone = X[which(A == 1),ind]
    szero = X[which(A == 0),ind]
    
    if (is.null(epsilon)){
      ss.test = test.surrogate(yone = yone, 
                                      yzero = yzero, 
                                      sone = sone, 
                                      szero = szero, 
                                      power.want.s = power_desired,
                                      mode = mode,
                                      test = test)
    } else {
      ss.test = test.surrogate(yone = yone, 
                                      yzero = yzero, 
                                      sone = sone, 
                                      szero = szero, 
                                      epsilon = epsilon,
                                      mode = mode,
                                      test = test)
    }
    
    # compute p-value
    # p = pnorm(ss.test$delta.estimate, ss.test$epsilon.used, ss.test$sd.delta)
    # return delta, sd, epsilon and p-value
    return(c(ss.test$delta.estimate, 
             ss.test$ci.delta,
             ss.test$sd.delta, 
             ss.test$epsilon.used, 
             ss.test$p.delta))
  }
  
  results = pbmclapply(c(1:P),
                       function(col_index) surrogate_test(X, col_index),
                       mc.cores = cores)
  
  results_vec = cbind(colnames(X), do.call(rbind, results)) %>% as.data.frame()
  results_vec[,c(2:7)] = results_vec[,c(2:7)] %>% sapply(as.numeric)
  colnames(results_vec) = c("marker", "delta","ci_lower", "ci_upper","sd", "epsilon", "p_unadjusted")
  results_vec = results_vec %>%
    mutate(
      p_adjusted = p.adjust(p_unadjusted, method = p_correction)
    ) %>%
    dplyr::select(marker, epsilon, delta, sd, ci_lower, ci_upper, p_unadjusted, p_adjusted)
  
  significant_markers = (results_vec %>% filter(p_adjusted < alpha))$marker
  
  weights = data.frame("marker" = significant_markers,
                       "weight" = (results_vec %>% filter(p_adjusted < alpha))$delta)
  
  results_screening = list("results" = results_vec, "significant_markers" = significant_markers, "weights" = weights)
  return(results_screening)
}


