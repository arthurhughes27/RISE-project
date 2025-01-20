#' Applies the Rank-based Identification of high-dimensional SurrogatE markers
#' screening process to identify a set of candidate surrogates for a given response
#' 
#' @param Y Numeric vector containing the primary response values
#' @param X Numeric matrix or dataframe with each column a surrogate candidate. 
#' The rows must be in the same order as in \code{Y}.
#' @param A Binary vector giving the treatment values.
#' The rows must be in the same order as in \code{Y}. 
#' @param reference String giving the reference level (i.e. untreated group) in \code{A}.
#' @param alpha Desired significance level for screening, defaults to 0.05. 
#' @param power_desired value between 0 and 1 giving the user's desired power for the individual markers 
#' to be able to detect a treatment effect. 
#' Either this or the \code{epsilon} argument must be supplied.
#' @param epsilon value between 0 and 1 giving the margin for the non-inferiority testing. 
#' Either this or the \code{power_desired} argument must be supplied.
#' @param p_correction Describes the multiple testing correction to be applied to the raw p-values. 
#' One of \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, 
#' \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.
#' @param cores number of cores to commit to parallel computation. Default is 1. 
#' 
#' @author Arthur Hughes
#' 
#' @import dplyr 
#' @import pbmcapply
#' @import SurrogateRank
#' 

rise_screen = function(Y, 
                       X,
                       A,
                       reference = NULL, 
                       alpha = 0.05,
                       power_desired = NULL,
                       epsilon = NULL,
                       p_correction = "none", 
                       cores = 1){

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
                             power.want.s = power_desired)
    } else {
    ss.test = test.surrogate(yone = yone, 
                             yzero = yzero, 
                             sone = sone, 
                             szero = szero, 
                             epsilon = epsilon)
    }
    
    # compute p-value
    p = pnorm(ss.test$delta.estimate, ss.test$epsilon.used, ss.test$sd.delta)
    # return delta, sd, epsilon and p-value
    return(c(ss.test$delta.estimate, 
             ss.test$sd.delta, 
             ss.test$epsilon.used, 
             p))
    }

  results = pbmclapply(c(1:P),
                       function(col_index) surrogate_test(X, col_index),
                       mc.cores = cores)
  
  results_vec = cbind(colnames(X), do.call(rbind, results)) %>% as.data.frame()
  results_vec[,c(2:5)] = results_vec[,c(2:5)] %>% sapply(as.numeric)
  colnames(results_vec) = c("marker", "delta","sd", "epsilon", "p_unadjusted")
  results_vec = results_vec %>%
    mutate(
      z_alpha = qnorm(1 - alpha),
      ci_upper = delta + z_alpha * sd,
      ci_lower = -1,
      p_adjusted = p.adjust(p_unadjusted, method = p_correction)
    ) %>%
    dplyr::select(marker, epsilon, delta, sd, ci_lower, ci_upper, p_unadjusted, p_adjusted)
  
  significant_markers = (results_vec %>% filter(p_adjusted < alpha))$marker
  
  weights = data.frame("marker" = significant_markers,
                       "weight" = (results_vec %>% filter(p_adjusted < alpha))$delta)
  
  results_screening = list("results" = results_vec, "significant_markers" = significant_markers, "weights" = weights)
  return(results_screening)
}
  
  
