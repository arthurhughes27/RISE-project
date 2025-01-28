#' Applies the evaluation stage of the RISE procedure to assess the surrogate strength 
#' of a weighted linear combination of standardized genes. The selected genes are 
#' typically the output of the \code{rise_screen()} function, which identifies strong surrogate candidates.
#' 
#' @param Y_evaluate Numeric vector containing the primary response values in the evaluation set.
#' @param X_evaluate Numeric matrix or dataframe where each column represents a surrogate candidate in the
#' evaluation set. Rows must correspond to those in \code{Y_evaluate}.
#' @param A_evaluate Binary vector indicating treatment values in the evaluation set.
#' Rows must correspond to those in \code{Y_evaluate}.
#' @param markers Character vector specifying the markers (column names in \code{X_evaluate}) to be
#' combined for evaluation as a single predictor.
#' @param power_desired Numeric value between 0 and 1 specifying the desired power for the combined
#' surrogate to detect a treatment effect. Either this or \code{epsilon} must be provided.
#' @param epsilon Numeric value between 0 and 1 specifying the margin for non-inferiority testing.
#' Either this or \code{power_desired} must be provided.
#' @param individual Logical: If \code{TRUE}, evaluates each marker individually in the evaluation set.
#' Default is \code{TRUE}.
#' @param weights Dataframe or matrix specifying weights for each marker for the combined predictor.
#' Must have two columns: \code{marker} (marker name) and \code{weight} (inverse weights from the screening stage).
#' If not provided, weights are not used. Typically derived from the \code{rise_screen()} function.
#' @param plot Logical: If \code{TRUE}, generates a plot of ranks of the evaluation set response
#' against the ranks of the new surrogate. Displays the Spearman correlation coefficient.
#' @param alpha Numeric value specifying the significance level. Default is 0.05.
#' 
#' @author Arthur Hughes
#' 
#' @import dplyr 
#' @import pbmcapply
#' @import SurrogateRank
#' @import ggplot2
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{X_evaluate_standardised}: Dataframe of standardised, weighted surrogate 
#'   candidates used to construct the univariate surrogate \code{gamma}, whose strength is evaluated.
#'   \item \code{evaluation_results}: Dataframe summarizing the surrogacy test for \code{gamma},
#'   including delta estimate, confidence intervals, standard deviation, non-inferiority margin (\code{epsilon}), and p-value.
#'   \item \code{individual_evaluation}: If \code{individual = TRUE}, a dataframe summarizing the 
#'   surrogacy test results for each individual marker.
#'   \item \code{plot}: If \code{plot = TRUE}, a ggplot2 object showing ranks of the primary response
#'   versus the univariate surrogate \code{gamma}, including the Spearman rank correlation coefficient.
#' }


rise_evaluate = function(Y_evaluate, 
                         X_evaluate,
                         A_evaluate,
                         markers,
                         power_desired = NULL,
                         epsilon = NULL,
                         individual = TRUE, 
                         weights = NULL,
                         plot = TRUE,
                         alpha = 0.05){
  
  # validity checks
  n = length(Y_evaluate)
  P = ncol(X_evaluate)
  
  if(nrow(X_evaluate) != n){
    stop("Number of rows in X_evaluate does not match length of Y_evaluate")
  }
  
  if(length(A_evaluate) != n){
    stop("Length of A_evaluate does not match length of Y_evaluate")
  }
  
  if(is.null(power_desired) & is.null(epsilon)){
    stop("either power_desired or epsilon argument must be specified")
  }
  
  # If individual = TRUE, evaluate markers individually
  if (individual == TRUE){
    
    X_evaluate = as.matrix(X_evaluate)
    Y_evaluate = as.matrix(Y_evaluate)
    A_evaluate = as.character(A_evaluate)
    A_evaluate = factor(A_evaluate)  
    P = length(markers)
    
    # function to perform surrogate test in parallel
    surrogate_test = function(df, ind){
      yone = Y_evaluate[which(A_evaluate == 1)]
      yzero = Y_evaluate[which(A_evaluate == 0)]
      sone = X_evaluate[which(A_evaluate == 1),ind]
      szero = X_evaluate[which(A_evaluate == 0),ind]
      
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
                         function(col_index) surrogate_test(X_evaluate[,markers], col_index),
                         mc.cores = 1)
    if (length(markers) > 1){
      results_vec = cbind(colnames(X_evaluate), do.call(rbind, results)) %>% as.data.frame()
    } else {
      results_vec = cbind(markers, do.call(rbind, results)) %>% as.data.frame()
    }
    results_vec[,c(2:5)] = results_vec[,c(2:5)] %>% sapply(as.numeric)
    colnames(results_vec) = c("marker", "delta","sd", "epsilon", "p_unadjusted")
    individual_evaluation = results_vec %>%
      mutate(
        z_alpha = qnorm(1 - alpha),
        ci_upper = delta + z_alpha * sd,
        ci_lower = -1,
        p_adjusted = p.adjust(p_unadjusted, method = "BH")
      ) %>%
      dplyr::select(marker, epsilon, delta, sd, ci_lower, ci_upper, p_unadjusted, p_adjusted)
  } else {
    individual_evaluation = NULL
  }
  
  
  X_standard = X_evaluate %>% scale()
  
  if (!is.null(weights)){
    weight_vector <- setNames(weights$weight, weights$marker)
    # Divide each column in X_standard by the corresponding weight
    if (length(markers) > 1){
      X_standard <- sweep(X_standard, 2, weight_vector[colnames(X_standard)], "/")
    } else {
      X_standard <- sweep(X_standard, 2, weight_vector[markers], "/")
    }
    
    gamma = X_standard %>% rowSums()
  } else {
    gamma = X_standard %>% rowSums()  
  }
  
  yone = Y_evaluate[which(A_evaluate == 1)]
  yzero = Y_evaluate[which(A_evaluate == 0)]
  sone = gamma[which(A_evaluate == 1)]
  szero = gamma[which(A_evaluate == 0)]
  
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
  res_evaluate = c("gamma_S",
                   ss.test$delta.estimate, 
                   ss.test$sd.delta, 
                   ss.test$epsilon.used,
                   pnorm(ss.test$delta.estimate, ss.test$epsilon.used, ss.test$sd.delta))
  
  results_vec = t(res_evaluate) %>% as.data.frame()
  results_vec[c(2:5)] = results_vec[c(2:5)] %>% sapply(as.numeric)
  colnames(results_vec) = c("marker", "delta","sd", "epsilon", "p_unadjusted")
  results_vec = results_vec %>%
    mutate(
      z_alpha = qnorm(1 - alpha),
      ci_upper = delta + z_alpha * sd,
      ci_lower = -1
    ) %>%
    dplyr::select(marker, epsilon, delta, sd, ci_lower, ci_upper, p_unadjusted)
  
  
  rank_df = data.frame("treatment" = A_evaluate,
                       "response_rank" = rank(Y_evaluate),
                       "gamma_rank" = rank(gamma))
  
  rank_plot = rank_df %>% ggplot(aes(x = response_rank, y = gamma_rank, col = as.factor(treatment))) + 
    geom_point(size = 5) +
    xlab("Primary Response (rank)") +
    ylab(expression(gamma[S]~"(rank)")) +
    ggtitle("Ranks of primary response vs new surrogate in evaluation data") +
    guides(col = guide_legend(title="Treatment")) +
    scale_color_manual(labels = levels(A_evaluate %>% as.factor()), values = c("#1D8A99","#c1121f")) +
    theme_minimal() +
    geom_abline(slope=1, col= "red", linetype = "dashed", size = 1.4) +
    coord_fixed(ratio = 1) +
    annotate(size = 15, 
             "text",
             x = Inf, 
             y = Inf, 
             label = bquote(rho == .(round(cor(as.numeric(rank_df[,2]), as.numeric(rank_df[,3])), 2))), 
             vjust = 2, 
             hjust = 3.3, 
             color = "red") +
    theme_minimal(base_size = 22) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.background = element_rect(fill = 'white', color = 'white'))
  
  results_evaluate = list("X_evaluate_standardised" = X_standard,
                          "evaluation_results" = results_vec,
                          "individual_evaluation" = individual_evaluation,
                          "plot" = rank_plot)
  return(results_evaluate)
}
