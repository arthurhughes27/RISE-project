#' Applies the validation stage of the RISE procedure to evaluate the surrogate strength 
#' of a weighted linear combination of standardised genes. The selected genes may typically be 
#' the output of the \code{rise_screen()} function, which screens for strong surrogate candidates.
#' 
#' @param Y_validate Numeric vector containing the primary response values in the validation set
#' @param X_validate Numeric matrix or dataframe with each column a surrogate candidate in the
#' validation set. The rows must be in the same order as in \code{Y}.
#' @param A_validate Binary vector giving the treatment values in the validation set.
#' The rows must be in the same order as in \code{Y}. 
#' @param markers character vector giving the markers (column names in X_validate) to be
#' combined for validation as a single predictor. 
#' @param power_desired value between 0 and 1 giving the user's desired power for the new combined
#' surrogate to be able to detect a treatment effect. 
#' Either this or the \code{epsilon} argument must be supplied.
#' @param epsilon value between 0 and 1 giving the margin for the non-inferiority testing. 
#' Either this or the \code{power_desired} argument must be supplied.
#' @param individual logical : if \code{TRUE}, outputs the marker-wise validation 
#' (i.e. each significant marker is re-tested individually on the validation set). 
#' Default is \code{TRUE}.
#' @param weights the weights for each marker for the combined predictor. 
#' The format is a dataframe or matrix with two columns : \code{marker} giving the marker name
#'and \code{weight} giving the inverse weights (i.e. delta from the screening stage).
#'If not provided, weights will not be used. 
#'Typically given as the \code{weights} argument from the \code{rise_screen()} function.
#'@param plot logical : plot the ranks of the validation set response against the ranks of the
#'new surrogate. The spearman correlation coefficient will be displayed. 
#' 
#' @author Arthur Hughes
#' 
#' @import dplyr 
#' @import pbmcapply
#' @import SurrogateRank
#' @import ggplot2
#' 

rise_validate = function(Y_validate, 
                         X_validate,
                         A_validate,
                         markers,
                         power_desired = NULL,
                         epsilon = NULL,
                         individual = TRUE, 
                         weights = NULL,
                         plot = TRUE){
  
  # validity checks
  n = length(Y_validate)
  P = ncol(X_validate)
  
  if(nrow(X_validate) != n){
    stop("Number of rows in X_validate does not match length of Y_validate")
  }
  
  if(length(A_validate) != n){
    stop("Length of A_validate does not match length of Y_validate")
  }
  
  if(is.null(power_desired) & is.null(epsilon)){
    stop("either power_desired or epsilon argument must be specified")
  }
  
  # If individual = TRUE, validate markers individually
  if (individual == TRUE){
    
    X_validate = as.matrix(X_validate)
    Y_validate = as.matrix(Y_validate)
    A_validate = as.character(A_validate)
    A_validate = factor(A_validate)  
    P = length(markers)
    
    # function to perform surrogate test in parallel
    surrogate_test = function(df, ind){
      yone = Y_validate[which(A_validate == 1)]
      yzero = Y_validate[which(A_validate == 0)]
      sone = X_validate[which(A_validate == 1),ind]
      szero = X_validate[which(A_validate == 0),ind]
      
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
                         function(col_index) surrogate_test(X_validate[,markers], col_index),
                         mc.cores = 1)
    if (length(markers) > 1){
      results_vec = cbind(colnames(X_validate[,markers]), do.call(rbind, results)) %>% as.data.frame()
    } else {
      results_vec = cbind(markers, do.call(rbind, results)) %>% as.data.frame()
    }
    results_vec[,c(2:5)] = results_vec[,c(2:5)] %>% sapply(as.numeric)
    colnames(results_vec) = c("marker", "delta","sd", "epsilon", "p_unadjusted")
    individual_validation = results_vec %>%
      mutate(
        z_alpha = qnorm(1 - alpha),
        ci_upper = delta + z_alpha * sd,
        ci_lower = -1,
        p_adjusted = p.adjust(p_unadjusted, method = "BH")
      ) %>%
      select(marker, epsilon, delta, sd, ci_lower, ci_upper, p_unadjusted, p_adjusted)
  } else {
    individual_validation = NULL
  }

  
  X_standard = X_validate[ ,markers] %>% scale()
  
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
  
  yone = Y_validate[which(A_validate == 1)]
  yzero = Y_validate[which(A_validate == 0)]
  sone = gamma[which(A_validate == 1)]
  szero = gamma[which(A_validate == 0)]
    
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
    res_validate = c("gamma_S",
                     ss.test$delta.estimate, 
                     ss.test$sd.delta, 
                     ss.test$epsilon.used,
                     pnorm(ss.test$delta.estimate, ss.test$epsilon.used, ss.test$sd.delta))
    
    results_vec = t(res_validate) %>% as.data.frame()
    results_vec[c(2:5)] = results_vec[c(2:5)] %>% sapply(as.numeric)
    colnames(results_vec) = c("marker", "delta","sd", "epsilon", "p_unadjusted")
    results_vec = results_vec %>%
      mutate(
        z_alpha = qnorm(1 - alpha),
        ci_upper = delta + z_alpha * sd,
        ci_lower = -1
      ) %>%
      select(marker, epsilon, delta, sd, ci_lower, ci_upper, p_unadjusted)
    
    
    rank_df = data.frame("treatment" = A_validate,
                         "response_rank" = rank(Y_validate),
                         "gamma_rank" = rank(gamma))
    
    rank_plot = rank_df %>% ggplot(aes(x = response_rank, y = gamma_rank, col = as.factor(treatment))) + 
      geom_point(size = 5) +
      xlab("Primary Response (rank)") +
      ylab(expression(gamma[S]~"(rank)")) +
      ggtitle("Ranks of primary response vs new surrogate in validation data") +
      guides(col = guide_legend(title="Treatment")) +
      scale_color_manual(labels = levels(A_validate %>% as.factor()), values = c("#1D8A99","#c1121f")) +
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
    
    results_validate = list("X_validate_standardised" = X_standard,
                            "validation_results" = results_vec,
                            "individual_validation" = individual_validation,
                            "plot" = rank_plot)
    return(results_validate)
    }