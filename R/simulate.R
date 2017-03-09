
#' Simulate different unimodal distributions.
#' 
#' @description Simulate different unimodal, skewed distributions based
#' on different mean and variance parameters.
#' 
#' @param n Size of the simulated distribution
#' @param x_mean Mean
#' @param x_var Variance
#' @param N_item Number of items simulated
#' @param seed (Optional) Set a random seed
#' @param name (Optional) Sample based on the names given by 
#' Keats & Lord (1962). 
#' Options are: "GANA", "MACAA", "TQS8", "WM8", "WMI".
#' 
#' @details All the inner working of this procedure can be seen 
#' in Keats & Lord (1962).
#' 
#' 
#' @return Simulated values
#' 
#' @export
sim_unimodal <- function(n, x_mean, x_var, N_item, seed=NULL, name=NULL){
  data_names <- c("GANA", "MACAA", "TQS8", "WM8", "WMI")
  
  # Data from Keats & Lord (??)
  kl_mean = c(27.06, 27.06, 32.93, 32.93, 25.82, 25.82, 6.76, 6.75, 23.75, 23.88)
  kl_var = c(8.19, 8.19, 8.04, 8.04, 7.28, 7.28, 5.12, 5.11, 5.59, 5.62)^2
  kl_N_item = c(rep(40,2),rep(50,4),rep(30,4))
  kl_n = c(2354, 2000, 6103, 6103, 2000, 1800, 1000, 1200, 1000, 1200)
  
  requireNamespace("emdbook")
  
  if(!is.null(seed))
    set.seed(seed)
  
  if(!is.null(name)){
    idx_name <- grep(name, data_names)
    idx_x <- idx_name + idx_name - 1
    
    return(list(X=sim_unimodal(kl_n[idx_x], kl_mean[idx_x], kl_var[idx_x], kl_N_item[idx_x]), 
                Y=sim_unimodal(kl_n[idx_x+1], kl_mean[idx_x+1], kl_var[idx_x+1], kl_N_item[idx_x+1])))
  }
  
  pi <- x_mean / N_item
  theta <- ((N_item^2) * pi * (1-pi) - (x_var * (n-1))/n) / ((x_var * (n-1))/n - N_item * pi * (1-pi))
  
  rbetabinom(n=n,prob=pi,size=N_item,theta=theta)
}


#' Contaminate a sample on different quantile based cutpoints
#' 
#' @param x Input sample
#' @param percent Percentage of contamination
#' @param location Options: "bq" - Both quantiles, "low" - Inferior quantile,
#' "high" - Superior quantile, "bm" - Both minimum and maximum, 
#' "min" - Minimum, "max" - Maximum
#' 
#' @return Contaminated sample
#' 
#' @export
contaminate_sample <- function(x, percent=0.04, location="bq"){
  n <- length(x)
  n_r <- floor(percent*n)
  
  if(location == "bq"){
    lo <- as.numeric(rep(quantile(x, probs=0.025), n_r/2))
    hi <- as.numeric(rep(quantile(x, probs=0.975), n_r/2))
    ret <- c(lo, x, hi)
  } else if(location == "low"){
    lo <- as.numeric(rep(quantile(x, probs=0.025), n_r))
    ret <- c(lo, x)
  } else if(location == "high"){
    hi <- as.numeric(rep(quantile(x, probs=0.975), n_r))
    ret <- c(x, hi)
  } else if(location == "bm"){
    lo <- as.numeric(rep(min(x), n_r/2))
    hi <- as.numeric(rep(max(x), n_r/2))
    ret <- c(lo, x, hi)
  } else if(location == "min"){
    lo <- as.numeric(rep(min(x), n_r))
    ret <- c(lo, x)
  } else if(location == "max"){
    hi <- as.numeric(rep(max(x), n_r))
    ret <- c(x, hi)
  }
  
  ret
}
