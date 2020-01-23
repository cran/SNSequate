sim_unimodal <- function(n, x_mean, x_var, N_item, seed=NULL, name=NULL){
  data_names <- c("GANA", "MACAA", "TQS8", "WM8", "WMI")
  
  # Data from Keats & Lord (??)
  kl_mean = c(27.06, 27.06, 32.93, 32.93, 25.82, 25.82, 6.76, 6.75, 23.75, 23.88)
  kl_var = c(8.19, 8.19, 8.04, 8.04, 7.28, 7.28, 5.12, 5.11, 5.59, 5.62)^2
  kl_N_item = c(rep(40,2),rep(50,4),rep(30,4))
  kl_n = c(2354, 2000, 6103, 6103, 2000, 1800, 1000, 1200, 1000, 1200)
  
#  require(emdbook)
  
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


