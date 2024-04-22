irt.eq <- function(n_items, param_x, param_y, theta_points=NULL, weights=NULL, n_points=10, w=1, 
                   A=NULL, B=NULL, link=NULL, method_link=NULL, common=NULL, method="TS", D=1.7)
UseMethod("irt.eq")

irt.eq.default <- function(n_items, param_x, param_y, theta_points=NULL, weights=NULL, n_points=10, w=1, 
                           A=NULL, B=NULL, link=NULL, method_link=NULL, common=NULL,  method="TS", D=1.7){
  if(method == "TS"){
    return(irteq.tse(n_items, param_x, param_y, D,
                      A, B, link, method_link, common))
  }
  else if(method == "OS"){
    return(irteq.ose(n_items, param_x, param_y, D, theta_points, weights, n_points, w, 
                      A, B, link, method_link, common))
  }
  
}

## Generates a funcion containing P for IRT
## given a fixed set of parameters
get_irt <- function(a,b,c,D){
  function(theta){
    t <- D * a * (theta - b)
    c + (1-c)*exp(t)/(1+exp(t))
  }
}


## Recursive formula from Lord-Wingersky (1984), using memoization.
LWmemo <- function(x, r, theta, irt){
  memo <- matrix(NA, nrow=r, ncol=r+1)
  
  f <- function(x, r, theta, irt){  
    if(!is.na(memo[r,x+1])) return(memo[r,x+1])
    # Base case for r
    if(r == 1){
      if(x == 0)
        memo[r,x+1] <<- 1 - irt[[r]](theta)
      
      if(x == 1)
        memo[r,x+1] <<- irt[[r]](theta)
      
      return(memo[r,x+1])
    }
    
    if(x == 0) {
      memo[r,x+1] <<- f(x, r-1, theta, irt) * (1 - irt[[r]](theta))
    } else if(x == r){
      memo[r,x+1] <<- f(x-1, r-1, theta, irt) * irt[[r]](theta)
    } else {
      memo[r,x+1] <<- f(x, r-1, theta, irt) * (1 - irt[[r]](theta)) + f(x-1, r-1, theta, irt) * irt[[r]](theta)
    }
    
    return(memo[r,x+1])
  }
  
  f(x, r, theta, irt)
}

irteq.tse <-function(n_items, param_x, param_y, D=1.7, A, B, link, method_link, common){
  # TODO: check input
  # Parameters
  a_x <- param_x$a
  b_x <- param_x$b
  c_x <- param_x$c
  D_x <- rep(D, length(a_x))
  
  a_y <- param_y$a
  b_y <- param_y$b
  c_y <- param_y$c
  D_y <- rep(D, length(a_y))
  
  # Number of items 
  scores <- 0:n_items
  
  
  ## Big ol' decision tree for scaling parameters
  ## If A or B are NULL
  if(is.null(A) || is.null(B)){
    
    ## If an irt.link object is given
    if(!is.null(link)){
#      if(class(link) != "irt.link"){
#        stop("Link must be an irt.link object")
#      }
    }
    else{
      ## Common items have not been set
      if(is.null(common)){
        common <- scores
      }
      
      param <- cbind(a_y, b_y, c_y, a_x, b_x, c_x)
      link <- irt.link(param, common=common, icc="logistic", model="3PL", D=D)
    }
    
    ## Finally, set the A and B values
    ## mean/sigma as default
    if(is.null(method_link)){
      A <- link$ms[1]
      B <- link$ms[2]
    } else {
      if(method_link == "mean/mean"){
        A <- link$mm[1]
        B <- link$mm[2]
      } 
      else if(method_link == "mean/sigma"){
        A <- link$ms[1]
        B <- link$ms[2]
      }
      else if(method_link == "Haebara"){
        A <- link$Haebara[1]
        B <- link$Haebara[2]
      }
      else if(method_link == "StockLord"){
        A <- link$StockLord[1]
        B <- link$StockLord[2]
      }
      else{
        stop("You must choose a valid linking method")
      }
    }
  }
  
  # Set of indexes indicating which scores 
  # can be equated
  admissible_scores <- logical(length(scores))
  
  # Admissible range to compute equating
  # \sum c_j < \tau < K, eq (6.19), pg 193
  c_sum <- sum(c_x)
  
  if(c_sum > 0){
    idx <- (ceiling(c_sum):(n_items-1) + 1)
  }else{
    idx <- (1:(n_items-1) + 1)
  }

  admissible_scores[idx] <- TRUE # +1 due to score 0
  
  # A list with all the P's and dP's for each parameter set
  irt_x <- lapply(1:length(a_x), function(j){ get_irt(a_x[j]/A,b_x[j]*A + B,c_x[j],D_x[j]) })
  
  # Theta equivalent vector of values
  theta_equivalent <- rep(NA, length(scores))
  
  # Function to find root for a given tau
  theta_find <- function(tau){
    fn <- function(theta){ 
      p_ij <- lapply(1:length(irt_x), function(j){ irt_x[[j]](theta) })
      return(tau - Reduce(sum, p_ij))  
    }
    
    # Fixed wide range to find root
    uniroot(fn, interval=c(-40,40))$root
  }
  
  # Find theta_equivalent for each admissible score
  theta_equivalent[admissible_scores] <- sapply(scores[admissible_scores], theta_find)
  theta_equivalent
  
  irt_y <- lapply(1:length(a_y), function(j){ get_irt(a_y[j],b_y[j],c_y[j],D_y[j]) })
  tau_y <- 0:n_items
  tau_y[admissible_scores] <- unlist(
    sapply(theta_equivalent[admissible_scores], 
           function(theta){
             Reduce(sum, lapply(1:length(irt_y), function(j){ irt_y[[j]](theta) }))
           }
    )
  )
  
  if(c_sum > 0)
    idx <- 2:ceiling(c_sum)
  else
    idx <- c(2)
  
  pc <- ifelse(c_sum > 0, sum(c_y)/sum(c_x), 1)
  tau_y[idx] <- tau_y[idx] * pc
  
  res = list(n_items=n_items, theta_equivalent=theta_equivalent, tau_y=tau_y, method="TS")
  class(res) <- "irt.eq"
  return(res)
}

irteq.ose <- function(n_items, param_x, param_y, D=1.7, theta_points, weights, n_points, w, 
                       A, B, link, method_link, common){
  
  a_x <- param_x$a
  b_x <- param_x$b
  c_x <- param_x$c
  D_x <- rep(D, length(a_x))
  
  a_y <- param_y$a
  b_y <- param_y$b
  c_y <- param_y$c
  D_y <- rep(D, length(a_y))
  
  # Number of items 
  scores <- 0:n_items
  
  # If theta or weights are not defined, calculate
  # them via Gaussian Quadrature
  if(is.null(theta_points) || is.null(weights)){
  #  require(statmod)
    gq <- gauss.quad.prob(n_points, dist="normal")
    theta_points <- gq$nodes
    weights <- gq$weights
  }
  
  ## Big ol' decision tree for scaling parameters
  ## If A or B are NULL
  if(is.null(A) || is.null(B)){
    
    ## If an irt.link object is given
    if(!is.null(link)){
#      if(class(link) != "irt.link"){
#        stop("Link must be an irt.link object")
#      }
    }
    else{
      ## Common items have not been set
      if(is.null(common)){
        common <- scores
      }
      
      param <- cbind(a_y, b_y, c_y, a_x, b_x, c_x)
      link <- irt.link(param, common=common, icc="logistic", model="3PL", D=D)
    }

    ## Finally, set the A and B values
    ## mean/sigma as default
    if(is.null(method_link)){
      A <- link$ms[1]
      B <- link$ms[2]
    } else {
      if(method_link == "mean/mean"){
        A <- link$mm[1]
        B <- link$mm[2]
      } 
      else if(method_link == "mean/sigma"){
        A <- link$ms[1]
        B <- link$ms[2]
      }
      else if(method_link == "Haebara"){
        A <- link$Haebara[1]
        B <- link$Haebara[2]
      }
      else if(method_link == "StockLord"){
        A <- link$StockLord[1]
        B <- link$StockLord[2]
      }
      else{
        stop("You must choose a valid linking method")
      }
    }
  }
  ## In the "else" of the first case above, the A and B are set so we do nothing
  
  
  #### Population 1
  # A list with all the P's and dP's for each parameter set
  irt_x <- lapply(1:length(a_x), function(j){ get_irt(a_x[j]/A,b_x[j]*A + B,c_x[j],D_x[j]) })
  irt_y <- lapply(1:length(a_y), function(j){ get_irt(a_y[j],b_y[j],c_y[j],D_y[j]) })
  
  xt_1 <- sapply(theta_points, 
               function(theta){ sapply(0:n_items, function(x){LWmemo(x, n_items, theta, irt_x)}) }
              )
  tmp <- LWmemo(0, n_items, theta_points[length(theta_points)], irt_y)
  yt_1 <- sapply(theta_points, 
               function(theta){ sapply(0:n_items, function(y){LWmemo(y, n_items, theta, irt_y)}) }
              )
  
  #### Population 2
  xt_2 <- sapply(theta_points*A + B, 
               function(theta){ sapply(0:n_items, function(x){LWmemo(x, n_items, theta, irt_x)}) }
  )
  tmp <- LWmemo(0, n_items, theta_points[length(theta_points)], irt_y)
  yt_2 <- sapply(theta_points*A + B, 
               function(theta){ sapply(0:n_items, function(y){LWmemo(y, n_items, theta, irt_y)}) }
  )
  
  if(is.null(weights))
    weights <- rep(1/length(theta_points), length(theta_points))
  
  f1 <- as.numeric(xt_1 %*% weights)
  f2 <- as.numeric(xt_2 %*% weights)
  
  g1 <- as.numeric(yt_1 %*% weights)
  g2 <- as.numeric(yt_2 %*% weights)
  
  # Synthetic population
  f_hat <- f1*w + f2*(1-w)
  g_hat <- g1*w + f2*(1-w)
  
  F_hat <- cumsum(f_hat)
  G_hat <- cumsum(g_hat)
  
  e_Y_x <- eqp_eq(0:n_items, 0:n_items, F_hat, G_hat)
  
  res = list(n_items=n_items, f_hat=f_hat, g_hat=g_hat, e_Y_x=e_Y_x, method="OS")
  class(res) <- "irt.eq"
  return(res)
}

print.irt.eq <- function(x,digits=3,...) {
  method <- ifelse(x$method == "TS", "True Score", "Observed Score")
  
  cat(paste("\nEquated values using IRT", method ,"method:\n"))	
  
  if(x$method=="TS"){
    df <- data.frame(Scores=0:x$n_items, theta_eqi=x$theta_equivalent, tau_y=x$tau_y)
  }
  else if(x$method=="OS"){
    df <- data.frame(Scores=0:x$n_items, f_hat=x$f_hat, g_hat=x$g_hat, eYx=x$e_Y_x)
  }
  
#  if(suppressMessages(suppressWarnings(require(knitr)))){
    print(kable(round(df, digits=digits), padding=2, format="pandoc"), row.names=F)
#  }
#  else{
#    print(round(df, digits=digits), row.names=F)
#  }
}	

