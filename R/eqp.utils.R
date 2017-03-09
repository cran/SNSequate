## These are functions called by eqp.eq() to perform equipercentile equating

F.x<-function(x,S){
  max(0,min(1,mean(S<=x)))
}
 
P.x<-function(x,S){
  x.s<-floor(x+.5)
  max(0,min(100, 100*(F.x(x.s-1,S)+(x-(x.s-.5))*(F.x(x.s,S)-F.x(x.s-1,S)))))
}

x.s.up<-function(p.s,S){
  val<-S[max(1,floor((p.s/100)*length(S)))]
  while (F.x(val,S)<=(p.s/100)){
    val<-val+1
}
  val
}

x.s.lo<-function(p.s,S){
  val<-S[min(length(S),ceiling((p.s/100)*length(S)))]
  while (F.x(val,S)>=(p.s/100)){
    val<-val-1
}
  val
}

P.inv<-function(p.s,S,Kx=max(S)){
  if (p.s==100){x.up.ps<-Kx+.5}
  else{
    x.s.u<-x.s.up(p.s,S)
    x.up.ps<-(p.s/100-F.x(x.s.u-1,S))/(F.x(x.s.u,S)-F.x(x.s.u-1,S))+(x.s.u-.5)
  }
  if (p.s==0){x.lo.ps<--.5}
  else{
    x.s.l<-x.s.lo(p.s,S)
    x.lo.ps<-(p.s/100-F.x(x.s.l,S))/(F.x(x.s.l+1,S)-F.x(x.s.l,S))+(x.s.l+.5)
  }
  mean(c(x.lo.ps,x.up.ps))
}

####################################
##  Reimplemented some methods
##  so it does not break the working 
##  dependant functions.
##  Needs posterior Refactoring
####################################
P_x <- function(x, F) {
  x_s <- floor(x + 0.5)
  F_x <- function(x){
            if(x < 0){
              return(0)
            } else if(x > length(F)){
              return(1)
            }

            return(F[x+1]) # Because index start at 1
         }
  
  t <- 100*(F_x(x_s - 1) + (x_s - (x_s - 0.5)) * (F_x(x_s) - F_x(x_s - 1)))
  
  return(max(0, min(100, t)))
}


P_inv_p <- function(p, F){
  K_x <- length(F) - 1
  F_x <- function(x){
            if(x < 0){
              return(0)
            } else if(x > length(F)){
              return(1)
            }
            
            return(F[x+1]) # Because index start at 1
          }
  
  
  if(p >= 100){ 
    x_u <- K_x + 0.5 
  } else { 
    x_u_s <- length(F) - 1
    if(length(idx <- which(F*100 > p)))
      x_u_s <- min(idx) - 1
    x_u <- (p/100 - F_x(x_u_s - 1)) / (F_x(x_u_s) - F_x(x_u_s - 1)) + (x_u_s - 0.5)
  }
  
  if(p <= 0){ 
    x_l <- -0.5 
  } else{
    x_l_s <- 0
    if(length(idx <- which(F*100 < p))) 
      x_l_s <- max(idx) - 1 
    x_l <- (p/100 - F_x(x_l_s)) / (F_x(x_l_s + 1) - F_x(x_l_s)) + (x_l_s + 0.5)
  }
  
  return(mean(c(x_u, x_l)))
}

eqp_eq <- function(x, y, F, G){
  p_x <- sapply(x, P_x, F)
  return(sapply(p_x, P_inv_p, G))
}
  
  
  




