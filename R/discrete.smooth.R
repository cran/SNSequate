###
### This is the discrete.smooth() function which performs presmoothing of 
### score distributions using nonparamenterc discrete kernels. 

### The function makes a call to other functions from the Ake package. 
### There is a problem with the Ake package which was already reported to the
### author. I believe the reason why Dirac is giving NAs is because the diracDU 
### function in file kef.R is returning no output. 
### In despite of this, in what follows I compile the discrete.smooth() funtion
### while waiting for the Ake authors to chek their package. 


discrete.smooth<-function(scores,kert,h,x=NULL) 
  UseMethod("discrete.smooth")

#library(Ake)
discrete.smooth.default<-function(scores,kert,h,x=NULL){
  
  cl <- match.call()
  
  sp.est<-kpmfe.fun(Vec=scores,h=h,type_data="discrete",ker=kert,x=x)$est.fn
  res<-list(call=cl,prob.est=sp.est)
  class(res)<-"discrete.smooth"
  return(res)
}


#### Print function

print.discrete.smooth<-function(x,...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nEstimated score probabilities:\n")
  cat("\n")
  print(data.frame(Est.Score.Prob.=x$prob.est))
  #cat("\nParameters:\n")
  #cat("\n")
  #print(x$parameters)
  
}
