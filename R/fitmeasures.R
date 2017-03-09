### fitmeasures.R                   
### A collection of measures to investigate goodness of fit
###
### Copyright: Jorge Gonzalez, 2012.
### Last modification: 25-05-2012.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### Author contact information:
###
###      Jorge Gonzalez B.
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3545467  URL  : http://www.mat.puc.cl/~jgonzale
###      Fax  : +56-2-3547729  Email: jgonzale@mat.puc.cl
###

gof <- function(obs, fit, methods=c("FT"), p.out=FALSE){
  retvals <- list()
  
  if("FT" %in% methods){
    ftres <- ft.res(obs, fit, p.out=p.out)
    retvals <- c(retvals, list(ft.res=ftres))
  }
  
  if("KL" %in% methods){
    kl <- kl.divergence(obs, fit) + kl.divergence(fit, obs)
    retvals <- c(retvals, list(kl.div=kl))
  }
  
  if("Chisq" %in% methods){
    chisq <- chisq.test(obs, fit)
    chisq <- list(statistic=chisq$stat, df=chisq$parameter, p.value=chisq$p.value)
    retvals <- c(retvals, list(chisq=chisq))
  }
  
  retvals$methods <- methods
  class(retvals) <- "snse.gof"
  retvals
}

print.snse.gof <- function(x, ...){
  if("FT" %in% x$methods){
    cat("Freeman-Tukey Residuals:\n")
    cat("------------------------\n")
    print(x$ft.res)
    cat("\n")
  }
  
  if("KL" %in% x$methods){
    cat("Symmetrised Kullback-Leibler divergence:\n")
    cat("----------------------------------------\n")
    print(x$kl.div)
    cat("\n")
  }
  
  if("Chisq" %in% x$methods){
    cat("Pearson's Chi-squared Test:\n")
    cat("---------------------------\n")
    cat("H0: The data is consistent with the specified distribution\n\n")
    cat(paste0("X-squared = ",x$chisq$stat, ", df = ",x$chisq$df, ", p-value = ",round(x$chisq$p.value, 3)))
    cat("\n")
  }
}


ft.res <- function(obs, fit, p.out=TRUE){
  FT <- as.numeric(sqrt(obs) + sqrt(obs+1) - sqrt(4*fit+1))
  
  if(p.out){
    par(mfrow=c(1,1)) # Reset the plot grid
    qqnorm(FT, main="Freeman-Tukey Residuals QQPlot")
    qqline(FT, col=2)
  }
  
  return(FT)
}


kl.divergence <- function(p, q){
  pn <- p / sum(p)
  qn <- q / sum(q)
  return( sum(pn*(log(pn) - log(qn)))  )
}
