\name{irt.eq}
\alias{irt.eq}
\alias{irt.eq.default}

\title{ IRT methods for Test Equating }

\description{ Implements methods to perform Test Equating over IRT models. }

\usage{
irt.eq(n_items, param_x, param_y, theta_points=NULL, weights=NULL, n_points=10, w=1, 
      A=NULL, B=NULL, link=NULL, method_link=NULL, common=NULL,  method="TS", D=1.7)
}

\arguments{  
  \item{n_items}{ Number of items of the test }
  
  \item{param_x}{ Estimated parameters for IRT model on test X. This list must have the following structure:
   list(a, b, c), where each parameter is a vector with the respective estimate for each subject. If you    want to perform other models (i.e. Rasch), replace according with a vector of zeros. 
  }
  
  \item{param_y}{ Estimated parameters for IRT model on test Y. This list must have the following structure:
   list(a, b, c), where each parameter is a vector with the respective estimate for each subject. If you    want to perform other models (i.e. Rasch), replace according with a vector of zeros. 
  }
  
  \item{method}{ A string, either "TS" or "OS". Each one stands for "True Score Equating" and "Observed score equating". Notice that OS requires the additional arguments "theta_points" and "weigths".  }
  
  \item{theta_points}{ For "OS" only. Points over a grid of possible values of  \eqn{\theta} to integrate out the ability term. }
  
  \item{weights}{ For "OS" only. Weigths for integrate out the ability term. If is NULL, the method assumes the distribution of ability is characterized by a finite number of abilities (Kolen and Brennan 2013, pg 199).}
  
  \item{n_points}{ In case theta_ponints is not provided, is the length of the grid for the gaussian quadrature.  }
  
  \item{A, B}{ Scaling parameters. In the case they are not provided, they will be calculated depending on the next described inputs.}
  
  \item{link}{ An \link{irt.link} object.  }
  
  \item{method_link}{ Method used to estimate A and B. Default is "mean/sigma". Others are "mean/mean", "Haebara" and "Stocklord". For more information see \link{irt.link} }
  
  \item{common}{ Common items to estimate A and B. Default asume all items are common.  }
  
  \item{w}{ Weight of the synthetic population.  }
  
  \item{D}{ Sclaing constant}
  
}

\details{

This function implements two methods to perform Test Equating over Item Response Theory models (Kolen and Brennan 2013). 

"True Score Equating" relate number-correct scores on Form X and Form Y. Assumes that the true score associated with each \eqn{\theta} is equivalent to the true score on another form associated with that \eqn{\theta}.

"Observed Score Equating" uses the IRT model to produce an estimated distribution of observed number-correct scores on each form. Using the compound binomial distribution (Lord and Wingersky 1984) to find the conditional distributions \eqn{f(x\mid\theta)}, and then integrate out the \eqn{\theta} parameter. Afterwards, an Equipercentile Equating process is done over the estimated distributions. 

}


\value{	

An object of the clas \code{irt.eq} is returned. Depending on the method used, the outputs are:

\describe{
    \item{True Score Equating}{ A list(n_items, theta_equivalent, tau_y) containing the number of items, the theta equivalent values on Form X to Form Y and the equivalent scores.}
    
    \item{Observed Score Equating}{A list(n_items, f_hat, g_hat, e_Y_x) containing the number of items, the estimated distributions and the equated values.}
}

}

\references{
    Kolen, M. J., and Brennan, R. L. (2014). \emph{Test Equating, Scaling, and Linking: Methods and Practices, Third Edition}. Springer Science & Business Media.
}

\author{Daniel Acuna Leon. \email{dnacuna@uc.cl}}

\seealso{ 
  \code{\link{irt.link}} 
}

\examples{
data(KB36_t)
dfo <- KB36_t

param_x <- list(a=dfo[,3],b=dfo[,4],c=dfo[,5])
param_y <- list(a=dfo[,7],b=dfo[,8],c=dfo[,9])

theta_points=c(-5.2086,-4.163,-3.1175,-2.072,-1.0269,0.0184,
               1.0635,2.109,3.1546,4.2001)
weights=c(0.000101,0.00276,0.03021,0.142,0.3149,0.3158,
         0.1542,0.03596,0.003925,0.000186)


irt.eq(36, param_x, param_y, method="TS", A=1, B=0)
irt.eq(36, param_x, param_y, theta_points, weights, method="OS", A=1, B=0)
}

\concept{kernel equating}
