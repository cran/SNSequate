\name{gof}
\alias{gof}
\title{Functions to assess model fitting.}

\description{This function contains various measures to assess the model's goodness of fit. }

\usage{
gof(obs, fit, methods=c("FT"), p.out=FALSE)
}

\arguments{
  \item{obs}{		A vector containing the observed values.}

  \item{fit}{		A vector containing the fitted values.}

  \item{methods}{		A character vector containing one or many of the following methods:
    \describe{
      \item{"FT"}{Freeman-Tukey Residuals. This is the default test.}
      \item{"Chisq"}{Pearson's Chi-squared test.}
      \item{"KL"}{Symmetrised Kullback-Leibler divergence [1].}
    }
  }
  
  \item{p.out}{  Boolean. Decides whether or not to display plots (on corresponding methods).}
}

\references{

        Gonzalez, J. (2014). SNSequate: Standard and Nonstandard Statistical Models and Methods for Test 
        Equating. \emph{Journal of Statistical Software, 59(7),} 1-30.

      	Kolen, M., and Brennan, R. (2004). \emph{Test Equating, Scaling and Linking}. 
      	New York, NY: Springer-Verlag.
      	
      	[1] \url{https://en.wikipedia.org/wiki/Kullback-Leibler_divergence}
}

\author{Daniel Leon A. \email{dnacuna@mat.uc.cl}}

\examples{
data(Math20EG)
mod <- ker.eq(scores=Math20EG,kert="gauss",degree=c(2,3),design="EG")

gof(Math20EG[,1], mod$rj*mod$nx, method=c("FT", "KL"))

}

\keyword{Traditional equating methods}
