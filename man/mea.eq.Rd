\name{mea.eq}
\alias{mea.eq}
\alias{mea.eq.default}	
\title{The mean method of equating}

\description{This function implements the mean method of test equating  as described in Kolen and Brennan (2004). 
}

\usage{
mea.eq(sx, sy, scale)
}

\arguments{
  \item{sx}{		A vector containing the observed scores of the sample taking test \eqn{X}{X}.}

  \item{sy}{		A vector containing the observed scores of the sample taking test \eqn{Y}{Y}.}

  \item{scale}{		Either an integer or vector containing the values on the scale to be equated.}

}

\details{The function implements the mean method of equating as described in Kolen and Brennan (2004). Given observed scores
\eqn{sx}{sx} and \eqn{sy}{sy}, the functions calculates 
	\deqn{\varphi(x;\mu_x,\mu_y)=x-\mu_x+\mu_y}{phi(x;mu_x,mu_y)=x-mu_x+mu_y}
	where \eqn{\mu_x}{mu_x} and \eqn{\mu_y}{mu_y} are the score means on test \eqn{X}{X} and \eqn{Y}{Y}, respectively.
}

\value{A two column matrix with the values of \eqn{\varphi()}{phi()} (second column) for each scale value \code{x} (first column)}

\references{

        Gonzalez, J. (2014). SNSequate: Standard and Nonstandard Statistical Models and Methods for Test 
        Equating. \emph{Journal of Statistical Software, 59(7),} 1-30.

	Kolen, M., and Brennan, R. (2004). \emph{Test Equating, Scaling and Linking}. 
	New York, NY: Springer-Verlag.
}

\author{Jorge Gonzalez \email{jorge.gonzalez@mat.uc.cl}}

\seealso{\code{\link{lin.eq}}, \code{\link{eqp.eq}}, \code{\link{ker.eq}}, \code{\link{le.eq}}
}

\examples{
#Artificial data for two two 100 item tests forms and 5 individuals in each group
x1<-c(67,70,77,79,65,74)
y1<-c(77,75,73,89,68,80)

#Score means
mean(x1); mean(y1)

#An equivalent form y1 score of 72 on form x1
mea.eq(x1,y1,72)

#Equivalent form y1 score for the whole scale range
mea.eq(x1,y1,0:100)
}

\concept{Traditional equating methods}
