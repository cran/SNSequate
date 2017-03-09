\name{eqp.eq}
\alias{eqp.eq}
\alias{eqp.eq.default}

\title{The equipercentile method of equating}

\description{This function implements the equipercentile method of test equating  as described in Kolen and Brennan (2004). 
}

\usage{
eqp.eq(sx, sy, X, Ky = max(sy))
}

\arguments{
  \item{sx}{		A vector containing the observed scores on test \eqn{X}{X}}

  \item{sy}{		A vector containing the observed scores on test \eqn{Y}{Y}}

  \item{X}{		Either an integer or vector containing the values on the scale 
			to be equated.}

  \item{Ky}{		The total number of items in test form \eqn{Y}{Y} to which form 
			\eqn{X}{X} scores will be equated}
}

\details{The function implements the equipercentile method of equating as described in Kolen and Brennan (2004). Given observed scores
\code{sx} and \code{sy}, the functions calculates 
	\deqn{\varphi(x)=G^{-1}(F(x))}{G^{-1}(F(x))}
	where \eqn{F}{F} and \eqn{G}{G} are the cdf of scores on test forms \eqn{X}{X} and \eqn{Y}{Y}, 
	respectively.
}

\value{A two column matrix with the values of \eqn{\varphi()}{phi()} (second column) for each scale value \code{x} (first column)}

\references{

        Gonzalez, J. (2014). SNSequate: Standard and Nonstandard Statistical Models and Methods for Test 
        Equating. \emph{Journal of Statistical Software, 59(7),} 1-30.

	Kolen, M., and Brennan, R. (2004). \emph{Test Equating, Scaling and Linking}. 
	New York, NY: Springer-Verlag.
}

\author{Jorge Gonzalez B. <jgonzale@mat.puc.cl>}

\seealso{\code{\link{mea.eq}}, \code{\link{lin.eq}}, \code{\link{ker.eq}}
}

\examples{
### Example from Kolen and Brennan (2004), pages 41-42:
### (scores distributions have been transformed to vectors of scores)

sx<-c(0,0,1,1,1,2,2,3,3,4)
sy<-c(0,1,1,2,2,3,3,3,4,4)
x<-2
eqp.eq(sx,sy,2)

# Whole scale range (Table 2.3 in KB)
eqp.eq(sx,sy,0:4)

}

\keyword{Traditional equating methods}

