\name{le.eq}
\alias{le.eq}
\alias{le.eq.default}

\title{Local equating methods}

\description{This function implements the local method of equating as descibed in van der Linden (2011). }

\usage{
le.eq(S.X, It.X, It.Y, Theta)
}

\arguments{
  \item{S.X}{		A vector containing the observed scores of the sample taking test \eqn{X}{X}.}

  \item{It.X}{		A matrix of item parameter estimates coming from an IRT model for test form \eqn{X}{X} 
			(difficulty, discrimation and guessing parameters are located in the first, second and 
			third column, respectively).}

  \item{It.Y}{		A matrix of item parameter estimates coming from an IRT model for test form \eqn{Y}{Y}.}

  \item{Theta}{		Either a number or vector of values representing the value of \code{theta} where to 
			condition on (see details)}
}

\details{The function implements the local equating method as described in van der Linden (2011). Based on 
Lord (1980) principle of equity, local equating methods utilizes the conditional on abilities distributions 
of scores to obtain the transformation \eqn{\varphi}{phi}. The method leads to a family of transformations 
of the form
\deqn{\varphi(x;\theta)=G_{Y\mid\theta}^{-1}(F_{X\mid\theta}(x)),\quad \theta\in\mathcal{R}}{phi(x;theta)=G_{Y|theta}^{-1}(F_{X|theta}(x)), theta in R}
The conditional distributions of \eqn{X}{X} and \eqn{Y}{Y} are obtained using the algorithm described by 
Lord and Wingersky (1984). Among other possibilities, a value for \eqn{\theta}{theta} can be a EAP, ML or MAP estimation of it, for and underlying 
IRT model (for example, using the \code{ltm} R package (Rizopoulos, 2006)). 
}


\value{A list containing the observed scores to be equated, the corresponding ability estimates where to condition on, and the equated values}

\references{

	van der Linden, W. (2011). \emph{Test Equating, Scaling and Linking}. 
	New York, NY: Springer-Verlag.

	Lord, F. (1980). \emph{Applications of Item Response Theory to Practical Testing Problems}. 
	Lawrence Erlbaum Associates, Hillsdale, NJ.

	Lord, F. and Wingersky, M. (1984). Comparison of IRT True-Score and Equipercentile Observed-Score Equatings. 
	\emph{Applied Psychological Measurement,8(4),} 453--461.

	Rizopoulos, D. (2006). ltm: An R package for latent variable modeling and item response theory analyses. 
	\emph{Journal of Statistical Software, 17(5),} 1--25.
}

\author{Jorge Gonzalez B. \email{jgonzale@mat.puc.cl}}

\seealso{\code{\link{mea.eq}}, \code{\link{eqp.eq}}, \code{\link{lin.eq}} \code{\link{ker.eq}}
}

\examples{
## Artificial data for two 5-items tests forms. Both forms are assumed
## being fitted by a 3PL model.

## Create (artificial) item parameters matrices for test form X and Y
ai<-c(1,0.8,1.2,1.1,0.9)
bi<-c(-2,-1,0,1,2)
ci<-c(0.1,0.15,0.05,0.1,0.2)
itx<-rbind(bi,ai,ci)
ai<-c(0.5,1.4,1.2,0.8,1)
bi<-c(-1,-0.5,1,1.5,0)
ci<-c(0.1,0.2,0.1,0.15,0.1)
ity<-rbind(bi,ai,ci)

#Two individuals with different ability (1 and 2) obtain the same score 2.
#Their corresponding equated scores values are:
le.eq(c(2,2),itx,ity,c(1,2))

}

\keyword{local equating}
