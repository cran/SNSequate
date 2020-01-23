\name{lin.eq}
\alias{lin.eq}
\alias{lin.eq.default}
\title{The linear method of equating}

\description{This function implements the linear method of test equating  as described in Kolen and Brennan (2004). 
}

\usage{
lin.eq(sx, sy, scale)
}

\arguments{
  \item{sx}{		A vector containing the observed scores of the sample taking test \eqn{X}{X}.}

  \item{sy}{		A vector containing the observed scores of the sample taking test \eqn{Y}{Y}.}

  \item{scale}{		Either an integer or vector containing the values on the scale to be equated.}
}

\details{The function implements the linear method of equating as described in Kolen and Brennan (2004). Given observed scores
\eqn{sx}{sx} and \eqn{sy}{sy}, the functions calculates 
	\deqn{\varphi(x;\mu_x,\mu_y,\sigma_x,\sigma_y)=\frac{\sigma_x}{\sigma_y}(x-\mu_x)+\mu_y}{phi(x;mu_x,mu_y,sigma_x,\sigma_y)=(sigma_x/sigma_y)*(x-mu_x)+mu_y}
	where \eqn{\mu_x,\mu_y,\sigma_x,\sigma_y}{mu_y,mu_y,sigma_x,sigma_y} are the score means and standard deviations on test \eqn{X}{X} and \eqn{Y}{Y}, 
	respectively.
}

\value{A two column matrix with the values of \eqn{\varphi()}{phi()} (second column) for each scale value \code{x} (first column)}

\references{

        Gonzalez, J. (2014). SNSequate: Standard and Nonstandard Statistical Models and Methods for Test 
        Equating. \emph{Journal of Statistical Software, 59(7),} 1-30.

	Kolen, M., and Brennan, R. (2004). \emph{Test Equating, Scaling and Linking}. 
	New York, NY: Springer-Verlag.
}

\author{Jorge Gonzalez \email{jorge.gonzalez@mat.uc.cl}}

\seealso{\code{\link{mea.eq}}, \code{\link{eqp.eq}}, \code{\link{ker.eq}}
}

\examples{
#Artificial data for two two 100 item tests forms and 5 individuals in each group
x1<-c(67,70,77,79,65,74)
y1<-c(77,75,73,89,68,80)

#Score means and sd
mean(x1); mean(y1)
sd(x1); sd(y1)

#An equivalent form y1 score of 72 on form x1
lin.eq(x1,y1,72)

#Equivalent form y1 score for the whole scale range
lin.eq(x1,y1,0:100)

#A plot comparing mean, linear and identity equating
plot(0:100,0:100, type='l', xlim=c(-20,100),ylim=c(0,100),lwd=2.0,lty=1,
ylab="Form Y raw score",xlab="Form X raw score")
abline(a=5,b=1,lwd=2,lty=2)
abline(a=mean(y1)-(sd(y1)/sd(x1))*mean(x1),b=sd(y1)/sd(x1),,lwd=2,lty=3)
arrows(72, 0, 72, 77,length = 0.15,code=2,angle=20)
arrows(72, 77, -20, 77,length = 0.15,code=2,angle=20)
abline(v=0,lty=2)
legend("bottomright",lty=c(1,2,3), c("Identity","Mean","Linear"),lwd=c(2,2,2))

}

\concept{Traditional equating methods}
