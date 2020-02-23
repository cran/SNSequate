\name{SEED}
\alias{SEED}
\alias{SEED.default}
\title{Standard error of equating difference
}
\description{This function calculates the standard error of equating diference (SEED) as described in 
Von Davier et al. (2004). 
}
\usage{
SEED(eq1, eq2)
}

\arguments{
\item{eq1}{	An object of class \code{ker.eq} which contains one of the two estimated equated functions to be used 
		for the SEED.}

\item{eq2}{	An object of class \code{ker.eq} which contains one of the two estimated equated functions to be used 
		for the SEED.}

}

\details{The SEED can be used as a measure to choose whether to support or not a certain equating function on another 
another one. For instance, when \eqn{h_X}{h_X} and \eqn{h_Y}{h_Y} tends to infinity, then the (gaussian kernel) 
	\eqn{\hat{e}_Y(x)}{\hat{e}_Y(x)} equating function tends to the linear equating function 
	(see Theorem 4.5 in Von Davier et al, 2004 for more details). Thus, one can calculate the measure
	\deqn{SEED_Y(x)=\sqrt{Var(\hat{e}_Y(x)-\widehat{Lin}_Y(x))}}{SEED_Y(x)={Var(\hat{e}_Y(x)-\widehat{Lin}_Y(x))}^{1/2}}
	to decide between \eqn{\hat{e}_Y(x)}{\hat{e}_Y(x)} and \eqn{\widehat{Lin}_Y(x)}{\widehat{Lin}_Y(x)}.
} 


\value{A two column matrix with the values of \code{SEEYx} for each \code{x} in the first column and the values of 
\code{SEEXy} for each \code{y} in the second column}

\references{

        Gonzalez, J. (2014). SNSequate: Standard and Nonstandard Statistical Models and Methods for Test 
        Equating. \emph{Journal of Statistical Software, 59(7),} 1-30.

	Von Davier, A., Holland, P., and Thayer, D. (2004). \emph{The Kernel Method of Test Equating}. 
	New York, NY: Springer-Verlag.
}
\author{Jorge Gonzalez \email{jorge.gonzalez@mat.uc.cl}
}

\seealso{\code{\link{ker.eq}}
}

\examples{
#Example: Figure7.7 in Von Davier et al, (2004)
data(Math20EG)

mod.gauss<-ker.eq(scores=Math20EG,kert="gauss", hx = NULL, hy = NULL,degree=c(2, 3),design="EG")
mod.linear<-ker.eq(scores=Math20EG,kert="gauss", hx = 20, hy = 20,degree=c(2, 3),design="EG")

Rx<-mod.gauss$eqYx-mod.linear$eqYx
seed<-SEED(mod.gauss,mod.linear)$SEEDYx

plot(0:20,Rx,ylim=c(-0.8,0.8),pch=15)
abline(h=0)
points(0:20,2*seed,pch=0)
points(0:20,-2*seed,pch=0)

#Example Figure 10.4 in Von Davier (2011)
mod.unif<-ker.eq(scores=Math20EG,kert="unif", hx = NULL, hy = NULL,degree=c(2, 3),design="EG")
mod.logis<-ker.eq(scores=Math20EG,kert="logis", hx = NULL, hy = NULL,degree=c(2, 3),design="EG")

Rx1<-mod.logis$eqYx-mod.gauss$eqYx
Rx2<-mod.unif$eqYx-mod.gauss$eqYx

seed1<-SEED(mod.logis,mod.gauss)$SEEDYx
seed2<-SEED(mod.unif,mod.gauss)$SEEDYx

plot(0:20,Rx1,ylim=c(-0.2,0.2),pch=15,main="LK vs GK",ylab="",xlab="Scores")
abline(h=0)
points(0:20,2*seed1,pch=0)
points(0:20,-2*seed1,pch=0)

plot(0:20,Rx2,ylim=c(-0.2,0.2),pch=15,main="UK vs GK",ylab="",xlab="Scores")
abline(h=0)
points(0:20,2*seed2,pch=0)
points(0:20,-2*seed2,pch=0)

}
\concept{kernel equating}

