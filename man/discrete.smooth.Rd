\name{discrete.smooth}
\alias{discrete.smooth}
\alias{discrete.smooth.default}

\title{Pre-smoothing using discrete kernels.}

\description{This function fits discrete kernels to score data and provides estimates of 
  the (vector of) score probabilities.
}

\usage{discrete.smooth(scores,kert,h,x)}

\arguments{
  \item{scores}{	Data.}
  
  \item{kert}{	kernel type.}
  
  \item{h}{	bandwidth.}
  
  \item{x}{ The points of the grid at which the density is to be estimated.}
  
}

\details{This function fits discrete kernels as described in XXXX, and XXXXX. 
  
  Particular cases of this general equation for each of the equating designs can be found in 
  Von Davier et al (2004) (e.g., Equations (7.1) and (7.2) for the "EG" design, Equation (8.1) for the "SG" design, 
                           Equations (9,1) and (9.2) for the "CB" design).
  
}

\value{
  \item{prob.est}{The estimated score probabilities}
 }
\references{
  Gonzalez, J. (2014). SNSequate: Standard and Nonstandard Statistical Models and Methods for Test 
  Equating. \emph{Journal of Statistical Software, 59(7),} 1-30.
  
  Holland, P. and Thayer, D. (1987). Notes on the use of loglinear models for fitting 
  discrete probability distributions. Research Report 87-31, Princeton NJ: Educational 
  Testing Service.
  
  Von Davier, A., Holland, P., and Thayer, D. (2004). \emph{The Kernel Method of Test Equating}. 
  New York, NY: Springer-Verlag.
  
  [1]	Moses, T. "Paper SA06_05 Using PROC GENMOD for Loglinear Smoothing Tim Moses and Alina A. von Davier, Educational Testing Service, Princeton, NJ".
  
}

\author{Jorge Gonzalez \email{jorge.gonzalez@mat.uc.cl}}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{\code{\link{glm}}, \code{\link{ker.eq}}}

\examples{
  data("SEPA", package = "SNSequate")
  
  # create score frequency distributions using freqtab from package equate
  library(equate)
  
  SEPAx<-freqtab(x=SEPA$xscores,scales=0:50)
  SEPAy<-freqtab(x=SEPA$yscores,scales=0:50)
  
  psxB<-discrete.smooth(scores=rep(0:50,SEPAx),kert="bino",h=0.25,x=0:50)
  psxT<-discrete.smooth(scores=rep(0:50,SEPAx),kert="triang",h=0.25,x=0:50)
  psxD<-discrete.smooth(scores=rep(0:50,SEPAx),kert="dirDU",h=0.0,x=0:50)

  plot(0:50,as.matrix(SEPAx)/sum(as.matrix(SEPAx)),lwd=2.0,xlab="Scores", 
  ylab="Relative    Frequency",type="h")
  points(0:50,psxB$prob.est,type="b",pch=0)
  points(0:50,psxT$prob.est,type="b",pch=1)

}

\concept{smoothing}

