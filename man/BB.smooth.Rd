\name{BB.smooth}
\alias{BB.smooth}
\alias{BB.smooth.default}

\title{Pre-smoothing using beta4 models.}

\description{This function fits beta models to score data and provides estimates of 
  the (vector of) score probabilities.
}

\usage{BB.smooth(x,nparm=4,rel)}

\arguments{
  \item{x}{	Data.}
  
  \item{nparm}{	parameters.}
  
  \item{rel}{	reliability.}
  
}

\details{This function fits beta models as described in XXXX, and XXXXX. 
  
  Particular cases of this general equation for each of the equating designs can be found in 
  Von Davier et al (2004) (e.g., Equations (7.1) and (7.2) for the "EG" design, Equation (8.1) for the "SG" design, 
                           Equations (9,1) and (9.2) for the "CB" design).
  
}

\value{
  \item{prob.est}{The estimated score probabilities}
  \item{freq.est}{The estimated score frequencies}
  \item{parameters}{The parameters estimates}
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
  
  beta4nx<-BB.smooth(SEPAx,nparm=4,rel=0) 
  beta4ny<-BB.smooth(SEPAy,nparm=4,rel=0) 
  
  plot(0:50,as.matrix(SEPAx)/sum(as.matrix(SEPAx)),type="b",pch=0, 
       ylim=c(0,0.06),ylab="Relative Frequency",xlab="Scores")
}

\concept{smoothing}

