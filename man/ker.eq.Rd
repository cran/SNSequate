\name{ker.eq}
\alias{ker.eq}
\alias{ker.eq.default}

\title{The Kernel method of test equating}

\description{This function implements the kernel method of test equating  as described in Holland and Thayer (1989), 
and Von Davier et al. (2004). Nonstandard kernels others than the gaussian are available. Associated standard error 
of equating are also provided. 
}

\usage{
ker.eq(scores, kert, hx = NULL, hy = NULL, degree, design, Kp = 1, scores2, 
degreeXA, degreeYA, J, K, L, wx, wy, w, gapsX, gapsY, gapsA, lumpX, lumpY, 
lumpA, alpha, h.adap)
}

\arguments{  
  {		Note that depending on the specified equating design, not all arguments are necessary 
			as detailed below.} 


  \item{scores}{	If the "EG" design is specified, a two column matrix containing the raw sample frequencies 
			coming from the two groups of scores to be equated. It is assumed that the data in the first 
			and second columns come from tests \eqn{X}{X} and \eqn{Y}{Y}, respectively.
			
			If the "SG" design is specified, a matrix containing the (joint) bivariate sample 
			frequencies for \eqn{X}{X} (raws) and \eqn{Y}{Y} (columns).
			
			If the "CB" design is specified, a two column matrix containing the observed scores 
			of the sample taking test \eqn{X}{X} first, followed by test \eqn{Y}{Y}. The \code{scores2} 
			argument is then used for the scores of the sample taking test Y first followed by test 
			\eqn{X}{X}.
			
			If either the "NEAT_CB" or "NEAT_PSE" design is selected, a two column matrix containing 
			the observed scores on	test \eqn{X}{X} (first column) and the observed scores on the anchor 
			test \eqn{A}{A} (second column). The \code{scores2} argument is then used for the observed 
			scores on test \eqn{Y}{Y}.}

  \item{kert}{		A character string giving the type of kernel to be used for continuization. 
			Current options include "\code{gauss}", "\code{logis}", "\code{uniform}", "\code{epan}" and "\code{adap}" 
			for the gaussian, logistic, uniform, Epanechnikov and Adaptative kernels, respectively}


  \item{hx}{		An integer indicating the value of the bandwidth parameter to be used for kernel continuization 
			of \eqn{F(x)}{F(x)}. If not provided (Default), this value is automatically calculated (see 
			details).}

  \item{hy}{		An integer indicating the value of the bandwidth parameter to be used for kernel continuization 
			of \eqn{G(y)}{G(y)}. If not provided (Default), this value is automatically calculated (see 
			details).}


  \item{degree}{	A vector indicating the number of power moments to be fitted to the marginal distributions 
			("EG" design), and/or the number or cross moments to be fitted to the joint distributions 
			(see Details).}

  \item{design}{	A character string indicating the equating design (one of "EG", "SG", "CB", "NEAT_CE", 
			"NEAT_PSE")}

  \item{Kp}{		A number which acts as a weight for the second term in the combined penalization function used 
			to obtain \code{h} (see details).}

  \item{scores2}{	Only used for the "CB", "NEAT_CE" and "NEAT_PSE" designs. See the description of
			\code{scores}.}

  \item{degreeXA}{	A vector indicating the	number of power moments to be fitted to the marginal distributions 
			\eqn{X}{X} and \eqn{A}{A}, and the number or cross moments to be fitted to the joint 
			distribution \eqn{(X,A)}{(X,A)}	(see details). Only used for the "NEAT_CE" and "NEAT_PSE" designs.}

  \item{degreeYA}{	Only used for the "NEAT_CE" and "NEAT_PSE" designs (see the description for
			\code{degreeXA})}

\item{J}{		The number of possible \eqn{X}{X} scores. Only needed for "CB", "NEAT_CB" and "NEAT_PSE" designs}

  \item{K}{		The number of possible \eqn{Y}{Y} scores. Only needed for "CB", "NEAT_CB" and "NEAT_PSE" designs}

  \item{L}{		The number of possible \eqn{A}{A} scores. Needed for "NEAT_CB" and "NEAT_PSE" designs}

  \item{wx}{		A number that satisfies \eqn{0\leq w_X\leq 1}{0<=w_x<=1} indicating the weight put on the data 
			that is not subject to order effects. Only used for the "CB" design.}

  \item{wy}{		A number that satisfies \eqn{0\leq w_Y\leq 1}{0<=w_y<=1} indicating the weight put on the data 
			that is not subject to order effects. Only used for the "CB" design.}

  \item{w}{		A number that satisfies \eqn{0\leq w\leq 1}{0<=w<=1} indicating the weight given to 
			population \eqn{P}{P}. Only used for the "NEAT" design.}
			
	\item{gapsX}{ A \emph{list} object containing:
	                \describe{
	                  \item{\code{index}}{A vector of indices between \eqn{0} and \eqn{J} to smooth "gaps", usually ocurring at regular intervals due to scores rounded to integer values and other methodological factors. }
	                  \item{\code{degree}}{An integer indicating the maximum degree of the moments fitted by the log-linear model.}
	                 }
	                  Only used for the "NEAT" design.
   }
   
   	\item{gapsY}{ A \emph{list} object containing:
	                \describe{
	                  \item{\code{index}}{A vector of indices between \eqn{0} and \eqn{K}.}
	                  \item{\code{degree}}{An integer indicating the maximum degree of the moments fitted.}
	                 }
	                  Only used for the "NEAT" design.
   }
   
   	\item{gapsA}{ A \emph{list} object containing:
	                \describe{
	                  \item{\code{index}}{A vector of indices between \eqn{0} and \eqn{L}. }
	                  \item{\code{degree}}{An integer indicating the maximum degree of the moments fitted.}
	                 }
	                  Only used for the "NEAT" design.
   }
   
   \item{lumpX}{An integer to represent the index where an artificial "lump" is created in the marginal distribution of frecuencies for \eqn{X} due to recording of negative rounded formulas or any other methodological artifact.}
   
   \item{lumpY}{An integer to represent the index where an artificial "lump" is created in the marginal distribution of frecuencies for \eqn{Y}.}
   
   \item{lumpA}{An integer to represent the index where an artificial "lump" is created in the marginal distribution of frecuencies for \eqn{A}.}
   
   \item{alpha}{Only for Adaptative Kernel. Sensitivity parameter.}
   
   \item{h.adap}{Only for Adaptative Kernel. A list(hx, hy) containing bandwidths for Adaptative kernel for each Form.}
}

\details{This is a generic function that implements the kernel method of test equating as described in Von Davier et al.
	(2004). Given test scores \eqn{X}{X} and \eqn{Y}{Y}, the functions calculates 
	\deqn{\hat{e}_Y(x)=G_{h_{Y}}^{-1}(F_{h_{X}}(x;\hat{r}),\hat{s})}{\hat{e}_Y(x)=G_{hy}^{-1}(F_{hx}(x;\hat{r}),\hat{s})}
	where \eqn{\hat{r}}{\hat{r}} and \eqn{\hat{s}}{\hat{s}} are estimated score probabilities obtained via loglinear
	smoothing (see \code{\link{loglin.smooth}}). The value of \eqn{h_X}{h_X} and \eqn{h_Y}{h_Y} can either be specified 
	by the user or left unspecified (default) in which case they are automatically calculated. For instance, one can
	specifies large values of \eqn{h_X}{h_X} and \eqn{h_Y}{h_Y}, so that the \eqn{\hat{e}_Y(x)}{\hat{e}_Y(x)} tends to the 
	linear equating function (see Theorem 4.5 in Von Davier et al, 2004 for more details). 
}

\value{	An object of class \code{ker.eq} representing the kernel equating process. Generic functions such as 
	\code{print}, and \code{summary} have methods to show the results of the equating. The results include 
	summary statistics, equated values, standard errors of equating, and others.  
  
	The function \code{\link{SEED}} can be used to obtain standard error of equating differences (SEED) of two 
	objects of class \code{ker.eq}. The function \code{\link{PREp}} can be used on a \code{ker.eq} object to 
	obtain the percentage relative error measure (see Von Davier et al, 2004). 


\item{Scores}{The possible values of \eqn{x_j}{xj} and \eqn{y_k}{yk}}
\item{eqYx }{The equated values of test \eqn{X}{X} in test \eqn{Y}{Y} scale}
\item{eqXy }{The equated values of test \eqn{Y}{Y} in test \eqn{X}{X} scale}
\item{SEEYx}{The standard error of equating for equating \eqn{X}{X} to \eqn{Y}{Y}}
\item{SEEXy}{The standard error of equating for equating \eqn{Y}{Y} to \eqn{X}{X}}
}

\references{

        Gonzalez, J. (2014). SNSequate: Standard and Nonstandard Statistical Models and Methods for Test 
        Equating. \emph{Journal of Statistical Software, 59(7),} 1-30.

	Holland, P. and Thayer, D. (1989). The kernel method of equating score distributions. 
	(Technical Report No 89-84). Princeton, NJ: Educational Testing Service.

	Holland, P., King, B. and Thayer, D. (1989). The standard error of equating for the kernel method 
	of equating score distributions (Tech. Rep. No. 89-83). Princeton, NJ: Educational Testing Service. 

	Von Davier, A., Holland, P., and Thayer, D. (2004). \emph{The Kernel Method of Test Equating}. 
	New York, NY: Springer-Verlag.
}

\author{Jorge Gonzalez \email{jorge.gonzalez@mat.uc.cl}}

\seealso{\code{\link{loglin.smooth}}, \code{\link{SEED}}, \code{\link{PREp}}
}

\examples{
#Kernel equating under the "EG" design
data(Math20EG)
mod<-ker.eq(scores=Math20EG,kert="gauss",hx=NULL,hy=NULL,degree=c(2,3),design="EG") 

summary(mod)

#Reproducing Table 7.6 in Von Davier et al, (2004)

scores<-0:20
SEEXy<-mod$SEEXy
SEEYx<-mod$SEEYx

Table7.6<-cbind(scores,SEEXy,SEEYx)
Table7.6

#Other nonstandard kernels. Table 10.3 in Von Davier (2011).

mod.logis<-ker.eq(scores=Math20EG,kert="logis",hx=NULL,hy=NULL,degree=c(2,3),design="EG") 
mod.unif<-ker.eq(scores=Math20EG,kert="unif",hx=NULL,hy=NULL,degree=c(2,3),design="EG") 
mod.gauss<-ker.eq(scores=Math20EG,kert="gauss",hx=NULL,hy=NULL,degree=c(2,3),design="EG") 

XtoY<-cbind(mod.logis$eqYx,mod.unif$eqYx,mod.gauss$eqYx)
YtoX<-cbind(mod.logis$eqXy,mod.unif$eqXy,mod.gauss$eqXy)

Table10.3<-cbind(XtoY,YtoX)
Table10.3

## Examples using Adaptive and Epanechnikov kernels
x_sim = c(1,2,3,4,5,6,7,8,9,10,11,10,9,8,7,6,5,4,3,2,1)
prob_sim = x_sim/sum(x_sim)
set.seed(1)
sim = rmultinom(1, p = prob_sim, size = 1000)

x_asimD = c(1,7,13,18,22,24,25,24,20,18,16,15,13,9,5,3,2.5,1.5,1.5,1,1)
probas_asimD = x_asimD/sum(x_asimD)
set.seed(1)
asim = rmultinom(1, p = probas_asimD, size = 1000)

scores = cbind(asim,sim)

mod.adap  = ker.eq(scores,degree=c(2,2),design="EG",kert="adap")
mod.epan  = ker.eq(scores,degree=c(2,2),design="EG",kert="epan")
}

\concept{kernel equating}


