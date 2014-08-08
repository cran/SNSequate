\name{bandwidth}
\alias{bandwidth}
\alias{bandwidth.default}

\title{Automatic selection of the bandwidth parameter \code{h}}

\description{This functions implements the minimization of the combined penalty function
	     described by Holland and Thayer (1989); Von Davier et al, (2004). It returns 
	     the optimal value of \code{h} for kernel continuization, according to the above 
	     mentioned criteria. Different types of kernels (others than the gaussian) are accepted.
}

\usage{
bandwidth(scores, kert, degree, design, Kp = 1, scores2, degreeXA, degreeYA, 
J, K, L, wx, wy, w, ...)
}

\arguments{

  \item{}{      	Note that depending on the specified equating design, not all arguments are necessary
			as detailed below.}

  \item{scores}{	If the "EG" design is specified, a vector containing the raw sample frequencies 
			coming from one group taking the test.
			
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
			Current options include "\code{gauss}", "\code{logis}", and "\code{uniform}" 
			for the gaussian, logistic and uniform kernels, respectively}


  \item{degree}{	Either a number or vector indicating the number of power moments to be fitted to the 
			marginal distributions, or the number or cross moments to be fitted to the joint 
			distributions, respectively. For the "EG" design it will be a number (see Details).}


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

  \item{\dots}{		Further arguments currently not used.}
}

\details{	To automatically select \code{h}, the function minimizes 
		\deqn{PEN_1(h)+K\times PEN_2(h)}{PEN_1(h)+K*PEN_2(h),}
		where \eqn{PEN_1(h)=\sum_j(\hat{r}_j-\hat{f}_h(x_j))^2}{PEN_1(h)=\sum_j(\hat{r}_j-\hat{f}_h(x_j))^2}, and 
		\eqn{PEN_2(h)=\sum_jA_j(1-B_j)}{PEN_2(h)=\sum_jA_j(1-B_j)}. The terms \eqn{A}{A} and \eqn{B}{B} 
		are such that \eqn{PEN_2}{PEN_2} acts as a smoothness penalty term that avoids rapid 
		fluctuations in the approximated density (see Chapter 10 in Von Davier, 2011 for more details). The 
		\eqn{K}{K} term corresponds to the \code{Kp} argument of the \code{bandwidth} function. The 
		\eqn{\hat{r}}{\hat{r}} values are assumed to be estimated by polynomial loglinear models of specific
		\code{degree}, which come from a call to \code{\link{loglin.smooth}}.
}

\value{		A number which is the optimal value of \code{h}.} 

\references{

        Gonzalez, J. (2014). SNSequate: Standard and Nonstandard Statistical Models and Methods for Test 
        Equating. \emph{Journal of Statistical Software, 59(7),} 1-30.

	Von Davier, A., Holland, P., and Thayer, D. (2004). \emph{The Kernel Method of Test Equating}. 
	New York, NY: Springer-Verlag.

	A. von Davier (Ed.) (2011). \emph{Statistical Models for Equating, Scaling, and Linking. New York: Springer}
}
\author{
Jorge Gonzalez B. \email{jgonzale@mat.puc.cl}
}


\seealso{\code{\link{loglin.smooth}}
}
\examples{
#Example: The "Standard" column and firsts two rows of Table 10.1 in 
#Chapter 10 of Von Davier 2011

data(Math20EG)

hx.logis<-bandwidth(scores=Math20EG[,1],kert="logis",degree=2,design="EG")$h
hx.unif<-bandwidth(scores=Math20EG[,1],kert="unif",degree=2,design="EG")$h 
hx.gauss<-bandwidth(scores=Math20EG[,1],kert="gauss",degree=2,design="EG")$h

hy.logis<-bandwidth(scores=Math20EG[,2],kert="logis",degree=3,design="EG")$h
hy.unif<-bandwidth(scores=Math20EG[,2],kert="unif",degree=3,design="EG")$h 
hy.gauss<-bandwidth(scores=Math20EG[,2],kert="gauss",degree=3,design="EG")$h

partialTable10.1<-rbind(c(hx.logis,hx.unif,hx.gauss),
				c(hy.logis,hy.unif,hy.gauss))

dimnames(partialTable10.1)<-list(c("h.x","h.y"),c("Logistic","Uniform","Gaussian"))
partialTable10.1

}

\keyword{kernel equating}
