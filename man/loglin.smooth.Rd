\name{loglin.smooth}
\alias{loglin.smooth}
\alias{loglin.smooth.default}

\title{Pre-smoothing using log-linear models.}

\description{This function fits log-linear models to score data and provides estimates of 
the (vector of) score probabilities as well as the \code{C} matrix decomposition of their 
covariance matrix, according to the specified equating design (see Details).
}

\usage{loglin.smooth(scores, degree, design, scores2, degreeXA, degreeYA, 
J, K, L, wx, wy, w, ...)}

\arguments{
  \item{}{		Note that depending on the specified equating design, not all arguments are necessary 
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

  \item{degree}{	Either a number or vector indicating the number of power moments to be fitted to the 
			marginal distributions, or the number or cross moments to be fitted to the joint 
			distributions, respectively. For the "EG" design it will be a number (see Details).}

  \item{design}{	A character string indicating the equating design (one of "EG", "SG", "CB", "NEAT_CE", 
			"NEAT_PSE")}

  \item{scores2}{       Only used for the "CB", "NEAT_CE" and "NEAT_PSE" designs. See the description of
			\code{scores}.}

  \item{degreeXA}{	A vector indicating the	number of power moments to be fitted to the marginal distributions \eqn{X}{X} and \eqn{A}{A},
			and the number or cross moments to be fitted to the joint distribution \eqn{(X,A)}{(X,A)}
			(see details). Only used for the "NEAT_CE" and "NEAT_PSE" designs.}

  \item{degreeYA}{	Only used for the "NEAT_CE" and "NEAT_PSE" designs (see the description for
			\code{degreeXA})}

  \item{J}{		The number of possible \eqn{X}{X} scores. Only needed for "CB", "NEAT_CB" and "NEAT_PSE" designs}

  \item{K}{		The number of possible \eqn{Y}{Y} scores. Only needed for "CB", "NEAT_CB" and "NEAT_PSE" designs}

  \item{L}{		The number of possible \eqn{A}{A} scores. Needed for "NEAT_CB" and "NEAT_PSE" designs}

  \item{wx}{		A number that satisfies \eqn{0\leq w_X\leq 1}{0<=w_x<=1} indicating the weight put on the 
			data that is not subject to order effects. Only used for the "CB" design.}

  \item{wy}{		A number that satisfies \eqn{0\leq w_Y\leq 1}{0<=w_y<=1} indicating the weight put on the 
			data that is not subject to order effects. Only used for the "CB" design.}

  \item{w}{		A number that satisfies \eqn{0\leq w\leq 1}{0<=w<=1} indicating the weight given to 
			population \eqn{P}{P}. Only used for the "NEAT" design.}

  \item{\dots}{		Further arguments currently not used.}
}

\details{This function fits loglinear models as described in Holland and Thayer (1987), and Von Davier 
	et al. (2004). The following general equation can be used to represent the models according to 
	the different designs used, in which the vector \eqn{o}{o} (or matrix) of (marginal or bivariate) 
	score probabilities satisfies the log-linear model: 
	\deqn{\log(o_{gh})=\alpha_m+Z_m(z_g)+W_m(w_h)+ZW_m(z_g,w_h)}{log(o_{gh})=alpha_m+Z_m(z_g)+W_m(w_h)+ZW_m(z_g,w_h) ,}
	where \eqn{Z_m(z_g)=\sum_{i=1}^{T_{Zm}}\beta_{zmi}(z_g)^i}{Z_m(z_g)=\sum_{i=1}^{T_{Zm}}\beta_{zmi}(z_g)^i}, 
	\eqn{W_m(w_h)=\sum_{i=1}^{T_{Wm}}\beta_{Wmi}(w_h)^i}{W_m(w_h)=\sum_{i=1}^{T_{Wm}}\beta_{Wmi}(w_h)^i}, and, 
	\eqn{ZW_m(z_g,w_h)=\sum_{i=1}^{I_{Zm}}\sum_{i'=1}^{I_{Wm}}\beta_{ZWmii'}(z_g)^i(w_h)^{i'}}{ZW_m(z_g,w_h)=\sum_{i=1}^{I_{Zm}}\sum_{i'=1}^{I_{Wm}}\beta_{ZWmii'}(z_g)^i(w_h)^{i'}}.
 
	The symbols will vary according to the different equating designs specified. Possible values are: 
	\eqn{o=p_{(12)}, p_{(21)}, p, q}{o=p_{(12)}, p_{(21)}, p, q}; \eqn{Z=X, Y}{Z=X, Y}; \eqn{W=Y, A}{W=Y, A}; 
	\eqn{z=x, y}{z=x, y}; \eqn{w=y, a}{w=y, a}; \eqn{m=(12), (21), P, Q}{m=(12), (21), P, Q}; \eqn{g=j, k}{g=j, k}; 
	\eqn{h=l, k}{h=l, k}.
	
	Particular cases of this general equation for each of the equating designs can be found in 
	Von Davier et al (2004) (e.g., Equations (7.1) and (7.2) for the "EG" design, Equation (8.1) for the "SG" design, 
	Equations (9,1) and (9.2) for the "CB" design).
	 
}
	
\value{
  \item{sp.est}{The estimated score probabilities}
  \item{C}{The C matrix which is so that \eqn{\Sigma=CC^t}{Sigma=CC^t}}}

\references{
    Holland, P. and Thayer, D. (1987). Notes on the use of loglinear models for fitting 
	discrete probability distributions. Research Report 87-31, Princeton NJ: Educational 
	Testing Service.

	Von Davier, A., Holland, P., and Thayer, D. (2004). \emph{The Kernel Method of Test Equating}. 
	New York, NY: Springer-Verlag.}

\author{Jorge Gonzalez B.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{glm}}, \code{\link{ker.eq}}}

\examples{
#Table 7.4 from Von Davier et al. (2004)
data(Math20EG)
rj<-loglin.smooth(scores=Math20EG[,1],degree=2,design="EG")$sp.est
sk<-loglin.smooth(scores=Math20EG[,2],degree=3,design="EG")$sp.est
score<-0:20
Table7.4<-cbind(score,rj,sk)
Table7.4
}

\keyword{presmoothing}
\keyword{log-linear}
