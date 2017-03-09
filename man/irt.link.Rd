\name{irt.link}
\alias{irt.link}
\alias{irt.link.default}

\title{IRT parameter linking methods}

\description{The function implements parameter linking methods to transform IRT scales. Mean-mean, 
mean-sigma, Haebara, and Stocking and Lord methods are available (see details). 

}
\usage{
irt.link(parm, common, model, icc, D, ...)
}

\arguments{
  \item{parm}{		A 6 column matrix containing item parameter estimates from an IRT model. The
			first three columns contains the parameters for the form \code{Y} fit, and the 
			last three those of form \code{X}. The order for item paramters in the matrix
			is discrimination, difficulty, and guessing. See details.}

  \item{common}{	A vector indicating the position where common items are located}

  \item{model}{		A character string indicating the underlying IRT model: "1PL", "2PL", "3PL".}

  \item{icc}{		A character string indicating the type of \code{icc} used in the characteristic
			curve methods (see details). Available options are "logistic" and "cloglog".}

  \item{D}{		A number indicating the value of the constant \code{D} 
			(see details)}

  \item{\dots}{		Further arguments currently not used.}

}

\details{The function implments various methods of IRT parameter linking (a.k.a, scale transformation 
methods). It calculates the linking constants \code{A} and \code{B} to tranform parameter estimates. 
When assuming a 1PL model, the matrix \code{parm} should contain a column of ones and a column of zeroes 
in the places where discrimination and guessing parameters are located, respectively. 

The characteristic curve methods (Haebara and Stocking and Lord) rely on the item characteristic curve 
\eqn{p_{ij}}{p_{ij}}assumed for the probability of a correct answer
\deqn{p_{ij}=P(Y_{ij}=1\mid\theta_i)=c_j+(1-c_j)\frac{\exp[Da_j(\theta_i-\beta_j)]}{1+\exp[Da_j(\theta_i-\beta_j)]}}{p_{ij}=P(Y_{ij}=1|theta_i)=cj+(1-cj){exp[Da_j(theta_i-beta_j)]}/{1+exp[Da_j(theta_i-beta_j)]}}

Besides the traditional logistic model, the \code{irt.link()} function allows the use of an asymetric 
cloglog ICC. See the help for \code{KB36.1PL} data set, where some details on how to fit a 1PL model with
cloglog link in \code{lmer} are given.

For more details on characteristic curve methods see Kolen and Brennan (2004).
}

\note{Currently, the cloglog ICC is only implmented for the 1PL model. A 1PL model with asymetric cloglog 
link can be fitted in R using the \code{lmer()} function in package \code{lme4}
}


\value{A list with the constants \code{A} and \code{B} calculated using the four different methods}

\references{

        Gonzalez, J. (2014). SNSequate: Standard and Nonstandard Statistical Models and Methods for Test 
        Equating. \emph{Journal of Statistical Software, 59(7),} 1-30.

	Kolen, M., and Brennan, R. (2004). \emph{Test Equating, Scaling and Linking}. 
	New York, NY: Springer-Verlag.

	Estay, G. (2012). \emph{Characteristic Curves Scale Transformation Methods Using 
	Asymmetric ICCs for IRT Equating}. Unpublished MSc. Thesis. Pontificia Universidad 
	Catolica de Chile
}

\author{Jorge Gonzalez B. \email{jgonzale@mat.puc.cl}}

\seealso{\code{\link{mea.eq}}, \code{\link{lin.eq}}, \code{\link{ker.eq}}
}

\examples{
#### Example. KB, Table 6.6
data(KB36)
parm.x = KB36$KBformX_par
parm.y = KB36$KBformY_par	
comitems = seq(3,36,3)
parm = as.data.frame(cbind(parm.y, parm.x))

# Table 6.6 KB
irt.link(parm,comitems,model="3PL",icc="logistic",D=1.7)


# Same data but assuming a 1PL model. The parameter estimates are load from 
# the KB36.1PL data set. See the help for KB36.1PL data for details on how these
# estimates were obtained using \code{lmer()} (see also Table 6.13 in KB)
 
data(KB36.1PL)

#preparing the input data matrices for irt.link() function
b.log.y<-KB36.1PL$b.logistic[,2]
b.log.x<-KB36.1PL$b.logistic[,1]
b.clog.y<-KB36.1PL$b.cloglog[,2]
b.clog.x<-KB36.1PL$b.cloglog[,1]

parm2 = as.data.frame(cbind(1,b.log.y,0, 1,b.log.x, 0))
parm3 = as.data.frame(cbind(1,b.clog.y,0, 1,b.clog.x,0))

#vector indicating common items
comitems = seq(3,36,3)

#Calculating the B constant under the logistic-link model
irt.link(parm2,comitems,model="1PL",icc="logistic",D=1.7)

#Calculating the B constant under the cloglog-link model
irt.link(parm3,comitems,model="1PL",icc="cloglog",D=1.7)
}

\keyword{IRT parameter linking}

