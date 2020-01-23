\name{KB36.1PL}
\alias{KB36.1PL}
\docType{data}

\title{Difficulty parameter estimates for KB36 data under a 1PL model}


\description{This data set contains the estimated item difficuty parameters for the 
\code{KB36} data, assuming a 1PL model. Two sets of parameters estimates for test forms 
\code{X} and \code{Y} are available: one that results from a fit assuming the traditional 
logistic link, and one which comes from the fit using a cloglog (asymmetric) link. 
}

\usage{data(KB36.1PL)}

\format{A list of 2 elements containing item (difficulty) parameters estimates for test 
forms \code{X} and \code{Y} under the logistic-link model (\code{b.logistic}), and under 
the cloglog-link model (\code{b.cloglog})
}

\details{This data set is used to illustrate the characteristic curve methods (Haebara and 
	 Stocking-Lord) which can use an asymmetric cloglog ICC for the calculations, as 
	 described in Estay (2012).

	 A 1PL model using both logistic and cloglog link can be fitted using the \code{lmer()} 
	 function in the \code{lme4} R package (see De Boeck et. al, 2011 for details). 
}

\source{The item parameter estimates for the 1PL model with logistic link are also shown in
Table 6.13 of Kolen and Brennan (2004).
}

\references{
	
	De Boeck, P., Bakker, M., Zwitser, R., Nivard, M., Hofman, A.,Tuerlinckx, F., Partchev, I. 
	(2011). The Estimation of Item Response Models with the \code{lmer} Function from the 
	\pkg{lme4} Package in \code{R}. \emph{Journal of Statistical Software, 39(12),} 1-28.

	Kolen, M., and Brennan, R. (2004). \emph{Test Equating, Scaling and Linking}. 
	New York, NY: Springer-Verlag.

	Estay, G. (2012). \emph{Characteristic Curves Scale Transformation Methods Using 
	Asymmetric ICCs for IRT Equating}. Unpublished MSc. Thesis. Pontificia Universidad 
	Catolica de Chile

}

\examples{
data(KB36.1PL)
## maybe str(KB36.1PL) ; plot(KB36.1PL) ...
}
\keyword{datasets}
