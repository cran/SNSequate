\name{PREp}
\alias{PREp}
\alias{PREp.default}
\title{Percent relative error}

\description{This function calculates the percent relative error as described in Von Davier et al. (2004). 
}
\usage{
PREp(eq, p)
}

\arguments{
  \item{eq}{	An object of class \code{ker.eq} previously obtained using \code{\link{ker.eq}}.

}
  \item{p}{	The number of moments to be calculated.

}
}

\details{PREp (when equating form \code{X} to \code{Y}) is calculated as 
\deqn{\mbox{PREp}=100\frac{\mu_p(e_Y(X))-\mu_p(Y)}{\mu_p(Y)}}{100*(\mu_p(e_Y(X))-\mu_p(Y))/\mu_p(Y)}
where \eqn{\mu_p(Y)=\sum_k(y_k)^ps_k}{\mu_p(Y)=\sum_k(y_k)^p*s_k} and 
\eqn{\mu_p(e_Y(X))=\sum_j(e_Y(x_j))^pr_j}{\mu_p(e_Y(X))=\sum_j(e_Y(x_j))^p*r_j}. Similar formulas can be found 
when equating from \code{Y} to \code{X}.
}

\value{ A matrix containing the PREp for both \code{X} to \code{Y} (first column) and \code{Y} to \code{X} 
(second column) cases.
}
\references{
	Von Davier, A., Holland, P., and Thayer, D. (2004). \emph{The Kernel Method of Test Equating}. 
	New York, NY: Springer-Verlag.
}
\author{Jorge Gonzalez B.
}

\seealso{\code{\link{ker.eq}}
}

\examples{
#Example: Table 7.5 in Von Davier et al. (2004)

data(Math20EG)
mod.gauss<-ker.eq(scores=Math20EG,kert="gauss", hx = NULL, hy = NULL,degree=c(2, 3),design="EG")
PREp(mod.gauss,10)

}
\keyword{PREp}

