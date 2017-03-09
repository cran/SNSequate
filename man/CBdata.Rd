\name{CBdata}
\alias{CBdata}
\docType{data}
\title{Observed (raw) score values for two different tests
}
\description{The data set is from a small field study from an international testing program. 
It contains the observed scores for two tests \eqn{X}{X} (with 75 items) and \eqn{Y}{Y} (with 76 items) 
administered to two independent, random samples of examinees from a single population \eqn{P}. 
For more details, see Chapter 9 in Von Davier et al, (2004) from where the data were obtained. 
}
\usage{data(CBdata)}
\format{
A list with elements containing the observed scores of the sample taking test X first,
followed by test Y (datX1Y2), and the scores of the sample taking test Y first followed by 
test X (datX2Y1).
}

\references{
	Von Davier, A., Holland, P., and Thayer, D. (2004). \emph{The Kernel Method of Test Equating}. 
	New York, NY: Springer-Verlag. 
}
\examples{
data(CBdata)
## maybe str(CBdata) ; ...
}
\keyword{datasets}
