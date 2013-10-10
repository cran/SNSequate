\name{ACTmKB}
\alias{ACTmKB}
\docType{data}
\title{Scores on two 40-items ACT mathematics test forms
}

\description{
The data set contains raw sample frequencies of number-right scores for two 
multiple choice 40-items mathematics tests forms. Form \code{X} 
was administered to 4329 examinees and form \code{Y} to 4152 examinees. This data has 
been described and analized by Kolen and Brennan (2004).
}

\usage{data(ACTmKB)}

\format{A 41x2 matrix containing raw sample frequencies (raws) for two tests (columns). 
}

\source{The data come with the distribution of the RAGE-RGEQUATE software which is freely available 
at \url{http://www.education.uiowa.edu/docs/default-source/casma-computerprograms/}
}

\references{

	Kolen, M., and Brennan, R. (2004). \emph{Test Equating, Scaling and Linking}. 
	New York, NY: Springer-Verlag.
}

\examples{
data(ACTmKB)
## maybe str(ACTmKB) ; plot(ACTmKB) ...
}
\keyword{datasets}
