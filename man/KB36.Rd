\name{KB36}
\alias{KB36}
\docType{data}

\title{Data on two 36-items test forms}

\description{
The data set contains both response patterns and item parameters estimates following a 3PL model 
for two 36-items tests forms. Form \code{X} was administered to 1655 examinees and form \code{Y} 
to 1638 examinees. Also, 12 out of the 36 items are common between both test forms (items 3, 6, 9, 
12, 15, 18, 21, 24, 27, 30, 33, 36). This data has been described and analized by Kolen and Brennan (2004).}

\usage{data(KB36)}

\format{A list with four elements containing binary data matrices of responses (\code{KBformX} and 
\code{KBformY}) and the corresponding parameter estimates which result from a 3PL fit to both data
matrices (\code{KBformX_par} and \code{KBformY_par}).
}

\source{The data come with the distribution of the CIPE software which is freely available 
at \url{https://education.uiowa.edu/casma/computer-programs}. 
The list of item parameters estimates can be found in Table 6.5 of Kolen and Brennan (2004).
}

\references{

	Kolen, M., and Brennan, R. (2004). \emph{Test Equating, Scaling and Linking}. 
	New York, NY: Springer-Verlag.
}

\examples{
data(KB36)
## maybe str(KB36) ; plot(KB36) ...
}
\keyword{datasets}
