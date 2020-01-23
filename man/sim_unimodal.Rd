\name{sim_unimodal}
\alias{sim_unimodal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate test scores.
}
\description{
Simulate test scores from a negative-hypergeometric (beta-binomial) distribution, according to
Keats & Lord (1962).
}
\usage{
sim_unimodal(n, x_mean, x_var, N_item, seed = NULL, name = NULL)
}

\arguments{
  \item{n}{
    Size of the resulting sample.
  }
  \item{x_mean}{
    Mean of the target distribution.
  }
  \item{x_var}{
    Variance of the target distribution.
  }
  \item{N_item}{
    Number of items in the test.
  }
  \item{seed}{
    Optional. Seed for the random number generator.
  }
  \item{name}{
    Optional. Generate X and Y scores from the data according 5 of the proposed distributions 
    in Keats & Lord (1967). Overrides any other previous parameter input set.
  }
}
\details{
Simulate test scores from a negative-hypergeometric (beta-binomial) distribution, according to
Keats & Lord (1962).
}
\value{
Simulated values.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Keats, J. A., & Lord, F. M. (1962). A theoretical distribution for mental test scores. Psychometrika, 27(1), 59-72.

}
\author{
Daniel Leon Acuna, \email{dnacuna@uc.cl}
}

\examples{

sim_unimodal(2354, 27.06, 8.19^2, 40)  # GANA
sim_unimodal(name="TQS8")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ Simulation }
\concept{ Scores }% __ONLY ONE__ keyword per line
