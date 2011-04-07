\name{nobs-me}
\alias{nobs.lme}
\alias{nobs.mer}
\alias{nobs.gls}
\encoding{utf-8}
\title{Extract the number of observations from a mixed effect fit}
\description{
Methods for extracting the number of sQuote(observations) from a mixed effect model fit.
}

\usage{
\method{nobs}{mer}(object, nall = FALSE, ...)
}

\arguments{
  \item{object}{a fitted model object}
  \item{nall}{logical, should actual number of observations be returned?}
  \item{\dots}{Further arguments to be passed to methods (currently ignored)}

}

\details{
These are additional methods to the generic function \code{\link{nobs}}, defined
in package \pkg{stats}, for classes \code{mer}, \code{lme}
and \code{gls}.

By default, the number of observations corrected by the number of fixed effects
is returned.

}

\author{Kamil \enc{Bartoń}{Barton}}

\seealso{
\code{\link{nobs}}, \code{\link[nlme]{logLik.lme}}
}

\keyword{models}
