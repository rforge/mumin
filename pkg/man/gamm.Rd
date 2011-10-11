\name{gamm}
\alias{gamm}
\alias{MuMIn-gamm}
\encoding{utf-8}

\title{GAMM wrapper for use with MuMIn}
\description{
Adapter function for use mixed Generalized Additive Models from packages
\pkg{mgcv} and \pkg{gamm} with \pkg{MuMIn}
}

\usage{
gamm(formula, random = NULL, ..., lme4 = inherits(random, "formula"))

}

\arguments{
  \item{formula, random, \dots}{ arguments passed to \code{gamm} or
	\code{gamm4} }
  \item{lme4}{ logical, if true \code{gamm4} is used rather than \code{gamm},
  If \code{TRUE}, the \code{random} argument must be provided as a
  \code{formula}.
  }
}

\value{
	Depending on the value of the 'lme4' switch, either a \code{gamm} or
	\code{gamm4} fitted model object. The only difference from the originals is
	the \code{call} component, allowing for \code{\link{update}}'ing of the
	object.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link[mgcv]{gamm}} and \code{\link[gamm4]{gamm4}}
}


\keyword{models}
