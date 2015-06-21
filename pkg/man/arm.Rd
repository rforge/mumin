\name{arm.glm}
\alias{arm.glm}
\encoding{utf-8}
\title{Adaptive Regression by Mixing}
\description{
Combine all-subsets GLMs using the ARM algorithm.
}

\usage{
arm.glm(object, nrepl = 250, weight.by = c("aic", "loglik"), trace = FALSE)
}

\arguments{
  \item{object}{a fitted \dQuote{global} \code{glm} object. }
  \item{nrepl}{number of permutations. }
  \item{weight.by}{indicates whether model weights should be calculated with AIC
    or log-likelihood. }
  \item{trace}{if \code{TRUE}, information is printed during the running of
    \code{arm.glm}. }
}

\details{
The data is randomly split into two parts. For all-subsets of the 
\dQuote{global} model, parameters are estimated using the first half of the data
and log-likelihood given the second half of the data is used to calculate 
AIC weights. Mean weights from \code{nrepl} permutations of the data are 
subsequently used to average all-subsets parameters estimated using complete 
data.
}

\note{
All-subsets respect marginality constraints. Number of parameters is limited to
\code{floor(nobs(object) / 2) - 1}.
}

\value{
An object of class \code{"averaging"} contaning only \dQuote{full} averaged 
coefficients. See \code{\link{model.avg}} for object description.
}

\references{
Yang Y. (2001) Adaptive Regression by Mixing. 
\emph{Journal of the American Statistical Association} 96: 574–588.

Yang Y. (2003) Regression with multiple candidate models: selecting or mixing? 
\emph{Statistica Sinica} 13: 783–810.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{model.avg}}, \code{\link{par.avg}}

Other implementation: \code{\link[MMIX]{arms}} in package \pkg{MMIX}.
}

\examples{
fm <- glm(y ~ X1 + X2 + X3 + X4, data=Cement)

summary(arm.glm(fm, nrepl = 25))
}

\keyword{models}
