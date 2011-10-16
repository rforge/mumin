\name{MuMIn-package}
\alias{MuMIn-package}
\alias{MuMIn}
\docType{package}
\encoding{utf-8}
\title{
Multi-model inference
}
\description{
The package \code{MuMIn} contains functions to streamline model selection and
perform model averaging based on information criteria (AIC, AICc and alike).
}

\details{
User level functions include:
\describe{
	\item{\code{\link{dredge}}}{performs automated model selection based on 
		subsets of the supplied \sQuote{global} model, and optional
		choices of other model properties (such as different link functions). 
		Model set may be generated either with \sQuote{all possible} 
		combinations, or tailored according to the conditions specified }
	\item{\code{\link{model.sel}}}{ creates a model selection table from 
		handpicked models }	
	\item{\code{\link{model.avg}}}{ calculates model averaged parameters }
	\item{\code{\link{AICc}}}{calculates second-order Akaike information 
		criterion for one or several fitted model objects.}
}

Model selection can be done according to any information criterion, including 
\code{\link{AIC}}, \code{\link{AICc}}, \code{\link{BIC}}, \code{\link{QAIC}}, 
\code{\link{ICOMP}} or \code{\link{Mallows' Cp}}.

}

\author{ Kamil Barto\enc{ń}{n} }

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed.
}


\seealso{
\code{\link{AIC}}, \code{\link{step}}

}
\examples{

fm1 <- lm(Fertility ~ . , data = swiss)

dd <- dredge(fm1)
top.models.1 <- get.models(dd, subset = delta < 4)
model.avg(top.models.1) # get averaged coefficients

top.models.2 <- get.models(dd, cumsum(weight) <= .95)
model.avg(top.models.2)

# Mixed models:
# modified example(lme)
data(Orthodont, package="nlme")
require(nlme)
fm2 <- lme(distance ~ age + Sex, data = Orthodont,
    random = ~ 1 | Subject, method="ML")
dredge(fm2)

}

\keyword{package}
\keyword{models}
