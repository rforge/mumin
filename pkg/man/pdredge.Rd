\name{pdredge}
\alias{pdredge}

\encoding{utf-8}
\title{Automated model selection using parallel computation}
\description{
This is a parallel version of \code{dredge}.
}

\usage{
pdredge(global.model, cluster = FALSE, beta = FALSE, evaluate = TRUE,
    rank = "AICc", fixed = NULL, m.max = NA, m.min = 0, subset,
    marg.ex = NULL, trace = FALSE, varying, extra, check = TRUE,
    ...)
}

\arguments{
	\item{global.model, beta, evaluate, rank, trace}{
		see \code{\link{dredge}}. }
	\item{fixed, m.max, m.min, subset, marg.ex, varying, extra, ...}{
		see \code{\link{dredge}}. }
	\item{cluster}{ either a valid cluster object, or \code{FALSE}. }
	\item{check}{ logical, whether to evaluate the \code{global.model} in the
	cluster	and compare with the original one. Only a simple checking is done. }

}

\details{
All the dependencies for fitting the \code{global.model} (including the data,
all required packages, and any object the modelling function may use) must be
exported/loaded into the cluster before passing it to \code{pdredge}.

This function is still experimental. 
Use of \code{pdredge} should be considered only with large datasets and complex 
	models, for which the standard version takes long time to complete.
	In other cases \code{pdredge} may even perform slower than \code{dredge}.
}


\value{
 See \code{\link{dredge}}.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
	\code{makeCluster} and other functions in packages
	\pkg{parallel} or \pkg{snow}.
}


\examples{


# Normally this should be simply "require(parallel) || require(snow)"
# but this is an ugly trick to avoid both R-check's complaining, and including 
# the dependency on these packages (it is just temporarily while 'pdredge' 
# is in experimental stage):
if (do.call("require", list("parallel", quietly = TRUE)) ||
	do.call("require", list("snow", quietly = TRUE))) {

# Example(Beetle)

data(Beetle)

Beetle100 <- Beetle[sample(nrow(Beetle), 100, replace = TRUE),]

fm1 <- glm(Prop ~ dose + I(dose^2) + log(dose) + I(log(dose)^2),
    data = Beetle100, family = binomial)

msubset <- expression(xor(dose, `log(dose)`) & (dose | !`I(dose^2)`)
    & (`log(dose)` | !`I(log(dose)^2)`))
varying.link <- list(family = alist(logit = binomial("logit"),
	probit = binomial("probit"), cloglog = binomial("cloglog") ))

# Set up the cluster
clust <- makeCluster(getOption("cl.cores", 2))
clusterExport(clust, "Beetle100")

# improvement only noticeable when data has about 3000 rows (on Windows 2-core 
# machine)
print(system.time(dredge(fm1, subset = msubset, varying = varying.link)))
print(system.time(pdredge(fm1, cluster = FALSE, subset = msubset, 
	varying = varying.link)))
print(system.time(pdredge(fm1, cluster = clust, subset = msubset, 
	varying = varying.link)))

stopCluster(clust)

}

}

\keyword{models}
