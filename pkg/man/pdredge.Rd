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
	\item{global.model, beta, evaluate rank, fixed, m.max, m.min, subset}{}
	\item{marg.ex, trace, varying, extra, ...}{see \code{\link{dredge}}.
	}

	\item{cluster}{ either a valid cluster object, or \code{FALSE}. }

	\item{check}{ logical, whether to evaluate the \code{global.model} in the
	cluster	and compare with the original one. Only a simple checking is done. }

}

\details{
All the dependencies for fitting the \code{global.model} (including the data,
all required packages, and any object the modelling function may use) must be
exported/loaded into the cluster before passing it to \code{pdredge}.
}


\value{
 See \code{\link{dredge}}.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
	\code{\link[parallel]{makeCluster}} and other functions in package
	\pkg{parallel}.
}


\examples{

if (require(parallel)) {

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

system.time(dredge(fm1, subset = msubset, varying = varying.link))
system.time(pdredge(fm1, cluster = FALSE, subset = msubset, varying = varying.link))
system.time(pdredge(fm1, cluster = clust, subset = msubset, varying = varying.link))

stopCluster(clust)
}
}

\keyword{models}
