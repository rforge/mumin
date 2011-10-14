\name{par.avg}
\alias{par.avg}
\encoding{utf-8}

\title{Parameter averaging}
\description{
Averages single model coefficient based on provided weights
}

\usage{
par.avg(x, se, weight, df = NULL, level = 1 - alpha, alpha = 0.05,
    revised.var = TRUE)
}

\arguments{
  \item{x}{vector of parameters}
  \item{se}{vector of standard errors}
  \item{weight}{vector of weights}
  \item{df}{(optional) vector of degrees of freedom}
  \item{alpha, level}{significance level for calculating confidence intervals}
  \item{revised.var}{logical, should the revised formula for standard errors be
	used? See \sQuote{Details}.}
}

\value{
  \code{par.avg} returns a vector with named elements:
  \item{Coefficient}{model coefficients}
  \item{SE}{unconditional standard error}
  \item{Adjusted SE}{adjusted standard error}
  \item{Lower CI, Upper CI}{unconditional confidence intervals}
}

\details{
Unconditional standard errors are square root of the variance estimator,
calculated either according to the original formula in
Burnham and Anderson (2002, p. 160, equation 4.7),
or a newer, revised formula from Burnham and Anderson (2004, equation 4)
(if \code{revised.var = TRUE}, this is the default).
If degrees of freedom are given, the confidence intervals are based on adjusted
standard error estimator (Burnham and Anderson 2002, page 164).
}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed.

Burnham, K. P. and Anderson, D. R. (2004). \emph{Multimodel inference -
understanding AIC and BIC in model selection}. Sociological Methods & Research
33(2): 261-304.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{model.avg}} for model averaging.
}

\examples{


# Compare estimates and variance, model averaged with different methods
# ("NA", "0", see ?model.avg)

data(Cement)
fm1 <- lm(y ~ ., data = Cement)

avgmod <- model.avg(dredge(fm1), cumsum(weight) <= .95, method="NA")
weight <- avgmod$summary$Weight

na2z <- function(x) ifelse(is.na(x), 0, x) #helper function (changes NA to 0)

ret <- sapply(colnames(avgmod$coefficients), function(i) {
	x <- avgmod$coefficients[, i]
	se <- avgmod$se[, i]
	df <- avgmod$dfs[, i]
	rbind(
	'coef(x)se(x)' = par.avg(x, se, weight, df),
	'coef(0)se(x)' = par.avg(na2z(x), se, weight, df),
	'coef(0)se(0)' = par.avg(na2z(x), na2z(se), weight, df)
	)
}, simplify=FALSE)


if(require(graphics)) {

# Plot the comparison: coefficients +- confidence intervals
y <- 1:3
cap <- .1
op <- par(mfrow=c(2,3), mar=c(3,5,1,3), oma=c(3,0,3,0))
lapply(names(ret), function(i) {
	x <- ret[[i]]
	plot(y, x[,1], xlim=c(.5,3.5), ylim=range(x[,c(1,4,5)]), pch=20, cex=2,
		ylab=i, xlab=NA, axes=FALSE, col=c(1,2,2))
	segments(y, x[,4], y, x[,5], lty=2, col=c(1,1,2))
	segments(y - cap, x[,4], y + cap, x[,4], col=c(1,1,2))
	segments(y - cap, x[,5], y + cap, x[,5], col=c(1,1,2))
	box(); axis(2); axis(1, at=y, labels=rownames(x))
	abline(h=0, lty=3)
})
title("Comparison of model averaging methods",
	sub="Cement hardening example", outer=TRUE, line=1)
par(op)

}




}


\keyword{models}
