\name{model.avg}
\alias{model.avg}
\alias{print.averaging}
\alias{coef.averaging}
\alias{predict.averaging}
\encoding{utf-8}

\title{Model averaging}
\description{
Model averaging based on an information criterion.
}

\usage{
model.avg(m1, ..., beta = FALSE, method = c("0", "NA"), rank = NULL,
	rank.args = NULL, alpha = 0.05)

\method{coef}{averaging}(object, ...)

\method{predict}{averaging}(object, newdata, se.fit = NULL, interval = NULL,
	type=NULL, ...)
}


\arguments{
	\item{m1}{A fitted model object or a list of such objects.}
	\item{beta}{Logical, should standardized coefficients be returned rather
		than normal ones?}
	\item{method}{If set to \dQuote{0} (default), terms missing in one model are
		assumed to be 0's, otherwise they are omitted from the weighted
		average.}
	\item{rank}{Custom rank function (information criterion) to use instead of
	\code{AICc}, e.g. \code{QAIC} or \code{BIC}, may be omitted if \code{m1} is
		a list returned by \code{dredge}. See \sQuote{Details}.}
	\item{rank.args}{Optional \code{list} of arguments for the \code{rank}
		function.}
	\item{alpha}{Significance level for calculatinq confidence intervals.}
	\item{object}{An object returned by \code{model.avg}.}
	\item{newdata}{ An optional data frame in which to look for variables with
		which to predict. If omitted, the fitted values are used.}
	\item{se.fit,interval}{Currently not used.}
	\item{type}{Ignored. Only predictions on the link scale are allowed.
		Warning is given if user tries something else here. }
	\item{\dots}{ for \code{model.avg} - more fitted model objects, for
		\code{predict} - arguments to be passed to respective \code{predict}
		method }

}

\value{
	\code{model.avg} returns a list with elements:

	\item{summary}{Model table with deviance, AICc, Delta and weight.}
	\item{coefficients}{the model coefficients}
	\item{variance}{variance of coefficients}
	\item{avg.model}{averaged model summary (\code{\link{data.frame}} with
	  columns: coef - averaged coefficients, var - unconditional variance estimator,
	  ase - adjusted standard error estimator, lci, uci - unconditional confidence
	  intervals)}
	\item{relative.importance}{relative variable importances}
	\item{variable.codes}{Variable names with numerical codes used in the summary}
	\item{relative.importance}{Relative importance of variables}
	\item{weights}{}
	\item{beta}{(logical) were standardized coefficients used?}
	\item{model}{the model matrix, analogical to one that would be used in a
	  single model.}
	\item{residuals}{ the residuals (response minus fitted values).}

}

\details{
\code{rank} is found by a call to \code{match.fun} and typically is specified as
a function or a symbol (e.g. a backquoted name) or a character string specifying
a function to be searched for from the environment of the call to lapply.\cr
Function \code{rank} must be able to accept model as a first argument and must
always return a scalar. \\
Apart from \code{predict} and \code{coef}, other default methods, such as
 \code{formula} and \code{residuals} may be used.
}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. }

\author{Kamil Bartoń}

\seealso{
\code{\link{dredge}}, \code{\link{get.models}}.
\code{\link{QAIC}} has examples of using custom rank function.
}

\examples{
# Example from Burnham and Anderson (2002), page 100:
data(Cement)
lm1 <- lm(y ~ ., data = Cement)
dd <- dredge(lm1)
dd

#models with delta.aicc < 4
model.avg(get.models(dd, subset = delta < 4)) # get averaged coefficients

#or as a 95\% confidence set:
top.models <- get.models(dd, cumsum(weight) <= .95)

model.avg(top.models) # get averaged coefficients

#topmost model:
top.models[[1]]

\dontrun{
# using BIC (Schwarz's Bayesian criterion) to rank the models
BIC <- function(x) AIC(x, k = log(length(residuals(x))))
mav <- model.avg(top.models, rank=BIC)
}


# Predicted values
nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE),
	length=len)

# New predictors: X1 along the range of original data, other variables held
#	constant at their means
newdata <- as.data.frame(lapply(lapply(Cement[1:5], mean), rep, 25))
newdata$X1 <- nseq(Cement$X1, 25)

# Predictions from each of the models in a set:
pred <- sapply(top.models, predict, newdata=newdata)
# Add predictions from the averaged model:
pred <- cbind(pred, averaged=predict(model.avg(top.models), newdata))

matplot(x = newdata$X1, y = pred, type="l", lwd=c(rep(1,ncol(pred)-1), 2))

legend("topleft", legend=c(lapply(top.models, formula), "Averaged model"),
	   col=seq(ncol(pred)), lty=1:5, lwd=c(rep(1,ncol(pred)-1), 2))

}

\keyword{models}
