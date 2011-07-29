\name{model.avg}
\alias{model.avg}
\alias{print.averaging}
\encoding{utf-8}

\title{Model averaging}
\description{
Model averaging based on an information criterion.
}

\usage{
model.avg(object, ..., beta = FALSE, method = c("0", "NA"), rank = NULL,
	rank.args = NULL, revised.var = TRUE)
}

\arguments{
	\item{object}{A fitted model object or a list of such objects. Alternatively 
		an object of class "model.selection". See \sQuote{Details}.}
		
	\item{\dots}{ more fitted model objects }

	\item{beta}{Logical, should standardized coefficients be returned?}
	
	\item{method}{The method of averaging parameter estimators that are not
		common for all the models. Either \code{"0"} (default) or \code{"NA"}.
		See \sQuote{Details}.}

	\item{rank}{Optional, custom rank function (information criterion) to use
		instead of \code{AICc}, e.g. \code{QAIC} or \code{BIC}, may be omitted
		if \code{m1} is a model list returned by \code{get.models}. See
		\sQuote{Details}.}

	\item{rank.args}{Optional \code{list} of arguments for the \code{rank}
		function. If one is an expression, an \code{x} within it is substituted
		with a current model.}

	\item{revised.var}{Logical, indicating whether to use revised formula for
		standard errors. See \code{\link{par.avg}}.}
}

\value{

An object of class \code{averaging} with following elements:

 \item{summary}{a \code{data.frame} with deviance, AICc, Delta and weights
	for the component models.}
 \item{coefficients, variance}{matrices of component models' coefficients and
	their variances}
 \item{variable.codes}{names of the variables with numerical codes used in the
	summary}
 \item{avg.model}{the averaged model summary, (\code{data.frame} containing
	averaged coefficients, unconditional standard error, adjusted SE,
	and confidence intervals)}
 \item{importance}{the relative importance of variables}
 \item{beta}{(logical) were standardized coefficients used?}
 \item{term.names}{character vector giving names of all terms in the model}
 \item{residuals}{the residuals (response minus fitted values).}
 \item{x, formula}{the model matrix and formula analogical to those that would
	be used in a single model.}
 \item{method}{how the missing terms were handled ("NA" or "0").}
 \item{call}{the matched call.}

}

\details{

\code{model.avg} has been tested to work with the following model classes:
\itemize{
	\item \code{lm}, \code{glm}
	\item \code{gam} (\pkg{mgcv})
	\item \code{lme}, \code{gls} (\pkg{nlme})
	\item \code{lmer} (\pkg{lme4})
	\item \code{rlm}, \code{glm.nb} (\pkg{MASS})
	\item \code{multinom} (\pkg{nnet})
	\item \code{sarlm}, \code{spautolm} (\pkg{spdep})
	\item \code{coxph} (\pkg{survival})
}

\code{model.avg} may be used with a list of models, or an object returned by 
	\code{dredge}). In the latter case, the models from the model selection table 
	are evaluated (with a call to \code{get.models}) prior to averaging.
	A warning is given if the subset argument is not provided, and the default 
	\code{delta <= 4} will be used.

Other model types are also likely to be supported, in particular those
	inheriting from one of the above classes. See \sQuote{Details} section of
	the \sQuote{\link{Miscellaneous}} page to see how to provide support for
	other types of models.

With \code{method = "0"} (default) all predictors are averaged as if they were
	present in all models in the set, and the value of parameter estimate is
	taken to be 0 if it is not present in a particular model. If \code{method =
	"NA"}, the predictors are averaged only over the models in which they
	appear.

\code{rank} is found by a call to \code{\link[base]{match.fun}} and typically is
	specified as a function or a symbol (e.g. a back-quoted name) or a character
	string specifying a function to be searched for from the environment of the
	call to lapply. \code{rank} must be a function able to accept model as a
	first argument and must always return a scalar.

Some generic methods such as \code{\link[MuMIn]{predict}}, \code{\link{coef}},
	\code{\link[stats]{formula}}, 
	\code{\link[stats]{residuals}} and \code{\link[stats]{vcov}} are supported. \cr
	\code{logLik} method returns a list of \code{\link[stats]{logLik}} objects
	for the component models.
}

\note{

From version 1.0.1, \code{print} method provides only a concise output
(similarly as for \code{lm}), to print a full summary of the results use
\code{summary} function. Confidence intervals can be obtained with
\code{\link{confint}}.

}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. }

\author{Kamil \enc{Bartoń}{Barton}}

\seealso{

See \code{\link{par.avg}} for details of averaged model calculation.

\code{\link{dredge}}, \code{\link{get.models}}.
\code{\link{QAIC}} has examples of using custom rank function and prediction
with confidence intervals. \cr
\code{\link{AICc}} has examples of averaging models fitted by REML.


\code{\link[AICcmodavg]{modavg}} in package \pkg{AICcmodavg}, and
\code{\link[glmulti]{coef.glmulti}}
in package \pkg{glmulti} also perform model averaging.
}

\examples{

require(graphics)

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
BIC <- function(x) AIC(x, k=log(length(residuals(x))))
mav <- model.avg(top.models, rank=BIC)
}


# Predicted values
nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE),
	max(x, na.rm=TRUE),	length=len)

# New predictors: X1 along the range of original data, other variables held
#	constant at their means
newdata <- as.data.frame(lapply(lapply(Cement[1:5], mean), rep, 25))
newdata$X1 <- nseq(Cement$X1, nrow(newdata))

# Predictions from each of the models in a set:
pred <- sapply(top.models, predict, newdata=newdata)
# Add predictions from the models averaged using two methods:
pred <- cbind(pred,
	averaged.0=predict(model.avg(top.models, method="0"), newdata),
	averaged.NA=predict(model.avg(top.models, method="NA"), newdata)
	)

matplot(x=newdata$X1, y=pred, type="l", lwd=c(rep(1,ncol(pred)-2), 2, 2),
	xlab="X1", ylab="y")

legend("topleft",
	legend=c(lapply(top.models, formula),
		paste("Averaged model (method=", c("0", "NA"), ")", sep="")),
	col=1:6, lty=1:5, lwd=c(rep(1,ncol(pred)-2), 2, 2), cex = .75
	)


\dontrun{
# Example with gam models (based on "example(gam)")
require(mgcv)
dat <- gamSim(1, n = 500, dist="poisson", scale=0.1)

gam1 <- gam(y ~ s(x0) + s(x1) + s(x2) +  s(x3) + (x1+x2+x3)^2,
	family = poisson, data = dat, method = "REML")

cat(dQuote(getAllTerms(gam1)), "\n")

# include only models with smooth OR linear term (but not both)
# for each variable:
dd <- dredge(gam1, subset=xor(`s(x1)`, x1) & xor(`s(x2)`, x2) & xor(`s(x3)`, x3))
# ...this may take a while.

subset(dd, cumsum(weight) < .95)

top.models <- get.models(dd, cumsum(weight) <= .95)

newdata <- as.data.frame(lapply(lapply(dat, mean), rep, 50))
newdata$x1 <- nseq(dat$x1, nrow(newdata))
pred <- cbind(
	sapply(top.models, predict, newdata=newdata),
	averaged=predict(model.avg(top.models), newdata)
)

matplot(x=newdata$x1, y=pred, type="l", lwd=c(rep(1,ncol(pred)-2), 2, 2),
	xlab="x1", ylab="y")
}


}

\keyword{models}
