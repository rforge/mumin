\name{model.avg}
\alias{model.avg}
\alias{print.averaging}
\encoding{utf-8}

\title{Model averaging}
\description{
Model averaging based on an information criterion.
}

\usage{
model.avg(object, ..., beta = FALSE, method = c("0", "NA"),
	method.var = c("NA", "0"), rank = NULL, rank.args = NULL,
	revised.var = TRUE)
}

\arguments{
    \item{object}{A fitted model object or a list of such objects. Alternatively
        an object of class \code{model.selection}. See \sQuote{Details}.}

    \item{\dots}{ more fitted model objects }

    \item{beta}{Logical, should standardized coefficients be returned?}

    \item{method, method.var}{The method of averaging estimators and
		their variance for parameters that are not common for all the models.
		Either \code{"0"} (default) or \code{"NA"}. See \sQuote{Details}.}

    \item{rank}{Optional, custom rank function (information criterion) to use
        instead of \code{AICc}, e.g. \code{BIC} or \code{QAIC}, may be omitted
        if \code{object} is a model list returned by \code{get.models}
        or a \code{model.selection} object.
        See \sQuote{Details}.}

    \item{rank.args}{Optional \code{list} of arguments for the \code{rank}
        function. If one is an expression, an \code{x} within it is substituted
        with a current model.}

    \item{revised.var}{Logical, indicating whether to use revised formula for
        standard errors. See \code{\link{par.avg}}.}
}

\value{

An object of class \code{averaging} is a list with components:

List of 13
\item{summary}{ a \code{data.frame} with deviance, AICc, Delta and weights
    for the component models }
\item{coefficients, se, dfs}{ matrices of component models' coefficients, their
	standard errors, and degrees of freedom }
\item{variable.codes}{ names of the variables with numerical codes used in 
    \code{summary} }
\item{avg.model}{ the model averaged parameters. A \code{data.frame} containing
    averaged coefficients, unconditional standard error, adjusted SE (if 
	\emph{df}s are available) and z-values (coefficient / SE) and 
	significance (assuming a normal error distribution) }
\item{importance}{ relative importance of the predictor variables,
    calculated as a sum of the Akaike weights over all of the models in
    which the parameter of interest appears }
\item{term.names}{ character vector giving names of all terms in the model }
\item{x, formula}{ the model matrix and formula corresponding to the one that 
	would be used in a single model }
\item{residuals}{ model averaged residuals (response minus fitted values) }
\item{call}{ the matched call. }

In addition, the object has following attributes:
\item{mList}{ a list of component model objects }
\item{method}{ a character vector of length 2 (for coefficient and variance), 
	describing how the missing parameters were handled ("NA" or "0").  }
\item{beta}{ logical, were standardized coefficients used? }
\item{revised.var}{ if TRUE, the standard errors were calculated with the 
	revised formula (See \code{\link{par.avg}}) }

}

\details{

\code{model.avg} has been tested to work with the following model classes:
	\code{lm}, \code{glm};
	\code{gam}, \code{gamm} (\pkg{mgcv}); 
	\code{gamm4} (\pkg{gamm4});
	\code{lme}, \code{gls} (\pkg{nlme});
	\code{lmer} (\pkg{lme4});
	\code{rlm}, \code{glm.nb}, \code{polr} (\pkg{MASS});
	\code{multinom} (\pkg{nnet});
	\code{sarlm}, \code{spautolm} (\pkg{spdep});
	\code{glmmML} (\pkg{glmmML});
	\code{coxph} (\pkg{survival});
	and \code{unmarkedFit} (\pkg{unmarked}).
Other model types are also likely to be supported, in particular those
    inheriting from one of the above classes. See package vignette 
	\sQuote{Extending \pkg{MuMIn}'s functionality} for a demonstration on how to
	provide support for other types of models.

\code{model.avg} may be used with a list of models, but also directly a
	\code{model.selection} object returned by \code{dredge}. In the latter case,
	the models from the model selection table are evaluated (with a call to 
	\code{get.models}) prior to averaging.
    A warning is given if the \code{subset} argument is not provided, and the 
	default \code{delta <= 4} will be used.

With \code{method = "0"} (default) all predictors are averaged as if they were
    present in all models in the set, and the value of parameter estimate is
    taken to be 0 if it is not present in a particular model.
	If \code{method = "NA"}, the predictors are averaged only over the models in
	which they appear, which biases them away from zero.
	Analogically, the argument \code{method.var} defines the method of
	calculating unconditional variance estimators (in this case the default
	is \code{"NA"}). 
	Example in \code{link{par.avg}} shows comparison of the two methods for both
	the coefficients and variance.

\code{rank} is found by a call to \code{\link[base]{match.fun}} and typically is
    specified as a function or a symbol (e.g. a back-quoted name) or a character
    string specifying a function to be searched for from the environment of the
    call to lapply. \code{rank} must be a function able to accept model as a
    first argument and must always return a scalar.

Some generic methods such as \code{\link{predict.averaging}}, \code{\link{coef}},
    \code{\link[stats]{formula}},
    \code{\link[stats]{residuals}} and \code{\link[stats]{vcov}} are supported. \cr
    \code{logLik} method returns a list of \code{\link[stats]{logLik}} objects
    for the component models.
}

\note{

From version 1.0.1, \code{print} method provides only a concise output
(similarly as for \code{lm}). To print a full summary of the results use
\code{summary} function. Confidence intervals can be obtained with
\code{\link{confint}}.

}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. }

\author{Kamil Barto\enc{ń}{n}}

\seealso{

See \code{\link{par.avg}} for details of model averaged parameter calculation.

\code{\link{dredge}}, \code{\link{get.models}}. \cr
\code{\link{AICc}} has examples of averaging models fitted by REML.

\code{\link[AICcmodavg]{modavg}} in package \pkg{AICcmodavg}, and
\code{\link[glmulti]{coef.glmulti}} in package \pkg{glmulti} also perform model
averaging.
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
\dontrun{
# The same result
model.avg(dd, cumsum(weight) <= .95)
}

\dontrun{
# using BIC (Schwarz's Bayesian criterion) to rank the models
BIC <- function(x) AIC(x, k=log(length(residuals(x))))
mav <- model.avg(top.models, rank=BIC)
}

}

\keyword{models}
