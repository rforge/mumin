\name{r.squaredLR}
\alias{r.squaredLR}
\alias{null.fit}

\title{likelihood-ratio based pseudo R-squared}
\description{
Calculate likelihood-ratio based coefficient of determination,
\eqn{R_{LR}^{2}}{pseudo R^2}
}

\usage{
r.squaredLR(x, null = null.fit(x, TRUE))

null.fit(x, evaluate = FALSE,
    envir = environment(as.formula(formula(x))))
}

\arguments{
  \item{x}{a fitted model object}
  \item{null}{a fitted \emph{null} model, if not provided, a \code{glm} with
  only intercept and appropriate family will be used }
  \item{evaluate}{If \code{TRUE} returns the fitted model object else return the
  call }
  \item{envir}{the environment in which the \emph{null} model is to be
  evaluated, by default an environment of the original model's formula. }
}

\value{
	\code{r.squaredLR} returns a value of \eqn{R_{LR}^{2}}{R^2}, with additional
	attribute \code{"adj.r.squared"} gives the modification proposed by
	Nagelkerke (1991). Note that this is not the same as the classical
	\sQuote{adjusted R squared}.

	\code{null.fit} returns the fitted \emph{null} model object (if evaluate =
	TRUE) or an unevaluated call to fit a \emph{null} model.

}

\details{

This is one of the several pseudo R-squared goodness-of-fit statistics that have
been proposed (mainly with application in logistic regression models). It is
based on an improvement from \emph{null} (intercept only) model to the fitted
model, and calculated as
\deqn{
R_{LR}^{2}=1-\exp(-\frac{2}{n}(\mathit{logLik}(x)-\mathit{logLik}(0)))
}{R^2 = 1 - exp(-2/n * logLik(x) - logLik(0))}
where \emph{logLik(x)} and \emph{logLik(0)} are the log likelihoods of the fitted
and the \emph{null} model respectively.


For ordinary linear models \eqn{R_{LR}^{2}}{R^2} is consistent with classical
\eqn{R^{2}}. For discrete models, which likelihood is a product of probabilities
(such as logistic regression), the maximum \eqn{R_{LR}^{2}}{R^2} is less than 1, and
equals
\deqn{
max(R_{LR}^{2}) = 1 - \exp(\frac{2}{n}\mathit{logLik}(\textrm{0}))
}{max(R^2) = 1 - exp(2 / n * \emph{logLik}(0))}

The modification proposed by Nagelkerke (1991) adjusts the \eqn{R_{LR}^{2}} to
achieve 1 at its maximum:
\eqn{
\bar{R}^{2} = R_{LR}^{2} / \max(R_{LR}^{2})
}{
Radj^2 = R^2 / max(R^2)
}

\code{null.fit} tries to guess an intercept only (\emph{null}) model (currently
only as a \code{glm}), given the provided fitted model object. It is primarily
intended for internal use.

}

\note{
Although this statistic is not related to information theory, it is provided
here because it may be often useful to compare it with the value of information
criterions.
}

\references{
Cox, D. R. and Snell, E. J. (1989) \emph{The analysis of binary Data}, 2nd ed.
London, Chapman and Hall.

Nagelkerke, N. J. D. (1991) A note on a general definition of the coefficient of
determination. \emph{Biometrika} 78: 691-692 (and references for original
formulation therein)
}

\seealso{
\code{\link{summary.lm}}
}

\keyword{models}
