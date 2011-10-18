\name{importance}
\alias{importance}
\encoding{utf-8}

\title{Relative variable importance}
\description{
Sum of \sQuote{Akaike weights} over all models including the explanatory
variable
}

\usage{
importance(x)
}

\arguments{
	\item{x}{Either a list of fitted model objects, or a \code{"model.selection"}
	or \code{"averaging"} object }

}

\value{
	a numeric vector of relative importance values, named as the predictor
	variables
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}}, \code{\link{model.avg}}, \code{\link{mod.sel}}
}

\examples{

# Generate some models
data(Cement)
lm1 <- lm(y ~ ., data = Cement)
dd <- dredge(lm1)

# Importance can be calculated/extracted from various objects:
importance(dd)
\dontrun{
importance(subset(mod.sel(dd), delta <= 4))
importance(model.avg(dd, subset = delta <= 4))
importance(subset(dd, delta <= 4))
importance(get.models(dd, delta <= 4))
}

# Re-evaluate the importances according to BIC
# note that re-ranking involves fitting the models again

# 'nobs' is not used here for backwards compatibility
lognobs <- log(length(resid(lm1)))

importance(subset(mod.sel(dd, rank=AIC, rank.args=list(k = lognobs)),
	cumsum(weight) <= .95))

# This gives a different result than previous command, because 'subset' is
# applied to the original selection table that is ranked with 'AICc'
importance(model.avg(dd, rank=AIC, rank.args=list(k = lognobs),
	subset=cumsum(weight) <= .95))

}

\keyword{models}
