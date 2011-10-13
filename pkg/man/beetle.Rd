\name{Beetle}
\alias{Beetle}
\encoding{utf-8}
\docType{data}
\title{Flour beetle mortality data}
\description{
Mortality of flour beetles (\emph{Tribolium confusu}) due to exposure to gaseous
carbon disulfide (CS_{2}), from Bliss (1935)
}

\usage{
data(Beetle)
}

\format{
\code{Beetle} is a data frame with 4 variables.

\describe{
   \item{dose}{The dose of CS_{2} in mg/L}
   \item{n.tested}{Number of beetles tested}
   \item{n.killed}{Number of beetles killed}
   \item{mortality}{Observed mortality rate}

}}

\source{
Bliss C. I. (1935) The calculation of the dosage-mortality curve.
\emph{Annals of Applied Biology}, 22: 134–167
}

\references{
Burnham, K. P. and Anderson, D. R. (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed.
}

\author{Kamil Barto\enc{ń}{n}}

\examples{

# This follows the "logistic regression example"
# from Burnham & Anderson (2002), section 4.11

data(Beetle)

# Fit a global model with all the considered variables
globmod <- glm(cbind(n.killed, n.tested- n.killed)~dose + I(dose^2) + log(dose)
	+ I(log(dose)^2), data=Beetle, family=binomial)


# A logical expression defining the subset of models to use:
# * either log(dose) or dose
# * the quadratic terms can appear only together with linear terms
msubset <- expression(xor(dose, `log(dose)`) & (dose | !`I(dose^2)`)
	& (`log(dose)` | !`I(log(dose)^2)`))

# Note the use of 'alist' rather than 'list' in order to keep the 'family'
# objects unevaluated
alt.family <- list(family=alist(binomial("logit"), binomial("probit"),
	binomial("cloglog")))
# or, provide the alternatives as a named list
alt.family <- list(family=list(logit=binomial("logit"),
	probit=binomial("probit"), cloglog=binomial("cloglog")))


# Table 4.6
(dd12 <- dredge(globmod, subset = msubset, varying = alt.family))

# Table 4.7 "models justifiable a priori"
(dd3 <- dredge(update(globmod, . ~ dose), m.min = 1, varying = alt.family))

mod3 <- get.models(dd3, 1:3)

# Figure 4.3
newdata <- list(dose = seq(min(Beetle$dose), max(Beetle$dose), length.out = 25))
matplot(newdata$dose, sapply(mod3, predict, newdata, type="response"),
	type="l", xlab=quote(list("Dose of"~ CS[2],(mg/L))),
	ylab="Mortality"
)



}

\keyword{datasets}
