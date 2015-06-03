\name{std.coef}
\alias{std.coef}
\alias{beta.weights}
\alias{partial.sd}

\encoding{utf-8}

\title{Standardized model coefficients}
\description{
Standardize model coefficients by Standard Deviation or Partial Standard Deviation.
}

\usage{

std.coef(x, partial.sd, ...)

partial.sd(x)

# Deprecated:
beta.weights(model)
}
\arguments{
  \item{x, model}{a fitted model object. }
  \item{partial.sd}{logical, if set to \code{TRUE}, model coefficients are multiplied by
	Partial \acronym{SD}, otherwise they are multiplied by the ratio of the standard
	deviations of the independent variable and dependent variable.
	}
 \item{\dots}{additional arguments passed to \code{\link{coefTable}}, 
    e.g. \code{dispersion}. }
 }

\details{
	Standardizing coefficients has the same effect as centering and scaling the
	input variables. \dQuote{Classical} standardized coefficients are calculated
	as: \eqn{\beta^{*}_i = \beta_i\frac{s_{X_{i}}}{s_{y}}}, where \eqn{\beta} is
	the unstandardized coefficient, \eqn{s_{X_{i}}} is the standard deviation of
	associated depenent variable \eqn{X_{i}} and \eqn{s_{y}} is \acronym{SD} of
	the response variable.
	
	The usual standardized coefficients are difficult to interpret if the variables
	are intercorrelated. The standard deviation of \eqn{X_i} used in computing the
	standardized coefficients \eqn{\beta_i^{*}} should be replaced by a partial standard
	deviation of \eqn{X_i} which is adjusted for the multiple correlation of
	\eqn{X_i} with the other \eqn{X} variables included in the regression equation.
	The partial standard deviation is calculated as
	\eqn{s_{X_{i}}^{*} = s_{X_{i}} \var{VIF}(X_i)^{-0.5} (\frac{n - 1}{n -
	p} )^{0.5}}, where \var{VIF} is the variance inflation factor, \var{n} is
	the number of observations and \var{p} number of predictors in the model.
	Coefficient is then transformed as \eqn{\beta^{*}_i = \beta_i s_{X_{i}}^{*}}.

}

\value{
   A matrix with at least two columns for standardized coefficient estimate and
   its standard error. Optionally, third column holds degrees of
   freedom associated with the coefficients.
}


\author{Kamil Barto\enc{ń}{n}

The functions include code for calculating variance-inflation factors based on
function \code{vif} from package \pkg{car} written by Henric Nilsson and John
Fox.
}

\references{

Cade, B.S. (in press) Model averaging and muddled multimodel inference. \emph{Ecology} \cr
http://dx.doi.org/10.1890/14-1639.1

Afifi A., May S., Clark V.A. (2011) \emph{Practical Multivariate Analysis},
Fifth Edition. CRC Press.

Bring, J. (1994). How to standardize regression coefficients. \emph{The American
Statistician} 48, 209-213.
}

\seealso{
\code{partial.sd} can be used with \code{\link{stdize}}.

\code{\link{coef}} or \code{\link{coeffs}} and \code{\link{coefTable}} for
unstandardized coefficients.
}
\examples{

# partial SD for the default formula: y ~ x1 + x2 + x3 + x4
psd <- partial.sd(lm(GPA))[-1] # remove first element for intercept

# first element of 'scale' is set to NA to ignore the first column 'y'
zGPA <- stdize(GPA, scale = c(NA, psd), center = TRUE)
            
fmz <- lm(c.y ~ z.x1 + z.x2 + z.x3 + z.x4, data = zGPA)
fm <- lm(y ~ x1 + x2 + x3 + x4, data = GPA)

# standardized data:
zapsmall(coefTable(fmz))
# standardized coefficients:
std.coef(fm, partial = TRUE)


# standardizing nonlinear models:
fam <- Gamma("inverse")
fmg <- glm(log(y) ~ x1 + x2 + x3 + x4, GPA, family = fam)

psdg <- partial.sd(fmg)
zGPA <- stdize(GPA, scale = c(NA, psdg[-1]), center = FALSE)
fmgz <- glm(log(y) ~ z.x1 + z.x2 + z.x3 + z.x4, zGPA, family = fam)

# standardized data:
coef(fmgz)
# (intercept is unchanged because the variables haven't been centered)
# standardized coefficients:
coef(fmg) * psdg

}


\keyword{models}
