\name{stdize}
\alias{stdize}
\alias{stdize.default}
\alias{stdize.logical}
\alias{stdize.formula}
\alias{stdize.data.frame}
\encoding{utf-8}

\title{Standardize data}
\description{
Standardize variables by centering and scaling.
}

\usage{
%stdize(x, ...)

\method{stdize}{default}(x, center = TRUE, scale = TRUE, ...)

# also for two-level factors
\method{stdize}{logical}(x, binary = c("center", "scale", "binary", "half", "omit"),
    center = TRUE, scale = FALSE, ...)

\method{stdize}{data.frame}(x, binary = c("center", "scale", "binary", "half", "omit"),
    center = TRUE, scale = TRUE, omit.cols = NULL, source = NULL,
    prefix = TRUE, ...)

\method{stdize}{formula}(x, data = NULL, response = FALSE,
    binary = c("center", "scale", "binary", "half", "omit"),
    center = TRUE, scale = TRUE,
	omit.cols = NULL, prefix = TRUE,  ...)
}

\arguments{
	\item{x}{a numeric or logical vector, factor, numeric matrix, 
	     \code{data.frame} or a formula. }
	\item{center, scale}{ either a logical value or vector, or a numeric vector 
	     of length equal to the number of columns of \code{x}. 
		 \code{scale} can be also a function to be used for scaling. 
		 }
	 \item{source}{a reference \code{data.frame}, being a result of previous  
	     \code{stdize}, from which scale and center values are taken. Column 
		 names are matched. This is used for scaling new data using statistics
		 of another data.  
		  }
	\item{omit.cols}{ column names or numeric indices of columns that should 
	     be left unaltered. }
	\item{prefix}{ either a logical value specifying whether the names of
		transformed columns should be prefixed, or a two-element character vector
	    giving the prefixes. The prefixes default to \dQuote{z.} for scaled and
	    \dQuote{c.} for centred variables. }
	\item{binary}{specifies how binary variables (logical or two-level factors)
		are scaled. Default is to \code{"center"} them; use \code{"scale"} to
		both centre and scale by \acronym{SD}, \code{"binary"} to center to 0 /
		1, \code{"half"} to center to -0.5 / 0.5, and \code{"omit"} to leave
		binary variables unmodified. 
		If \code{center} or \code{scale} are given, this argument is ignored. 
		}
	\item{response}{logical, stating whether the response be standardized. By
	    default only variables on the right-hand side of formula are standardized. }
	\item{data}{ an object coercible to \code{data.frame}, containing the
	    variables in \code{formula}. Passed to, and used by \code{\link{model.frame}}. }
	\item{\dots}{ for the \code{formula} method, additional arguments passed to
	    \code{\link{model.frame}}. For other methods, silently ignored. }
}

\value{
Returns a vector or object of the same dimensions as \code{x}, where the values
are centred and/or scaled. Transformation is carried out column-wise in
\code{data.frame}s and matrices. 

If \code{center} or \code{scale} are logical scalars or vectors of length equal
to the number of columns of \code{x}, the centring is done by subtracting the
mean (if \code{center} corresponding to the column is \code{TRUE}), and scaling
is done by dividing the (centred) value by standard deviation (if corresponding
\code{scale} is \code{TRUE}). 
If \code{center} or \code{scale} are numeric vectors with length equal 
to the number of columns of \code{x} (or numeric scalars for vector methods), then 
these are used instead.

Binary variables, logical or factors with two levels, are are converted to
numeric variables and transformed according to the argument \code{binary},
unless \code{center} or \code{scale} are explicitly given.

The returned value is compatible with that of \code{\link{scale}} in that the
numeric centering and scalings used are stored in attributes attributes
\code{"scaled:center"} and \code{"scaled:scale"} (these can be \code{NA} if no
centering or scaling has been done).
}

\details{
\code{stdize} is very similar to \code{\link{scale}}, but uses special rules
for factors, similarly to \code{\link[arm]{standardize}} in package \pkg{arm}.

This function differs from \code{\link[arm]{standardize}} in that it is used on
data rather than on the fitted model object. The scaled data should afterwards
be passed to the modelling function, instead of the original data.

Unlike \code{standardize}, it applies special \sQuote{binary} scaling only to
two-level \code{factor}s and logical variables, rather than to any variable with
two unique values.

By default, \code{stdize} scales by dividing by \acronym{SD} rather than twice
the \acronym{SD} as \code{standardize} does. Unlike \code{scale}, scaling by
\acronym{SD} is used also on uncentered values, which is different from where
root-mean-square centering is done.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
Compare with \code{\link{scale}} or \code{\link[arm]{standardize}} and
\code{\link[arm]{rescale}} (the latter two in package \pkg{arm}).

\code{\link{apply}} and \code{\link{sweep}} for arbitrary transformations. 
}

\examples{
# compare "stdize" and "scale"
nmat <- matrix(runif(15, 0, 10), ncol = 3)

stdize(nmat)
scale(nmat)

rootmeansq <- function(v) {
    v <- v[!is.na(v)]
    sqrt(sum(v^2) / max(1, length(v) - 1L))
}

scale(nmat, center = FALSE)
stdize(nmat, center = FALSE, scale = rootmeansq)

if(require(lme4)) {
# define scale function as twice the SD for compatibility with "standardize"
twosd <- function(v) 2 * sd(v, na.rm = TRUE)

# standardize data (scaled variables are prefixed with "z.")
z.CO2 <- stdize(uptake ~ conc + Plant, data = CO2, omit = "Plant", scale = twosd)
summary(z.CO2)

fmz <- lmer(uptake ~ z.conc + I(z.conc^2) + (1 | Plant), data = z.CO2)

# create new data as a sequence along "conc"
newdata <-  data.frame(conc = seq(min(CO2$conc), max(CO2$conc), length = 10))

# scale new data using scale and center of the original data: 
z.newdata <- stdize(newdata, z = z.CO2)

\dontshow{ if(require(graphics)) \{ }
# plot predictions against "conc" on real scale:
plot(newdata$conc, predict(fmz, z.newdata, re.form = NA))
\dontshow{ \} }

# compare with "arm::standardize"
\dontrun{
library(arm)
    fms <- standardize(lmer(uptake ~ conc + I(conc^2) + (1 | Plant), data = CO2))
    plot(newdata$conc, predict(fms, z.newdata, re.form = NA))
}
}
}

\keyword{models}
