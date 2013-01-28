\name{updateable}
\alias{updateable}
\encoding{utf-8}

\title{Make a function returning updateable result}
\description{
Creates a function wrapper that stores a call in the values returned by its
argument \code{FUN}. }

\usage{
updateable(FUN) 
}

\arguments{
    \item{FUN}{function to be modified, found via \code{\link{match.fun}}. }
}

\details{
	Most model fitting functions in \R returns an object that can be updated or
	re-fitted via \code{\link{update}}. This is thanks to the \code{call} stored
	in the object, which can be used (possibly modified) later on. It is also
	utilised by \code{dredge} to generate sub-models.
	
	Some functions (such as \code{gamm} or \code{MCMCglmm}) do not provide their
	result with the \code{call} element. In this case \code{updateable} can be
	used on that function to add it. The resulting wrapper should be used in
	exactly the same way as the original function.
	
}

\value{
   A function with the same arguments as \code{FUN}, wrapping a call to
   \code{FUN} and adding an element named \code{call} to its result if possible,
   or an attribute \code{"call"} (if the returned value is atomic or a formal S4
   object).
}

\note{
  \pkg{MuMIn} has a replacement default method for \code{getCall} (original
  being from the \pkg{stats} package) that can recognize the call as an
  \code{\link[=attr]{attribute}} (rather than an element of the object). 
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{update}}, \code{\link{getCall}}, \code{\link{getElement}},
\code{\link{attributes}}

\code{\link{gamm.wrapper}} for an updateable wrapper for \code{gamm} models. 

}

\examples{

# Simple (and not very useful) example with cor.test:

# From example(cor.test)
x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
y <- c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)

ct1 <- cor.test(x, y, method = "kendall", alternative = "greater")

uCor.test <- updateable(cor.test)

ct2 <- uCor.test(x, y, method = "kendall", alternative = "greater")

getCall(ct1) # --> NULL
getCall(ct2)

#update(ct1, method = "pearson") --> Error
update(ct2, method = "pearson")
update(ct2, alternative = "two.sided")



\dontrun{
library(caper)
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)

fm1 <- crunch(Egg.Mass ~ F.Mass * M.Mass, data = shorebird)

uCrunch <- updateable(crunch)

fm2 <- uCrunch(Egg.Mass ~ F.Mass * M.Mass, data = shorebird)

getCall(fm1)
getCall(fm2)
update(fm2) # Error with 'fm1'
dredge(fm2)

}



}

\keyword{utils}
