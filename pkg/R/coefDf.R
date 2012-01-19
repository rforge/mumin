`coefDf` <- function(x) UseMethod("coefDf")
`coefDf.lme` <- function(x) x$fixDF$X
`coefDf.mer` <- function(x) rep(NA, x@dims[["p"]])
`coefDf.gls` <- function(x) rep(x$dims$N - x$dims$p, x$dims$p)
`coefDf.default` <- function(x) {
	dfres <- tryCatch(df.residual(x), error = function(e) NA)
	if(is.null(dfres)) dfres <- NA_integer_
	rep(dfres, length(coeffs(x)))
}
