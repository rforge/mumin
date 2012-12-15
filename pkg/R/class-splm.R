`logLik.splm` <- function (object, ...) {
	ret <- object$logLik
	if(is.null(ret)) return(NA)
	attr(ret, "nobs") <- length(resid(object))
	attr(ret, "df") <- length(object$coefficients) + length(object$errcomp) +
		length(object$arcoef) + 1L
	class(ret) <- "logLik"
	ret
}

`nobs.splm` <- function (object, ...) length(resid(object))
`coeffs.splm` <- function (model) {
	c(model$coefficients, model$arcoef,
	if(is.matrix(model$errcomp)) model$errcomp[, 1L] else model$errcomp)
}

`coefTable.splm` <- function (model, ...) {
	cf <- sapply(c("coefficients", "arcoef", "errcomp"), function(i)
		if(is.matrix(model[[i]])) model[[i]][, 1L] else model[[i]],
		simplify = FALSE)
	
	ncf <- sapply(cf, length)
	vcovlab <- c(coefficients = "vcov", arcoef = "vcov.arcoef", errcomp = "vcov.errcomp")
	se <- sqrt(unlist(lapply(names(vcovlab), function(i) {
		vcv2 <- diag(model[[vcovlab[i]]])
		c(vcv2, rep(NA_real_, ncf[[i]] - length(vcv2)))
	})))
	
	.makeCoefTable(unlist(cf, use.names = FALSE), se,
		coefNames = unlist(lapply(cf, names), use.names = FALSE))
}