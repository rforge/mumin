# cbind list of data.frames omitting duplicated column (names)
`cbindDataFrameList` <-
function(x) {
	dfnames <- unlist(lapply(x, colnames))
	uq <- !duplicated(dfnames)
	res <- do.call("cbind", x)[,uq]
	colnames(res) <- dfnames[uq]
	return(res)
}

# same for rbind, check colnames and add NA's when any are missing
`rbindDataFrameList` <-
function(x) {
	all.colnames <- unique(unlist(lapply(x, colnames)))
	x <- lapply(x, function(y) {
		y[all.colnames[!(all.colnames %in% colnames(y))]] <- NA
		return(y[all.colnames])
	})
	return(do.call("rbind", x))
}

# test for marginality constraints
`formulaAllowed` <-
function(frm, except=NULL) {
	if(isTRUE(except)) return(TRUE)
	factors <- attr(terms(frm), "factors")
	if(length(factors) == 0) return(TRUE)
	if(is.character(except))
		factors <- factors[!(rownames(factors) %in% except), ]
	return(all(factors < 2))
}


# Calculate Akaike weights
`Weights` <-
function(aic, ...) {
	delta <- aic - min(aic)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	return (weight)
}


#sorts alphabetically interaction components in model term names
`fixCoefNames` <-
function(x) {
	if(!is.character(x)) return(x)
	return(sapply(lapply(strsplit(x, ":"), sort), paste, collapse=":"))
}


# logLik for survival::coxph model
# https://stat.ethz.ch/pipermail/r-help/2006-December/122118.html
# originally by Charles C. Berry, mod. by KB: correction for the null model
`logLik.coxph` <- function(object,...) {
# Thx to Mathieu Basille:
    y <-  object$loglik[length(object$loglik)]
	#y <-  if(length(object$loglik) > 1)
	#	-1 * (object$loglik[1] - object$loglik[2]) else object$loglik
    class(y) <- "logLik"
    attr(y,'df') <- if(is.null(object$coef)) 0 else sum(!is.na(object$coef))
    return(y)
}

if (!existsFunction("nobs"))
`nobs` <- function(object, ...) UseMethod("nobs")
`nobs.mer` <- function(object, ...) object@dims[["n"]]

`nobs.gls` <- function(object, nall = FALSE, ...) {
	p <- object$dims$p
	N <- object$dims$N
	if (nall) return (N)
	REML <- object$method == "REML"
	N - REML * p
}

`nobs.lme` <- function(object, nall = FALSE, ...) {
    p <- object$dims$ncol[object$dims$Q + 1]
	N <- object$dims$N
	if (nall) return (N)
	REML <- object$method == "REML"
	N - REML * p
}





`nobs.glmmML` <- function(object, ...) length(object$coefficients) + object$cluster.null.df
`nobs.default` <- function(object, ...) NROW(resid(object, ...))

`coefDf` <- function(x) UseMethod("coefDf")
`coefDf.lme` <- function(x) x$fixDF$X
`coefDf.mer` <- function(x) rep(NA, x@dims[["p"]])
`coefDf.gls` <- function(x) rep(x$dims$N - x$dims$p, x$dims$p)
`coefDf.default` <- function(x) rep(df.residual(x), length(coef(x)))
