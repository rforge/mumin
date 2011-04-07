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

# logLik for survival::coxph model
# https://stat.ethz.ch/pipermail/r-help/2006-December/122118.html
# originally by Charles C. Berry, mod. by KB: correction for the null model
`logLik.coxph` <- function(object,...) {
# Thx to Mathieu Basille:
    y <-  object$loglik[length(object$loglik)]
	#y <-  if(length(object$loglik) > 1)
	#	-1 * (object$loglik[1] - object$loglik[2]) else object$loglik
    class(y) <- "logLik"
	#attr(y,"nall") <-
	#attr(y,"nobs") <-
    attr(y,'df') <- if(is.null(object$coef)) 0 else sum(!is.na(object$coef))
    return(y)
}

if (!existsFunction("nobs"))
`nobs` <- function(object, ...) UseMethod("nobs")

`nobs.mer` <- function(object, nall = FALSE, ...) {
	N <- object@dims[["n"]]
	p <- object@dims[["p"]]
	if (nall) return (N)
	REML <- object@dims[['REML']]
	N - REML * p
}

`nobs.gls` <- function(object, nall = FALSE, ...) {
	p <- object$dims$p
	N <- object$dims$N
	if (nall) return (N)
	REML <- object$method == "REML"
	N - REML * p
	# p - the number of coefficients in the linear model.
}


`nobs.lme` <- function(object, nall = FALSE, ...) {
    p <- object$dims$ncol[object$dims$Q + 1]
	N <- object$dims$N
	if (nall) return (N)
	REML <- object$method == "REML"
	N - REML * p
	#N - the number of observations in the data,
	#Q - the number of grouping levels
	#ncol - the number of columns in the model matrix for each level of grouping from innermost to outermost
	#  (last two values are equal to the number of fixed effects and one).
}

`nobs.glmmML` <- function(object, ...) length(object$coefficients) + object$cluster.null.df
`nobs.default` <- function(object, ...) NROW(resid(object, ...))

`coefDf` <- function(x) UseMethod("coefDf")
`coefDf.lme` <- function(x) x$fixDF$X
`coefDf.mer` <- function(x) rep(NA, x@dims[["p"]])
`coefDf.gls` <- function(x) rep(x$dims$N - x$dims$p, x$dims$p)
`coefDf.default` <- function(x) rep(df.residual(x), length(coef(x)))


# Hidden functions

`.getLogLik` <- function()
	if ("stats4" %in% loadedNamespaces())
        stats4:::logLik else
		logLik

`.getCall` <- function(x) {
	if(mode(x) == "S4") {
		if ("call" %in% slotNames(x)) slot(x, "call") else
			NULL
	} else {
		if(!is.null(x$call)) {
			x$call
		} else if(!is.null(attr(x, "call"))) {
			attr(x, "call")
		} else
			NULL
	}
}

`.isREMLFit` <- function(x) {
	if (inherits(x, "mer")) return (x@dims[["REML"]] != 0)
	if (inherits(x, c("lme", "gls", "gam")) && !is.null(x$method))
		return(x$method %in% c("lme.REML", "REML"))
	if (any(inherits(x, c("lmer", "glmer"))))
		return(x@status["REML"] != 0)
	return(NA)
}


#sorts alphabetically interaction components in model term names
`fixCoefNames` <-
function(x) {
	if(!is.character(x)) return(x)
	return(sapply(lapply(strsplit(x, ":"), sort), paste, collapse=":"))
}
