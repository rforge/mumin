`DIC` <-
function (object, ...) {
	if (length(list(...))) {
		lls <- sapply(list(object, ...), function(x) {
			c(extractDIC(x), attr(logLik(x), "df"))
		})
		val <- data.frame(df = lls[2L, ], DIC = lls[1L, ])
		Call <- match.call()
		row.names(val) <- make.unique(as.character(Call[-1L]))
		val
	} else extractDIC(object)
}

if(!exists("extractDIC", mode = "function")) {
	extractDIC <- function (fit, ...) UseMethod("extractDIC")
}

# from package 'arm'
`extractDIC.mer` <- function (fit, ...) {
	dev <- deviance(fit, REML = fit@dims["REML"])
    devML <- deviance(fit, REML = FALSE)
    as.vector(2 * devML - dev)
}

`extractDIC.merMod` <- function (fit, ...) {
	dev <- deviance(fit, REML = isREML(fit))
    devML <- deviance(fit, REML = FALSE)
    as.vector(2 * devML - dev)
}

`extractDIC.MCMCglmm` <- function (fit, ...) fit$DIC

`extractDIC.lme` <- function (fit, ...) {
	ll <- as.vector(logLik(fit, REML = fit$method == "REML"))
    llML <- as.vector(logLik(fit, REML = FALSE))
    2 * ll - 4 * llML
}


`formula.MCMCglmm` <-
function (x, ...) x$Fixed$formula

`nobs.MCMCglmm` <-
function (object, ...) object$Residual$nrl

`family.MCMCglmm` <-
function (object, ...) object$family

`logLik.MCMCglmm` <-
function (object, ...)
	structure(NA, df = object$Fixed$nfl + object$Random$nfl,
			  nobs = object$Residual$nrl, class = "logLik")

`coeffs.MCMCglmm` <-
function (model) summary(model)$solutions[, 1L]

`coefTable.MCMCglmm` <-
function (model, ...) {
	cf <- coeffs(model)
	.makeCoefTable(cf, se = rep(NA_real_, length.out = length(cf)))
}

`getAllTerms.MCMCglmm` <- function (x, ...) {
	res <- MuMIn:::getAllTerms.default(x, ...)
	attr(res, "random") <- .formulaEnv(.~., environment(formula(x)))
	attr(res, "random.terms") <- deparse(x$Random$formula, control = NULL)[1]
	res
}

