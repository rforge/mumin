`deviance.mark` <- function(object, ...) object$results[['deviance']]
`nobs.mark` <- function (object, ...) object$results[['n']]

`coeffs.mark` <- function(model) {
	cf <- model$results$beta[, 1L]
	names(cf) <- gsub("^([a-zA-Z]+):(.*)$", "\\1(\\2)",
		rownames(model$results$beta), perl = TRUE)
	cf
}

`coefTable.mark` <- function (model, orig.names = FALSE, ...) {
    dfs <- model$results[['n']] - model$results[['npar']]
    beta <- model$results[['beta']]
    MuMIn:::.makeCoefTable(beta[, 1L], beta[, 2L], dfs,
		coefNames = if(orig.names) rownames(beta) else
			gsub("^([a-zA-Z]+):(.*)$", "\\1(\\2)", rownames(beta), perl = TRUE))
}


`logLik.mark` <- function (object, adjust = TRUE, ...) {
	res <- -0.5 * object$results$lnl
	attr(res, "df") <- object$results[[if(!adjust && !is.null(object$results$npar.unadjusted))
		'npar.unadjusted' else 'npar']]
	attr(res, "nobs") <- object$results$n 
	class(res) <- "logLik"
	res
}

`confint.mark` <- function (object, parm, level = 0.95, ...) {
    cf <- object$results$beta[, 1L]
	nm <- names(cf) <- rownames(object$results$beta)
	df.residual <- object$results$n - object$results$npar
	vcv <- object$results$beta.vcv
	dimnames(vcv) <- list(nm, nm)
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    fac <- qt(a, df.residual)
    pct <- stats:::format.perc(a, 3L)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
        pct))
    ses <- sqrt(diag(vcv))[parm]
    ci[] <- cf[parm] + ses %o% fac
    ci
}

`formula.mark` <- function (x, expand = TRUE, ...) {
	param <- if(is.null(x$model.parameters)) x$parameters else  x$model.parameters
	f <- lapply(param, "[[", 'formula')
	f <- f[!vapply(f, is.null, logical(1L))]
	
	npty <- length(f)
	z <- vector(npty, mode = "list")
	pty <- names(f)
	
	if(expand) {
		for(i in seq_len(npty)) z[[i]] <- paste(pty[i], "(",
				getAllTerms(f[[i]], intercept = TRUE), ")", sep = "")
		res <- reformulate(gsub("((Intercept))", "(1)", unlist(z), fixed = TRUE))
	} else {
		for(i in seq_len(npty)) z[[i]] <- call(pty[i], f[[i]][[2L]])
		res <- z[[1L]]
		if(npty > 1L) for(i in seq(2L, npty)) res <- call("+", res, z[[i]])
		res <- eval(call("~", res))		
	}
	environment(res) <- environment(f[[1L]])
	res
}


`getAllTerms.mark` <- function (x, intercept = FALSE, ...) {
	
	f <- formula(x, expand = FALSE)[[2L]]
	ret <- list()
	while(length(f) == 3L && f[[1L]] == "+") {
		ret <- append(f[[3L]], ret)
		f <- f[[2L]]
	}
	ret <- append(f, ret)
	res <- lapply(ret, function(x) {
		func <- deparse(x[[1L]], control = NULL)
		tt <- terms(eval(call("~", x[[2L]])))
		tlab <- attr(tt, "term.labels")
		torder <- attr(tt, "order")
		if(attr(tt, "intercept")) {
			tlab <- append("(Intercept)", tlab)
			torder <- c(0L, torder)
		}
		res1 <- lapply(fixCoefNames(tlab), function(z) paste(func, "(", z, ")", sep = ""))
		attr(res1, "order") <- torder
		res1
	})
	
	ord <- order(rep(seq_along(res), sapply(res, length)),
		unlist(lapply(res, attr, "order")))
	res <- unlist(res, recursive = TRUE)[ord]
	ints <- grep("((Intercept))", res, fixed = TRUE)
	attr(res, "intercept") <- as.numeric(ints != 0L)
	attr(res, "interceptLabel") <- res[ints]
	if(!intercept) {
		res <- do.call("structure", c(list(res[-ints]), attributes(res)))
		attr(res, "order") <- order(ord[-ints])
	} else {
		attr(res, "order") <- order(ord)
	}
	
	res
}

`makeArgs.mark` <- function(obj, termNames, comb, opt, ...) {
	
	interceptLabel <- "(Intercept)"
	termNames <- sub(interceptLabel, "1", termNames, fixed = TRUE)

	rxres <- regexpr("^([a-zA-Z]+)\\((.*)\\)$", termNames, perl = TRUE)
	cs <- attr(rxres, "capture.start")
	cl <- attr(rxres, "capture.length")
	parname <- substring(termNames, cs[, 1L], cs[, 1L] + cl[,1L] - 1L)
	parval <- substring(termNames, cs[, 2L], cs[, 2L] + cl[,2L] - 1L)
	
	formulaList <- lapply(split(parval, parname), function(x) {
		int <- x == "1"
		x <- x[!int]
		res <- if(!length(x))
				if(int) ~ 1 else ~ 0 else 
			reformulate(x, intercept = any(int))
		environment(res) <- opt$gmFormulaEnv
		res
	})
	
	
	
	mpar <- if(is.null(obj$model.parameters))
		eval(opt$gmCall$model$parameters) else
		obj$model.parameters
	for(i in names(mpar)) mpar[[i]]$formula <- formulaList[[i]]
	#ret <- list(model.parameters = mpar)
	
	if(opt$gmCall[[1L]] == "run.mark.model") {
		arg.model <- opt$gmCall$model
		arg.model$parameters <- mpar
		ret <- list(model = arg.model)
	} else {
		ret <- list(model.parameters = mpar)
	}
	
	attr(ret, "formulaList") <- formulaList
	ret
}
