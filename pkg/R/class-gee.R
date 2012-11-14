

`coefTable.gee` <-
`coefTable.geeglm` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$coefficients
	type <- match.arg(type)
	j <- if(type == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, 2L], coefNames = rownames(cf))
}

`coefTable.geese` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$mean
	type <- match.arg(type)
	j <- if(type == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, 2L], coefNames = rownames(cf))
}

`coef.geese` <-
function(object, ...)  object$beta


.qlik <- function(y, mu, fam) {
	ret <- switch(fam,
		   gaussian = -sum((y - mu)^2)/2,
		   binomial = sum(y * log(mu/(1 - mu)) + log(1 - mu)),
		   #binomial.sqvar = sum(((2 * y - 1) * log(mu /(1 - mu))) - (y / mu) - ((1 - y)/(1 - mu))),
		   poisson = sum(y * log(mu) - mu),
		   Gamma = -sum(y/mu + log(mu)),
		   inverse.gaussian = sum(-y/(2 * mu^2) + 1/mu),
		   stop("do not know how to calculate quasi-likelihood for family ",
				dQuote(fam))
		   )
	ret
}

quasiLik <- function (object, ...) UseMethod("quasiLik")

`quasiLik.geeglm` <-
`quasiLik.gee` <-
function(object, ...) {
	ret <- .qlik(object$y, object$fitted.values, family(object)$family)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(object$y)
	class(ret) <- "quasiLik"
	ret
}


print.quasiLik <- function (x, digits = getOption("digits"), ...) {
    cat("'quasi Lik.' ", paste(format(c(x), digits = digits), collapse = ", "), 
        "\n", sep = "")
    invisible(x)
}

`QIC` <- function(object, ...) UseMethod("QIC", object)

.qic <- function(mu, vbeta, i.vbeta.naiv, qlik) {
	AIinv <- solve(i.vbeta.naiv) # solve via indenity
	tr <- sum(diag(AIinv %*% vbeta))
	px <- length(mu) # number non-redunant columns in design matrix
	# QIC
	ret <- -2 * qlik + 2 * tr
	QICu <- -2 * qlik + 2 * px    # Approximation assuming model structured correctly
	attr(ret, "QICu") <- QICu
	ret
}

getQIC <- function(x) UseMethod("getQIC")
	
getQIC.gee <- function(x) {
	capture.output(suppressMessages(xi <- update(x, corstr = "independence",
		silent = TRUE)))
	mu <- x$fitted.values 
	c(.qic(mu, x$robust.variance, xi$naive.variance, .qlik(x$y, mu,
		family(x)$family)), length(x$y))
}

getQIC.geeglm <- function(x) {
	xi <- update(x, corstr = "independence")
	mu <- x$fitted.values 
	c(.qic(mu, x$geese$vbeta, xi$geese$vbeta.naiv,
		.qlik(x$y, mu, family(x)$family)), length(x$y))
}

getQIC.default <- function(x) evalq(.NotYetImplemented(), parent.frame())


getCall.yagsResult <- function(x, ...) x@Call

`coefTable.yagsResult` <-
function(model, ..., type = c("naive", "robust")) {
	type <- match.arg(type)
	vcv <- slot(model, if(type == "naive") "naive.parmvar" else "robust.parmvar")
	.makeCoefTable(model@coefficients, vcv, coefNames = model@varnames)
}




family.default <- function (object, ...) {
    cl <- getCall(object)
    if (is.null(cl)) 
        return(NULL)
    fam <- cl$family
    if (is.null(fam)) 
        fam <- formals(match.fun(cl[[1L]]))$family
    if (is.null(fam)) {
        return(NA)
    }
    switch(mode(fam), call = eval(fam), name = , character = match.fun(fam)())
}


#getQIC.yagsResult <- function(x) x@pan.aic
getQIC.yagsResult <- function(x) {
	xi <- update(x, corstruct = "independence")
	##
	#cl <- match.call(call = getCall(yags1), yags::yags)
	#cl[[1L]] <- as.name("model.frame.default")
	#cl$formula[[3L]] <- 1L
	#cl <- cl[c(TRUE, (names(cl)[-1L] %in% c("formula", "data", "subset")))]
	#y <- eval(cl, parent.frame())[, 1L]
	
	mu <- x@fitted.values
	y <- mu + x@residuals
	c(.qic(mu, x@robust.parmvar, xi@naive.parmvar,
		.qlik(y, mu, family(x)$family)), length(y))
}


`QIC` <- function (object, ...) {
	if (length(list(...))) {
		res <- sapply(list(object, ...), getQIC)
		val <- data.frame(QIC = res[1L, ])
		Call <- match.call()
		row.names(val) <- as.character(Call[-1L])
		val
	} else getQIC(object)[1L]
}





