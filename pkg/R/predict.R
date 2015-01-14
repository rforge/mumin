##XXX
#require(MuMIn)
#predict.merMod <- MuMIn:::predict.merMod
#predict.lme <- MuMIn:::predict.lme
##XXX

std_predict <-
function (object, ...) UseMethod("std_predict")

invtransform <-
function (object, y) UseMethod("invtransform")

terms.averaging <-
function (x, ...) {
	terms(formula(x))
}

## Caution! {linkinv = "log"} will transform the prediction on link scale using
## exponential function (through make.link("log"), whereas {linkinv = log} will
## log-transform them.

## use.lincomb currently does not compute predictions properly for models using
## fitting weights and offset.

## component model must have 'model.matrix' method

avgpred <-
function(x, newdata, type = c("response", "invlink", "link", "terms"),
		 se.fit = FALSE, full =  TRUE, linkinv = NULL,
		 use.lincomb = FALSE, ...) {
	type <- match.arg(type)
	noModelList <- is.null(attr(x, "modelList"))
	
	if((type != "terms" && !full) || noModelList || use.lincomb) {
		if(se.fit && !full) stop("cannot calculate \"conditional\" predictions with standard errors")
		if(se.fit) stop("cannot calculate averaged prediction std. errors without model list")
		if(type %in% c("response", "terms")) {
			if(!full) cry(NA, "cannot calculate \"conditional\" predictions of type \"%s\"", type)
			cry(NA, "cannot calculate predictions of type \"%s\" without model list", type)
		}
		message("Computing prediction using averaged coefficients.")
		
		y <- predict_lincomb(x, betas = coef(x, full = full), X = model.matrix(x), trms = terms(x))
		if(type == "invlink") {
			if(is.null(linkinv)) stop("argument 'linkinv' is missing")
			if(inherits(linkinv, c("family", "link-glm",  "character")))
				linkinv <- glm.link(linkinv)$linkinv
			y <- linkinv(y)
		}
	} else {
		stopifnot(!is.null(attr(x, "modelList")))
		models <- attr(x, "modelList")
		xtype <- if(type == "invlink") "link" else type
		
		yall <- vector("list", length(models))
		for(i in seq_along(models))
			yall[[i]] <- std_predict(models[[i]], newdata = newdata,
				type = xtype, se.fit = se.fit, ...)
		
		y1 <- yall[[1L]]
		mod1 <- models[[1L]]
		
		if(!all(vapply(models, inherits, FALSE, class(mod1))))
			stop("cannot average predictions from different model classes (yet)")
			
		if(type == "invlink") {
			if(!is.null(linkinv)) warning("argument 'linkinv' ignored")
			
			links <- tryCatch(vapply(models, function(m) glm.link(m)$name, ""),
				error = function(e) NULL)
			if (is.null(links)) stop("link functions of at least some models not recognized")
			if(any(links[1L] != links[-1L]))
				stop("cannot inverse-transform averaged prediction of models with different link functions")
		}
		dfs <- numeric(0L)
		fn <- if(is.matrix(y1)) {
			#if(all(attr(terms(mod1),"term.labels") == colnames(y1)))
			if(type == "terms") {
				dfs <- tryCatch(vapply(models, df.residual, 1), error = function(e) NULL)
				avgterms
			} else avgmat
		} else if (is.list(y1) && all(c("fit", "se.fit") %in% names(y1))) {
			if(is.matrix(y1$fit)) {
			#if(type == "terms" && is.matrix(y1$fit)) {
				dfs <- tryCatch(vapply(models, df.residual, 1), error = function(e) NULL)
				avgterms
			} else avgsefit
		} else if(is.numeric(y1)) {
			avgvec
		} else stop("'predict' on component model returned unsupported object")
		y <- fn(yall, Weights(x), revised.var = attr(x, "revised.var"),
				full = full, dfs = dfs)
		if(type == "invlink") y <- invtransform(mod1, y)
	}
	y
}

predict_lincomb <-
function (object, newdata = NULL, betas = coef(object),
		  X = model.matrix(object), trms = terms(object),
		  xlev = .getXlevels(trms, model.frame(trms, data = newdata)),
		  contrasts.arg = attr(X, "contrasts")) {
    if (!missing(newdata) && !is.null(newdata))
        X <- model.matrix(trms, data = newdata, contrasts.arg = contrasts.arg, xlev = xlev)
	(X %*% betas)[, 1L]
	
#		if (se.fit) {
#		covmx <- solve(t(X) %*% X)
#		se <- sqrt(diag((Xnew %*% covmx) %*% t(Xnew))) * sqrt(scale) ## TODO: use matmult
#		return(list(fit = y, se.fit = se))
#	} else y
}

avgvec <-
function(yall, w, ...) {
	ret <- yall[[1L]]
	ydim <- c(length(ret), length(yall))
	yarr <- array(unlist(yall), dim = ydim)
	for(i in seq_len(ydim[1L]))
		ret[i] <- weighted.mean(yarr[i, ], w = w)
	ret
}

## helper function
toarray <- function(x)
	array(unlist(x), dim = c(dim(x[[1L]]), length(x)),
	 dimnames = c(dimnames(x[[1L]]),  list(names(x))))


## averages prediction of each model term (predict.[g]lm):
avgterms <-
function(yall, w, revised.var, full, dfs, ...) {
	fitonly <- vapply(yall, is.atomic, FALSE)
	stopifnot(!any(fitonly[1L] != fitonly[-1L])) ## inconsistent prediction form
	fitonly <- fitonly[1L]

	if(fitonly) {
		all.coefnames <- unique(fixCoefNames(unlist(lapply(yall, colnames))))
		fits <- toarray(lapply(yall, function(x) x[, match(all.coefnames, colnames(x)), drop = FALSE]))
		if(full) fits[is.na(fits)] <- 0
		sefits <- array(0, dim = dim(fits))
	} else {
		all.coefnames <- unique(fixCoefNames(unlist(lapply(yall, function(x) colnames(x$fit)))))
		yall <- lapply(yall, function(x) {	
			j <- match(all.coefnames, colnames(x$fit))
			x$fit <- x$fit[, j]
			x$se.fit <- x$se.fit[, j]
			x
		})
		fits <- toarray(lapply(yall, "[[", "fit"))
		sefits <- toarray(lapply(yall, "[[", "se.fit"))
		
		if(full) {
			sefits[is.na(sefits)] <- 0
			fits[is.na(fits)] <- 0
		}
	}
	
	hasDfs <- !is.null(dfs) && !any(is.na(dfs))
	k <- if(hasDfs) 3L else 2L
	d <- dim(fits)
	avgfit <- avgsefit <- array(NA_real_, dim = d[1L:2L], dimnames = dimnames(fits)[1L:2L])
	for (i in 1L:d[1L]) for (j in 1L:d[2L]) {
		res <- par.avg(fits[i, j, ], sefits[i, j, ], weight = w, df = dfs) 
		avgfit[i, j] <- res[1L]
		avgsefit[i, j] <- res[2L]
	}
	colnames(avgfit) <- colnames(avgsefit) <- all.coefnames
	if(fitonly) avgfit else list(fit = avgfit, se.fit = avgsefit)
}

avgsefit <-
function(yall, w, revised.var, full, dfs,  ...) {
	fit <- do.call("cbind", lapply(yall, "[[", "fit"))
	se.fit <- do.call("cbind", lapply(yall, "[[", "se.fit"))
	n <- nrow(fit)
	y <- matrix(0, ncol = 2L, nrow = n)
	hasDfs <- !is.null(dfs) && !any(is.na(dfs))
	k <- c(1L, if(hasDfs) 3L else 2L) 
	for(i in 1L:n) y[i, ] <-
		par.avg(fit[i, ], se.fit[i, ], weight = w,
			df = dfs, revised.var = revised.var)[k]
	list(fit = y[, 1L], se.fit = y[, 2L])
}

avgmat <-
function(yall, w, ...) {
	ret <- yall[[1L]]
	ydim <- c(dim(ret), length(yall))
	yarr <- array(unlist(yall), dim = ydim)
	kseq <- seq_len(ydim[2L])
	for(i in seq_len(ydim[1L]))
		for(k in kseq) ret[i, k] <- weighted.mean(yarr[i, k, ], w = w)
	ret
}

std_predict.default <-
function(object, newdata, ...) {
	cl <- match.call()
	if(missing(newdata) || is.null(newdata)) cl$newdata <- NULL
	cl[[1L]] <- as.name("predict")
	eval.parent(cl)
}

std_predict.multinom <-
std_predict.polr <-
function(object, newdata, type = c("link", "response"), ...) {
	type <- match.arg(type)
	cl <- match.call()
	cl[[1L]] <- as.name("predict")
	cl[['type']] <- "probs"
	rval <- eval.parent(cl)
	if(type == "link") glm.link(object)$linkfun(rval) else rval
}
		
		
std_predict.merMod <-
function(object, newdata = NULL, type = c("link", "response"),
		 se.fit = FALSE, re.form = NA, ...) {
	if(missing(newdata)) newdata <- NULL
	predict.merMod(object, newdata, type, se.fit, re.form = re.form, ...)
}

std_predict.lme <-
function(object, newdata = NULL, type = NA, se.fit = FALSE, ...) {
	predict.lme(object, newdata, se.fit = se.fit, ...)
}

std_predict.lm <-
function(object, newdata = NULL, type = c("link", "response", "terms"), se.fit = FALSE, ...) {
	type <- match.arg(type)
	if(!inherits(object, "glm") && type == "link") type <- "response"
	std_predict.default(object, newdata = newdata, type = type, se.fit = se.fit, ...)
}

invtransform.gls <-
invtransform.lme <-
invtransform.lm <-
function(object, y) y

invtransform.glm <-
invtransform.default <-
function(object, y) {
	link <- glm.link(object)
	if(is.list(y) && all(c("fit", "se.fit") %in% names(y))) {
		y$se.fit <- y$se.fit * abs(link$mu.eta(y$fit))
		y$fit <- link$linkinv(y$fit)
		y
	} else {
		link$linkinv(y)
	}
}

std_predict.zeroinfl <-
function(object, newdata, type = c("link", "response", "prob"), ...) {
	type <- match.arg(type)
	cl <- match.call()
	cl[[1L]] <- as.name("predict")
	if(is.null(cl$newdata)) cl$newdata <- NULL
	if(type == "link") {
		cl$type <- "zero"
		ret <- family(object)$linkfun(eval.parent(cl))
		cl$type <- "count"
		ret <- cbind(ret, log(eval.parent(cl)))
		colnames(ret) <- c("linkphi", "logmu")
		ret
	} else {
		eval.parent(cl)
	}
}

std_predict.hurdle <-
function(object, newdata = NULL, type = c("link", "response", "prob"), ...) {
	type <- match.arg(type)
	cl <- match.call()
	cl[[1L]] <- as.name("predict")
	if(is.null(cl$newdata)) cl$newdata <- NULL
	if(type == "link") {
		cl$type <- "zero"
		ret <- log(eval.parent(cl))
		cl$type <- "count"
		ret <- cbind(ret, log(eval.parent(cl)))
		colnames(ret) <- c("logphi", "logmu")
		ret
	} else {
		eval.parent(cl)
	}
}

std_predict.unmarkedFit <-
function (object, newdata, type = NA, ...) {
	if(!is.na(type) && type != "response")
		stop(gettextf("prediction of type '%s' is not implemented", type))
	cl <- match.call()
	cl[[1L]] <- as.name("predict")
	if(missing(newdata)) cl$newdata <- NULL
	types <- names(object@estimates)
	ntypes <- length(types)
	y <- vector("list", ntypes)
	names(y) <- types <- types
	for(i in seq_len(ntypes)) {
		cl$type <- types[i]
		y[[i]] <- eval.parent(cl)
	}
	rval <- vector("list", 2L)
	for(i in 1L:2L) rval[[i]] <- do.call("cbind", lapply(y, "[", , i))
	names(rval) <- c("fit", "se.fit")
	rval
}

invtransform.zeroinfl <-
function(object, y) {
	(1 - object$linkinv(y[, "linkphi"])) * exp(y[, "logmu"])
}

invtransform.hurdle <-
function(object, y) {
	exp(y[, "logphi"] + y[, "logmu"])
}
