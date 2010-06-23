`model.avg` <-
function(m1, ..., beta = FALSE, method = c("0", "NA"), rank = NULL, rank.args = NULL,
		alpha = 0.05
) {
	method <- match.arg(method)

	if (!is.null(rank)) {
	   	rankFn <- match.fun(rank)
		#browser()
		rank.call <- as.call(c(as.name(substitute(rank)), NA, rank.args))
		rank <- substitute(rank)

	} else if (!is.null(attr(m1, "rank.call"))) {
		rank.call <- attr(m1, "rank.call")
		rank.args <- as.list(attr(m1, "rank.call"))[-(1:2)]
		rankFn <- match.fun(rank.call[[1]])
		rank <- as.character(rank.call[[1]])
	}

	if (inherits(m1, "list")) {
		models <- m1
		m1 <- models[[1]]
	} else {
		models <- list(m1, ...)
	}

	if (!is.null(rank)) {
		# test the rank function
  		rank.call[[2]] <- quote(m1)
  		x <- eval(rank.call)
  		if (!is.numeric(x) || length(x) != 1) {
			stop(sQuote("rank"), " should return numeric vector of length 1")
		}

	}

	if (length(models) == 1) {
		stop("Only one model supplied. Nothing to do")
	}
	if (is.null(rank)) {
		aicc <- sapply (models, AICc)
	} else {
		cl <- as.call(c(as.name("sapply"), quote(models), quote(rankFn), rank.args))
		aicc <- eval(cl)
		#aicc <- do.call("sapply", list(models, rankFn, rank.args), quote = T)

	}

	if (!is.null(tryCatch(deviance(models[[1]]), error = function(...)
	NULL)))
		dev <- sapply (models, deviance)
	else
		dev <- NULL

	delta <- aicc - min(aicc)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))

	model.order <- order(weight, decreasing=TRUE)
	aicc <- aicc[model.order]
	delta <- delta[model.order]
	weight <- weight[model.order]
	models <- models[model.order]
	dev <- dev[model.order]

	selection.table <- data.frame(AICc = aicc, Delta = delta, Weight = weight)

	if (!is.null(dev)) {
          selection.table <- cbind(Deviance = dev, selection.table)
	}

	all.par <- unique(unlist(lapply(models, function(m) names(coeffs(m)))))

	all.terms <- unique(unlist(lapply(models, getAllTerms)))
	all.terms <- all.terms[order(sapply(gregexpr(":", all.terms),
										function(x) if(x[1] == -1) 0 else length(x)), all.terms
								 )]

	all.par <- all.par[order(sapply(gregexpr(":", all.par),
									function(x) if(x[1] == -1) 0 else length(x)), all.par
							 )]

	all.coef <- all.var <- all.df <- numeric(0)

	ac <- rep(0, length = length(all.par))

 	for (m in models) {
		m.tTable <- tTable(m)
		n <- NROW(residuals(m))
		m.coef <- m.tTable[,1]
		m.var <- m.tTable[,2]
 		m.df <- n - length(m.coef)

		if (beta) {
			response.sd <- sd(model.frame(m1)[, attr(terms(m1) ,"response")])
			m.vars.sd <- sd(model.matrix(m))
			bx <- m.vars.sd / response.sd
			m.coef <- m.coef * bx
			m.var <- m.var * bx
		}

		m.vars <- match(all.par, rownames(m.tTable))

		all.coef <- rbind(all.coef, model = c(m.coef[m.vars]))
		all.var <- rbind(all.var, model = c(m.var[m.vars]))
		all.df <- append(all.df, m.df)

	}

	all.model.names <- sapply(models, function(x) paste(match(getAllTerms(x), all.terms), collapse="+"))

	importance <- apply(weight * t(sapply(models, function(x) all.terms %in% getAllTerms(x))), 2, sum)

	names(importance) <- all.terms

	# check if models are unique:
	dup <- duplicated(all.model.names)
	if (any(dup)) {
  		dup <- table(all.model.names)
		dup <- seq(all.model.names)[all.model.names %in% names(dup[dup > 1])]
		stop("Models are not unique. Duplicates: ", paste(dup, collapse=", "))
	}
	##

	rownames(all.var) <- rownames(all.coef) <- rownames(selection.table) <- all.model.names

	if (method == "0") {
		all.coef[is.na(all.coef)] <- 0
		all.var[is.na(all.var)] <- 0
	}

	avg.model <- t(sapply(seq_along(all.par),
		function(i) par.avg(all.coef[,i], all.var[,i], all.df, weight, alpha)))

	all.coef[all.coef == 0] <- NA
	all.var[all.var == 0] <- NA

	importance <- sort(importance, decreasing=T)
	colnames(all.coef) <- colnames(all.var) <- rownames(avg.model) <-  all.par

    names(all.terms) <- seq_along(all.terms)

	if (!is.null(rank)) {
		colnames(selection.table)[2] <- as.character(rank)
	}


	mmxs <- lapply(models, model.matrix)
	mx <- mmxs[[1]];
	for (i in mmxs[-1])
		mx <- cbind(mx, i[,!(colnames(i) %in% colnames(mx)), drop=FALSE])

	# residuals averaged with brute force
	residuals <- apply(sapply(models, residuals), 1, weighted.mean, w=weight)

	trm <- terms(models[[1]])
	frm <- reformulate(all.terms,
				response = attr(trm, "variables")[-1][[attr(trm, "response")]])

	ret <- list(
		summary = selection.table,
		coefficients = all.coef,
		variable.codes = all.terms,
		variance = all.var,
		avg.model = avg.model,
		relative.importance = importance,
		weights = weight,
		beta = beta,
		terms = all.par,
		model.matrix = mx,
		residuals = residuals,
		formula = frm
	)

	attr(ret, "mList") <- models

	class(ret) <- "averaging"
	return(ret)
}

`coef.averaging` <-
function(object, ...) object$avg.model[,1]


`predict.averaging` <-
function(object, newdata = NULL, se.fit = NULL, interval = NULL, type = NULL, ...) {

	#if(("type" %in% names(match.call())) && type != "link") {
	if(!missing("type") && type != "link") {
		warning("Only predictions on the link scale are allowed. Argument ",
				dQuote("type"), " ignored")
	}
	if (!missing(se.fit)) .NotYetUsed("se.fit", error = FALSE)
	if (!missing(interval)) .NotYetUsed("interval", error = FALSE)

	models <- attr(object, "mList")

	# If all models inherit from lm:
	if (all(sapply(models, inherits, what="lm"))) {
		coeff <- coef(object)
		frm <- formula(object)

		tt <- delete.response(terms(frm))
		X <- object$model.matrix

		if (missing(newdata) || is.null(newdata)) {
			Xnew <- X
		} else {
			Xnew <- model.matrix(tt, data = newdata)
		}
		Xnew <- Xnew[, match( names(coeff),colnames(Xnew), nomatch = 0)]
		ny <- (Xnew %*% coeff)[, 1]

		#if (se.fit) {
		#	covmx <- solve(t(X) %*% X)
		#	se <- sqrt(diag(Xnew %*% covmx %*% t(Xnew))) * sqrt(scale)
		#	return(list(fit = y, se.fit = se))
		#}
	} else {
		# otherwise, use brute force:
		ny <- sapply(models, predict, newdata = newdata, ...)
		ny <- apply(ny, 1, weighted.mean, w = object$weight)
	}

	return(ny)
}


`summary.averaging` <-
function (object, ...) print.averaging(object)


`print.averaging` <-
function(x, ...) {
	cat("\nModel summary:\n")
	print(round(x$summary,1))

	cat("\nVariables:\n")
	print(x$variable.codes, quote= F)

	cat("\nAveraged model parameters:\n")
	print(signif(x$avg.model, 3))

	cat("\nRelative variable importance:\n")
	print(round(x$relative.importance, 2))
}
