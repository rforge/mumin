`AICc` <-
function(object, ..., k = 2) {
	if(length(list(...))) {
		object <- list(object, ...)
		val <- as.data.frame(t(sapply(object,
 							function(el) {
					     		z <- getAICc(el, k = k)
					     		c(attr(z, "df"), attr(z, "AIC"), z)
						  	}
						  )))
		names(val) <- c("df", "AIC", "AICc")
		   Call <- match.call()
		   Call$k <- NULL
		row.names(val) <- as.character(Call[-1])
		return(val)
	} else {
		return(getAICc(object, k = k))
	}
}

`getAICc` <-
function(model, k = 2) {
	if (any(inherits(model,  c("lmer", "glmer")))) {
		mLogLik <- logLik(model, model@status["REML"])
		N <- NROW(model@frame)
	} else {
		mLogLik <- logLik(model)
		N <- length(resid(model))
	}

	mK <- attr(mLogLik, "df")
	mAIC <- -2 * c(mLogLik) + k * mK
	ret <- mAIC + 2 * mK * (mK + 1)/(N - mK - 1)
	attr(ret, "df") <- mK
	attr(ret, "AIC") <- mAIC
	return (ret)
}


`QAIC` <- function(object, ..., chat) {
	#chat <- summary(gm)$dispersion
	if(length(list(...))) {
		object <- list(object, ...)
		val <- data.frame(QAIC=sapply(object, getQAIC, chat = chat))
		   Call <- match.call()
		   Call$chat <- NULL
		row.names(val) <- as.character(Call[-1])
		return(val)
	} else {
		return(getQAIC(object, chat = chat))
	}
}


`getQAIC` <-
function(model, chat) {
	if (any(inherits(model,  c("lmer", "glmer")))) {
		mLogLik <- logLik(model, model@status["REML"])
		N <- NROW(model@frame)
	} else {
		mLogLik <- logLik(model)
		N <- length(resid(model))
	}

	k <- attr(mLogLik, "df") + 1
 	ret <- (deviance(model) / chat) + 2 * k
	return (ret)
}


`beta.weights` <-
function(model) {

	summ <- summary(model)
	m.coef <- summ$coefficients[,1]
	m.se <- summ$coefficients[,2]

	response.sd <- sd(eval(attr(model$terms ,"variables"), envir=model$model)[[attr(model$terms ,"response")]])
	m.terms.sd <- sd(model.matrix(model))
	bx <- m.terms.sd / response.sd

	m.b.coef <- m.coef * bx
	m.b.se <- m.se * bx

	ret <- data.frame(m.coef, m.se, m.b.coef, m.b.se)

	colnames(ret) <- c("Estimate", "Std. Err.", "Beta", "Std. Err. Beta")
	rownames(ret) <- names(model$coefficients)

	ret <- as.matrix(ret)
	return (ret)
}

`coeffs` <-
function (model) UseMethod("coeffs")

`coeffs.lme` <-
function(model) { model$coefficients$fixed}


`coeffs.mer` <-
 function(model) { return(model@fixef)}


`coeffs.glmer` <-
`coeffs.lmer` <-
function(model) { ret <- model@fixef; names(ret) <- model@cnames$.fixed; ret}

`coeffs.spautolm` <-
function(model) { model$fit$coefficients}

`coeffs.gls` <-
function (model) return(summary(model)$coefficients)

`coeffs.default` <-
function(model) { model$coefficients}



`dredge` <-
function(global.model, beta = FALSE, eval = TRUE, rank = "AICc", ...) {

	rankFn <- match.fun(rank)
	if (is.function(rank)) {
  		rank <- deparse(substitute(rank))
	}

	if (rank != "AICc") {
		arg <- list(...)
		rankFnCall <- as.call(c(as.name("rankFn"), substitute(global.model), arg))

		# test the rank function
  		x <- eval(rankFnCall)
  		if (!is.numeric(x) || length(x) != 1) {
			stop(sQuote("rank"), " should return numeric vector of length 1")
		}

	}
	intercept <- "(Intercept)"

	all.terms <- getAllTerms(global.model)
	has.int <- attr(all.terms, "intercept")

	n.vars <- length(all.terms)
	ms.tbl <- numeric(0)
	formulas <- character(0)

	is.glm <- inherits(global.model, "glm")
	is.lm <- !is.glm & inherits(global.model, "lm")

	if (
			(inherits(global.model, c("mer")) && ("REML" %in% names(deviance(global.model))))
		|| 	(inherits(global.model, c("lme")) && global.model$method == "REML")
		||   (any(inherits(global.model, c("lmer", "glmer"))) && global.model@status["REML"] != 0)
	) {
			warning("Comparing models with different fixed effects fitted by REML")
	}

	if (!is.lm && beta) {
		warning("Cannot calculate beta weigths (yet) for ", class(global.model)[1])
          beta <- FALSE
	}

	has.rsq <- "r.squared" %in% names(summary(global.model))
	has.dev <- !is.null(deviance(global.model))

	all.comb <- expand.grid(split(c(1:n.vars, rep(0, n.vars)), all.terms))

	formulas <- apply(all.comb, 1, function(.x) as.formula(paste(". ~", paste(c(1, all.terms[.x]), sep=" ", collapse=" + "))))

	ss <- sapply(formulas, formulaAllowed)

	all.comb <- as.matrix(all.comb[ss, ])
	formulas <- formulas[ss]
	names(formulas) <- seq(formulas)

	if (any(inherits(global.model,  c("mer", "lmer", "glmer")))) {
          formulas <- lapply(formulas, update, attr(all.terms, "random"))
	}

	if (!eval) {
		return(formulas)
	}


	###
	for(b in seq(NROW(all.comb))) {
          cterms <- all.terms[unlist(all.comb[b,])]

		frm <- formulas[[b]]

		c.row <- rep(NA, n.vars)
		c.row[match(cterms, all.terms)] <- rep(1, length(cterms))


		cl <- call("update", substitute(global.model), frm)

		cmod <- try(eval(cl, parent.frame()))
		if (inherits(cmod, "try-error")) {
			print(formulas[[as.character(b)]])
			formulas[[as.character(b)]] <- NA

			next;
		}
		mod.coef <- c(na.omit(match(all.terms, names(coeffs(cmod)))))

          icept <- if (attr(all.terms, "intercept")) coeffs(cmod)["(Intercept)"] else NA

	     cmod.all.coef <- if (beta) beta.weights(cmod)[,3] else coeffs(cmod)

		mod.coef <- cmod.all.coef[mod.coef]
		mod.coef.names <- names(mod.coef)
		mod.coef.index <- match(mod.coef.names, all.terms)
		c.row[match(mod.coef.names, all.terms)] <- mod.coef

		aicc <- AICc(cmod)
		aic <- attr(aicc, "AIC")

		c.row <- c(icept, c.row, k=attr(aicc,"df"))
		if (has.rsq) {
			cmod.summary <- summary(cmod)
			c.row <- c(c.row, r.squared=cmod.summary$r.squared, adj.r.squared=cmod.summary$adj.r.squared)
		}
		if (has.dev)
			c.row <- c(c.row, deviance(cmod))


		# TODO:
		if (rank != "AICc") {
			rankFnCall[[2]] <- cmod
			ic <- eval(rankFnCall)
			c.row <- c(c.row, IC=ic)

			#c.row <- c(c.row, IC=rankFn(cmod, ...))
		} else {
		     c.row <- c(c.row, AIC=aic, AICc=aicc)
		}


		ms.tbl <- rbind(ms.tbl, c.row)

	}

	formulas[is.na(formulas)] <- NULL
	ms.tbl <- data.frame(ms.tbl, row.names=1:NROW(ms.tbl))


	cnames <- c("(int.)", all.terms, "k")
	if (has.rsq) {
		cnames <- append(cnames, c("R.sq", "Adj.R.sq"))
	}
	if (has.dev) {
		cnames <- append(cnames, ifelse (is.lm, "RSS", "Dev."))
	}

	if (rank == "AICc") {
		cnames <- append(cnames, c("AIC", "AICc"))
	} else {
		cnames <- append(cnames, rank)
	}

	colnames(ms.tbl) <- cnames

	o <- order(ms.tbl[, rank], decreasing = FALSE)

	ms.tbl <- ms.tbl[o,]
	ms.tbl$delta <- ms.tbl[, rank] - min(ms.tbl[, rank])
	ms.tbl$weight <- exp(-ms.tbl$delta / 2) / sum(exp(-ms.tbl$delta / 2))

	class(ms.tbl) = c("model.selection", "data.frame")

	attr(ms.tbl, "formulas") <- formulas[o]
	attr(ms.tbl, "global") <- global.model
	attr(ms.tbl, "terms") <- c("(Intercept)", all.terms)

	if (rank != "AICc") {
		rankFnCall[[1]] <- as.name(rank)
		rankFnCall[[2]] <- substitute(global.model)
		attr(ms.tbl, "rank.call") <- rankFnCall
	}


	if (!is.null(attr(all.terms, "random.terms")))
		attr(ms.tbl, "random.terms") <- attr(all.terms, "random.terms")

	return(ms.tbl)
}

`formulaAllowed` <-
function(frm) {
	factors <- attr(terms(frm), "factors")
	return(all(factors < 2))
}

`get.models` <-
function(dd, subset = delta <= 4, ...) {

	subset <- eval(substitute(subset), envir = dd, enclos = parent.frame())
	gmod <- attr(dd, "global")
	frm <- attr(dd, "formulas")[subset]

	sgmod <- substitute(gmod)
	models <- lapply(frm, function(.x) eval(call("update", sgmod, .x), sys.parent(3)))

	if (!is.null(attr(dd, "rank.call"))) {
  		attr(models, "rank.call") <- attr(dd, "rank.call")
	}

	return(models)
}


`getAllTerms.default` <-
function(x, ...) {
	return(getAllTerms(as.formula(formula(x))))
}

`getAllTerms.formula` <-
function(x, ...) {
	mTerms <- terms(x)
	ret <- attr(terms(x),"term.labels")
	if (length(ret) > 0) {
		ret <- ret[order(ret)]
		i <- grep(" ", ret)
		ret[i] <- paste("(", ret[i] , ")")

		mTerms <- terms(as.formula(paste(". ~", paste(ret, sep=" ", collapse=" + "))))
		ret <- attr(mTerms, "term.labels")
	}

	attr(ret, "intercept") <- attr(mTerms, "intercept")
	ret
}

`getAllTerms.lme` <-
function(x, ...) {
	getAllTerms(as.formula(formula(x)))
}


`getAllTerms.mer` <-
`getAllTerms.glmer` <-
`getAllTerms.lmer` <-
function(x, ...) {
     ret <- getAllTerms(as.formula(formula(x)))
     i <- grep(" \\| ", ret)

     intercept <- attr(ret, "intercept")

	rnd <- ret[i]
     ret <- ret[-i]
     attr(ret, "random.terms") <- rnd

     rnd.formula <- paste("(", rnd, ")", sep="", collapse=" + ")
     rnd.formula <- as.formula(paste(". ~ .", rnd.formula, sep="+"))

	attr(ret, "random") <- rnd.formula
	attr(ret, "intercept") <- intercept

     ret
}

`getAllTerms` <-
function (x, ...) UseMethod("getAllTerms")


Weights <- function(ic, ...) {


	delta <- ic - min(ic)
     weight <- exp(-delta / 2) / sum(exp(-delta / 2))
}

Weights.default <- function(ic) {
	delta <- ic - min(ic)
     weight <- exp(-delta / 2) / sum(exp(-delta / 2))
}

"AICc"


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

		#cat("rank call:", rank.call, "\n")

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
	if (!is.null(deviance(models[[1]])))
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

	selection.table <- data.frame (AICc = aicc, Delta = delta, Weight = weight)

	if (!is.null(dev)) {
          selection.table <- cbind(Deviance = dev, selection.table)
	}

	all.par <- unique(unlist(lapply(models, function(m) names(coeffs(m)))))

	all.terms <- unique(unlist(lapply(models, getAllTerms)))
	all.terms <- all.terms[order(sapply(gregexpr(":", all.terms), function(x) if(x[1] == -1) 0 else length(x)), all.terms)]

	all.par <- all.par[order(sapply(gregexpr(":", all.par), function(x) if(x[1] == -1) 0 else length(x)), all.par)]

	all.coef <- all.var <- all.df <- numeric(0)

	ac <- rep(0, length = length(all.par))

 	for (m in models) {
		m.tTable <- tTable(m)
		n <- length(resid(m))
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

	avg.model <- t(sapply(seq_along(all.par), function(i) par.avg(all.coef[,i], all.var[,i], all.df, weight)))

	all.coef[all.coef == 0] <- NA
	all.var[all.var == 0] <- NA

	importance <- sort(importance, decreasing=T)
	colnames(all.coef) <- colnames(all.var) <- rownames(avg.model) <-  all.par

     names(all.terms) <- seq_along(all.terms)


	if (!is.null(rank)) {
		colnames(selection.table)[2] <- as.character(rank)
	}

	ret <- list(
		summary = selection.table,
		coefficients = all.coef,
		variable.codes = all.terms,
		variance = all.var,
		avg.model = avg.model,
		relative.importance = importance,
		weights = weight,
		beta = beta,
		terms = all.par

	)

	class(ret) <- "averaging"
	return(ret)
}

par.avg <- function(x, se, npar, weight, alpha = 0.05) {

	wx <- weighted.mean(x, weight)
	x.sqdiff <- (x - wx)^2

	xvar <- se^2
	avar <- weighted.mean(xvar + x.sqdiff, weight)
	ase <- weighted.mean(sqrt(xvar + x.sqdiff), weight)

	z <- c((qt(1 - (alpha / 2), npar) / qnorm(1 - (alpha / 2)))^2)

	use <- weighted.mean(sqrt((xvar * z) + x.sqdiff), weight)

	ci <- qnorm(1 - (alpha / 2)) * use

	return(c(`Coefficient` = wx, `Variance` = avar,  `SE` = ase, `Unconditional SE` = use, `Lower CI` = wx - ci, `Upper CI` = wx + ci))
}



# x.sqdiff = (x - WM(x))^2
# avar = WM( sqrt(se^2 + x.sqdiff) )
# use = WM(sqrt((xvar * z) + x.sqdiff) )


`print.averaging` <-
function(x, ...) {
	cat("\nModel summary:\n")
	print(signif(x$summary,3))

	cat("\nVariables:\n")
	print(x$variable.codes, quote= F)

	cat("\nAveraged model parameters:\n")
	print(signif(x$avg.model, 3))

	cat("\nRelative variable importance:\n")
	print(round(x$relative.importance, 2))
}

`print.model.selection` <-
function(x, ...) {
	x$weight <- round(x$weight, 3)

	nn <- attr(x, "terms")

	names(x)[seq(along=nn)] <- sapply( strsplit(nn, ":"), function(xx) paste(sapply(xx, abbreviate, 6 )  , collapse=":") )

	cat ("Model selection table", "\n")
	print(signif(as.matrix(x), digits=4), na.print="")
	if (!is.null(attr(x, "random.terms"))) {
		cat("Random terms:", paste(attr(x, "random.terms"), collapse=", "), "\n")
 	}
}

`tTable` <-
function (model, ...) {
	UseMethod("tTable")
}

`tTable.default` <-
function(model, ...) {
	return(summary(model)$coefficients)
}

`tTable.lme` <-
function(model, ...) {
	return(summary(model)$tTable)
}

`tTable.mer` <- `tTable.glmer` <- `tTable.lmer` <-
function(model, ...) {
	sm <- eval(expression(summary), environment(lmer))
	return (sm(model)@coefs)
	#return((lme4::summary(model))@coefs)
}

`tTable.spautolm` <-
function(model, ...) {
return(summary(model)$Coef)
}

`tTable.gam` <-
function(model, ...) {
	cf <- model$coefficients
	se <- summary(model)$se
     return(cbind(`Estimate`=cf, `Std. Error` = se))
}

`tTable.gls` <-
function (model, ...) return(summary(model)$tTable)
