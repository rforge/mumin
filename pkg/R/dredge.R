`dredge` <-
function(global.model, beta = FALSE, eval = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, subset, ...) {

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
			(inherits(global.model, c("mer")) && ("REML" %in%
				names(deviance(global.model))))
		|| 	(inherits(global.model, c("lme")) && global.model$method == "REML")
		||  (any(inherits(global.model, c("lmer", "glmer")))
			  && global.model@status["REML"] != 0)
	) {
		warning("Comparing models with different fixed effects fitted by REML")
	}

	if (!is.lm && beta) {
		warning("Cannot calculate beta weigths (yet) for ", class(global.model)[1])
          beta <- FALSE
	}

	has.rsq <- "r.squared" %in% names(summary(global.model))
	has.dev <- !is.null(deviance(global.model))

	m.max <- if (missing(m.max)) n.vars else min(n.vars, m.max)

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1]] != "~" || length(fixed) != 2)
				warning(sQuote("fixed"), " formula should be of form ",
						dQuote("~ a + b + c"))
			fixed <- c(getAllTerms(fixed))
		} else if (!is.character(fixed)) {
			stop (sQuote("fixed"), " should be either a character vector with"
				  + " names of variables or a one-sided formula")
		}
		if (!all(fixed %in% all.terms)) {
			warning("Not all terms of ", sQuote("fixed"), " exist in ",
					sQuote("global.model"))
			fixed <- fixed[fixed %in% all.terms]
		}
		m.max <- m.max - length(fixed)
	}

	if (m.max > 0) {
		num.opt.vars <- 1:n.vars
		if (!is.null(fixed))
			num.opt.vars <- num.opt.vars[!(all.terms %in% fixed)]

		all.comb <- lapply(seq(m.max), combn, x = num.opt.vars)
		all.comb <- unlist(lapply(all.comb, function(.x) split(.x, col(.x))),
						   recursive = FALSE)
		all.comb <- c(`0` = list(0), all.comb)

	} else {
		all.comb <- list (0)
	}

	if (!is.null(fixed))
		all.comb <- lapply(all.comb, append, (1:n.vars)[all.terms %in% fixed])

	formulas <- lapply(all.comb, function(.x) reformulate(c("1", all.terms[.x]),
														  response = "." ))

	formulas <- lapply(formulas, `attr<-`, ".Environment",
				 attr(formula(global.model), ".Environment"))

	ss <- sapply(formulas, formulaAllowed)

	all.comb <- all.comb[ss]
	formulas <- formulas[ss]

	## Apply subset:
	if (!missing(subset)) {
		xtable <- as.data.frame(!is.na(do.call("rbind", lapply(all.comb, match,
			x = seq.int(n.vars)))))
		colnames(xtable) <- all.terms
		rownames(xtable) <- NULL

		ss <- eval(substitute(subset), envir = xtable)
		formulas <- formulas[ss]
		all.comb <- all.comb[ss]
		# 10 time slower! #xtable <- !is.na(t(sapply(all.comb, match, x=seq.int(n.vars))))
	}

	names(formulas) <- seq(formulas)
	if (any(inherits(global.model, c("mer", "lmer", "glmer")))) {
          formulas <- lapply(formulas, update, attr(all.terms, "random"))
	}

	if (!eval) {
		attr(formulas, "x") <- xtable
		return(formulas)
	}

	rank.custom <- rank != "AICc"

	###
	for(b in seq(length(all.comb))) {
		# print(all.comb[[b]])
        cterms <- all.terms[all.comb[[b]]]
		frm <- formulas[[b]]

		c.row <- rep(NA, n.vars)
		c.row[match(cterms, all.terms)] <- rep(1, length(cterms))

		cl <- call("update", substitute(global.model), frm)

		cmod <- tryCatch(eval(cl, parent.frame()), error=function(e) NULL)

		if (is.null(cmod)) {
			formulas[[as.character(b)]] <- NA
			next;
		}
		mod.coef <- c(na.omit(match(all.terms, names(coeffs(cmod)))))

		icept <- if (attr(all.terms, "intercept")) coeffs(cmod)[intercept]
			else NA

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
			c.row <- c(c.row, r.squared=cmod.summary$r.squared,
				adj.r.squared=cmod.summary$adj.r.squared)
		}
		if (has.dev)
			c.row <- c(c.row, deviance(cmod))

		if (rank.custom) {
			rankFnCall[[2]] <- cmod
			ic <- eval(rankFnCall)
			c.row <- c(c.row, IC=ic)
		} else {
		     c.row <- c(c.row, AIC=aic, AICc=aicc)
		}
		ms.tbl <- rbind(ms.tbl, c.row)
	}

	formulas[is.na(formulas)] <- NULL
	ms.tbl <- data.frame(ms.tbl, row.names=1:NROW(ms.tbl))

	# Convert columns with presence/absence of terms to factors
	tfac <- which(c(FALSE, !(all.terms %in% names(coeffs(global.model)))))
	ms.tbl[, tfac] <- lapply(ms.tbl[,tfac], factor, levels=1, labels="+")

	colnames(ms.tbl) <- c("(int.)", all.terms, "k",
		if (has.rsq) c("R.sq", "Adj.R.sq"),
		if (has.dev) if (is.lm) "RSS" else "Dev.",
		if (rank.custom) rank else c("AIC", "AICc")
	)

	o <- order(ms.tbl[, rank], decreasing = FALSE)
	ms.tbl <- ms.tbl[o,]
	ms.tbl$delta <- ms.tbl[, rank] - min(ms.tbl[, rank])
	ms.tbl$weight <- exp(-ms.tbl$delta / 2) / sum(exp(-ms.tbl$delta / 2))

	class(ms.tbl) = c("model.selection", "data.frame")

	attr(ms.tbl, "formulas") <- formulas[o]
	attr(ms.tbl, "global") <- global.model
	attr(ms.tbl, "terms") <- c(intercept, all.terms)

	if (rank.custom) {
		rankFnCall[[1]] <- as.name(rank)
		rankFnCall[[2]] <- substitute(global.model)
		attr(ms.tbl, "rank.call") <- rankFnCall
	}

	if (!is.null(attr(all.terms, "random.terms")))
		attr(ms.tbl, "random.terms") <- attr(all.terms, "random.terms")

	return(ms.tbl)
}

`subset.model.selection` <-
function(x, subset, select, ...) {
	cl <- match.call(expand.dots = FALSE)
    cl <- cl[c(1L, match(names(formals("subset.data.frame")), names(cl), 0L))]
    cl[[1L]] <- as.name("subset.data.frame")
    res <- eval(cl, parent.frame(2))

	if (missing(select)) {
		e <- substitute(subset)
		r <- eval(e, x, parent.frame(2))
		attr(res, "formulas") <- attr(x, "formulas")[r]
		attr(res, "global") <- attr(x, "global")
		attr(res, "terms") <- attr(x, "terms")
		class(res) <- class(x)
	} else {
		class(res) <- "data.frame"
	}
	return(res)
}

`[.model.selection` <-
function (x, i, j, ...) {
	res <- `[.data.frame`(x, i, j, ...)
	if (missing(j)) {
		attr(res, "global") <- attr(x, "global")
		attr(res, "terms") <- attr(x, "terms")
		attr(res, "formulas") <- attr(x, "formulas")[i]
	} else {
		class(res) <- "data.frame"
	}
	return(res)
}

`print.model.selection` <-
function(x, ...) {

	if(!is.null(x$weight))
		x$weight <- round(x$weight, 3)

	nn <- attr(x, "terms")

	if(is.null(nn)) {
		print.data.frame(x, ...)
	} else {
		names(x)[seq(along=nn)] <- sapply( strsplit(nn, ":"), function(xx) paste(sapply(xx, abbreviate, 6 )  , collapse=":") )

		cat ("Model selection table", "\n")
		i <- sapply(x, is.numeric)
		x[,i] <- signif(x[,i], 4)
		print.default(as.matrix(x), na.print="", quote=FALSE)
		if (!is.null(attr(x, "random.terms"))) {
			cat("Random terms:", paste(attr(x, "random.terms"), collapse=", "), "\n")
		}
	}
}
