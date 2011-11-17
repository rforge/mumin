`pdredge` <-
function(global.model, cluster = NA, beta = FALSE, evaluate = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL,
		 trace = FALSE, varying, extra, check = FALSE,
		  ...) {
#FIXME: m.max cannot be 0 - e.g. for intercept only model

	qlen <- 25L
	# Imports: clusterCall, clusterApply
	doParallel <- inherits(cluster, "cluster")
	if(doParallel) {
		.parallelPkgCheck() # XXX: workaround to avoid importing from 'parallel'
		clusterCall <- get("clusterCall")
		clusterApply <- get("clusterApply")
		clusterCall(cluster, "require", "MuMIn", character.only = TRUE)
		clusterCall(cluster, "eval", expression(assign("parGetMsRow", 
			MuMIn:::parGetMsRow, envir = .GlobalEnv), NULL))
		.getRow <- function(X) clusterApply(cluster, X, fun = "parGetMsRow")
	} else {
		.getRow <- function(X) lapply(X, parGetMsRow, parCommonProps)
		clusterCall <- function(...) NULL
		message("Not using cluster.")
	}

	# *** Rank ***
	rank.custom <- !missing(rank)
	rankArgs <- list(...)
	IC <- .getRank(rank, rankArgs)
	ICName <- as.character(attr(IC, "call")[[1L]])

	allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE)

	# Intercept(s)
	interceptLabel <- attr(allTerms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"
	nInts <- sum(attr(allTerms, "intercept"))

	if(length(grep(":", all.vars(reformulate(allTerms))) > 0L))
		stop("variable names in the formula cannot contain \":\"")

	gmEnv <- parent.frame()
	gmCall <- .getCall(global.model)
	if (is.null(gmCall)) {
		gmCall <- substitute(global.model)
		if(!is.call(gmCall)) {
			if(inherits(global.model, c("gamm", "gamm4")))
				message("to use gamm models with 'dredge', use 'MuMIn::gamm' wrapper")
			stop("could not retrieve the call to 'global.model'")
		}
		#"For objects without a 'call' component the call to the fitting function \n",
		#" must be used directly as an argument to 'dredge'.")
		# NB: this is unlikely to happen:
		if(!exists(as.character(gmCall[[1L]]), parent.frame(), mode="function"))
			 stop("could not find function '", gmCall[[1L]], "'")
	} else {
		# if 'update' method does not expand dots, we have a problem
		# with expressions like ..1, ..2 in the call.
		# So, try to replace them with respective arguments in the original call
		is.dotted <- grep("^\\.\\.", sapply(as.list(gmCall), deparse))
		if(length(is.dotted) > 0L) {
			substGmCall <- substitute(global.model)
			if(is.name(substGmCall)) {
				stop("call to 'global.model' contains '...' arguments and cannot be updated: ",
					deparse(gmCall, control = NULL))
			} else {
				gmCall[is.dotted] <- substitute(global.model)[names(gmCall[is.dotted])]
			}
		}
	}
	logLik <- .getLogLik()

	# parallel: check whether the models would be identical:
	if(doParallel) testUpdatedObj(cluster, global.model, gmCall, do.eval = check)

	gmCoefNames0 <- names(coeffs(global.model))

	# Check for na.omit
	if (!is.null(gmCall$na.action) &&
		as.character(gmCall$na.action) %in% c("na.omit", "na.exclude")) {
		stop("'global.model' should not use 'na.action' = ", gmCall$na.action)
	}

	if(names(gmCall)[2L] == "") names(gmCall)[2L] <- 
		names(formals(deparse(gmCall[[1L]]))[1L])

	gmCoefNames <- fixCoefNames(gmCoefNames0)
	n.vars <- length(allTerms)

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		warning("comparing models with different fixed effects fitted by REML")

	if (beta && is.null(tryCatch(beta.weights(global.model), error=function(e) NULL,
		warning = function(e) NULL))) {
		warning("do not know how to calculate beta weights for ",
				class(global.model)[1L], ", argument 'beta' ignored")
		beta <- FALSE
	}

	m.max <- if (missing(m.max)) (n.vars - nInts) else min(n.vars - nInts, m.max)

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1]] != "~" || length(fixed) != 2L)
				warning("'fixed' should be a one-sided formula")
			fixed <- c(getAllTerms(fixed))
		} else if (!is.character(fixed)) {
			stop ("'fixed' should be either a character vector with"
				  + " names of variables or a one-sided formula")
		}
		if (!all(fixed %in% allTerms)) {
			warning("not all terms in 'fixed' exist in 'global.model'")
			fixed <- fixed[fixed %in% allTerms]
		}
	}
	fixed <- c(fixed, allTerms[allTerms %in% interceptLabel])
	n.fixed <- length(fixed)
	termsOrder <- order(allTerms %in% fixed)
	ordAllTerms <- allTerms[termsOrder]
	allTerms <- ordAllTerms
	isMER <- any(inherits(global.model, c("mer", "lmer", "glmer")))
	gmFormula <- as.formula(formula(global.model))
	gmFormulaEnv <- attr(gmFormula, ".Environment")

	### BEGIN:
	## varying BEGIN
	if(!missing(varying) && !is.null(varying)) {
		variantsIdx <- expand.grid(lapply(varying, seq_along))
		seq.variants <- seq.int(nrow(variantsIdx))
		nvarying <- length(varying)
		varying.names <- names(varying)
	} else {
		variantsIdx <- NULL
		seq.variants <- 1L
		nvarying <- 0L
		varying.names <- character(0L)
	}
	nvariants <- length(seq.variants)
	## varying END

	## extra BEGIN
	if(!missing(extra) && length(extra) != 0L) {
		extraNames <- sapply(extra, function(x) switch(mode(x),
			call = deparse(x[[1]]), name = deparse(x), character = , x))
		if(!is.null(names(extra)))
			extraNames <- ifelse(names(extra) != "", names(extra), extraNames)

		extra <- structure(as.list(unique(extra)), names = extraNames)
		if("R^2" %in% extra) {
			null.fit <- null.fit(global.model, TRUE, gmFormulaEnv)
			extra[extra == "R^2"][[1L]] <- function(x) r.squaredLR(x, null.fit)
		}
		extra <- sapply(extra, match.fun, simplify = FALSE)
		applyExtras <- function(x) unlist(lapply(extra, function(f) f(x)))
		extraResult <- applyExtras(global.model)
		if(!is.numeric(extraResult))
			stop("function in 'extra' returned non-numeric result")

		nextra <- length(extraResult)
		extraNames <- names(extraResult)
	} else {
		nextra <- 0L
		extraNames <- character(0L)
	}
	## extra END

	nov <- as.integer(n.vars - n.fixed)
	ncomb <- 2L ^ nov
	if(nov > 31L) stop(gettextf("maximum number of predictors is 31, but %d is given", nov))
	if(nov > 10L) warning(gettextf("%d predictors will generate up to %.0f possible combinations", nov, ncomb))
	nmax <- ncomb * nvariants
	if(evaluate) {
		ret.nchunk <- 25L
		ret.ncol <- n.vars + nvarying + 3L + nextra
		ret <- matrix(NA_real_, ncol = ret.ncol, nrow = ret.nchunk)
	} else {
		ret.nchunk <- nmax
	}

	calls <- vector(mode = "list", length = ret.nchunk)

	if(hasSubset <- !missing(subset))  {
		if(!tryCatch(is.language(subset), error = function(e) FALSE))
			subset <- substitute(subset)
		if(inherits(subset, "formula")) {
			if (subset[[1]] != "~" || length(subset) != 2L)
				stop("'subset' should be a one-sided formula")
			subset <- subset[[2L]]
		}
		if(!all(all.vars(subset) %in% allTerms))
			warning("not all terms in 'subset' exist in 'global.model'")
	}

	comb.sfx <- rep(TRUE, n.fixed)
	comb.seq <- if(nov != 0L) seq_len(nov) else 0L
	k <- 0L
	ord <- extraResult1 <- integer(0L)

	argsOptions <- list(
		response = attr(allTerms0, "response"),
		intercept = nInts,
		interceptLabel = interceptLabel,
		random = attr(allTerms0, "random"),
		gmCall = gmCall,
		gmEnv = gmEnv,
		allTerms = allTerms0,
		gmCoefNames = gmCoefNames,
		gmDataHead = if(!is.null(gmCall$data)) {
			if(eval(call("is.data.frame", gmCall$data), gmEnv))
				eval(call("head", gmCall$data, 1L), gmEnv) else gmCall$data
			} else NULL,
		gmFormulaEnv = gmFormulaEnv
		)

	# BEGIN parallel
	qi <- 0L
	queued <- vector(qlen, mode = "list")
	parCommonProps <- list(gmEnv = gmEnv, IC = IC, beta = beta,
		allTerms = allTerms, nextra = nextra)
	if(nextra) {
		parCommonProps$applyExtras <- applyExtras
		parCommonProps$extraResult <- extraResult
	}
	if(doParallel) clusterVExport(cluster, clustDredgeProps = parCommonProps)
	# END parallel
	
	retColIdx <- if(nvarying) -n.vars - seq_len(nvarying) else TRUE

	for(iComb in seq.int(ncomb)) {
		comb <- c(as.logical(intToBits(iComb - 1L)[comb.seq]), comb.sfx)

		nvar <- sum(comb) - nInts
		if(!(nvar > m.max || nvar < m.min) && (!hasSubset || eval(subset, 
			structure(as.list(comb), names = allTerms)))) {
		
			newArgs <- makeArgs(global.model, allTerms[comb], comb, argsOptions)
			formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
				attr(newArgs, "formulaList")
			if(all(vapply(formulaList, formulaAllowed, logical(1L), marg.ex))) {
				if(!is.null(attr(newArgs, "problems"))) {
					print.warnings(structure(vector(mode = "list",
						length = length(attr(newArgs, "problems"))),
							names = attr(newArgs, "problems")))
				} # end if <problems>

				cl <- gmCall
				cl[names(newArgs)] <- newArgs

				for (v in seq.variants) { ## --- Variants ---------------------------
					modelId <- ((iComb - 1L) * nvariants) + v
					clVariant <- cl
					if(nvarying) {
						newVaryingArgs <- sapply(varying.names, function(x)
							varying[[x]][[variantsIdx[v, x]]], simplify = FALSE)
						clVariant[varying.names] <- newVaryingArgs
					}
					
					if(trace) {
						cat(modelId, ": "); print(clVariant);
						utils::flush.console()
						}
					if(evaluate) {
						qi <- qi + 1L
						queued[[(qi)]] <- list(call = clVariant, id = modelId)
					} else { # if !evaluate
						k <- k + 1L # all OK, add model to table
						calls[[k]] <- clVariant
					}
				}
			} # end if <formulaAllowed>
		} # end if <subset, m.max >= nvar >= m.min>
		
		if(evaluate && qi && (qi + nvariants > qlen || iComb == ncomb)) {
			#DebugPrint(paste(qi, nvariants, qlen, iComb, ncomb))
			qseq <- seq_len(qi)
			qresult <- .getRow(queued[qseq])
			cat(sprintf("queue done: %d\n", qi)) # DEBUG
			withoutProblems <- which(!sapply(qresult, inherits, "condition"))
			if(!length(withoutProblems)) withoutProblems <- seq_along(qresult)
			#for (i in qresult[!withoutProblems]) warning(i)
			lapply(qresult[-withoutProblems], warning)
			qresult <- qresult[withoutProblems]
			qresultLen <- length(qresult)
			retNrow <- nrow(ret)
			if(k + qresultLen > retNrow) {
				nadd <- min(ret.nchunk, (ncomb * nvariants) - retNrow)
				ret <- rbind(ret, matrix(NA, ncol = ret.ncol, nrow = nadd), deparse.level = 0L)
				calls <- c(calls, vector("list", nadd))
				ord <- c(ord, integer(nadd))
			}
			qseq <- seq_len(qresultLen)
			for(m in qseq) ret[k + m, retColIdx] <- qresult[[m]]
			ord[k + qseq] <- vapply(queued[withoutProblems], "[[", 1L, "id")
			calls[k + qseq] <- lapply(queued[withoutProblems], "[[", "call")
			k <- k + qresultLen
			qi <- 0L
		}
	} ### for (iComb ...)

	if(k == 0L) stop("the result is empty")
	if(!evaluate) return(calls[seq_len(k)])

	if(k < nrow(ret)) {
		i <- seq_len(k)
		ret <- ret[i, , drop = FALSE]
		ord <- ord[i]
		calls <- calls[i]
	}
	
	if(nvarying) {
		varlev <- ord %% nvariants; varlev[varlev == 0L] <- nvariants
		ret[, n.vars + seq_len(nvarying) ] <- as.matrix(variantsIdx)[varlev, ]
	}
	
	ret <- as.data.frame(ret)
	row.names(ret) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(!(allTerms %in% gmCoefNames))
	ret[tfac] <- lapply(ret[tfac], factor, levels = NaN, labels = "+")

	i <- seq_along(allTerms)
	v <- order(termsOrder)
	ret[, i] <- ret[, v]
	allTerms <- allTerms[v]
	colnames(ret) <- c(allTerms, varying.names, extraNames, "df", "logLik", ICName)

	if(nvarying) {
		variant.names <- lapply(varying, function(x) 
			make.unique(if(is.null(names(x))) as.character(x) else names(x)))
		for (i in varying.names) ret[, i] <-
			factor(ret[, i], levels = seq_along(variant.names[[i]]),
				labels = variant.names[[i]])
	}

	o <- order(ret[, ICName], decreasing = FALSE)
	ret <- ret[o, ]
	ret$delta <- ret[, ICName] - min(ret[, ICName])
	ret$weight <- exp(-ret$delta / 2) / sum(exp(-ret$delta / 2))

	ret <- structure(ret,
		class = c("model.selection", "data.frame"),
		calls = calls[o],
		global = global.model,
		global.call = gmCall,
		terms = allTerms,
		rank = IC,
		rank.call = attr(IC,"call"),
		call = match.call(expand.dots = TRUE)
	)

	if (!is.null(attr(allTerms0, "random.terms")))
		attr(ret, "random.terms") <- attr(allTerms0, "random.terms")

	if(doParallel) clusterCall(cluster, "rm",
		list = c("parGetMsRow", "clustDredgeProps"), envir = .GlobalEnv)
	return(ret)
} ######

`parGetMsRow` <- function(modelProps, commonProps = get("clustDredgeProps", .GlobalEnv)) {
	### modelProps <- list(call = clVariant, id = modelId)

	fit1 <- tryCatch(eval(modelProps$call, commonProps$gmEnv), error = function(err) {
		err$message <- paste(conditionMessage(err), "(model",
			modelProps$id, "skipped)", collapse = "")
		class(err) <- c("simpleError", "warning", "condition")
		return(err)
	})
	if (inherits(fit1, "condition")) return(fit1)

	if(commonProps$nextra != 0L) {
		extraResult1 <- commonProps$applyExtras(fit1)
		if(length(extraResult1) < commonProps$nextra) {
			tmp <- rep(NA_real_, commonProps$nextra)
			tmp[match(names(extraResult1), names(commonProps$extraResult))] <-
				extraResult1
			extraResult1 <- tmp
		}
	} else extraResult1 <- NULL
	ll <- MuMIn:::.getLogLik()(fit1)
	c(	MuMIn:::matchCoef(fit1, all.terms = commonProps$allTerms,
			beta = commonProps$beta), # [commonProps$allTerms]
		extraResult1, df = attr(ll, "df"),
		ll = ll,
		ic = commonProps$IC(fit1)
	)
}

.test_pdredge <- function(dd) {
	cl <- attr(dd, "call")
	cl$cluster <- cl$check <- NULL
	cl[[1]] <- as.name("dredge")
	if(!identical(c(dd), c(eval(cl)))) stop("buuu...") 
	dd
}