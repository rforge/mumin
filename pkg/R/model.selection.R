`coefTable.model.selection` <-
function (model, ...) {
	#structure(attr(model, "coefTables"), names = rownames(model))
	ret <- attr(model, "coefTables")
	names(ret) <- rownames(model)
	ret
}

`coef.model.selection` <-
function (object, ...) {
	ct <- attr(object, "coefTables")
	n <- length(ct)
	allcf <- unique(unlist(lapply(ct, rownames)))
	ret <- matrix(NA_real_, nrow = n, ncol = length(allcf),
		dimnames = list(rownames(object), allcf))
	for(i in seq_len(n))
		ret[i, match(rownames(ct[[i]]), allcf)] <- ct[[i]][, 1L]
	ret
}

`coeffs.model.selection` <-
function (model) coef.model.selection(model)

`coefArray` <- function(object) {
	coefNames <- fixCoefNames(unique(unlist(lapply(object, rownames),
		use.names = FALSE)))
	nCoef <- length(coefNames)
	nModels <- length(object)
	ret <- array(NA_real_, dim = c(nModels, 3L, nCoef),
		dimnames = list(names(object), c("Estimate", "Std. Error", "df"), coefNames))
	for(i in seq_along(object)) {
		z <- object[[i]]
		ret[i, 1:3, ]
		ret[i, seq_len(ncol(z)), ] <- t(z[match(coefNames, fixCoefNames(rownames(z))), ])
	}
	ret
}

`getCall.model.selection`  <-
function (x, i = NULL, ...) {
	if(is.null(i)) return(attr(x, "call", exact = TRUE))
	if(length(i) == 1L) return(attr(x, "model.calls", exact = TRUE)[[i]])
	return(attr(x, "model.calls", exact = TRUE)[i])
}

evalSubsetExpr <-
function(ss, dfr) {
	ss <- .exprapply(.exprapply(.exprapply(
		ss,
		"dc", .sub_dc_has, as.name(".subset_vdc")),
		c("{", "Term"), .sub_Term),
		"has", .sub_has)
	ss <- subst(ss, . = dfr)
	DebugPrint(ss)
	eval(ss, dfr)
}

#TODO: update 'subset' to new '['
`subset.model.selection` <-
function(x, subset, select, recalc.weights = TRUE, recalc.delta = FALSE, ...) {
	
	if (missing(select)) {
		if(missing(subset)) return(x)
		i <- evalSubsetExpr(substitute(subset), x)
		return(`[.model.selection`(x, i, recalc.weights = recalc.weights, 
			recalc.delta = recalc.delta, ...))
	} else {
		cl <- match.call(expand.dots = FALSE)
		if(!missing(subset)) cl$subset <- evalSubsetExpr(substitute(subset), x)
		
	    cl <- cl[c(1L, match(names(formals("subset.data.frame")), names(cl), 0L))]
	    cl[[1L]] <- as.name("subset.data.frame")
		DebugPrint(cl)
		ret <- eval.parent(cl)
		if(recalc.weights && ("weight" %in% colnames(ret)))
			ret[, 'weight'] <- ret[, 'weight'] / sum(ret[, 'weight'])
		if(recalc.delta && ("delta" %in% colnames(ret)))
			ret[, 'delta'] <- ret[, 'delta'] - min(ret[, 'delta'])
	    return(ret)
	}
}



getModelClass <-
function(x) {
	if(inherits(x, "model.selection")) {
		if(!is.null(attr(x, "global"))) return(class(attr(x, "global"))[1L])
		if("class" %in% colnames(x)) return(as.character(x[, "class"]))
		if(!is.null(attr(x, "model.class"))) return(attr(x, "model.class"))
	}
	return(NULL)
}

`merge.model.selection` <-
function (x, y, suffixes = c(".x", ".y"), ...) {
	rval <- rbind(x, y, make.row.names = FALSE)
	if (!is.null(suffixes)) row.names(rval) <-
		c(paste0(row.names(x), suffixes[1L]), 
            paste0(row.names(y), suffixes[2L]))
	rval
}

`rbind.model.selection` <- 
function (..., deparse.level = 1, make.row.names = TRUE) {
	allargs <- list(...)
	n <- length(allargs) 
	if(n == 1L) return(allargs[[1L]])

	if(!all(vapply(allargs, inherits, FALSE, "model.selection")))
		stop("need all \"model.selection\" objects")

	allargs <- lapply(allargs, "class<-", "data.frame") ### XXX: NO COPYING - modifies
											            ### original objects!!!
	on.exit({
		lapply(allargs, "class<-", c("model.selection", "data.frame"))
	})
	
	allitemsidentical <- function(x) all(vapply(x[-1L], identical, FALSE, x[[1L]]))
	
	if(!allitemsidentical(lapply(allargs, "attr", "rank.call")))
		stop("tables are not ranked by the same IC")
	if(!allitemsidentical(lapply(allargs, "attr", "nobs")))
		stop("models are fitted to different number of observations")
	

	.combine <-
	function(x, y, pos, len = length(y)) {
		if(is.factor(x) || is.factor(y)) {
			if(is.factor(x)) {
				if(!is.factor(y)) y <- factor(y)
			} else if(is.factor(y)) x <- factor(x)
			alllev <- unique(c(levels(x), levels(y)))
			x <- factor(x, levels = alllev, labels = alllev)
		}
		x[pos:(pos + len - 1L)] <- y
		x
	}
	
	ct <- unname(lapply(allargs, attr, "column.types"))
	vct <- unlist(ct)
	vct <- vct[order(as.integer(unlist(ct)), unlist(lapply(ct, seq_along)))]
	vct <- vct[!duplicated(names(vct))]
	# TODO: check mismatch in column.types
	nm <- names(vct)

	rval <- as.data.frame(array(NA, dim = c(sum(sapply(allargs, nrow)), length(nm)),
								dimnames = list(NULL, nm)))
	row1 <- 1L
	for(z in allargs) {
		n <- nrow(z)
		nmz <- nm[nm %in% names(z)]
		for(j in nmz) rval[, j] <- .combine(rval[, j], z[, j], row1, n)
		row1 <- row1 + n
	}

	newattr <- list(column.types = vct)
	for(i in c("model.calls", "coefTables"))
		newattr[[i]] <- unlist(lapply(allargs, attr, i), recursive = FALSE, use.names = FALSE)
	k <- c("rank", "rank.call", "nobs")
	newattr[k] <- attributes(allargs[[1L]])[k]
	
	tmp <- lapply(allargs, attr, "terms")
	newattr[["terms"]] <- structure(unique(unlist(tmp, recursive = FALSE, use.names = FALSE)),
			  interceptLabel = unique(unlist(lapply(tmp, attr, "interceptLabel"))))
	
	for(i in names(newattr)) attr(rval, i) <- newattr[[i]]
	class(rval) <- c("model.selection", "data.frame")
	if(make.row.names) {
		rn1 <- rep(names(allargs), sapply(allargs, nrow))
		rn1[i] <- paste0(rn1[i <- rn1 != ""], ".")
		rlabs <- paste0(rn1, unlist(lapply(allargs, rownames)))
		if(anyDuplicated(rlabs))
			rlabs <- make.unique(as.character(rlabs), sep = "")
	} else {
		rlabs <- as.character(1L:nrow(rval))
	}
	rownames(rval) <- rlabs	

	o <- order(rval[, names(vct)[vct == "ic"]])
	rval <- rval[o, recalc.delta = TRUE]
	rval
}

`print.model.selection` <-
function(x, abbrev.names = TRUE, warnings = getOption("warn") != -1L, ...) {
	#cat("new\n")
	origx <- x
	class(x) <- "data.frame"
	xterms <- attr(origx, "terms")
	if(is.null(xterms) || !all(xterms %in% colnames(x)[seq_along(xterms)])) {
		print.data.frame(x, ...)
	} else {
		if(abbrev.names) xterms <- abbreviateTerms(xterms, 6L, 3L, deflate = TRUE)
		colnames(x)[seq_along(xterms)] <- xterms
		globcl <- attr(origx, "global.call")
		if(!is.null(globcl)) {
			cat("Global model call: ")
			print(globcl)
			cat("---\n")
			random.terms <- attr(getAllTerms(attr(origx, "global")), "random.terms")
			if(!is.null(random.terms)) random.terms <- list(random.terms)
		} else random.terms <- attr(origx, "random.terms")
		cat("Model selection table \n")
		
		dig <- c(varying = NA, extra = NA, df = 0L, loglik = 3L, ic = 1L, delta = 2L,
				 weight = 3L, terms = NA)
		column.types <- attr(origx, "column.types")
		decprint <- dig[column.types[colnames(x)]]

		i <- vapply(x, is.numeric, FALSE) & is.na(decprint)
		x[, i] <- signif(x[, i], 4L)
		k <- which(!is.na(decprint))
		for(i in k) x[, i] <- round(x[, i], digits = decprint[i])
			
		vLegend <- NULL
		if(abbrev.names) {
			
			vCols <- names(column.types)[column.types == "varying"]
			vCols <- vCols[(vCols %in% colnames(x)) & !(vCols %in% c("class"))]
			vlen <- nchar(vCols)
			vLegend <- vector(length(vCols), mode = "list")
			names(vLegend) <- vCols
			if(!is.null(vCols)) {
				for(i in vCols) {
					if(!is.factor(x[, i])) next
					lev <- levels(x[, i])
					lev <- lev[!(lev %in% c("", "NULL"))]
					shlev <- abbreviateTerms(lev, nchar(i), deflate = TRUE)
					x[, i] <- factor(x[, i], levels = lev, labels = shlev)
					if(any(j <- shlev != lev)) vLegend[[i]] <-
						paste(shlev[j], "=", sQuote(lev[j]))
				}
				vLegend <- vLegend[!vapply(vLegend, is.null, TRUE)]
			}
			
		}

		uqran <- unique(unlist(random.terms, use.names = FALSE))
		abbran <- abbreviateTerms(gsub("1 | ", "", uqran, fixed = TRUE), 1L,
			deflate = TRUE)
		colran <- vapply(random.terms, function(s) paste(abbran[match(s, uqran)],
			collapse = "+"), "")

		if(addrandcol <- length(unique(colran)) > 1L) {
			k <- which(colnames(x) == "df")[1L]
			x <- cbind(x[, 1L:(k - 1L)], random = colran, x[, k:ncol(x)])
		}

		print.default(as.matrix(x)[, !vapply(x, function(y) all(is.na(y)), FALSE),
			drop = FALSE], na.print = "", quote = FALSE, right = TRUE)

		if(abbrev.names && length(vLegend)) {
			cat("Abbreviations:", sep = "\n")
			for(i in names(vLegend)) {
				cat(vLegend[[i]], sep = ", ", fill = TRUE, labels =
					c(paste0(i, ":"), rep(paste(rep(" ", nchar(i) + 1L),
					collapse = ""), length(vLegend[[i]]) - 1L)))
			}
		}
		
		cat("Models ranked by", asChar(attr(attr(origx, 'rank'), "call")), "\n")

		if(!is.null(random.terms)) {
			if(addrandcol) {
				cat("Random terms: \n")
				cat(paste(abbran, "=", sQuote(uqran)), sep = "\n")
			} else {
				cat("Random terms (all models): \n")
				cat(paste(sQuote(uqran)), sep = ", ")
				cat("\n")
			}

		}
		if (warnings && !is.null(attr(origx, "warnings"))) {
			cat("\n"); print.warnings(attr(origx, "warnings"))
		}
	}
	invisible(origx)
}

`update.model.selection` <- function (object, global.model, ..., evaluate = TRUE) {
    cl <- attr(object, "call")
    if (is.null(cl)) stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...

	if(!missing(global.model))
		extras <- c(list(global.model = substitute(global.model)), extras)

    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(cl)))
        for (a in names(extras)[existing]) cl[a] <- extras[a]
        if (any(!existing)) {
            cl <- c(as.list(cl), extras[!existing])
            cl <- as.call(cl)
        }
    }
    return(if (evaluate) eval.parent(cl) else cl)
}

`logLik.model.selection` <- function (object, ...) {
	nobs <- attr(object, "nobs")
	n <- nrow(object)
	ret <- vector(n, mode = "list")
	for(i in 1:n) ret[[i]] <-
		structure(object[i, "logLik"], df = object[i, "df"], nobs = nobs,
			class = "logLik")
	ret
}

`family.model.selection` <-
function (object, ...) {
	if(!is.null(attr(object, "global"))) {
		model.calls <- attr(object, "model.calls")
		if(!is.null(model.calls[[1L]][["family"]])) {
			fam <- lapply(model.calls, "[[", "family")
			fam1 <- unique(fam)
			ret <- lapply(unique(fam), eval)[
				as.integer(as.factor(vapply(fam, asChar, "")))
				]
			names(ret) <- rownames(object)
			#index <- split(seq_along(fam), vapply(fam, asChar, ""))
			#for(i in seq_along(fam1)) fam1[[i]] <- list(family = eval(fam1[[i]]), index = index[[i]])
			#fam <- family(dd1)
			#index <- lapply(fam, "[[", "index")
			#ret <- rep(lapply(fam, "[[", "family"), vapply(index, length, 1L))[order(unlist(index))]
			return(ret)
		} else return(family(attr(object, "global")))
	} else {
		return(attr(object, "model.family"))
	}
}

#### XXX
.argTable <-
function(cl, family = NULL, class = NULL,
		 args.omit = NULL, different.only = FALSE) {
	haveNoCall <-  vapply(cl, is.null, FALSE)
	cl[haveNoCall] <- lapply(cl[haveNoCall], function(x) call("<unknown>", formula = NA))
	arg <- lapply(cl, function(x) sapply(x, function(argval)
		switch(mode(argval), character = , logical = argval,
		numeric = signif(argval, 3L), asChar(argval))))
	arg <- rbindDataFrameList(lapply(lapply(arg, t), as.data.frame))
	if(!is.null(args.omit)) arg <- arg[, !(colnames(arg) %in% args.omit)]

	arg[] <- lapply(arg, as.factor)
	
	if(!is.null(family)) {
		.getFam <- function(x) unlist(x[c("family",	"link")])
		fam <-  if(inherits(family, "family"))
			matrix(.getFam(family), dimnames = list(c("family", "link"), NULL)) else
			sapply(family, .getFam)
		f <- fam[1L, ]
		f[is.na(f)] <- ""
		f <- vapply(strsplit(f, "(", fixed = TRUE), "[", "", 1L)
		f[f == "Negative Binomial"] <- "negative.binomial"
		fam[2L, fam[2L, ] == vapply(unique(f), function(x) if(is.na(x))
									NA_character_ else formals(get(x))$link,
									FUN.VALUE = "")[f]] <- NA_character_
		j <- !is.na(fam[2L,])
		famname <- fam[1L, j]
		famname <- ifelse(substring(famname, nchar(famname)) != ")",
			paste0(famname, "("), paste0(substring(famname, 1L, nchar(famname) - 1L),
				", "))
		fam[1L, j] <- paste0(famname, fam[2L, j], ")")
		arg <- cbind(arg, t(fam))
	}
	if(!is.null(class)) arg[, "class"] <- rep(class, length.out = nrow(arg))

	#arg <- as.matrix(arg)
	#arg[is.na(arg) | arg == "NULL"] <- ""
	colnames(arg)[1L] <- "FUN"
	for (i in seq_len(ncol(arg))) {
		v <- arg[, i]
		if(any(j <- is.na(v) | v == "NULL" | v == "")) {
			levels(v) <- c(levels(v), "")
			v[j] <- ""
			arg[, i] <- v
		}
	}
	
	if(different.only)
		arg <- arg[, vapply(arg, nlevels, 1L) != 1L, drop = FALSE]

	#if(ncol(arg) != 0L) arg <- gsub("([\"'\\s]+|\\w+ *=)","", arg, perl = TRUE)
	arg
}

`nobs.model.selection` <-
function (object, ...)
attr(object, "nobs")

`row.names<-.model.selection` <-
function (x, value)  {
	oldnames <- dimnames(x)[[1L]]
	x <- NextMethod()
	newnames <- dimnames(x)[[1L]]
	rowattrib <- c("model.calls", "coefTables", "random.terms", "order",
		if(!is.null(attr(x, "modelList")))"modelList")
	for(i in rowattrib) if(!is.null(attr(x, i))) names(attr(x, i)) <- newnames
	x
}

`names<-.model.selection` <-
function (x, value) {
	oldnames <- names(x)
	if(any(attr(x, "column.types")[oldnames[oldnames != value]] %in%
	   c('df', 'loglik', 'ic', 'delta', 'weight', 'terms'))) {
		class(x) <- "data.frame"
		message("rename protected column ->data.frame")
	}
	NextMethod()
}

verify_model_selection_object <-
function(x, attrib, modif = NULL, rowchange = TRUE) {
	#message(sys.call()[[1L]])
	protectedcoltypes <- c('df', 'loglik', 'ic', 'delta', 'weight', 'terms')
	excludeattr <- c("names", "row.names", "class")
	column.types <- attrib[["column.types"]]
	keepattr <- names(attrib)[!(names(attrib) %in% excludeattr)]
	.setattr <- function(x, newattr = NULL, which = keepattr) {
		#for(i in which) attr(x, i) <- newattr[[i]]
		attributes(x)[which] <- if(is.null(newattr)) NULL else newattr[which]
		x
	}
	if(inherits(x, "model.selection")) {
		#message("is model.selection")
		if(!is.null(modif) && any(modif %in% names(column.types)[
			column.types %in% protectedcoltypes])) {
			#message("modifying protected column ", sQuote(modif), ", ->data.frame")
			class(x) <- "data.frame"
			return(.setattr(x))
		} else {
			s <- dimnames(x)[[2L]]
			k <- match(names(column.types), colnames(x), nomatch = 0L)
			if(any(column.types[k == 0L] %in% protectedcoltypes)) {
				#message("protected columns missing ", paste(names(column.types)[k == 0], collapse = ", "),
						#" ->data.frame")
				class(x) <- "data.frame"
				return(.setattr(x))
			} else {
				if(any(column.types[k == 0L] %in% c("varying", "extra"))) {
					#message("varying/extra columns missing")
					column.types <- column.types[k != 0L]
					attrib[["column.types"]] <- column.types
				}
			}
		}
		
		oldrownames <- attrib[['row.names']]
		newrownames <- dimnames(x)[[1L]]
		if(rowchange && (length(oldrownames) != length(newrownames) ||
						 any(oldrownames != newrownames))) {
			rowattrib <- c("model.calls", "coefTables", "random.terms", "order",
			   if(!is.null(attr(x, "modelList")))"modelList")
			k <- match(newrownames, oldrownames)
			#message("row order changed: ", k)
			attrib[rowattrib] <- lapply(attrib[rowattrib], `[`, k)
		}# else {
			#message("rows not changed")
		#}
		
		x <- .setattr(x, attrib)
		if(!is.null(warningList <- attrib$warnings))
			attr(x, "warnings") <- warningList[sapply(warningList, attr, "id")
											   %in% newrownames]

	} else {
		#message("not model.selection")
		return(.setattr(x))
	}
	x
}

`[<-.model.selection` <-
function (x, i, j, value)  {
	verify_model_selection_object(NextMethod("[<-"),
		attributes(x), if(is.character(j)) j else colnames(x)[j])
}

`$<-.model.selection` <-
function (x, name, value) {
	verify_model_selection_object(NextMethod("$<-"), attributes(x), name, rowchange = FALSE)
}

`[.model.selection` <-
function (x, i, j, recalc.weights = TRUE, recalc.delta = FALSE, ...) {
	x <- verify_model_selection_object(`[.data.frame`(x, i, j, ...), attributes(x))
	if(inherits(x, "model.selection")) {
		column.types <- attr(x, "column.types")
		nct <- names(column.types)
		ic <- `[.data.frame`(x, , nct[column.types == "ic"])
		if(recalc.weights) x <- `[<-.data.frame`(x, , nct[column.types == "weight"], Weights(ic))
		if(recalc.delta) x <- `[<-.data.frame`(x, , nct[column.types == "delta"], ic - min(ic))
	}
	x
}