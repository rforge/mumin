`getAllTerms.default` <-
function(x, ...) {
	return(getAllTerms(as.formula(formula(x))))
}

`getAllTerms.formula` <-
function(x, ...) {
	mTerms <- terms(x)
	intercept <- attr(mTerms, "intercept")
	if (!is.null(attr(mTerms, "offset"))){
		offs <-
		sapply(as.list(attr(mTerms,"variables")[-1])[attr(mTerms,"offset")],
		deparse)
	} else {
		offs <- NULL
	}
	ret <- attr(mTerms, "term.labels")

	if (length(ret) > 0) {

		#TODO: handle terms with ":" in the name
		if(length(grep(":", sapply(attr(mTerms, "variables")[-1],
		"as.character"))) > 0) {
			stop("Variable names cannot contain \":\"")
		}

		ia <- attr(mTerms, "order") > 1
		ret[ia] <- unlist(sapply(lapply(strsplit(ret[ia], ":", fixed = TRUE),
			sort), paste, collapse=":"))

		ret <- ret[order(attr(mTerms, "order"), ret)]

		#WTF?
		#i <- grep(" ", ret)
		#ret[i] <- paste("(", ret[i], ")")
	}

	attr(ret, "intercept") <- intercept

	if (!is.null(offs[1]))
		attr(ret, "offset") <- offs
	return(ret)
}

`getAllTerms.lme` <-
function(x, ...) {
	ret <- getAllTerms.formula(x)
	attr(ret, "random") <- . ~ .
	attr(ret, "random.terms") <- deparse(x$call$random)
	return(ret)

	#x <- fm2
	#reStruct <- x$modelStruct$reStruct
	#nobj <- length(reStruct)
	#if (is.null(namx <- names(reStruct)))
	#	names(reStruct) <- nobj:1
	#aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	#aux[lower.tri(aux)] <- ""
	#reStruct[] <- rev(reStruct)
	#aux <- t(array(rep(names(rs), nobj), c(nobj, nobj)))
	#aux[lower.tri(aux)] <- ""
	#rev(apply(aux, 1, function(x) paste(x[x != ""], collapse = " %in% ")))
	#paste(lapply(reStruct, attr, "formula"), "|",
	#	  rev(apply(aux, 1, function(x) paste(x[x != ""], collapse = " %in% "))))
}


`getAllTerms.glmer` <- # For backwards compatibility
`getAllTerms.lmer` <-  # with older versions of lme4
`getAllTerms.mer` <-
function(x, ...) {
	#fixed <- fixCoefNames(attr(terms(x), "term.labels"))

	frm <- as.formula(formula(x))
	#tt <- terms(frm)
	#allc <- attr(tt, "term.labels")
	allc <- getAllTerms.formula(frm)

	j <- grep(" | ", allc, fixed = TRUE)
	rnd <- allc[j]
	rnd.formula <- reformulate(c(".", paste("(", rnd, ")", sep = "")),
							   response = ".")

	ret <- fixCoefNames(allc[-j])
	attr(ret, "intercept") <- attr(allc, "intercept")
	attr(ret, "random") <- rnd.formula
	attr(ret, "random.terms") <- rnd
	return(ret)
}


`getAllTerms` <-
function(x, ...) UseMethod("getAllTerms")
