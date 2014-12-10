isFALSE <- function(x) identical(FALSE, x)

stdize <-
function(x, ...) UseMethod("stdize")

stdize.default <-
stdize.numeric <-
function(x, center = TRUE, scale = TRUE, ...) {
	if(is.function(scale)) {
		scaleFunc <- scale
		scale <- TRUE
	} else scaleFunc <- function(x) sd(x, na.rm = TRUE)
	
	#if(length(list(...))) warning("additional arguments ignored")
	if(length(scale) != 1L) warning("only first element of 'center' is used")
	if(length(center) != 1L) warning("only first element of 'scale' is used")
	scale <- scale[1L]
	center <- center[1L]
	if(is.logical(scale)) scale <- if(scale) scaleFunc(x) else 1
	if(is.logical(center)) center <- if(center) mean(x, na.rm = TRUE) else 0
	x <- (x - center) / scale
	attr(x, "scaled:center") <-  center
	attr(x, "scaled:scale") <- scale
	x
}

stdize.matrix <-
function(x, center = TRUE, scale = TRUE, ...) {
	if(!is.numeric(x)) return(x)
	#if(length(list(...))) warning("additional arguments ignored")
	if(is.function(scale)) {
		scaleFunc <- scale
		scale <- TRUE
	} else scaleFunc <- function(x) sd(x, na.rm = TRUE)
	if(is.logical(scale)) scale <- if(scale) apply(x, 2L, scaleFunc) else 1
	if(is.logical(center)) center <- if(center) colMeans(x, na.rm = TRUE) else 0
	nc <- ncol(x)
	center <- rep(center, length.out = nc)
	scale <- rep(scale, length.out = nc)
	for(i in 1L:nc) x[, i] <- (x[, i] - center[i]) / scale[i]
	attr(x, "scaled:center") <- center
	attr(x, "scaled:scale") <- scale
	x
}

stdize.factor <-
function(x, binary = c("center", "scale", "binary", "half", "omit"),
center = TRUE, scale = FALSE, 
...) {
	#if(length(list(...))) warning("additional arguments ignored")
	if(nlevels(x) == 2L) {
		cl <- match.call()
		cl$x <- call("-", call("as.numeric", cl$x), 1)
		cl[[1L]] <- as.name("stdize.logical")
		eval(cl, envir = parent.frame(), enclos = asNamespace("MuMIn"))
	} else x
}

rootmeansq <- function(v) {
	v <- as.numeric(v[!is.na(v)])
	sqrt(sum(v^2) / max(1, length(v) - 1L))
}

stdize.logical <-
function(x, binary = c("center", "scale", "binary", "half", "omit"),
center = TRUE, scale = FALSE, 
...) {
	#if(length(list(...))) warning("additional arguments ignored")

	if(!missing(center) || !missing(scale)) {
		if(!missing(binary)) warning("argument 'binary' ignored when 'center' or 'scale' is given")
		stdize.numeric(as.numeric(x), center = center, scale = scale)
	} else {
		switch(match.arg(binary),
			   center = stdize.numeric(x, center = TRUE, scale = 1),
			   scale = stdize.numeric(x, center = TRUE, scale = TRUE),
			   half = stdize.numeric(x, center = 0.5, scale = 1),
			   binary = stdize.numeric(x, center = 0, scale = 1),
			   omit = x, NA)
	}
}

stdize.data.frame <-
function(x, binary = c("center", "scale", "binary", "half", "omit"),
	center = TRUE, scale = TRUE,
	omit.cols = NULL,
	source = NULL, prefix = TRUE, ...) {
	if(is.function(scale)) {
		scaleFunc <- scale
		scale <- TRUE
	} else scaleFunc <- function(x) sd(x, na.rm = TRUE)

	if(!is.null(source)) {
		if(!missing(center) || !missing(scale))
			warning("arguments 'center' and 'scale' ignored if 'source' is given")
	
		j <- match(colnames(x), attr(source, "orig.names"))
		if(any(is.na(j))) stop("some columns in 'x' are missing from 'source'")
		center <- attr(source, "scaled:center")[j]
		scale <- attr(source, "scaled:scale")[j]
		if(is.null(center) || is.null(scale)) stop("invalid 'source' object")
	}
	binary <- match.arg(binary)
    nc <- ncol(x)

	dataClasses <- vapply(x, function(x) {
		if (is.logical(x)) return("logical")
		if (is.factor(x)) if(nlevels(x) == 2L) return("factor2") else 
			return("other")
		if (is.matrix(x) && is.numeric(x)) return("nmatrix")
		if (is.numeric(x)) return("numeric")
		 return("other")
	}, "")
	
	if(is.character(omit.cols))
		dataClasses[colnames(x) %in% omit.cols] <- "omit"
	else if(is.numeric(omit.cols))
		dataClasses[omit.cols] <- "omit"
	
	numData <- dataClasses == "numeric"
	binaryData <- if(binary == "omit") FALSE else
		dataClasses == "factor2" | dataClasses == "logical"
	if(binary != "omit") for (i in which(binaryData)) x[, i] <- as.numeric(x[, i]) -
		if(dataClasses[i] == "factor2") 1 else 0
	
    if (is.logical(center)) {
		if(length(center) != 1 && length(center) != nc)
			stop("length of 'center' must equal one or the number of columns of 'x'")
		center <- rep(center, length.out = nc)
		binCenter <- switch(binary, center = TRUE, scale = TRUE, half = .5, 
			binary = 0, omit = FALSE, NA)
		jCenter <- (center & numData) | (isTRUE(binCenter) & binaryData)
		center[jCenter] <- colMeans(x[, jCenter, drop = FALSE], na.rm = TRUE)
		if(is.numeric(binCenter)) {
			jCenter[binaryData] <- TRUE
			center[binaryData] <- binCenter
		}
    } else if (is.numeric(center)) {
		if(length(center) != nc) 
			stop("length of 'center' must equal the number of columns of 'x'")
		jCenter <- !is.na(center) & (numData | binaryData)
	} else stop("invalid 'center' argument")
	if (is.logical(scale)) {
		if(length(scale) != 1 && length(scale) != nc)
			stop("length of 'scale' must equal one or the number of columns of 'x'")
		scale <- rep(scale, length.out = nc)
		binScale <- switch(binary, center = 1, scale = TRUE, half = 1, 
			binary = 1, omit = FALSE, NA)
		jScale <- (scale & numData) | (isTRUE(binScale) & binaryData)
		scale <- numeric(nc)
		for (i in which(jScale)) scale[i] <- scaleFunc(as.numeric(x[, i]))
		if(is.numeric(binScale)) {
			jScale[binaryData] <- TRUE
			scale[binaryData] <- binScale
		}
    } else if (is.numeric(scale)) {
		if(length(scale) != nc) 
			stop("length of 'scale' must equal the number of columns of 'x'")
		jScale <- !is.na(scale) & (numData | binaryData)
	} else stop("invalid 'scale' argument")

	jTransformed <- jScale | jCenter
	center[jTransformed & !jCenter] <- 0
	scale[jTransformed & !jScale] <- 1
	for (i in which(jTransformed)) x[, i] <- (x[, i] - center[i]) / scale[i]
	
	attr(x, "scaled:center") <- ifelse(jCenter, center, NA)
	attr(x, "scaled:scale") <- ifelse(jScale, scale, NA)
	attr(x, "orig.names") <- colnames(x)
	doprefix <-  FALSE
	if(is.character(prefix) ||
	   (doprefix <- (is.logical(prefix) && isTRUE(prefix)))) {
		prefix <- if(doprefix) c("z.", "c.") else rep(prefix, length.out = 2L)
		colnames(x)[jTransformed] <- 	
			paste0(prefix[jTransformed + (jCenter & !jScale)], 
				colnames(x)[jTransformed])
	}
	return(x)
}

stdize.formula <-
function(x, data = NULL, response = FALSE,
binary = c("center", "scale", "binary", "half", "omit"),
center = TRUE, scale = TRUE,
omit.cols = NULL,
prefix = TRUE, ...) {
	mf <- model.frame(x, data = data, drop.unused.levels = TRUE, ...)
	if(!is.null(omit.cols)) 
		omit.cols <- if(is.character(omit.cols)) 
			which(colnames(mf) %in% omit.cols) else
			stop("'omit.cols' must be a character vector")
	if(!response) omit.cols <- unique(c(omit.cols, 1L))
	attr(mf, "terms") <- NULL
	mf <- stdize.data.frame(mf, center = center, scale = scale, 
		omit.cols = omit.cols, binary = binary, prefix = prefix)
	mf
}
