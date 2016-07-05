# Calculate Akaike weights
`Weights` <-
function(x)  UseMethod("Weights")

`Weights.model.selection` <-
function(x) {
	i <- type2col(x, "weight")
	structure(item(x, i) / sum(item(x, i)),	names = row.names(x))
}

`Weights.averaging` <-
function(x) {
	x$msTable[, ncol(x$msTable)]
}

`Weights.data.frame` <-
function(x) {
	if(ncol(x) == 2L && colnames(x)[1L] == "df"	&& is.numeric(x[, 2L]))
		return(Weights(x[, 2L]))
	if(ncol(x) == 1L && is.numeric(x[, 1L]))
		return(Weights(x[, 1L]))
	return(NA)
}

`Weights.numeric` <-
function(x) {
	x <- x - min(x)
	d <- exp(-x / 2)
	d / sum(d)
}

`Weights.default` <-
function(x) {
    cry(, "cannot use \"%s\" as 'x'", class(x)[1L])
}

`Weights<-` <-
function(x, value)  UseMethod("Weights<-")


`Weights<-.default` <-
function(x, value) {
	stop("can assign weights only to an \"averaging\" object")
}

`Weights<-.averaging` <-
function(x, value) {

	wi <- ncol(x$msTable)
	if(is.null(value)) {
		wts <- Weights(x$msTable[, wi - 1L])
		x$msTable[, wi] <- wts
		colnames(x$msTable)[wi] <- "weight"
	} else {
		x$msTable[, wi] <- value
		wts <- x$msTable[, wi]
		wts <- wts / sum(wts)
		x$msTable[, wi] <- wts
		colnames(x$msTable)[wi] <- "[weight]"
	}

	rv <-  attr(x, "revised.var")
	for(i in 1L:nrow(x$coefficients)) {
		  full <- rownames(x$coefficients)[i] == "full"
		  x$coefficients[i, ] <- .coefarr.avg(x$coefArray, wts, full = full, alpha = 0.05, revised.var = rv)[, 1L]
	}
	
	o <- order(wts, decreasing = TRUE)
	x$msTable <- x$msTable[o, ]
	x$coefArray <- x$coefArray[o,,]
	if(!is.null(attr(x, "modelList"))) attr(x, "modelList") <- attr(x, "modelList")[o]
	x
}


