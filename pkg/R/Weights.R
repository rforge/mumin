# Calculate Akaike weights
`Weights` <-
function(x)  UseMethod("Weights")

`Weights.model.selection` <-
function(x) structure(x[, "weight"] / sum(x[, "weight"]),
	names = rownames(x))

`Weights.averaging` <-
function(x) {
	x$msTable[, ncol(x$msTable)]
}

`Weights.data.frame` <-
function(x) {
	if(ncol(x) == 2L && colnames(x)[1L] == "df"	&& is.numeric(x[, 2L]))
		return(Weights.default(x[, 2L]))
	if(ncol(x) == 1L && is.numeric(x[, 1L]))
		return(Weights.default(x[, 1L]))
	return(NA)
}

`Weights.default` <-
function(x) {
	d <- exp(-x / 2)
	d / sum(d)
}