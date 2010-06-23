`Weights.default` <-
function(ic) {
	delta <- ic - min(ic)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	return (weight)
}

`Weights` <-
function(ic, ...) {
	delta <- ic - min(ic)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	return (weight)
}
