##=============================================================================
## Classes: gee & geeglm
##=============================================================================


`coefTable.gee` <-
`coefTable.geeglm` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$coefficients
	type <- match.arg(type)
	j <- if(type == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, 2L], coefNames = rownames(cf))
}


`coefTable.geese` <-
function(model, ..., type = c("naive", "robust")) {
	cf <- summary(model, ...)$mean
	type <- match.arg(type)
	j <- if(type == "naive") 2L else 4L
	.makeCoefTable(cf[, 1L], cf[, 2L], coefNames = rownames(cf))
}

`coef.geese` <- 
function (object, ...) object$beta

##=============================================================================
## Class: yags
##=============================================================================


`coefTable.yagsResult` <-
function(model, ..., type = c("naive", "robust")) {
	type <- match.arg(type)
		ql <- logLik(if(type == "I") indep else object)
	ql <- if(family(object)$family == "gaussian") {
}






