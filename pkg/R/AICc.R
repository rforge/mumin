`AICc` <-
function(object, ..., k = 2, REML = NULL) {
	if(length(list(...))) {
		object <- list(object, ...)
		val <- as.data.frame(t(sapply(object, function(x) {
			z <- getAICc(x, k = k, REML = REML)
			c(AICc=z, unlist(attributes(z)))
			})))

		Call <- match.call()
		Call$k <- Call$REML <- NULL
		row.names(val) <- make.unique(as.character(Call[-1]))
		return(val)
	} else {
		return(getAICc(object, k = k, REML = REML))
	}
}

`getAICc` <-
function(object, k = 2, REML = NULL) {
	ll <- .getLogLik()

	if(me <- inherits(object, c("mer", "lme", "gls"))) {
		if (is.null(REML)) {
			if (inherits(object, c("lme", "gls"))) {
				REML <- object$method == "REML"
			} else {
				REML <- object@dims['REML'] != 0
			}
		}
		mLogLik <- ll(object, REML=REML)
	} else
		mLogLik <- ll(object)

	N <- nobs(object)
	mK <- attr(mLogLik, "df")
	mAIC <- -2 * c(mLogLik) + k * mK
	ret <- mAIC + 2 * mK * (mK + 1)/(N - mK - 1)
	attr(ret, "df") <- mK
	attr(ret, "AIC") <- mAIC
	if(me) attr(ret, "REML") <- REML
	return (ret)
}
