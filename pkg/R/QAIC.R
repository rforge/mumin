`QAIC` <-
function(object, ..., chat) {
	#chat <- summary(gm)$dispersion

	`getQAIC` <- function(model, chat) {
		if (any(inherits(model,  c("lmer", "glmer")))) {
			mLogLik <- logLik(model, model@status["REML"])
			N <- NROW(model@frame)
		} else {
			mLogLik <- logLik(model)
			N <- length(resid(model))
		}

		k <- attr(mLogLik, "df") + 1
		ret <- (deviance(model) / chat) + 2 * k
		return (ret)
	}

	if(length(list(...))) {
		object <- list(object, ...)
		val <- data.frame(QAIC=sapply(object, getQAIC, chat = chat))
		   Call <- match.call()
		   Call$chat <- NULL
		row.names(val) <- as.character(Call[-1])
		return(val)
	} else {
		return(getQAIC(object, chat = chat))
	}
}
