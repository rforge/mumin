#' @title Leave-one-out cross-validation 
#' @aliases loo loo.default loo.lm
#' @rdname loo
#' @keywords models
#' @encoding utf-8
#' @description
#' Computes the RMSE/log-likelihood based on leave-one-out cross-validation.
#' @param object a fitted object model, currently only \code{lm}/\code{glm} is
#' accepted.
#' @param method the criterion to use, given as a character string,
#' either \code{"rmse"} for Root-Mean-Square Error or \code{"loglik"}
#' for log-likelihood.
#' @param ... other arguments are currently ignored.
#' @references
#' Dormann, C. et al. (\emph{in prep.}) Model averaging in ecology: evidence,
#' approach examples.
#' @author Kamil Barto\enc{ń}{n}, based on code by Carsten Dormann
#' @details
#' Leave-one-out cross validation is \var{K}-fold cross validation, with \var{K}
#' equal to the number of data points in the set \var{N}. For all \var{i} from 1 
#' to \var{N}, the model is fitted to all the data except for point \var{i} and a 
#' prediction is made for that point. The average error is computed and used to 
#' evaluate the model. 
#' @examples
#' fm <- lm(y ~ X1 + X2 + X3 + X4, Cement)
#' loo(fm, method = "l")
#' loo(fm, method = "r")
#'
#' ## Compare LOO_RMSE and AIC/c
#' options(na.action = na.fail)
#' dd <- dredge(fm, rank = loo, extra = list(AIC, AICc), method = "rmse")
#' plot(loo ~ AIC, dd, ylab = expression(LOO[RMSE]), xlab = "AIC/c")
#' points(loo ~ AICc, data = dd, pch = 19)
#' legend("topleft", legend = c("AIC", "AICc"), pch = c(1, 19))
#' @return \code{loo} returns a single numeric value of RMSE or mean log-likelihood.
#'
#' @export 
loo <-
function(object, method = c("loglik", "rmse"), ...)
	UseMethod("loo")

#' @method loo default
#' @export
loo.default <-
function(object, method = c("loglik", "rmse"), ...)
	.NotYetImplemented()
	# for other types of models use manipulated call - will be SLOW


#  #' @rdname loo
#  #' @method loo lm
#  #' @param start,etastart,mustart,control,intercept optional arguments to be
#  #' passed to \code{\link{glm.fit}}, but are currently ignored, with a warning.
#' @export
loo.lm <- # lm, glm, mgcv::gam
function(object, method = c("loglik", "rmse"), start, etastart, mustart,
		 control, intercept, ...) {
	
	if(!inherits(object, "lm"))
		stop("'object' must be a \"glm\" or \"lm\" object")
	
	## TODO: pass other arguments: 
	if(!missing(start) || !missing(etastart)  || !missing(mustart)  ||
	   !missing(control) || !missing(intercept)) {
		warning("arguments 'start', 'etastart', 'mustart', 'control', 'intercept' are ignored")
	}
	## binary response: object$y is always a vector and weights(object) give size

	method <- match.arg(method)
	tt <- terms(object)
	intercept <- attr(tt, "intercept") # == 0 or 1
	beta <- coef(object)

	fam <- family(object)
	x <- model.matrix(object)
	if(inherits(object, "glm")) {
		y <- object$y
	} else {
		y <- model.frame(object)[, asChar(attr(tt, "variables")[-1L][[attr(tt, "response")]])]
		if(is.matrix(y)) y <- y[, 1L] / rowSums(y) # binomial
	}
	n <- length(y)

	wt <- weights(object, "prior")
	if(is.null(wt)) wt <- rep(1, n)
	
	offset <- object$offset
	if(is.null(offset)) offset <- numeric(n)
	
	dev <- deviance(object)
	aicfun <- fam$aic
	aic_n <- if(fam$family == "binomial") wt else rep(1, n)
	yt <- fam$linkfun(y)

	func <- switch(method,
		loglik = function(fit, i) {
			py1 <- predict_glm_fit(fit$coefficients, x[i, , drop = FALSE],
				offset[i], fam)[, 1L]
			aicfun(y[i], aic_n[i], py1, wt[i], deviance(fit)) / 2 # == logLik
		},
		rmse = function(fit, i) { # alternatively: MSE on transformed[?] data
			py1 <- predict_glm_fit(fit$coefficients, x[i, , drop = FALSE],
				offset[i])[, 1L]
				# prediction on link scale
			yt[i] - py1 # inefficient to '^2' here
		})

	rval <- numeric(n)
	for (i in seq.int(n)) {
		fm1 <- glm.fit(y = y[-i], x = x[-i, , drop = FALSE],
				family = fam, offset = offset[-i], weights = wt[-i])
		rval[i] <- func(fm1, i)
	}

	if(method == "rmse") {
		sqrt(mean(rval^2))
	} else mean(rval) # XXX: shoult it be a mean of logliks ?
}
