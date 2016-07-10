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
	# for other types of models use manipulated call - SLOW

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
	beta <- coef(object)
	fam <- family(object)
	X <- model.matrix(object)
	y0 <- model.frame(object)[, asChar(attr(tt, "variables")[-1L][[attr(tt, "response")]])]
	nobs <- NROW(y0)
	wt <- weights(object, "prior")
	if(is.null(wt)) wt <- rep(1, nobs)
	
	if(NCOL(y0) == 2L) { # binomial
		n <- rowSums(y0)
		y <- y0[, 1L] / n
		wt0 <- wt / n
	} else  {
		n <- rep(1, nobs)
		y <- y0
		y0 <- as.matrix(y0)
		wt0 <- wt
	}


	offset <- object$offset
	if(is.null(offset)) offset <- numeric(nobs)
	
	aic <- fam$aic
	eta <- fam$linkfun(y)
	
	dev <- function(y, mu, wt, fam)  sum(fam$dev.resids(y, mu, wt))
	llik <- function(y, X, beta, fam, n, wt = 1, off = NULL) {
		# wt : fit$prior.weights
		no <- NROW(y)
		wt <- rep(wt, length.out = no)
		mu <- predict_glm_fit(beta, X, off, fam)[, 1L]
		z <- if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian")) 1 else 0
		(fam$aic(y, n, mu, wt, dev(y, mu, wt, fam)) / 2) + z # should be negative?
	}
	
	func <- switch(method,
	loglik = function(fit, i) {
		llik(y[i], X[i, , drop = FALSE], fit$coefficients, fit$family, n[i], wt[i], offset[i])
	},
	rmse = function(fit, i) { # alternatively: MSE on transformed[?] data
		py <- predict_glm_fit(fit$coefficients, X[i, , drop = FALSE],
			offset[i])[, 1L]
			# prediction on link scale
		eta[i] - py # inefficient to '^2' here
	})
	
	
	if (isTRUE(getOption("debug.MuMIn"))) {
		message("running test 1...")
		# XXX: DEBUG test
		testLL1 <- llik(y, X, object$coefficients, object$family, n, wt, offset)
		#print(testLL1)
		#print(logLik(object))
		#print(testLL1 - logLik(object))
		#print(rbind(n, wt, offset))
		stopifnot(all.equal(testLL1, c(logLik(object)), tolerance = 1e-5))
		message("OK")
		message("running test 2...")
		testFm <- glm.fit(y = y0, x = X, family = fam, offset = offset, weights = wt0)
		#print(rbind(testFm$coefficients, object$coefficients))
		stopifnot(all.equal(testFm$coefficients, object$coefficients))
		message("OK")
		message("running test 3...")
		testLL2 <- llik(y, X, testFm$coefficients, testFm$family, n, wt, offset)
		#print(c(testLL2,  logLik(object)))
		stopifnot(all.equal(testLL2, c(logLik(object)), tolerance = 1e-5))
		message("OK")
		#print(testLL2)
		#print(logLik(object))
	}

	rval <- numeric(nobs)
	for (i in seq.int(nobs)) {
		fm1 <- glm.fit(y = y0[-i, , drop = FALSE], x = X[-i, , drop = FALSE],
				family = fam, offset = offset[-i], weights = wt[-i])
		rval[i] <- func(fm1, i)
	}
	if(method == "rmse") {
		sqrt(mean(rval^2))
	} else mean(rval) # XXX: shoult it be a mean of logliks ?
}
