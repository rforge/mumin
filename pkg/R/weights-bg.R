#' @title Bates-Granger model weights
#' @rdname BGWeights
#' @aliases BGWeights
#' @keywords models
#' @encoding utf-8
#' @description
#' Computes empirical weights based on out of sample forecast variances,
#'    following Bates and Granger (1969).
#' @details
#' Bates-Granger model weights are calculated using prediction covariance. To
#' get the estimate of prediction covariance, the models are fitted to
#' randomly selected half of \code{data} and prediction is done on the
#' remaining half.
#' These predictions are then used to compute the variance-covariance between
#' models, \eqn{\Sigma}. Model weights are then calculated as
#' \eqn{w_{BG} = (1'\Sigma^{-1}1)^{-1} 1\Sigma^{-1} },
#' where \eqn{1} a vector of 1-s.
#' 
#' Bates-Granger model weights may be outside of the \eqn{[0,1]} range, which
#' may cause the averaged variances to be negative. Apparently this method
#' works best when data is large.
#' 
#' @note For matrix inversion, \code{\link[MASS]{ginv}} from package
#' \pkg{MASS} is more stable near singularities than \code{\link{solve}}. It
#' will be used as a fallback if \code{solve} fails and \pkg{MASS} is
#' available. 
#' 
#' @author Carsten Dormann, Kamil Barto\enc{ń}{n}
#' 
#' @param object,\dots two or more fitted \code{\link{glm}} objects, or a
#' \code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.
#' @param data a data frame containing the variables in the model.
#' @param seed optionally, the random seed, see \code{\link{set.seed}} .
#' @param force.update if \code{TRUE}, the much less efficient method of
#'    updating \code{glm} function will be  used rather than directly \emph{via}
#'    \code{\link{glm.fit}}. This only applies to \code{glm}s, in 
#' case of other model types \code{update} is always used.
#' @return The function returns a numeric vector of model weights.
#' @seealso \code{\link{Weights}}, \code{\link{model.avg}}
#' @references
#' Bates, J. M. & Granger, C. W. J. (1969) The combination of forecasts.
#'    \emph{Journal of the Operational Research Society}, 20: 451-468. 
#' @family model weights
#' @examples
#' fm <- glm(Prop ~ mortality + dose, family = binomial, Beetle, na.action = na.fail)
#' models <- lapply(dredge(fm, evaluate = FALSE), eval)
#' ma <- model.avg(models)
#' 
#' # this produces warnings because of negative variances:
#' Weights(ma) <- BGWeights(ma, data = Beetle , seed = 10)
#' coefTable(ma, full = TRUE)
#' 
#' # SE for prediction is not reliable if some or none of coefficient's SE
#' # are available
#' predict(ma, data = test.data, se.fit = TRUE)
#' coefTable(ma, full = TRUE)
#' 
#' @export
BGWeights <-
function(object, ..., data, seed = NULL, force.update = FALSE) {
     
	models <- getModelArgs()
    m <- length(models)
	if(m < 2) stop("need more than one model")
	.checkModels(models)

    set.seed(seed)

    n <- nrow(data)
    k <- sample.int(n, floor(n / 2))
    dat_train <- data[k, ]
	dat_test <- data[-k, ]
	# TODO: allow user to specify offset and weights  
    offset <- rep(0, n) 
    weights <- rep(1, n)
    weights_train <- weights[k]
    weights_test <- weights[-k]
    offset_train <- offset[k]
    offset_test <-  offset[-k]
    
    if(!force.update && all(vapply(models, inherits, FALSE, "glm"))) { # XXX: what about lm
        # XXX: here 'offset_train' DOES include offset specified in 'formula'
        
		py_test <- array(dim = c(nrow(dat_test), m))
        for(i in seq.int(length(models))) {
            fit <- models[[i]]
            tf <- terms(fit)
            fit_train <- update_glm_fit(fit, dat_train, weights_train, offset_train)
            py_test[, i] <- predict_glm_fit(fit_train$coefficients, model.matrix(tf, dat_test),
                offset = offset_test, family = family(fit))[, 1L]
        }
       
    } else { # for non-glm models, use update (2x slower)
		
		if(any(!vapply(models, function(x) is.null(get_call(x)$offset), FALSE)))
			stop("use of \"offset\" argument in model calls is not allowed. ",
				 "Specify 'offset()' term in the formula instead")
			
		
        # needed for `update` to work, otherwise it cannot find 'weights' and 'offset'
        pf <- parent.frame()
		
		`*tmp_dat*` <- NULL # to confuse Rcheck
        assign("*tmp_dat*", list(dat_train, weights_train), pf)
        on.exit(rm(list = "*tmp_dat*", envir = pf, inherits = FALSE))
        
        # XXX: here 'offset_train' DOES NOT include offset specified in 'formula'
        py_test <- array(dim = c(nrow(dat_test), m))
        for(i in seq.int(length(models)))
            py_test[, i] <- predict(update(
                models[[i]], data = `*tmp_dat*`[[1L]], weights = `*tmp_dat*`[[2L]],
				# XXX: weird behaviour of predict when offset= is given.
				#      prediction always has length of the offset. a bug in predict.glm?
				#offset = `*tmp_dat*`[[3L]]),
                 newdata = dat_test,
				 type = "response"))
         
        # NOTE: No speed gain with *apply - more memory needed to store train_fits
        ## train_fits <- lapply(models, function(x) update(x, data = `*tmp_dat`[[1L]], weights = `*tmp_dat`[[2L]], offset = `*tmp_dat`[[3L]]))
        ## py_test <- vapply(train_fits, predict, numeric(nrow(dat_test)), newdata = dat_test, type = "response")
    }
    
    # XXX: this assumes all models share the same response, do some checking
    y_test <- get.response(models[[1L]], data = dat_test)
    if(is.matrix(y_test)) y_test <- y_test[, 1L] / rowSums(y_test) # binomial
    

	Sigma <- cov(y_test - py_test)
	ones <- rep(1, m)
    # XXX: I want do avoid dependency on MASS
    #ginv <- if(use.MASS) getFrom("MASS", "ginv") else solve
    
   	ones <- rep(1, m)
	
	fn1 <- function(ones, Sigma, ginv) ginv(t(ones) %*% ginv(Sigma) %*% ones) %*% ones %*% ginv(Sigma)
	

    rval <- tryCatch(fn1(ones, Sigma, solve), error = function(e) {
		if(length(find.package("MASS", quiet = TRUE)) == 1L)
			fn1(ones, Sigma, getFrom("MASS", "ginv")) else
			stop(e)
		})[1L, ]
	
	structure(rval, name = "Bates-Granger", class = c("model.weights", class(rval)))
}
