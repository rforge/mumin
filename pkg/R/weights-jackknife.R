#' @title Jackknifed model weights
#' @rdname jackknifeWeights
#' @aliases jackknifeWeights
#' @keywords models
#' @encoding utf-8
#' @description
#' Computes model weights optimized for jackknifed model fits.
#' @details
#'     Model weights are chosen (using \code{\link{optim}}) to minimise RMSE of
#'     the prediction for data point \var{i}, of a model fitted omitting that
#'     data point \var{i}. The jackknife procedure is therefore run for all
#'     provided models and for all data points.
#' @author Carsten Dormann, Kamil Barto\enc{ń}{n}
#' @param object,\dots two or more fitted \code{\link{glm}} objects, or a
#'    \code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.
#' @param data a data frame in which to look for variables for use with
#'    \link[=predict]{prediction}.
#' @param method the function to minimize. Either \code{"RMSE"} or
#'    \code{"likelihood"} (yet only the former is currently implemented).
#' @param weights prior model weights.
#' @param optim.method optional, optimisation method, passed to \code{\link{optim}}.
#' @param maxit optional, the maximum number of iterations, passed to
#'    \code{\link{optim}}.
#' @param seed optionally, the random seed, see \code{\link{set.seed}}.
#' @param optim.args optional list of other arguments passed to \code{\link{optim}}.
#' @param force.update for \code{glm}, the \code{glm.fit} function is used for
#'    fitting models to the train data, which is much more efficient. Set to
#'    \code{TRUE} to use \code{update} instead.
#' @return The function returns a numeric vector of model weights.
#' @seealso \code{\link{Weights}}
#' @references
#' Hansen, B. E. & Racine, J. S. (2012) Jackknife model averaging.
#'     \emph{Journal of Econometrics}, 979: 38–46
#' @family model weights
#' @examples
#' fm <- glm(Prop ~ mortality * dose, binomial(), Beetle, na.action = na.fail)
#' 
#' fits <- lapply(dredge(fm, eval = FALSE), eval)
#' 
#' amJk <- amAICc <- model.avg(fits)
#' Weights(amJk) <- jackknifeWeights(fits, data = Beetle)
#' 
#' coef(amJk)
#' coef(amAICc)
#'
#' @export
jackknifeWeights <- 
function(object, ...,
		 data,
		 method = c("RMSE", "likelihood"),
		 weights = rep(1, n),
		 optim.method = "BFGS",
		 maxit = 1000,
		 seed = NULL,
		 optim.args = list(),
		 force.update = FALSE
		) {
	
	#if(is.listOfCalls(object)) {
	#	stop("'object' given as a list of calls cannot be handled yet")
	#} else
	models <- getModelArgs()
	m <- length(models)

	if(m < 2) stop("need more than one model")

	.checkModels(models)
	
	method <- match.arg(method)
	
	# compute RMSE/likelihood for a value of w, given J:
	.weightsFn <- switch(method,
		likelihood = stop("method \"likelihood\" is not yet implemented"),
		RMSE = function(ww, pymat, y) {
			w <- c(1, exp(ww))
			w <- w / sum(w)
			py <- pymat %*% w
			return(sqrt(mean((py - y)^2)))
		})
	
	set.seed(seed)

	n <- nrow(data)
	pymat <- array(dim = c(n, m))
	
	if(!force.update && all(vapply(models, inherits, FALSE, "glm"))) {

		for(j in seq.int(m)) {
			fit <- models[[j]]
			tf <- terms(fit)
			fam <- family(fit)
			off <- fit$offset
			for(i in 1L:n) {
				coef1 <- do_glm_fit(tf, data[-i, , drop = FALSE], fam,
					weights[-i], off[-i])$coefficients
				pymat[i, j] <- predict_glm_fit(coef1,
					model.matrix(tf, data[i, , drop = FALSE]), offset = off[i],
					family = fam)[, 1L]
			}
		}
		
	} else {
		if(any(!vapply(models, function(x) is.null(get_call(x)$offset), FALSE)))
			stop("use of \"offset\" argument in model calls is not allowed. ",
				 "Specify 'offset()' term in the formula instead")

	  
		# TODO: check for offset argument and stop if found
		cll <- lapply(models, function(x) {
			update(x, evaluate = FALSE,
				data = `*tmp_dat*`[[1L]],
				weights = `*tmp_dat*`[[2L]] # TODO: make 'weights' optional
				)})

		pf <- parent.frame()
		`*tmp_dat*` <- NULL # confuse Rcheck
		on.exit(rm(list = "*tmp_dat*", envir = pf, inherits = FALSE))
		
		for(i in 1L:n) {
			assign("*tmp_dat*", list(data[-i, , drop = FALSE], weights[-i]), pf)
			pymat[i, ] <- vapply(cll,
				function(x, newdata, envir) predict(eval(x, envir), newdata = newdata, type = "response"),
				newdata = data[i, , drop = FALSE],
				envir = pf,
				FUN.VALUE = numeric(1L))
		}
	}
	
	y <- get.response(models[[1L]], data)
	if(is.matrix(y)) y <- y[, 1L] / rowSums(y) # binomial
		
	if("control" %in% names(optim.args)) {
		optim.control <- optim.args$control
		if(is.numeric(optim.args$control$maxit) && missing(maxit))
			maxit <- optim.args$control$maxit
		optim.args$control <- NULL
	} else optim.control <- list()
	optim.control$maxit <- maxit
	
	#.weightsFn(runif(ncol(pymat) - 1L), pymat = pymat, y = y)
	
	optres <- eval(as.call(c(alist(optim, .weightsFn, pymat = pymat,
			y = y,
			par = runif(ncol(pymat) - 1L),
			method = optim.method, control = optim.control),
		optim.args)))
	if (optres$convergence != 0) stop("not converged. 'optim' error code [", optres$convergence, "]")

	wts <- c(1, optres$par) / (sum(optres$par) + 1)
	
	structure(wts, name = "jackknife", class = c("model.weights", class(wts)))
}

	
is.listOfCalls <-
function(x)  is.list(x) && all(vapply(x, is.call, FALSE))
