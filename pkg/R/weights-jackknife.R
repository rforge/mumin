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
#' @param object,\dots two or more fitted \code{\link{glm}} objects, or a
#'    \code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.
#' @param data a data frame containing the variables in the model. It is
#'    optional if all models are \code{glm}.
#' @param method a character string specifying the function to minimize. Either
#'    \code{"RMSE"} or \code{"likelihood"} (yet only the former is currently
#'    implemented).
#' @param weights an optional vector of \sQuote{\link[=weights.glm]{prior
#'    weights}} to be used in the fitting process. Should be \code{NULL} or
#'    a numeric vector.
#' @param optim.method optional, optimisation method, passed to \code{\link{optim}}.
#' @param maxit optional, the maximum number of iterations, passed to
#'    \code{\link{optim}}.
#' @param seed optionally, the random seed, see \code{\link{set.seed}}.
#' @param optim.args optional list of other arguments passed to \code{\link{optim}}.
#' @param force.update for \code{glm}, the \code{glm.fit} function is used for
#'    fitting models to the train data, which is much more efficient. Set to
#'    \code{TRUE} to use \code{update} instead.
#' @return The function returns a numeric vector of model weights.
#' @seealso \code{\link{Weights}}, \code{\link{model.avg}}
#' @author Kamil Barto\enc{ń}{n}. Carsten Dormann
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
		 weights = NULL,
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
	.weightedPred <- function(ww, pymat) {
		w <- c(1, exp(ww))
		w <- w / sum(w)
		pymat %*% w
	}
	
	.weightsFn <- switch(method,
		likelihood = {
			stop("method \"likelihood\" is not yet implemented")
			
			famstr <- tryCatch(vapply(models, function(m) {f <- family(m); c(f[["link"]],
					f[["family"]]) }, character(2L)), error = function(e) NULL)
				if(is.null(famstr)) stop("cannot get 'family' function for some models")
				if(!all(famstr[1L, ] == famstr[1L, 1L]) && all(famstr[2L, ] == famstr[2L, 1L]))
					stop("\"likelihood\" method requires all models use the same 'family' function")

			fam1 <- family(models[[1L]])
			n <- NA
			wt <- NA
			
			llik2 <-
				function(y, mu, fam, n, wt) {
					wt <- rep(wt, length.out = NROW(y)) # wt : fit$prior.weights
					z <- if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian")) 1 else 0
					(fam$aic(y, n, mu, wt, sum(fam$dev.resids(y, mu, wt))) / 2) + z # +LL
				}
			
			function(ww, pymat, y) {
				py <- .weightedPred(ww, pymat)
				llik2(y, py, fam1, n, wt = 1)
			}
		},
		RMSE = function(ww, pymat, y) {
			py <- .weightedPred(ww, pymat)
			return(sqrt(mean((py - y)^2)))
		})
	
	set.seed(seed)
	
	if(!force.update && all(vapply(models, inherits, FALSE, "glm"))) {
		
		if (isTRUE(getOption("debug.MuMIn")))
			message("using glm.fit")
		
		mf <- mergeMF(models, check = FALSE)
		tf <- terms(mf)
		if(!missing(data))
			mf <- model.frame(tf, data = data)
	
		X <- model.matrix(tf, mf)
		y <- as.matrix(get.response(mf))
		colnames(X) <- fixCoefNames(colnames(X))
		xil <- lapply(models, function(fit, allcn) {
			match(fixCoefNames(names(coef(fit))), allcn)
			}, colnames(X))
		no <- NROW(y)
		
		useWeightsArg <- !missing(weights)
		prwts <- if(is.null(weights)) rep(1, no) else rep(weights, length.out = no)

		pymat <- array(dim = c(no, m))
		
		xseq <- seq.int(no)
		for(j in seq.int(m)) {
			fit <- models[[j]]
			tf <- terms(fit)
			fam <- family(fit)
			off <- fit$offset
			
			if(!useWeightsArg) {
				prwts <- fit$prior.weights
				if(fam$family == "binomial") prwts <- prwts / rowSums(y)
			}
			if(is.null(off)) off <- rep(0, NROW(y))
				
			if (isTRUE(getOption("debug.MuMIn"))) {
				message("testing glm.fit #",  j)
				cf1 <- glm.fit(X[, xil[[j]], drop = FALSE], y, family = family(fit),
					weights = prwts, offset = off)$coefficients
				cf2 <- models[[j]]$coefficients
				names(cf2) <- fixCoefNames(names(cf2))
				stopifnot(all.equal(cf1[names(cf2)], cf2))
			}
				
			for(i in xseq) {
				coef1 <- glm.fit(X[-i, xil[[j]], drop = FALSE], y[-i, , drop = FALSE], family = family(fit),
					weights = prwts[-i], offset = off[-i])$coefficients
				pymat[i, j] <- predict_glm_fit(coef1, X[i, xil[[j]], drop = FALSE],
					offset = off[i], family = fam)[, 1L]
			}
		}
		
	} else {
		
		if(isTRUE(getOption("debug.MuMIn")))
			message("using update")
		
		if(any(!vapply(models, function(x) is.null(get_call(x)$offset), FALSE)))
			stop("use of \"offset\" argument in model calls is not allowed. ",
				 "Specify 'offset' term in the model formula instead")
			
		y <- get.response(models[[1L]], data)
		no <- NROW(y)
		pymat <- array(dim = c(no, m))
		
		prwts <- if(is.null(weights)) rep(1, no) else rep(weights, length.out = no)

		# TODO: check for offset argument and stop if found
		cll <- lapply(models, function(x) {
			update(x, evaluate = FALSE,
				data = `*tmp_dat*`[[1L]],
				weights = `*tmp_dat*`[[2L]] # TODO: make 'weights' optional
				)})

		pf <- parent.frame()
		`*tmp_dat*` <- NULL # confuse Rcheck
		on.exit(rm(list = "*tmp_dat*", envir = pf, inherits = FALSE))
		
		for(i in 1L:no) {
			assign("*tmp_dat*", list(data[-i, , drop = FALSE], prwts[-i]), pf)
			pymat[i, ] <- vapply(cll,
				function(x, newdata, envir) predict(eval(x, envir), newdata = newdata, type = "response"),
				newdata = data[i, , drop = FALSE],
				envir = pf,
				FUN.VALUE = numeric(1L))
		}
		
	}
	
	if(NCOL(y) == 2L) y <- y[, 1L] / rowSums(y) # binomial
		
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
	if (optres$convergence != 0) stop("not converged. 'optim' gave error code [", optres$convergence, "]")

	wts <- c(1, optres$par) / (sum(optres$par) + 1)
	
	structure(wts, name = "jackknife", class = c("model.weights", class(wts)))
}

	
#is.listOfCalls <-
#function(x)  is.list(x) && all(vapply(x, is.call, FALSE))
