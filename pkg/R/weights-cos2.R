#' @title Cos-squared model weights
#' @rdname cos2Weights
#' @aliases cos2Weights
#' @keywords models
#' @encoding utf-8
#' @description Calculates cos-squared model weights, following the algorithm
#'     outlined in the appendix of Garthwaite & Mubwandarikwa (2010).
#' @param object,\dots two or more fitted \code{\link{glm}} objects, or a
#' \code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.
#'     Currently only \code{lm} and \code{glm} objects are accepted.
#' @param data a test data frame in which to look for variables
#' for use with \link[=predict]{prediction}. If omitted, the fitted linear
#'     predictors are used.
#' @param eps tolerance for determining convergence.
#' @param maxit maximum number of iterations.
#' @param predict.args optionally, a \code{list} of additional arguments to be
#'     passed to \code{predict}.
#' @return The function returns a numeric vector of model weights.
#' @seealso \code{\link{Weights}}, \code{\link{model.avg}}
#' @family model weights
#' @references
#' Garthwaite, P. H. and Mubwandarikwa, E. (2010) Selection of weights for
#'     weighted model averaging. \emph{Australian & New Zealand Journal of
#'     Statistics}, 52: 363–382.
#'
#' Dormann, C. et al. (\emph{in prep.}) Model averaging in ecology:
#'     evidence, approach examples.
#' @author Carsten Dormann, adapted by Kamil Barto\enc{ń}{n}
#' @examples
#' \dontshow{
#' if(length(find.package("expm", quiet = TRUE)) == 1) \{
#' }
#' fm <- lm(y ~ X1 + X2 + X3 + X4, Cement, na.action = na.fail)
#' # most efficient way to produce a list of all-subsets models
#' models <- lapply(dredge(fm, evaluate = FALSE), eval)
#' ma <- model.avg(models)
#' 
#' test.data <- Cement
#' Weights(ma) <- cos2Weights(models, data = test.data)
#' predict(ma, data = test.data)
#' \dontshow{
#' \} else message("Package 'expm' is needed to run this example")
#' }
#' @export
cos2Weights <- 
function(object, ..., data, eps = 1E-6, maxit = 100, predict.args = list()) {
	models <- getModelArgs()
	n <- length(models)
	if(n < 2) stop("need more than one model")
 
    if(!all(vapply(models, inherits, TRUE, "lm")))
	   stop("'models' must inherit from \"lm\" class")
  
  	#py <- sapply(models, predict, newdata = data, type = "response", ...)
	cl <- as.call(c(as.name("predict"), alist(models[[i]], newdata = data,
		type = "response"), predict.args))
	
	# TODO: glm.fit version
	
	i <- 1L
	py1 <- eval(cl)
	if(is.array(py1)) {
		stop(">1-dimensional predictions cannot be handled yet")
		# assuming prediction is a matrix: (add some checking for it)
		py <- array(dim = c(dim(py1), n))
		py[, , 1L] <- py1
		for(i in seq.int(2L, n)) py[, , i] <- eval(cl)
	} else {
		py <- array(dim = c(length(py1), n))
		py[, 1L] <- py1
		for(i in seq.int(2L, n)) py[, i] <- eval(cl)
	}
	
	# if one model is constant:
	if (any(g <- apply(py, 2L, "sd") == 0))
		py[, g] <- py[, g] + rnorm(NROW(py))
		
	sqrtm <- getFrom("expm", "sqrtm")
	
	R <- cor(py)
	nR <- NCOL(R)
	D1 <- diag(rep(2, nR))
	D2 <- diag(nR)
	counter <- 0L
	while (any(abs(diag(D1) - diag(D2)) > eps)) {
		ED <- eigen(D1 %*% R %*% D1)
		Q <- ED$vectors
		Lambda <- diag(ED$values)
		## test:
		#Q %*% Lambda %*% solve(Q) # fine
		Lambda12 <- sqrtm(Lambda)
		E <- solve(D1) %*% Q %*% Lambda12 %*% solve(Q)
		D2 <- D1
		D1 <- diag(diag(Re(E)))
		counter <- counter + 1L
		if (counter >= maxit) {
		  warning("maximum number of iterations reached without convergence")
		  break
		}
	}
	wts <- diag(D2)^2 / sum(diag(D2)^2)
	
	structure(wts, name = "cos-squared", class = c("model.weights", class(wts)))
}
