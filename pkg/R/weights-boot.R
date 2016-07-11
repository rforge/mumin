#' @title Bootstrap model weights
#' @rdname bootWeights
#' @aliases bootWeights
#' @keywords models
#' @encoding utf-8
#' @description
#' Computes model weights using bootstrap.
#' @details
#' The models are fitted repeatedly to a resampled data set and ranked
#' using AIC-type criterion. The model weights represent the proportion of
#' replicates when a model has the lowest IC value.
#' @author Kamil Barto\enc{ń}{n}, Carsten Dormann 
#' @param object,\dots two or more fitted \code{\link{glm}} objects, or a
#'    \code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.
#' @param R the number of replicates.
#' @param rank a character string, specifying the information criterion to use
#'     for model ranking. Defaults to \code{\link{AICc}}.
#' @param seed optionally, the random seed, see \code{\link{set.seed}}.
#' @return The function returns a numeric vector of model weights.
#' @seealso \code{\link{Weights}}, \code{\link{model.avg}}
#' @family model weights
#' @examples
#' # To speed up the bootstrap, use 'x = TRUE' so that model matrix is included
#' #     in the returned object
#'fm <- glm(Prop ~ mortality + dose, family = binomial, data = Beetle, 
#'    na.action = na.fail, x = TRUE)
#'
#'fml <- lapply(dredge(fm, eval = FALSE), eval)
#'am <- model.avg(fml)
#'
#'Weights(am) <- bootWeights(am, data = Beetle, R = 25)
#'
#'summary(am)
#'
#' @export

##~ environment(bootWeights) <- asNamespace("MuMIn") 
##~ fm <- glm(formula = Prop ~ mortality + dose, family = binomial, data = Beetle, 
##~     na.action = na.fail, x = FALSE)
##~ fml <- lapply(dredge(fm, eval = FALSE), eval)
##~ am <- model.avg(fml)
##~ 
##~ fm <- glm(formula = Prop ~ mortality + dose, family = binomial, data = Beetle, 
##~     na.action = na.fail, x = TRUE)
##~ fml <- lapply(dredge(fm, eval = FALSE), eval)
##~ amx <- model.avg(fml)
##~ 
##~ system.time(bootWeights(amx, R = 250))
##~ system.time(bootWeights(am, R = 250))
##~ 
##~ #glm.fit(model.matrix(fm), get.response(fm), fit$

bootWeights <- 
function(object, ...,
		 R,
		 rank = c("AICc", "AIC", "BIC"),
		 seed = NULL
		) {
	
	models <- getModelArgs()
	m <- length(models)

	if(m < 2) stop("need more than one model")

	.checkModels(models)
	
	for(fm in models) {
	  if(is.na(match("x", names(fm)))) {
		warning("for efficiency of the bootstrap procedure, 'glm' should be called with 'x = TRUE'")
		break
	  }
	}
	
	rank <- match.arg(rank)
	
	ic <- switch(rank, AICc = AICc, AIC = AIC, BIC = BIC)
	
	set.seed(seed)

	mseq <- seq.int(m)
	
	if(all(vapply(models, inherits, FALSE, "glm"))) { ## !force.update &&
	  	best <- integer(R)
		ics <- numeric(m)
		n <- nobs(models[[1L]]) # assuming nobs is the same across models
		r <- 1L
		while(r <= R) {
			g <- sample.int(n, replace = TRUE)

			for(j in mseq) {
				fit <- models[[j]]
				fit1 <- glm.fit(model.matrix(fit)[g, , drop = FALSE], fit$y[g], fit$prior.weights[g],
				  offset = fit$offset[g], family = fit$family)
				ics[j] <- ic(loglik_glm_fit(fit1))
			}
			best[r] <- which.min(ics)
			r <- r + 1L
		}
	} else stop("all model objects must be of \"glm\" class")
  
	wts <- tabulate(best, m)
	wts <- wts  / sum(wts)
	structure(wts, name = "bootstrap", class = c("model.weights", class(wts)))
}
