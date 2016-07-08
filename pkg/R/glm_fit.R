

# helper function: prediction from matrix, coefficients and inverse link
predict_glm_fit <-
function(beta, x, offset, family = NULL) {
	if(is.null(offset)) offset <- 0
	if(inherits(family, "family")) return(family$linkinv(offset + (x %*% beta)))
	return(offset + (x %*% beta))
}


update_glm_fit <-
function(fit, data, weights, offset, nobs = nrow(data), y = NULL) {
    tf <- terms(fit)
    x <- model.matrix(tf, data = data)
    if(is.null(y)) {
        y <- get.response(tf, data)
        if(is.matrix(y)) {
            wts <- rowSums(y)
            y <- y[, 1L] / wts
        } else wts <- rep(1, nobs)
        weights <- weights * wts
    }
    glm.fit(x, y, weights, offset = offset, family = family(fit))
}


do_glm_fit <-
function(tf, data, family, weights, offset, nobs = nrow(data), y = NULL) {
	x <- model.matrix(tf, data = data)
    if(is.null(y)) {
		y <- get.response(tf, data)
        if(is.matrix(y)) {
            wts <- rowSums(y)
            y <- y[, 1L] / wts
        } else wts <- rep(1, nobs)
        weights <- weights * wts
    }
    glm.fit(x, y, weights, offset = offset, family = family)
}

aicloglik_glm_fit <-
function(object, y, x, wt, offset = NULL) {
    fam <- object$family
    nobs <- NROW(x)
    n <- if (NCOL(y) == 1) 
        rep.int(1, nobs) else rowSums(y)
    #mu <- fam$linkinv((x %*% object$coefficients)[, 1L])
	mu <- predict_glm_fit(object$coefficients, x, offset, fam)[, 1L]
    dev <- sum(fam$dev.resids(y, mu, wt))
    aic <- fam$aic(y, n, mu, wt, dev) + 2 * object$rank
    p <- object$rank
    if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian")) 
        p <- p + 1
    ll <- p - aic/2
    # c(aic = aic, loglik = ll, nobs = nobs, df = p)
    c(aic, ll, p)
}

# list(coefficients =, family =, rank=)