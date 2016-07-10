#' @title Stacking model weights
#' @rdname stackingWeights
#' @keywords models
#' @encoding utf-8
#' @param object,\dots two or more fitted \code{\link{glm}} objects, or a
#' \code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.
#' @param data a data frame containing the variables in the model, used for
#'    fitting and prediction.
#' @param R the number of replicates.
#' @param p the proportion of the \code{data} to be used as training set.
#'      Defaults to 0.5.
#' @param seed optionally, the random seed. See \code{\link{set.seed}}.
#' @return \code{stackingWeights} returns a matrix with two rows, holding model weights
#'    calculated using \code{mean} and \code{median}.
#' @description Computes model weights based on a cross-validation-like procedure.
#' @details
#' Each model in a set is fitted to the training data: a subset of \code{p * N}
#' observations in \code{data}. From these models a prediction is produced on
#' the remaining part of \code{data} (the test
#' or hold-out data). These hold-out predictions are fitted to the hold-out
#' observations, by optimising the weights by which the models are combined. This
#' process is repeated \code{R} times, yielding a distribution of weights for each
#' model (which Smyth & Wolpert (1998) referred to as an \sQuote{empirical Bayesian
#' estimate of posterior model probability}). A mean or median of model weights for
#' each model is taken and re-scaled to sum to one.
#' @seealso \code{\link{Weights}}, \code{\link{model.avg}}
#' @family model weights
#' @note
#' This approach requires a sample size of at least \eqn{2\times}{2x} the number of models. 
#' @references
#' Wolpert, D. H. (1992) Stacked generalization. \emph{Neural Networks}, 5: 241-259.
#'
#' Smyth, P. & Wolpert, D. (1998) \emph{An Evaluation of Linearly Combining
#'    Density Estimators via Stacking. Technical Report No. 98-25.} Information
#'    and Computer Science Department, University of California, Irvine, CA.
#' @author Carsten Dormann, Kamil Barto\enc{ń}{n}
#' @examples
#' # global model fitted to training data:
#' fm <- glm(y ~ X1 + X2 + X3 + X4, data = Cement, na.action = na.fail)
#' # generate a list of *some* subsets of the global model
#' models <- lapply(dredge(fm, evaluate = FALSE, fixed = "X1", m.lim = c(1, 3)), eval)
#' 
#' wts <- stackingWeights(models, data = Cement, R = 10)
#' 
#' ma <- model.avg(models)
#' Weights(ma) <- wts["mean", ]
#' 
#' predict(ma)
#' 
##~ # generate predictions on test data from all sub-models fitted to train data
##~ py.test <- vapply(fmmc, function(x, newdata) {
##~        predict(eval.parent(x), newdata = newdata, type = "response")
##~    },
##~    predict(fm, newdata = dat$test), # vapply is used for efficiency, so we need to
##~                                      # pass it a returned value's template
##~    newdata = dat$test)
##~  XXX: in the Examples, can use `sapply` for simplicity
##~ 
##~ get response variable using model's "terms" object and data.
##~ y.test <- get.response(fm, dat$test)
##~  
##~  
##~  # Fit models for the whole dataset: [?]
##~  fits <- lapply(dredge(update(fm, data = Cement), evaluate = FALSE,
##~      fixed = "X1", m.lim = c(1, 3)), eval)
#' @export
stackingWeights <-
function(object, ..., data, R, p = 0.5, seed = NULL) {
    
    models <- getModelArgs()
    m <- length(models)
    # TODO: check for lm class  / add support for other models
    if(m < 2) stop("need more than one model")

    R <- as.integer(R[1L])
    if(R <= 1) stop("'R' must be positive")
    if(p <= 0 || p >= 1) stop("'p' must be in range <0,1>")

    n <- nrow(data)
    nt <- round(n * p)

    wmat <- array(dim = c(R, m))
    r <- counter <- 1L
    counterLimit <- R * 2L
    mode(R) <- mode(counterLimit) <- "integer"
    while(counter < counterLimit && r <= R) {
        counter <- counter + 1L
        k <- sample.int(n, size = nt)
        pymat <- array(dim = c(n - nt, m))
        for(j in 1L:m) {
            fit <- models[[j]]
            tf <- terms(fit)
            fam <- family(fit)
            off <- fit$offset
            wts <- fit$weights
            coef1 <- do_glm_fit(tf, data[k, , drop = FALSE], fam, wts[k],
                off[k])$coefficients
            pymat[, j] <- predict_glm_fit(coef1,
                model.matrix(tf, data[-k, , drop = FALSE]), offset = off[-k],
                    family = fam)[, 1L]
        }
        y.test <- get.response(fit, data[-k, , drop = FALSE])
        sw1 <- tryCatch(.stacking(pymat, y.test), error = function(...) NULL)
        if(!is.null(sw1)) {
            wmat[r, ] <- sw1
            r <- r + 1L
        }
    }
    wts <- rbind(colMeans(wmat), apply(wmat, 2L, median),
                 deparse.level = 0)
    dimnames(wts) <- list(c("mean", "median"), names(models))
	structure(wts / rowSums(wts), name = "stacking", class = c("model.weights", class(wts)))
}


.stacking <-
function(predicted, observed) {
   
    if(!is.matrix(predicted))
        stop("\"predicted\" must be a matrix")
    if(nrow(predicted) != length(observed))
        stop("number of rows in \"predicted\" is not equal to length of \"observed\"")

    if (NCOL(predicted) >= length(observed))
        stop("more models than test points. ",
             "Increase the test set or reduce the number of models")
        # TODO: make the error message more specific.

    # now do an internal splitting into "folds" data sets:
    weightsopt <-
    function(ww) {
        # function to compute RMSE on test data
        w <- c(1, exp(ww))
        w <- w / sum(w) 
        ## w all in (0,1) SIMON; set weight1 always to 1, other weights are
        ## scaled accordingly (this leads to a tiny dependence of optimal
        ## weights on whether model1 is any good or utter rubbish; see by moving
        ## the 1 to the end instead -> 3rd digit changes)
        pred <- as.vector(predicted %*% w)
        return(sqrt(mean((pred - observed)^2)))
    }

    ops <- optim(par = runif(NCOL(predicted) - 1L), weightsopt, method = "BFGS")
    if (ops$convergence != 0) stop("optimization not converged")
    round(c(1, exp(ops$par)) / sum(c(1, exp(ops$par))), 4L)
}


##~ #' @rdname stacking
##~ #' @description \code{splitTrainTest} is a simple utility function to split
##~ #'     \code{data.frame} into training and test parts.
##~ #' @param data  a \code{data.frame} to be split.
##~ #' @param p proportion of the \code{data} to be used as training set. Defaults
##~ #'      to 0.5.
##~ #' @return \code{splitTrainTest} returns a \code{list} with two elements:
##~ #'     \code{"train"} and \code{"test"}, each containing a part of \code{data}.
##~ #' @export
##~ splitTrainTest <-
##~ function(data, p = 0.5) {
##~     if(p <= 0 || p >= 1) stop("'p' must be in range <0,1>")
##~     n <- nrow(data)
##~     k <- sample.int(n, size = round(n * p))
##~     list(train = data[k, , drop = FALSE], test = data[-k, , drop = FALSE])
##~ }
