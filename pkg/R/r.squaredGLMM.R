
`r.squaredGLMM` <-
function(x, nullfx = NULL, ...)
	UseMethod("r.squaredGLMM")

`r.squaredGLMM.default` <-
function(x, nullfx = NULL, ...) 
	.NotYetImplemented()
	
`r.squaredGLMM.lme` <-
function(x, ...) {
	VarFx <- var(fitted(x, level = 0L))
		
	### x$data is the original data.frame, not subset'ted nor na.omit'ted
	cl <- x[['call']]
	mfArgs <- list(formula = asOneFormula(c(formula(x),
		lapply(x$modelStruct$reStruct, attr, 'formula'))),
		data = x$data, na.action = x$call$na.action)
	if (!is.null(cl$subset))  mfArgs[["subset"]] <-
		asOneSidedFormula(cl[["subset"]])[[2L]]
	mfArgs$drop.unused.levels <- TRUE
	dataMix <- do.call("model.frame", mfArgs)
	
	mMfull <- model.matrix(x$modelStruct$reStruct, data = dataMix)
	n <- nrow(mMfull)

	varRan <- sum(sapply(x$modelStruct$reStruct, function(z) {
		sig <- pdMatrix(z) * x$sigma^2
		Mm1 <-  mMfull[, rownames(sig)]
		sum(diag(Mm1 %*% sig %*% t(Mm1))) / n
	}))
	
	varAll <- sum(VarFx, varRan)
	res <- c(VarFx, varAll) / (varAll + x$sigma^2)
	names(res) <- c("R2m", "R2c")
	res
}


`r.squaredGLMM.merMod` <-
`r.squaredGLMM.mer` <-
function(x, nullfx = NULL) {
	cl <- getCall(x)
	envir <- environment(formula(x))
	## this replaces all '(x | y)' to '(x)' 
	frm <- 	.substFun4Fun(formula(x), "|", function(e, ...) e[[2L]])
	mmAll <- model.matrix(frm, data = model.frame(x), contrasts.arg = eval(cl$contrasts, envir = envir))

	vc <- VarCorr(x)
	n <- nrow(mmAll)
	fx <- fixef(x) # fixed effect estimates
	varRan <- sum(sapply(vc, function(sig) {
		mm1 <-  mmAll[, rownames(sig)]
		sum(diag(mm1 %*% sig %*% t(mm1))) / n
	}))
		
	.rsqGLMM(x, fam = family(x),
		varFx = var(as.vector(model.matrix(x) %*% fx)),
		varRan = varRan,
		resVar = attr(vc, "sc")^2,
		fxNullCoef = fixef(if(is.null(nullfx)) 
				null.fit(x, RE.keep = TRUE, evaluate = TRUE) else nullfx
				)
		)
}

`r.squaredGLMM.glmmML` <-
function(x, nullfx = NULL, ...) {
	if(is.null(x$x))
		stop("glmmML must be fitted with 'x = TRUE'")

	.rsqGLMM(x, family(x),
			 varFx = var(as.vector(x$x %*% coef(x))),
			 varRan = x$sigma^2, resVar = NULL,
			 fxNullCoef = coef(if(is.null(nullfx)) 
				null.fit(x, RE.keep = TRUE, evaluate = TRUE) else nullfx
				)
			)
}


`r.squaredGLMM.lm` <-
function(x, ...)
	.rsqGLMM(x, family(x),
		 varFx = var(as.vector(coef(x) %*% t(model.matrix(x)))),
		 #varFx = var(fitted(x)),
		 varRan = 0,
		 resVar = sum(if(is.null(x$weights)) resid(x)^2 else
					   resid(x)^2 * x$weights) / df.residual(x),
		 fxNullCoef = coef(null.fit(x, evaluate = TRUE)))


`.rsqGLMM` <-
function(x, fam, varFx, varRan, resVar, fxNullCoef) {
	v <- switch(paste(fam$family, fam$link, sep = "."), 
		gaussian.identity = resVar,
		binomial.logit = 3.28986813369645, #  = pi^2 / 3
		binomial.probit = 1,
		poisson.log = log(1 / exp(fxNullCoef) + 1),
		poisson.sqrt = 0.25,
		stop("do not know how to calculate variance for this family/link combination")
	)
	varAll <- sum(varFx, varRan)
	res <- c(varFx, varAll) / (varAll + v)
	if(fam$family == "poisson") {
		res[2L] <- NA_real_
		warning(simpleWarning("conditional statistic for 'poisson' family cannot be (yet) calculated", 
			sys.call(1L)))
	}
	names(res) <- c("R2m", "R2c")
	res
}

