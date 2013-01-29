`updateable` <-
function (FUN) {
	Fname <- substitute(FUN)
	FUN <- match.fun(FUN)
	FUNV <- function() {
		ocl <- cl <- match.call()
		cl[[1L]] <- Fname
		res <- eval(cl, parent.frame())
		if(!isS4(res) && is.list(res))
			res$call <- ocl else
			attr(res, "call") <- ocl
		res
	}
    formals(FUNV) <- formals(FUN)
    FUNV
}

`updateable2` <-
function (FUN, Class) {
	Fname <- substitute(FUN)
	FUN <- match.fun(FUN)
	FUNV <- function() {
		ocl <- cl <- match.call()
		cl[[1L]] <- Fname
		res <- eval(cl, parent.frame())
		if(!isS4(res) && is.list(res))
			res$call <- ocl else
			attr(res, "call") <- ocl
		class(res) <- Class
		res
	}
	if(missing(Class)) body(FUNV)[[6L]] <- NULL
    formals(FUNV) <- formals(FUN)
	rm(FUN)
    FUNV
}


`.getCall` <- function(x) {
	if(isS4(x)) {
		if(any(i <- (sln <- c("call", "CALL", "Call")) %in% slotNames(x)))
			slot(x, sln[i][1L]) else
			if(!is.null(attr(x, "call")))
				attr(x, "call") else NULL
	} else {
		if(!is.atomic(x) && !is.null(x$call)) {
			x$call
		} else if(!is.null(attr(x, "call"))) {
			attr(x, "call")
		} else
			NULL
	}
}





`getCall.default` <-
function (x, ...)
.getCall(x)



##==============================================================================

`updGamm` <-
function(...) {
	ocl <- cl <- match.call(definition = Fun <- get("gamm", asNamespace("mgcv")))
	cl[[1L]] <- call("get", "gamm", asNamespace("mgcv"))
	res <- eval(cl, parent.frame())
	res$call <- ocl
	class(res) <- c("gamm", "list")
	res
}

`updGamm4` <-
function(...) {
	ocl <- cl <- match.call(definition = Fun <- get("gamm4", asNamespace("gamm4")))
	cl[[1L]] <- call("get", "gamm4", asNamespace("gamm4"))
	res <- eval(cl, parent.frame())
	res$call <- ocl
	class(res) <- c("gamm4", "gamm", "list")
	res
}

