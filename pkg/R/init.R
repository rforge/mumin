.onLoad <- function(libname, pkgName) {
	# cat(sprintf("onLoad(%s, %s) \n", libname, pkgName))

	asNeeded <- function(pkgName, Fun) {
		if(paste("package", pkgName, sep = ":") %in% search()) Fun() else 
			setHook(packageEvent(pkgName, "attach"), Fun)
	}
	
	# ugly tricks to put MuMIn's replacement methods on top: 
	asNeeded("unmarked", function(...)
		do.call("setMethod", list("logLik", "unmarkedFit", logLik.unmarkedFit),
				envir = .GlobalEnv))

	asNeeded("lme4", function(...)
		do.call("registerS3method", list("predict", "merMod", predict.merMod),
				envir = .GlobalEnv))
	
	## XXX: nlme:::predict.gls, predict.lme
}
