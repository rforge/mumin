.onLoad <- function(libname, pkgName) {
	# cat(sprintf("onLoad(%s, %s) \n", libname, pkgName))
	
	# Do ugly tricks to put MuMIn's replacement methods on top:
	
	asNeeded <- function(pkgName, fun) {
		if(paste("package", pkgName, sep = ":") %in% search()) fun() else 
			setHook(packageEvent(pkgName, "attach"), fun)
	}
	
	regmethod <- function(funname, classname, s4 = FALSE,
		fun = get(paste0(funname, ".", classname), getNamespace(.packageName)),
		envir = .GlobalEnv)
		do.call(if(s4) "setMethod" else "registerS3method", list(funname, classname, fun),
			envir = envir)	
	
	
	asNeeded("unmarked", function(...) regmethod("logLik", "unmarkedFit", TRUE))
	asNeeded("lme4", function(...) regmethod("predict", "merMod"))
	asNeeded("nlme", function(...) regmethod("predict", "lme"))
	asNeeded("nlme", function(...) regmethod("predict", "gls"))
	
}
