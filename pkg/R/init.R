.onLoad <- function(libname, pkgName) {
	# cat(sprintf("onLoad(%s, %s) \n", libname, pkgName))

	asNeeded <- function(pkgName, Fun) {
		if(paste("package", pkgName, sep = ":") %in% search()) Fun() else 
			setHook(packageEvent(pkgName, "attach"), Fun)
	}
	
	asNeeded("unmarked", function(...) do.call("setMethod", list("logLik", "unmarkedFit", 
		logLik.unmarkedFit), envir = .GlobalEnv))

	asNeeded("mgcv", function(...) assign("uGamm", updateable2(mgcv::gamm, c("gamm", "list")), 1L))
	asNeeded("gamm4", function(...) assign("uGamm4", updateable2(gamm4::gamm4, c("gamm4", "gamm", 
		"list")), 1L))

}
