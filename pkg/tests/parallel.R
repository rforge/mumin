if (do.call("require", list("parallel", quietly = TRUE)) ||
	do.call("require", list("snow", quietly = TRUE))) {
	
	library(MuMIn)

	library(lme4)
	data(Orthodont, package = "nlme")

	# Orthodont$rand <- runif(nrow(Orthodont))
	# fm2 <- lmer(log(distance) ~ rand*Sex*age + (1|Subject) + (1|Sex), 
	#	data = Orthodont, REML=FALSE)
	fm2 <- lmer(log(distance) ~ Sex*age + (1|Subject) + (1|Sex),
		data = Orthodont, REML = FALSE)

	
	clust <- makeCluster(getOption("cl.cores", 2))
	clusterExport(clust, "Orthodont")
	#clusterEvalQ(clust, library(lme4))
	clusterCall(clust, "library", "lme4", character.only = TRUE)

	print(system.time(pddc <- pdredge(fm2, cluster = clust)))
	print(system.time(pdd1 <- pdredge(fm2, cluster = F)))
	print(system.time(dd1 <- dredge(fm2)))
	
	
	
	
# library(unmarked)
# library(MuMIn)
# data(mallard)
# mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site, obsCovs = mallard.obs)
# (fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest, mallardUMF, K=30))\
# (fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest, mallardUMF, K=30))
# do.call("require", list("snow", quietly = TRUE))
# clust <- makeCluster(getOption("cl.cores", 2))
# clusterCall(clust, "library", "unmarked", character.only = TRUE)
# clusterExport(clust, "mallardUMF")
## fm.mallard
# system.time(pdd <- pdredge(fm.mallard))
# system.time(pdd <- pdredge(fm.mallard, clust))
# system.time(pdd1 <- pdredge(fm.mallard, NA))
# pdd
# stopCluster(clust)
	
	

	stopCluster(clust)

	stopifnot(identical(c(pddc), c(pdd1)) && identical(c(pdd1), c(dd1)))

	
# suppressPackageStartupMessages(library(spdep))

# suppressMessages(example(NY_data, echo = FALSE))

# esar1f <- spautolm(Z ~ PEXPOSURE * PCTAGE65P + PCTOWNHOME,
 # data=nydata, listw=listw_NY, family="SAR", method="full", verbose=FALSE)
 
# clusterCall(clust, "library", "spdep", character.only = TRUE)
# clusterExport(clust, "listw_NY", "nydata")

# options(warn=1)

# varying <- list(family = list("CAR", "SAR"), method=list("Matrix_J", "full"))

# dd <- dredge(esar1f, m.max=1,  fixed=~PEXPOSURE, varying = varying, trace=FALSE)

}

#system.time(pdredge(fm2, cluster = clust))
#system.time(pdredge(fm2, cluster = F))
#system.time(dredge(fm2))

#traceback()
