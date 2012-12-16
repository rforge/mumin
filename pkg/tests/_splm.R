if(length(.find.package(c("splm", "boot"), quiet = TRUE)) == 2L) {

library(splm)
library(MuMIn)

data(Produc)
data(usaww)
Produc <- Produc[Produc$year<1975, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp

respaterr <- spml(fm,
	data = Produc, listw = mat2listw(usaww), model = "random",
	spatial.error = "b")

dd <- dredge(respaterr, m.max = 3)

mod <- get.models(dd[1:3])
summary(model.avg(mod))
summary(model.avg(dd[1:3]))

respatlag <- spml(fm, data = Produc, listw = mat2listw(usaww), model="random", spatial.error="none", lag=TRUE)

model.sel(respatlag,  respaterr)
	
dd2 <- dredge(respatlag,
	fixed = "log(emp)",
	#fixed = c("log(emp)", "log(pc)", "log(pcap)"),
	varying = list(
		spatial.error = list("b","kkp","none"),
		#effect = c("individual","time","twoways"),
		model = c("random", "pooling"),
		lag = c(TRUE, FALSE)
		), eval = T,
	subset = ~ (`*nvar*` > 3) && (V(model) == "pooling") && (V(lag) != TRUE) &&
		(spatial.error != "b")
	)

summary(model.avg(dd2))
summary(model.avg(get.models(dd2)))

}

