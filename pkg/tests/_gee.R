###===============================================================

if(length(.find.package("geepack", quiet = TRUE))) {

library(geepack)
#library(gee)
#library(yags)
library(MuMIn)

#example(geeglm)
data(dietox)



gee2 <- geeglm(Weight ~ Cu * Time,
	family = Gamma("inverse"), data = dietox, 
    id = Pig, corstr = "ar1")

#gee1 <- gee(Weight ~ Cu * Time,
#	family = Gamma("inverse"), data = dietox, 
#    id = Pig, corstr = "independence")


#gee3 <- yags(Weight ~ Cu * Time, id = Pig, 
#	family = Gamma, data = dietox, 
#    corstr = "independence", alphainit=0)

#QIC(gee1, gee2, gee3)

QIC(gee2)

coefTable(gee2, type = "robust")
coefTable(gee2, type = "naive")

dde <- dredge(gee2, rank = QIC)


dredge(gee2, rank = QIC, m.max = 3, ct.arg = list(type = "robust"))

#traceback()
ddr <- dredge(gee2, rank = QIC, m.max = 1, ct.arg = list(type = "robust"))
ddn <- dredge(gee2, rank = QIC, m.max = 1, ct.arg = list(type = "naive"))

attr(dde, "coefTables")

#dredge(gee2, rank = QIC)

mod <- get.models(ddn, 1:3)
modr <- get.models(ddr, 1:3)


coefTable(model.avg(mod))
coefTable(model.avg(ddn[1:3]))
coefTable(model.avg(mod, ct.args = list(type = "r")))
coefTable(model.avg(ddr[1:3]))

#traceback()

#qic <- sapply(mod, QIC)
#Weights(qic)

summary(model.avg(dde))

}
	