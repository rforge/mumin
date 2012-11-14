###===============================================================

if(length(.find.package("geepack", quiet = TRUE))) {

library(geepack)
library(MuMIn)

#example(geeglm)
data(dietox)



gee2 <- geeglm(Weight ~ Cu * Time,
	family = Gamma("inverse"), data = dietox, 
    id = Pig, corstr = "ar1")

#QIC(gee2)

dde <- dredge(gee2, rank = QIC)

dredge(gee2, rank = QIC, type = "I")
dredge(gee2, rank = QIC, type = "R")

mod <- get.models(dde, 1:3)

#qic <- sapply(mod, QIC)
#Weights(qic)


summary(model.avg(dde))

}
	