library(mgcv)
library(gamm4)
library(MuMIn)


RNGkind("Mersenne")
set.seed(0) # 16
dat <- gamSim(6, n=100, scale=5, dist="normal")

fmgs2 <- gamm(y ~s(x0)+ s(x3) + s(x2), family=gaussian, data=dat, random = list(fac=~1))
dd <- dredge(fmgs2)
dd
summary(model.avg(dd, delta <= 4))

fmg4s2 <- gamm(y ~s(x0)+ s(x3) + s(x2), family=gaussian, data=dat, random = ~ (1|fac))
dd4 <- dredge(fmg4s2)
dd4
summary(avg <- model.avg(dd4, delta <= 4))

confint(avg)
