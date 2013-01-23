if(length(.find.package("logistf", quiet = TRUE)) == 1) {

library(logistf)
library(MuMIn)

data(sex2)

fm <- logistf(formula = case ~ oc + age, data = sex2)
fms <- update(fm, firth = FALSE)
fmg <- glm(case ~ oc + age, data = sex2, family = binomial())

dredge(fmg, m.min = 1)
dredge(fms, m.min = 1)
dredge(fm, m.min = 1)

}