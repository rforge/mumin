library(ordinal)
library(MuMIn)

# example(clmm)

## Get test data:
data(soup)

## Cumulative link mixed model with two random terms:
mm1 <- clmm(SURENESS ~ PROD + (1|RESP) + (1|RESP:PROD), data = soup,
            link = "probit", threshold = "equidistant")
## test random effect:
mm3 <- clmm(formula = SURENESS ~ PROD + GENDER + (1 | RESP) + (1 | RESP:PROD),
	data = soup, link = "probit", threshold = "equidistant")


model.sel(mm1, mm3)

dd <- dredge(mm3, m.min = 1)
summary(model.avg(dd))


install.packages(c('gamm4', 'lme4', 'aod', 'coxme', 'glmmML', 'MCMCglmm', 'pscl', 'unmarked'))


