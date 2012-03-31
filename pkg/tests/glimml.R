if(length(.find.package("glmmML", quiet = TRUE))) {
library(MuMIn)
library(aod)
library(glmmML)
library(lme4)
library(MASS)

data(orob2)
data(salmonella)
#_______________________________________________________________________________

# Fit some models
fmbb <- betabin(cbind(y, n - y) ~ root, ~ seed, data = orob2)
fmb1 <- betabin(cbind(y, n - y) ~ root, ~ 1, data = orob2)
fmgm <- glmmML(cbind(y, n - y) ~ root, cluster = seed, data = orob2)
fmbr <- lmer(cbind(y, n - y) ~ root + (1 | seed), family = binomial("logit"),
	data = orob2)
fmgl <- glm(cbind(y, n - y) ~ root, family = binomial("logit"), data = orob2)
fmgl0 <- glm(cbind(y, n - y) ~ 1, family = binomial("logit"), data = orob2)

models <- list(fmbb, fmb1, fmgm, fmbr, fmgl)
ms <- model.sel(models)
summary(am <- model.avg(ms))
coef(am)
coefTable(am)
confint(am)

salmonella$r1 <- factor(sample(3, nrow(salmonella), replace = T))

fmng <-  negbin(y ~ log(dose + 10) + dose*r1, ~ 1, salmonella)
#fmngr <-  negbin(y ~ log(dose + 10) + dose, ~ r1, salmonella)
fmps <-  glm(y ~ log(dose + 10) + dose*r1, poisson, salmonella)
fm.nb <- glm.nb(y ~ log(dose + 10) + dose*r1, link = "log", data = salmonella)

models <- list(fmng, fm.nb, fmps)
(ms <- model.sel(models))
summary(model.avg(ms))

#ms$link <- NULL
#ms$init.theta <- NULL
ms
fmnbo <- negbin(y ~ group + log(trisk) + offset(log(trisk)), ~ group, dja)
#fmnbo <- negbin(y ~ group + offset(log(trisk)), ~ group, dja)
dredge(fmnbo)

#_______________________________________________________________________________

orob2$rand <- runif(nrow(orob2), 0, 100)
orob2$rand2 <- rnorm(nrow(orob2), 0, 10)

#fm1 <- betabin(cbind(y, n - y) ~ seed, ~ 1, data = orob2)
#fm2 <- betabin(cbind(y, n - y) ~ seed + root, ~ 1, data = orob2)
#fm3 <- betabin(cbind(y, n - y) ~ seed * root, ~ 1, data = orob2)
fm3r <- betabin(cbind(y, n - y) ~ root + rand * rand2, ~ seed, data = orob2)

#options(warn=1)

(dd <- dredge(fm3r, trace = T))
summary(model.avg(dd))
summary(model.avg(dd, delta < 4))

}
