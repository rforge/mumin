# Test mixed models
if(length(find.package(c("nlme", "lme4"), quiet = TRUE)) == 2) {

library(MuMIn)
library(nlme)
library(lme4)

# example(corGaus)
fmlme1 <- lme(weight ~ Time * Diet, data = BodyWeight, random = ~ Time, method="ML")

varying <- list(
	correlation = alist(corExp(form = ~ Time), corGaus(form = ~ Time)),
	weights = alist(NULL, varPower())
)

dd <- dredge(fmlme1, m.max = 2, fixed = ~Time, varying = varying)
models <- get.models(dd, 1:4)
ma <- model.avg(models, revised = TRUE)
summary(ma1 <- model.avg(models[1:4]))
summary(ma2 <- model.avg(dd[1:4]))

summary(model.avg(models))
stopifnot(isTRUE(all.equal(coefTable(ma2), coefTable(ma1))))

logLik(dd)
mod.sel(models)
summary(ma1)
confint(ma1)
predict(ma1, BodyWeight[1:10, ], se.fit = T, level = 0)

rm(list=ls())

# TEST nlme --------------------------------------------------------------------
data(Orthodont, package = "nlme")

#:: Model-averaging mixed models :::::::::::::::::::::::::::::::::::::::::::::::
# Fitting by REML
fm2 <- lme(distance ~ Sex*age + age*Sex, data = Orthodont,
		   random = ~ 1|Subject / Sex, method = "REML")

# Model selection: ranking by AICc which uses ML
dd <- dredge(fm2, rank = "AICc", REML = FALSE)
# Get models (which are fitted by REML, like the global model)
gm <- get.models(dd, 1:4)

summary(ma <- model.avg(gm, revised = T))
confint(ma)
newdata <- data.frame(Sex = "Male", Subject = "M01", age = 8:12)

predict(ma, newdata, se.fit = T, level = 0)

detach(package:nlme); rm(list=ls())


set.seed(1)
data(Orthodont, package = "nlme")

Orthodont$rand <- runif(nrow(Orthodont))
fm2 <- lmer(log(distance) ~ rand*Sex*age + (1|Subject), data = Orthodont, REML = FALSE)

r.squaredGLMM(fm2)

#packageVersion("MuMIn")
#packageVersion("lme4")
#lme4:::predict.merMod(fm2, REform = NA) - predict.merMod(fm2)

#MuMIn:::predict.merMod(fm2, se.fit = T)
#predict.merMod(fm2, Orthodont)

#######
# summary(fm2)

# library(MCMCglmm)
# fm2 <- glmer(log(distance) ~ rand*Sex*age + (1|Subject), data = Orthodont, REML = F)
# fm2a <- MCMCglmm(log(distance) ~ rand*Sex*age, random = ~Subject, data = Orthodont)
# fm2a$call <- quote(MCMCglmm(log(distance) ~ rand*Sex*age, random = ~Subject, data = Orthodont))
# DIC(fm2a)
# DIC(fm2)

# dd1 <- dredge(fm2, rank = DIC)
# dd1a <- dredge(fm2a, rank = DIC)

# plot(dd1a)
# plot(dd1)

# summary(fm2a)

####
dd <- dredge(fm2, trace=F)
gm <- get.models(dd, 1:6)
(ma <- model.avg(gm))

# update.mer does not expand dots, so here we have a call:
# lmer(formula = distance ~ Sex + (1 | Subject), data = Orthodont,
#    REML = ..2, model = ..3)
dd <- dredge(update(fm2, REML = FALSE, model = FALSE), trace = TRUE)

fm1 <- lmer(log(distance) ~ Sex * age + (1|Sex), data = Orthodont)
# fm2 <- lmer(log(distance) ~ rand * age + (1|Subject), data = Orthodont)
fm3 <- lmer(log(distance) ~ age + (1|Subject), data = Orthodont)
fm4 <- lmer(log(distance) ~ Sex + (1|Subject), data = Orthodont)
fm5 <- lmer(log(distance) ~ age + (1|Sex) + (1|Subject), data = Orthodont)
fm6 <- lm(log(distance) ~ age, data = Orthodont)

models <- list(fm1, fm2, fm3, fm4, fm5, fm6)
(dd2 <- model.sel(models, rank = AICc, rank.args = list(REML = FALSE)))

# Comparing model.avg on model list and applied directly:
#ma0 <- model.avg(models)
ma0 <- model.avg(get.models(dd2))
ma1 <- model.avg(dd2)

summary(ma1)
stopifnot(isTRUE(all.equal(coefTable(ma1), coefTable(ma0))))

# Comparing re-ranked model.sel on model list and applied directly:
msBIC <- model.sel(models, rank = "BIC")
msAIC <- model.sel(models, rank = "AIC")
msBIC2 <- model.sel(get.models(msAIC), rank = "BIC")
msBIC3 <- model.sel(msAIC, rank = "BIC")
msAIC2 <- model.sel(msBIC3, rank = "AIC")
msAIC3 <- model.sel(get.models(msAIC2), rank = "AIC")
all.equal(msAIC2, msAIC)
# !all.equal(msAIC3, msAIC)
# !all.equal(msBIC2, msBIC)

}
#===============================================================================

#library(MASS)
#library(lme4)
#data(Insurance, package = "MASS")
#
#zo <- lmer(Claims ~ District + Age + (1 | Group) + offset(log(Holders)),
#    data = Insurance, family = poisson)
#z <- lmer(Claims ~ District + Age + (1 | Group),
#    data = Insurance, family = poisson)
#
#y1 <- predict(zo, type = "l", se.fit = T)
#y1 <- predict(zo, Insurance, type = "r", se.fit = T)
#
#fmr3 <- lmer(log(distance) ~ age + (1|Subject), data = Orthodont, REML = F)
#fmr3reml <- lmer(log(distance) ~ age + (1|Subject), data = Orthodont, REML = T)
#fme3reml <- lme(log(distance) ~ age, ~ 1 | Subject, data = Orthodont, method = "REML")
#fme3 <- lme(log(distance) ~ age, ~ 1 | Subject, data = Orthodont, method = "ML")
#z <- lmer(Claims ~ District + Age + (1 | Group), data = Insurance, family = poisson)

# library(MCMCglmm)
# # library(glmmADMB)
# library(lme4)
# data(Insurance, package = "MASS")

# mcfm <- MCMCglmm(Claims ~ District * Age, random =~ Group, data = Insurance, family = "poisson")
# mlfm <- lmer(Claims ~ District * Age + (1 | Group), data = Insurance, family = "poisson", REML = F)
# remlfm <- lmer(Claims ~ District * Age + (1 | Group), data = Insurance, family = "poisson", REML = T)

# adfm <- glmmadmb(Claims ~ District * Age + (1 | Group), data = Insurance, family = "poisson")

# DIC(mcfm)
# DIC(mlfm)
# DIC(remlfm)

# mcfm$call <- quote(MCMCglmm(Claims ~ District * Age, random =~ Group, data = Insurance, family = "poisson"))

# ddmc <- dredge(mcfm, rank = DIC)
# ddmer <- dredge(mlfm, rank = DIC)
