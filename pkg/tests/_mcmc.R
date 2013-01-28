if(length(.find.package("MCMCglmm", quiet = TRUE))) {

library(MuMIn)
library(nlme)
library(MCMCglmm)
library(lme4)

uMCMCglmm <- updateable(MCMCglmm)

fmmcmc1 <- uMCMCglmm(weight ~ Time * Diet, data = BodyWeight, random = ~ Time)
fmmcmc2 <- uMCMCglmm(weight ~ Time, data = BodyWeight, random = ~ Time)
fmlmer1 <- glmer(weight ~ Time * Diet + (1|Time), data = BodyWeight)
fmlme1 <- lme(weight ~ Time * Diet, random = ~1|Time, data = BodyWeight)
fmlm1 <- lm(weight ~ Time * Diet, data = BodyWeight)
fmlm0 <- lm(weight ~ Time, data = BodyWeight)


#fmlmer1ML <- update(fmlmer1, REML = FALSE)
#fmlme1ML <- update(fmlme1, method = "ML")

model.sel(fmmcmc1, fmmcmc2, fmlmer1, fmlme1, rank = DIC)

ddmcmc <- dredge(fmmcmc1, rank = DIC, trace = T)
ddmcmc
model.sel(ddmcmc[1:4])
# ddlme <- dredge(fmlme1)


}
