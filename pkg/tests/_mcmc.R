if(length(.find.package("MCMCglmm", quiet = TRUE))) {

library(MuMIn)
library(nlme)
library(MCMCglmm)
library(lme4)

#detach(package:MuMIn, unload = T)

fmmcmc1 <- MCMCglmm(weight ~ Time * Diet, data = BodyWeight, random = ~ Time)
fmmcmc2 <- MCMCglmm(weight ~ Time, data = BodyWeight, random = ~ Time)
fmlmer1 <- glmer(weight ~ Time * Diet + (1|Time), data = BodyWeight)
fmlme1 <- lme(weight ~ Time * Diet, random = ~1|Time, data = BodyWeight)
fmlm1 <- lm(weight ~ Time * Diet, data = BodyWeight)
fmlm0 <- lm(weight ~ Time, data = BodyWeight)


#fmlmer1ML <- update(fmlmer1, REML = FALSE)
#fmlme1ML <- update(fmlme1, method = "ML")

model.sel(fmmcmc1, fmmcmc2, fmlmer1, fmlme1, rank = DIC)

ddmcmc <- dredge(MCMCglmm(weight ~ Time * Diet, data = BodyWeight, random = ~ Time), rank = DIC, trace = T)
model.sel(ddmcmc[1:3])
# ddlme <- dredge(fmlme1)


}
