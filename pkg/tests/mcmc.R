if(length(.find.package("MCMCglmm", quiet = TRUE))) {

library(MuMIn)
library(nlme)
#do.call(library, alist(MCMCglmm))
library(MCMCglmm)

# fmlme1 <- lme(weight ~ Time * Diet, data = BodyWeight, random = ~ Time, method="ML")
fmmcmc1 <- MCMCglmm(weight ~ Time * Diet, data = BodyWeight, random = ~ Time)
ddmcmc <- dredge(MCMCglmm(weight ~ Time * Diet, data = BodyWeight, random = ~ Time), rank = DIC, trace = T)
model.sel(ddmcmc[1:3])
# ddlme <- dredge(fmlme1)
}
