if(length(.find.package("caper", quiet = TRUE))) {

library(caper)
library(MuMIn)

data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv = TRUE, vcv.dim = 3)

mod1 <- pgls(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird, lambda = 'ML')

print(AICc(mod1))

print(dd <- dredge(mod1))
print(summary(model.avg(dd)))

shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)

#fm1 <- crunch(Egg.Mass ~ F.Mass * M.Mass, data=shorebird)

uCrunch <- updateable(crunch)
fm2 <- uCrunch(Egg.Mass ~ F.Mass * M.Mass, data=shorebird)

update(fm2) # Error with 'fm1'
print(dredge(fm2))


}

