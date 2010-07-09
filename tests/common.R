
# Test lm:

data(Cement)
lm1 <- lm(y ~ ., data = Cement)
dd.lm <- dredge(lm1)

# Examples of using 'subset':
dredge(lm1, subset = !X1 | !X2)
subset(dd.lm, is.na(X1) | is.na(X2))

# keep only models with X3
dredge(lm1, subset = X3)
# the same, but more effective:
dredge(lm1, fixed = ~ X3)

gm.lm <- get.models(dd.lm, subset = delta < 4)

ma.lm <- model.avg(gm.lm)

predict(ma.lm)





# Test lmer

require(lme4)
fm3 <- lmer(distance ~ age + Sex + (1 | Subject), data = Orthodont, REML=FALSE)
dd3 <- dredge(fm3)

# Get top-most models, but fitted by REML:
(top.models.3 <- get.models(dd3, subset = delta < 4, REML=TRUE))
# use: method = "REML" for older versions of lme4
