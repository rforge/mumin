if (MuMIn:::testStart("glmmML")) {

    set.seed(100)
    dat <- data.frame(y = rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1), 
        x = rnorm(100), x2 = rnorm(100), id = factor(rep(1:20, rep(5, 20))))
    
    fm1 <- glmmML(y ~ x * x2, data = dat, cluster = id, x = TRUE)
    dd <- dredge(fm1)
    # mod <- get.models(dd, subset = delta <= 4)
    summary(ma <- model.avg(dd, subset = delta <= 4))
    # vcov(ma)
    coefTable(ma)
    
}
