glm.link <- function(x) {
	UseMethod("glm.link", x)
}

glm.link.family <- function(x) {
	rval <- structure(x[c("linkfun", "linkinv", "mu.eta", "valideta", "link")], class = 'link-glm')
	names(rval)[5L] <- "name"
	rval
}

glm.link.character <- function(x) {
	switch(x[1L],
		loglog = .logloglink(),
		logistic = make.link("logit"),
		make.link(x[1L]),
		NA
		)
}

`glm.link.link-glm` <- function(x) x

glm.link.default <- function(x) {
	fam <- family(x)
	if(inherits(fam, "family")) make.link(family(x)[['link']])
		else NA
}

glm.link.glmmadmb <- function(x) {
	make.link(x$link)
}

glm.link.betareg <- function(x) {
	x$link$mean
}

glm.link.polr <-
function (x, ...) {
	glm.link(x$method)
}

.logloglink <- function() {
	ret <- structure(list(linkfun = function (mu) -log(-log(mu)),
		 linkinv = function (eta) pmax(pmin(exp(-exp(-eta)),
			1 - .Machine$double.eps), .Machine$double.eps),
		 mu.eta = function (eta) {
			eta <- pmin(eta, 700)
			pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
		},
		valideta = function (eta) TRUE),
		class = "link-glm")
	ns <- asNamespace("stats")
	for(i in names(ret)) environment(ret[[i]]) <- ns
	ret$name <- "loglog"
	ret
}

	