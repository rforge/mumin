\name{MuMIn-works-with}
\alias{MuMIn-works-with}
\alias{MuMin-models}
\encoding{utf-8}
\title{List of supported models}
\description{
List of model classes supported by \code{model.avg},
\code{model.sel}, and \code{dredge}.
}

\details{
Fitted model objects that can be used as a \code{global.model} in \code{dredge} and
for building model selection tables with \code{model.sel} and for model averaging
 include ones returned by
\itemize{
\item \code{lm}, \code{glm} (package \pkg{stats});
\item \code{gam}, \code{gamm} (\pkg{mgcv});
\item \code{gamm4} (\pkg{gamm4});
\item \code{lme}, \code{gls} (\pkg{nlme});
\item \code{lmer} (\pkg{lme4});
\item \code{rlm}, \code{glm.nb}, \code{polr} (\pkg{MASS});
\item \code{multinom} (\pkg{nnet});
\item \code{sarlm}, \code{spautolm} (\pkg{spdep});
\item \code{glmmML} (\pkg{glmmML});
\item \code{coxph}, \code{survreg} (\pkg{survival});
\item \code{coxme}, \code{lmekin} (\pkg{coxme});
\item \code{hurdle}, \code{zeroinfl} (\pkg{pscl});
\item \code{rq} (\pkg{quantreg});
\item model classes from package \pkg{unmarked};
}
}

\seealso{
\code{\link{model.avg}}, \code{\link{model.sel}} and \code{\link{dredge}}.
}

\keyword{package}
