\name{r.squaredGLMM}
\alias{r.squaredGLMM}

\encoding{utf-8}

\newcommand{\Rsqx}{\ifelse{latex}{\eqn{R_{#1}^{2}}{R_#1^2}}{\ifelse{html}{\enc{R_#1²}{R_#1^2}}{R_#1^2}}}
\newcommand{\Rsq}{\ifelse{latex}{\eqn{R^{2}}{R^2}}{\ifelse{html}{\enc{R²}{R^2}}{R^2}}}

\newcommand{\mydequation}{\ifelse{latex}{\deqn{#1}{#3}}{\ifelse{html}{\enc{\deqn{#2}}{#3}}{#3}}}
\newcommand{\myequation}{\ifelse{latex}{\eqn{#1}{#3}}{\ifelse{html}{\enc{\eqn{#2}}{#3}}{#3}}}

\title{Pseudo-R-squared for Generalized Mixed-Effect models}
\description{
Calculate conditional and marginal coefficient of determination for 
Generalized mixed-effect models (\Rsqx{GLMM}).
}

\usage{
r.squaredGLMM(object, null, ...)
\method{r.squaredGLMM}{merMod}(object, null, pj2014 = FALSE, ...)
}

\arguments{
  \item{object}{a fitted linear model object. }
  \item{null}{optionally, a null model, including only random effects. See \sQuote{Details}. }
  \item{pj2014}{logical, if \code{TRUE}, the result will include \Rsqx{GLMM}
  using original formulation of Johnson (2014) (i.e. \Rsqx{GLMM} extended for
  random slopes). This requires fitting \code{object} with an observation-level
  random effect term added. }
}

\value{
	\code{r.squaredGLMM} returns a two-column numeric matrix, each (possibly
	named) row holding values for marginal and conditional \Rsqx{GLMM}
	calculated with different methods, such as \dQuote{delta},
	\dQuote{log-normal}, \dQuote{trigamma}, or \dQuote{theoretical} for models
	of \code{binomial} family. 
}

\details{
For mixed-effects models, \Rsqx{GLMM} comes in two types: marginal and conditional.

\bold{Marginal \Rsqx{GLMM}} represents the variance explained by the fixed 
effects, and is defined as:

\mydequation{R_{GLMM(m)}^{2}= \frac{\sigma_f^2}{\sigma_f^2 + \sigma_{\alpha}^2 + \sigma_{\varepsilon }^2}
}{R_GLMM(m)² = (\sigma_f²) / (\sigma_f² + \sigma_\alpha² + \sigma_\epsilon²)
}{R_GLMM(m)^2 = (sigma_f^2) / (sigma_f^2 + sigma_alpha^2 + sigma_epsilon^2)
}

\bold{Conditional \Rsqx{GLMM}} is interpreted as a variance explained by the
entire model, including both fixed and random effects, and is calculated
according to the equation:

\mydequation{R_{GLMM(c)}^{2}= \frac{\sigma_f^2 + \sigma_{\alpha}^2}{\sigma_f^2 + \sigma_{\alpha}^2 + \sigma_{\varepsilon }^2}
}{R_GLMM(c)² = (\sigma_f² + \sigma_\alpha²) / (\sigma_f² + \sigma_\alpha² + \sigma_\epsilon²)
}{R_GLMM(c)^2 = (sigma_f^2 + sigma_alpha^2) / (sigma_f^2 + sigma_alpha^2 + sigma_epsilon^2)
}

where \myequation{\sigma_f^2}{\sigma_f²}{sigma_f^2}
is the variance of the fixed effect components,
\myequation{\sigma_{\alpha}}{\sigma_\alpha²}{sigma_alpha^2} 
is the variance of the random effects, and  
\myequation{\sigma_\epsilon^2}{\sigma_\epsilon²}{sigma_epsilon^2}
is the \dQuote{observation-level} variance.

Three different methods are available for deriving the observation-level variance
\eqn{\sigma_\varepsilon}: the delta method, lognormal approximation and using the
trigamma function.
%% \eqn{\Psi}.
The delta method can be used with for all distributions and link functions,
while lognormal approximation and trigamma function are limited to distributions
with logarithmic link.
Additionally, for binomial distributions, theoretical variances exist
specific for each link function distribution.

Calculation of the observation-level variance often involves fitting a
model containing no fixed effects other than intercept, otherwise identical to
the original model (including all the random effects). When comparing several 
models differing only in fixed effects, to avoid redundant calculations, such 
null model object can be provided \emph{via} the \code{null} argument.

}

\note{
\Rsqx{GLMM} can be calculated also for fixed-effect models. In
the simpliest case of \acronym{OLS} it reduces to \code{var(fitted) /
(var(fitted) + deviance / 2)}. Unlike likelihood-ratio based \Rsq for 
\acronym{OLS}, value of this statistic differs from that of
the classical \Rsq.

Currently methods exist for classes: \code{merMod}, \code{lme},
\code{glmmTMB}, \code{glmmADMB}, \code{glmmPQL}, \code{cpglm}(\code{m}) and
(\code{g})\code{lm}.

See note in \code{\link{r.squaredLR}} help page for comment on using \Rsq in 
model selection. 

}

\author{
This implementation is based on \R code from \sQuote{Supporting Information} for
Nakagawa et al (2014), and (the extension for random-slopes) Johnson (2014).
}

\references{
Nakagawa, S., Schielzeth, H. (2013) A general and simple method for obtaining
\Rsq from Generalized Linear Mixed-effects Models. \emph{Methods in
Ecology and Evolution} 4: 133–142

Johnson, P.C.D. (2014) Extension of Nakagawa & Schielzeth’s \Rsqx{GLMM} to random 
slopes models. \emph{Methods in Ecology and Evolution} 5: 44-946

Nakagawa, S., Johnson, P.C.D., Schielzeth, H. (2017) The coefficient of
determination \Rsq and intra-class correlation coefficient from generalized
linear mixed-effects models revisited and expanded.
\emph{J. R. Soc. Interface} 14: 20170213.

}

\seealso{
\code{\link{summary.lm}}, \code{\link{r.squaredLR}}
}

\examples{

\dontshow{ if(require(nlme)) \{   }
data(Orthodont, package = "nlme")

fm1 <- lme(distance ~ Sex * age, ~ 1 | Subject, data = Orthodont)

r.squaredGLMM(fm1)
r.squaredLR(fm1)
r.squaredLR(fm1, null.RE = TRUE)
\dontshow{  \}  }
}

\keyword{models}
