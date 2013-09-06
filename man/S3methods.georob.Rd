% 2013-09-06 A. Papritz
% R CMD Rdconv -t html -o bla.html S3methods.georob.Rd ; open bla.html; R CMD Rd2pdf --force S3methods.georob.Rd; 

\encoding{macintosh}
\name{georob-S3methods}
\alias{deviance.georob}
\alias{fixed.effects}
\alias{fixed.effects.georob}
\alias{fixef}
\alias{fixef.georob}
\alias{georobMethods}
\alias{logLik.georob}
\alias{model.frame.georob}
\alias{model.matrix.georob}
\alias{nobs.georob}
\alias{print.georob}
\alias{print.summary.georob}
\alias{random.effects}
\alias{random.effects.georob}
\alias{ranef}
\alias{ranef.georob}
\alias{resid.georob}
\alias{residuals.georob}
\alias{rstandard.georob}
\alias{rstudent.georob}
\alias{summary.georob}
\alias{vcov.georob}
\alias{waldtest}
\alias{waldtest.georob}

\title{Common S3 Methods for Class \code{georob}}

\description{This page documents the methods \code{fixef},
  \code{fixed.effects}, \code{logLik}, \code{model.frame},
  \code{model.matrix}, \code{nobs}, \code{print}, \code{ranef},
  \code{random.effects}, \code{resid}, \code{residuals}, \code{rstandard},
  \code{rstudent}, \code{summary}, \code{vcov} and \code{waldtest} for the
  class \code{georob}.  }

\usage{
\method{deviance}{georob}(object, ...)

\method{fixef}{georob}(object, ...)

\method{fixed.effects}{georob}(object, ...)

\method{logLik}{georob}(object, REML = FALSE, ...)

\method{model.frame}{georob}(formula, ...)

\method{model.matrix}{georob}(object, ...)

\method{nobs}{georob}(object, ...)

\method{print}{georob}(x, digits = max(3, getOption("digits") - 3), ...)

\method{ranef}{georob}(object, standard = FALSE, ...)

\method{random.effects}{georob}(object, standard = FALSE, ...)

\method{resid}{georob}(object, 
    type = c("working", "response", "deviance", "pearson", "partial" ), 
    terms = NULL,
    level = 1, ... )

\method{residuals}{georob}(object, 
    type = c("working", "response", "deviance", "pearson", "partial" ), 
    terms  = NULL,
    level = 1, ... )
    
\method{rstandard}{georob}(model, level = 1, ...)

\method{rstudent}{georob}(model, ...)

\method{summary}{georob}(object, correlation = FALSE, signif = 0.95, ...)

\method{vcov}{georob}(object, ...)

\method{waldtest}{georob}(object, ..., vcov = NULL, test = c("Chisq", "F"), 
    name = NULL, fixed = TRUE, verbose = 1)



}

\arguments{

  \item{object, model, x}{an object of class \code{georob}, see \code{\link{georobObject}}.}
  
  \item{formula}{a model \code{\link[stats]{formula}} or
  \code{\link[stats]{terms}} object or an object of class \code{georob}, see
  \code{\link{georobObject}}.}
  
  \item{correlation}{logical controlling whether the correlation matrix of
  the estimated regression coefficients and of the fitted variogram
  parameters (only for non-robust fits) is computed (default \code{FALSE}).}
  
  \item{digits}{positive integer indicating the number of decimal digits to print.}
  
  \item{fixed}{logical controlling whether the variogram parameters are
  \emph{not} adjusted when fitting alternative regression models for Wald
  tests (default \code{TRUE}).}
  
  \item{level}{an optional integer giving the level for extracting the
  residuals from \code{object}.  \code{level = 0} extracts the regression
  residuals \eqn{\widehat{B}(\mbox{\boldmath$s$\unboldmath}) +
  \widehat{\varepsilon}(\mbox{\boldmath$s$\unboldmath})}{hatB(s) +
  hat\epsilon(s)} and \code{level = 1} (default) only the estimated errors
  \eqn{\widehat{\varepsilon}(\mbox{\boldmath$s$\unboldmath})}{hat\epsilon(s)}.}
  
  \item{name}{a function for extracting a suitable name/description from a
  fitted model object.  By default the name is queried by calling
  \code{\link[stats]{formula}}, see
  \code{\link[lmtest]{waldtest}}.}
  
  \item{REML}{logical controlling whether the restricted loglikelihood
  should be extracted (default \code{TRUE}).}
  
  \item{signif}{confidence level for computing confidence intervals for
  variogram parameters (default \code{0.95}).}
  
  \item{standard}{logical controlling whether the spatial random effects
  \eqn{\mbox{\boldmath$B$\unboldmath}}{B} should be standardized (default
  \code{FALSE}).}
  
  \item{test}{character keyword specifying whether to compute the large
  sample Chi-squared statistic (with asymptotic Chi-squared distribution)
  or the finite sample F statistic (with approximate F distribution), see
  \code{\link[lmtest]{waldtest}}.}
  
  \item{type}{character keyword indicating the type of residuals to compute, see
  \code{\link[stats]{residuals.lm}}.}
  
  \item{terms}{If \code{type = "terms"}, which terms (default is all terms).}
  
  \item{verbose}{positive integer controlling logging of diagnostic
  messages to the console during model fitting by \code{\link{georob}}.}
  
  \item{vcov}{a function for estimating the covariance matrix of the
  regression coefficients, see \code{\link[lmtest]{waldtest}}.}
  
  \item{\ldots}{additional arguments passed to methods.}
  
}
\details{

  The methods \code{model.frame}, \code{model.matrix} and \code{nobs}
  extract the model frame, model matrix and the number of observations, see
  help pages of respective generic functions.
  
  The methods \code{residuals} (and \code{resid}) extract either the
  estimated independent errors
  \eqn{\widehat{\varepsilon}(\mbox{\boldmath$s$\unboldmath})}{hat\epsilon(s)}
  or the sum of the latter quantities and the spatial random effects
  \eqn{\widehat{B}(\mbox{\boldmath$s$\unboldmath})}{hatB(s)}.
  \code{rstandard} does the same but standardizes the residuals to unit
  variance.  \code{ranef} (\code{random.effects}) extracts the spatial
  random effects with the option to standardize them as well, and
  \code{fixef} (\code{fixed.effects}) extracts the fitted regression
  coefficients, which may of course also be obtained by \code{coef}.
  
  \code{deviance} and \code{logLik} return the deviance and the maximized
  (restricted) loglikelihood (non-robust fits only).
  
  Besides, the default methods of the generic functions
  \code{\link[stats]{coef}}, \code{\link[stats]{confint}},
  \code{\link[stats]{df.residual}}, \code{\link[stats]{fitted}},
  \code{\link[stats]{formula}}, \code{\link[stats]{termplot}} and
  \code{\link[stats]{update}} can be used for objects of class
  \code{georob}.
  
  In addition, the functions of the \R package \pkg{multcomp} can be used
  to test general linear hypotheses about the fixed effects of the model.
  
}


\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{georobObject}} for a description of the class \code{georob}.
}

\examples{
\dontrun{
  
data(meuse)

## Gaussian REML fit
r.logzn.reml <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "exponential",
    param = c( variance = 0.15, nugget = 0.05, scale = 200 ),
    tuning.psi = 1000,
    control = georob.control(cov.bhat = TRUE, cov.ehat.p.bhat = TRUE))
summary(r.logzn.reml, correlation = TRUE)
logLik(r.logzn.reml)

waldtest(r.logzn.reml, .~. + ffreq, verbose = 2)
waldtest(r.logzn.reml, .~. + ffreq, fixed = FALSE, verbose = 2)

## robust REML fit 
r.logzn.rob <- update(r.logzn.reml, tuning.psi = 1)
    
summary(r.logzn.rob, correlation = TRUE)

## residual diagnostics
old.par <- par(mfrow = c(2,3))

plot(fitted(r.logzn.reml), rstandard(r.logzn.reml))
abline(h = 0, lty = "dotted")
qqnorm(rstandard(r.logzn.reml))
abline(0, 1)
qqnorm(ranef(r.logzn.reml, standard = TRUE))
abline(0, 1)
plot(fitted(r.logzn.rob), rstandard(r.logzn.rob))
abline(h = 0, lty = "dotted")
qqnorm(rstandard(r.logzn.rob))
abline(0, 1)
qqnorm(ranef(r.logzn.rob, standard = TRUE))
abline(0, 1)

par(old.par)
}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
