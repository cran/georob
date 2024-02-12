% 2024-02-09 A. Papritz
% R CMD Rdconv -t html -o bla.html georobS3methods.Rd ; open bla.html; R CMD Rd2pdf --force georobS3methods.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  georobS3methods.Rd

\encoding{utf8}
\name{georobS3methods}
\alias{georobS3methods}
\alias{georobMethods}
\alias{coef.georob}
\alias{fixed.effects}
\alias{fixed.effects.georob}
\alias{fixef}
\alias{fixef.georob}
\alias{model.frame.georob}
\alias{model.matrix.georob}
\alias{nobs.georob}
\alias{print.coef.georob}
\alias{print.georob}
\alias{print.summary.georob}
\alias{random.effects}
\alias{random.effects.georob}
\alias{ranef}
\alias{ranef.georob}
\alias{resid.georob}
\alias{residuals.georob}
\alias{rstandard.georob}
% \alias{rstudent.georob}
\alias{summary.georob}
\alias{vcov.georob}

\title{Common S3 Methods for Class \code{georob}}

\description{This page documents the methods \code{coef}, \code{fixef},
  \code{fixed.effects}, \code{model.frame}, \code{model.matrix},
  \code{nobs}, \code{print}, \code{ranef}, \code{random.effects},
  \code{resid}, \code{residuals}, \code{rstandard}, % \code{rstudent},
  \code{summary} and \code{vcov} for the class \code{georob} which extract
  the respective components or summarize a \code{georob} object.}

\usage{
\method{coef}{georob}(object, what = c("trend", "variogram"), ...)

\method{fixef}{georob}(object, ...)

\method{fixed.effects}{georob}(object, ...)

\method{model.frame}{georob}(formula, ...)

\method{model.matrix}{georob}(object, ...)

\method{nobs}{georob}(object, ...)

\method{print}{georob}(x, digits = max(3, getOption("digits") - 3), ...)

\method{ranef}{georob}(object, standard = FALSE, ...)

\method{random.effects}{georob}(object, standard = FALSE, ...)

\method{resid}{georob}(object,
    type = c("working", "response", "deviance", "pearson", "partial"),
    terms = NULL,
    level = 1, ...)

\method{residuals}{georob}(object,
    type = c("working", "response", "deviance", "pearson", "partial"),
    terms  = NULL,
    level = 1, ...)

\method{rstandard}{georob}(model, level = 1, ...)

\method{summary}{georob}(object, correlation = FALSE, signif = 0.95, ...)

\method{vcov}{georob}(object, ...)

}

\arguments{

  \item{object, model, x}{an object of class \code{georob}, see
  \code{\link{georobObject}}.}

  \item{formula}{a model \code{\link[stats]{formula}} or
  \code{\link[stats]{terms}} object or an object of class \code{georob}, see
  \code{\link{georobObject}}.}

  \item{correlation}{a logical scalar controlling whether the correlation
  matrix of the estimated regression coefficients and of the fitted
  variogram parameters (only for non-robust fits) is computed (default
  \code{FALSE}).}

  \item{digits}{a positive integer indicating the number of decimal digits to print.}

  \item{level}{an optional integer giving the level for extracting the
  residuals from \code{object}.  \code{level = 0} extracts the regression
  residuals \eqn{\widehat{B}(\boldsymbol{s}) +
  \widehat{\varepsilon}(\boldsymbol{s})}{hatB(s) +
  hat\epsilon(s)} and \code{level = 1} (default) only the estimated errors
  \eqn{\widehat{\varepsilon}(\boldsymbol{s})}{hat\epsilon(s)}.}

  \item{signif}{a numeric with the confidence level for computing
  confidence intervals for variogram parameters (default \code{0.95}).}

  \item{standard}{a logical scalar controlling whether the spatial random effects
  \eqn{\boldsymbol{B}}{B} should be standardized (default
  \code{FALSE}).}

  \item{type}{a character keyword indicating the type of residuals to
  compute, see \code{\link[stats]{residuals.lm}}.
  \code{type = "huber"}
  computes `huberized' residuals \eqn{\widehat{\sigma} /
  \gamma_1\psi(\widehat{\varepsilon}(\boldsymbol{s}) /
  \widehat{\sigma})}{hat\sigma/\gamma_1\psi(hat\epsilon(s)/hat\sigma)}.
  }

  \item{terms}{If \code{type = "terms"}, which terms (default is all terms).}

  \item{what}{If \code{what = "trend"} (default) the function \code{coef}
  extracts the coefficients of the trend model and for \code{what =
  "variogram"} the variogram parameters.}

  \item{...}{additional arguments passed to methods.}

}
\details{

  For robust REML fits \code{deviance} returns (possibly with a warning)
  the deviance of the Gaussian REML fit of the equivalent Gaussian spatial
  linear model with heteroscedastic nugget.

  The methods \code{model.frame}, \code{model.matrix} and \code{nobs}
  extract the model frame, model matrix and the number of observations, see
  help pages of respective generic functions.

  The methods \code{residuals} (and \code{resid}) extract either the
  estimated independent errors
  \eqn{\widehat{\varepsilon}(\boldsymbol{s})}{hat\epsilon(s)}
  or the sum of the latter quantities and the spatial random effects
  \eqn{\widehat{B}(\boldsymbol{s})}{hatB(s)}.
  \code{rstandard} does the same but standardizes the residuals to unit
  variance.  \code{ranef} (\code{random.effects}) extracts the spatial
  random effects with the option to standardize them as well, and
  \code{fixef} (\code{fixed.effects}) extracts the fitted fixed-effects
  regression coefficients, which may of course also be obtained by
  \code{coef}.

  For Gaussian REML the method \code{summary} computes confidence intervals
  of the estimated variogram and anisotropy parameters from the Hessian
  matrix of the (restricted) log-likelihood (= observed Fisher
  information), based on the asymptotic normal distribution of (RE)ML
  estimates.  Note that the Hessian matrix with respect to the
  \emph{transformed} variogram and anisotropy parameters is used for this.
  Hence the inverse Hessian matrix is the covariance matrix of the
  transformed parameters, confidence intervals are first computed for the
  transformed parameters and the limits of these intervals are transformed
  back to the orginal scale of the parameters.  Optionally, \code{summary}
  reports the correlation matrix of the \emph{transformed} parameters, also
  computed from the Hessian matrix.

  Note that the methods \code{coef} and \code{summary} generate objects of
  class \code{coef.georob} and \code{summary.georob}, respectively, for
  which only \code{print} methods are available.

  Besides, the default methods of the generic functions
  \code{\link[stats]{confint}},
  \code{\link[stats]{df.residual}}, \code{\link[stats]{fitted}},
  \code{\link[stats]{formula}}, \code{\link[stats]{termplot}} and
  \code{\link[stats]{update}} can be used for objects of class
  \code{georob}.

}

\value{

  The methods \code{fixef.georob} and \code{fixed.effects.georob} return
  the numeric vector of estimated fixed-effects regression coefficients, and
  \code{vcov.georob} returns the covariance matrix of the estimated
  regression coefficients.

  The method \code{coef.georob} returns an object of class
  \code{coef.georob} which is a numeric vector with estimated fixed-effects
  regression coefficients or variogram and anisotropy parameters.  There is
  a \code{print} method for objects of class \code{coef.georob} which
  returns invisibly the object unchanged.

  The methods \code{resid.georob}, \code{residuals.georob} and
  \code{rstandard.georob} return numeric vectors of (standardized)
  residuals, and \code{ranef.georob} and \code{random.effects.georob} the
  numeric vector of (standardized) spatial random effects, see
  \emph{Details}.

  The methods \code{model.frame.georob} and \code{model.matrix.georob}
  return a model frame and the fixed-effects model matrix, respectively,
  and \code{nobs.georob} returns the number of observations used to fit a
  spatial linear model.

  The method \code{summary.georob} generates an object of class
  \code{summary.georob} which is a list with components extracted directly
  from \code{object} (\code{call}, \code{residuals}, \code{bhat},
  \code{rweights}, \code{converged}, \code{convergence.code}, \code{iter},
  \code{loglik}, \code{variogram.object}, \code{gradient},
  \code{tuning.psi}, \code{df.residual}, \code{control}, \code{terms})
  and complemented by the following components:

  \describe{

    \item{\code{scale}}{the square root of the estimated nugget effect
    \eqn{\tau^2}.}

    \item{\code{coefficients}}{a 4-column matrix with estimated regression
    coefficients, their standard errors, \var{t}-statistics and
    corresponding (two-sided) \var{p}-values.}

    \item{\code{correlation}}{an optional \code{\link{compress}}ed
    lower-triagonal matrix with the Pearson correlation coefficients of the
    estimated regression coefficients.}

    \item{\code{param.aniso}}{either a vector (robust REML) or a 3-column
    matrix (Gaussian REML) with estimated variogram and anisotropy
    parameters, complemented for Gaussian REML with confidence limits, see
    \emph{Details}.}

    \item{\code{cor.tf.param}}{an optional \code{\link{compress}}ed
    lower-triagonal matrix with the Pearson correlation coefficients of
    estimated transformed variogram and anisotropy parameters, see
    \emph{Details}.}

    \item{\code{se.residuals}}{a vector with the standard errors of the
    estimated \eqn{\varepsilon}{\epsilon}.}

  }
  There is a \code{print} methods for class \code{summary.georob} which
  invisibly returns the object unchanged.

  The method \code{print.georob} invisibly returns the object unchanged.

}


\author{
   Andreas Papritz \email{papritz@retired.ethz.ch}.
}

\seealso{
  \code{\link{georobPackage}} for a description of the model and a brief summary of the algorithms;

  \code{\link{georob}} for (robust) fitting of spatial linear models;

  \code{\link{georobObject}} for a description of the class \code{georob};

  \code{\link{profilelogLik}} for computing profiles of Gaussian likelihoods;

  \code{\link{plot.georob}} for display of RE(ML) variogram estimates;

  \code{\link{control.georob}} for controlling the behaviour of \code{georob};

  \code{\link{georobModelBuilding}} for stepwise building models of class \code{georob};

  \code{\link{cv.georob}} for assessing the goodness of a fit by \code{georob};
%  \code{\link{georobMethods}} for further methods for the class \code{georob};

  \code{\link{predict.georob}} for computing robust Kriging predictions;
%  \code{\link{validate.predictions}} for validating Kriging predictions;

  \code{\link{lgnpp}} for unbiased back-transformation of Kriging prediction
  of log-transformed data;

  \code{\link{georobSimulation}} for simulating realizations of a Gaussian process
  from model fitted by \code{georob}; and finally

  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
  for robust estimation and modelling of sample variograms.
}

\examples{
data(meuse)

## Gaussian REML fit
r.logzn.reml <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "RMexp",
    param = c(variance = 0.15, nugget = 0.05, scale = 200),
    tuning.psi = 1000)
summary(r.logzn.reml, correlation = TRUE)

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

\keyword{models}
\keyword{spatial}
\keyword{robust}
\keyword{methods}
\keyword{print}
