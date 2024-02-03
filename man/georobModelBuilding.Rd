% 2024-02-01 A. Papritz
% R CMD Rdconv -t html -o bla.html georobModelBuilding.Rd ; open bla.html;
% R CMD Rd2pdf --force georobModelBuilding.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  georobModelBuilding.Rd
\encoding{utf8}
\name{georobModelBuilding}
\alias{georobModelBuilding}
\alias{add1.georob}
\alias{deviance.georob}
\alias{drop1.georob}
\alias{extractAIC.georob}
\alias{logLik.georob}
\alias{step}
\alias{step.default}
\alias{step.georob}
\alias{waldtest}
\alias{waldtest.georob}

\title{S3 Methods for Stepwise Building Fixed-Effects Models for Class \code{georob}}

\description{This page documents the methods \code{deviance},
\code{logLik}, \code{extractAIC}, \code{add1}, \code{drop1},
  \code{step} and \code{waldtest} for the class \code{georob}.  The package
  \code{georob} provides a generic \code{step} function and a default
  method which is identical with the (non-generic) function
  \code{\link[stats]{step}}.}

\usage{

\method{deviance}{georob}(object, warn = TRUE, REML = FALSE, ...)

\method{logLik}{georob}(object, warn = TRUE, REML = FALSE, ...)

\method{extractAIC}{georob}(fit, scale = 0, k = 2, ...)

\method{add1}{georob}(object, scope, scale = 0, test = c("none", "Chisq"), k = 2,
    trace = FALSE, fixed = TRUE, use.fitted.param = TRUE, verbose = 0,
    ncores = 1, ...)

\method{drop1}{georob}(object, scope, scale = 0, test = c("none", "Chisq"), k = 2,
    trace = FALSE, fixed = TRUE, use.fitted.param = TRUE, verbose = 0,
    ncores = 1, ...)

step(object, ...)

\method{step}{default}(object, scope, scale = 0,
    direction = c("both", "backward", "forward"), trace = 1,
    keep = NULL, steps = 1000, k = 2, ...)

\method{step}{georob}(object, scope, scale = 0,
    direction = c("both", "backward", "forward"), trace = 1,
    keep = NULL, steps = 1000, k = 2,
    fixed.add1.drop1 = TRUE, fixed.step = fixed.add1.drop1,
    use.fitted.param = TRUE, verbose = 0, ncores = 1, ...)

\method{waldtest}{georob}(object, ..., vcov = NULL, test = c("F", "Chisq"),
    name = NULL)

}

\arguments{

  \item{object, fit}{an object of class \code{georob}, see
    \code{\link{georobObject}}.}

  \item{direction}{a character keyword with the mode of stepwise search,
  see \code{\link[stats]{step}}.  }

  \item{fixed,fixed.add1.drop1}{a logical scalar controlling whether the
  variogram parameters are \emph{not} adjusted when \code{add}ing or
  \code{drop}ping model terms by \code{add1} and \code{drop1} (default
  \code{TRUE}), see \emph{Details}.}

  \item{fixed.step}{a logical scalar controlling whether the variogram
  parameters are \emph{not} adjusted after having called \code{add1} and
  \code{drop1} in \code{step} (default \code{TRUE}), see \emph{Details}.}

  \item{k}{a numeric specifying the 'weight' of the equivalent degrees of
    freedom (=: edf) part in the AIC formula, see
    \code{\link[stats]{extractAIC}}.}

  \item{keep}{a filter function whose input is a fitted model object and the
    associated \code{AIC} statistic, and whose output is arbitrary,
    see \code{\link[stats]{step}}.}

  \item{name}{a function for extracting a suitable name/description from a
    fitted model object.  By default the name is queried by calling
    \code{\link[stats]{formula}}, see
    \code{\link[lmtest]{waldtest}}.}

    \item{ncores}{an integer specifying the number of cores used for
    parallelized execution of \code{add1} and \code{drop1}.  If larger than
    one then the minimum of \code{ncores}, \code{parallel::detectCores()}
    and the number of terms to be added or dropped determines the number of
    cores that is actually used.}

  \item{REML}{a logical scalar controlling whether the restricted log-likelihood
    should be extracted (default \code{TRUE}).}

  \item{scale}{a numeric, currently not used, see
    \code{\link[stats]{extractAIC}}.}

  \item{scope}{defines the range of models examined in the stepwise search.
    This should be either a single formula, or a list containing
    components \code{upper} and \code{lower}, both formulae,
    see \code{\link[stats]{step}} for details.}

  \item{steps}{a numeric with the maximum number of steps to be considered
  (default is 1000), see \code{\link[stats]{step}}.}

  \item{test}{a character keyword specifying whether to compute the large
    sample Chi-squared statistic (with asymptotic Chi-squared distribution)
    or the finite sample F statistic (with approximate F distribution), see
    \code{\link[lmtest]{waldtest}}.}

  \item{trace}{a numeric.  If positive, information is printed during the
  running of \code{step}, see \code{\link[stats]{step}}.}

  \item{use.fitted.param}{a logical scalar controlling whether fitted
    values of \code{param} (and \code{aniso} are used as initial values
    when variogram parameters are fitted afresh for \code{add}ing and
    \code{drop}ping terms from the model (default \code{TRUE}).}

  \item{vcov}{a function for estimating the covariance matrix of the
    regression coefficients, see \code{\link[lmtest]{waldtest}}.}

  \item{verbose}{a positive integer controlling logging of diagnostic
    messages to the console during model fitting, see \code{\link{georob}}
    (default \code{0}).}

  \item{warn}{a logical scalar controlling whether warnings should be
  suppressed.}


  \item{...}{additional arguments passed to methods (see in particular
    \code{waldtest.default}).}

}
\details{

  For a non-robust fit the  function \code{deviance} returns the residual deviance
  \deqn{(\boldsymbol{Y} - \boldsymbol{X} \widehat{\boldsymbol{\beta}})^{\mathrm{T}}
    (\widehat{\tau}^2 \boldsymbol{I} +
    \boldsymbol{\Gamma}_{\widehat{\theta}})^{-1}
    (\boldsymbol{Y} - \boldsymbol{X} \widehat{\boldsymbol{\beta}})
    }{(Y - X hat\beta)^T(\tau^2 I + \Gamma_hat\theta)^{-1}(Y - X hat\beta)
  }
  (see \code{\link{georobPackage}} for an explanation of the notation).
  For a robust fit the deviance is not defined. The function then computes with a warning
  the deviance of an equivalent Gaussian model with heteroscedastic nugget
  \eqn{\tau^2/\boldsymbol{w}}{\tau/w} where \eqn{\boldsymbol{w}}{w} are
  the \dQuote{robustness weights} \code{rweights}, see \code{\link{georobObject}}.

  \code{logLik} returns the maximized (restricted) log-likelihood.  For
  a robust fit, the log-likelihood is not defined.  The function then
  computes the (restricted) log-likelihood of an equivalent Gaussian model with
  heteroscedastic nugget (see above).

  The methods \code{extractAIC}, \code{add1}, \code{drop1} and \code{step}
  are used for stepwise model building. \cr If \code{fixed=TRUE} or
  \code{fixed.add1.drop1=TRUE} (default) then the variogram parameters are
  kept fixed at the values of \code{object}.  For
  \code{fixed=FALSE} or \code{fixed.add1.drop1=FALSE} the variogram
  parameters are fitted afresh for each model tested by \code{add1} and
  \code{drop1}.  Then either the variogram parameters in
  \code{object$initial.objects} (\code{use.fitted.param=FALSE}) or the
  fitted parameters of \code{object} (\code{use.fitted.param=TRUE}) are
  used as initial values.  For \code{fixed.step=TRUE} the variogram
  parameters are \emph{not} fitted afresh by \code{step} after the calls to
  \code{drop1} and \code{add1} have been completed, unlike for
  \code{fixed.step=FALSE} where the parameters are estimated afresh for
  the new model that minimized AIC (BIC) in the previous step.


  In addition, the functions of the \R package \pkg{multcomp} can be used
  to test general linear hypotheses about the fixed effects of the model.

}

\value{

  The method \code{deviance.georob} returns the deviance of the fitted
  spatial linear model with the attributes \code{log.det.covmat}
  containing the logarithm of the determinant of the covariance matrix
  \eqn{
    \tau^2 \boldsymbol{I} + \boldsymbol{\Gamma}_\theta
  }{ \tau^2 I + \Gamma_\theta}
  of the observations and optionally the attribute
  \code{log.det.xticovmatx} with the logarithm of the determinant of
  \eqn{
    \boldsymbol{X}^\mathrm{T} (\tau^2
    \boldsymbol{I} + \boldsymbol{\Gamma}_\theta)^{-1} \boldsymbol{X}
  }{ X^T(\tau^2 I + \Gamma_\theta)^{-1} X}, when \code{REML = true},
  see \emph{Details} above.

  The method \code{logLik.georob} returns an object of class \code{logLik}
  with the maximized (restricted) log-likelihood, see \emph{Details} above
  and \code{\link[stats]{logLik}}.

  The method \code{extractAIC.georob} returns a numeric vector of length 2
  with the first and second elements giving the equivalent degrees of
  freedom and the (generalized) Akaike Information Criterion for the fitted
  model \code{fit}.

  The methods \code{add1.georob} and \code{drop1.georob} return objects of
  class \code{anova} which are \code{\link{data.frame}}s summarizing the
  differences in fit between models.  In addition to the customary
  variables \code{Df} and \code{AIC} the output contains a logical variable
  \code{Converged} which signals (non-)convergence when fitting the
  respective sub-model.

  The generic \code{step} returns the stepwise selected model plus
  optionally some additional attributes depending on the method.

  The methods \code{step.default} and \code{step.georob} return the
  stepwise-selected models with up to two additional components
  (\code{anova}, \code{keep}), see \code{\link[stats]{step}} for details.

  The method \code{waldtest.georob} returns an object of class \code{anova}
  which contains the residual degrees of freedom, the difference in degrees
  of freedom, Wald statistic (either "Chisq" or "F") and corresponding
  \var{p}-value.

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

%  \code{\link{georobModelBuilding}} for stepwise building models of class \code{georob};
  \code{\link{cv.georob}} for assessing the goodness of a fit by \code{georob};

  \code{\link{georobMethods}} for further methods for the class \code{georob};

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
\donttest{

data(meuse)

## Gaussian REML fit
r.logzn.reml <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "RMexp",
    param = c(variance = 0.15, nugget = 0.05, scale = 200),
    tuning.psi = 1000)
summary(r.logzn.reml, correlation = TRUE)

deviance(r.logzn.reml)
logLik(r.logzn.reml)

waldtest(r.logzn.reml, .~. + ffreq)

step(r.logzn.reml, ~ sqrt(dist) + ffreq + soil)

## robust REML fit
r.logzn.rob <- update(r.logzn.reml, tuning.psi = 1)

deviance(r.logzn.rob)
logLik(r.logzn.rob)
logLik(r.logzn.rob, REML=TRUE)

step(r.logzn.rob, ~ sqrt(dist) + ffreq + soil, fixed.step=FALSE, trace=2)
}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
