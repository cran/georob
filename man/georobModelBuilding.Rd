% 2014-03-12 A. Papritz
% R CMD Rdconv -t html -o bla.html georobModelBuilding.Rd ; open bla.html; R CMD Rd2pdf --force georobModelBuilding.Rd; 

\encoding{macintosh}
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
  method which is identical with the (non-generic) function \code{\link[stats]{step}}.}

\usage{

\method{deviance}{georob}(object, warn = TRUE, ...)

\method{logLik}{georob}(object, REML = FALSE, ...)

\method{extractAIC}{georob}(fit, scale = 0, k = 2, ...)

\method{add1}{georob}(object, scope, scale = 0, test = c("none", "Chisq"), k = 2,
    trace = FALSE, data = NULL, fixed = TRUE, use.fitted.param = TRUE, verbose = 0, 
    ncores = 1, ...)
    
\method{drop1}{georob}(object, scope, scale = 0, test = c("none", "Chisq"), k = 2,
    trace = FALSE, data = NULL, fixed = TRUE, use.fitted.param = TRUE, verbose = 0, 
    ncores = 1, ...)
    
step(object, ...)

\method{step}{default}(object, scope, scale = 0, 
    direction = c("both", "backward", "forward"), trace = 1, 
    keep = NULL, steps = 1000, k = 2, ...)

\method{step}{georob}(object, scope, scale = 0, 
    direction = c("both", "backward", "forward"), trace = 1, 
    keep = NULL, steps = 1000, k = 2, data = NULL, fixed = TRUE, 
    use.fitted.param = TRUE, verbose = 0, ncores = 1, ...)

\method{waldtest}{georob}(object, ..., vcov = NULL, test = c("Chisq", "F"), 
    name = NULL, fixed = TRUE)

}

\arguments{

  \item{object, fit}{an object of class \code{georob}, see
    \code{\link{georobObject}}.}
    
  \item{data}{an optional data frame.}
    
  \item{direction}{the mode of stepwise search, see
    \code{\link[stats]{step}}.  }
  
  \item{fixed}{logical controlling whether the variogram parameters are
    \emph{not} adjusted when fitting alternative regression models for Wald
    test and for \code{add}ing or \code{drop}ping model terms (default
    \code{TRUE}).}
  
  \item{k}{numeric specifying the 'weight' of the equivalent degrees of
    freedom (=: edf) part in the AIC formula, see
    \code{\link[stats]{extractAIC}}.}
  
  \item{keep}{a filter function whose input is a fitted model object and the
    associated \code{AIC} statistic, and whose output is arbitrary, 
    see \code{\link[stats]{step}}.}
  
  \item{name}{a function for extracting a suitable name/description from a
    fitted model object.  By default the name is queried by calling
    \code{\link[stats]{formula}}, see
    \code{\link[lmtest]{waldtest}}.}
  
  \item{ncores}{integer specifying the number of cores used for
    parallelized execution of \code{add1} and \code{drop1}.  If larger than
    one then the minimum of \code{ncores}, detectCores() and the number of
    terms to be added or dropped determines the number of cores that is
    actually used.}
  
  \item{REML}{logical controlling whether the restricted loglikelihood
    should be extracted (default \code{TRUE}).}
  
  \item{scale}{numeric, currently not used, see
    \code{\link[stats]{extractAIC}}.}
  
  \item{scope}{defines the range of models examined in the stepwise search.
    This should be either a single formula, or a list containing
    components \code{upper} and \code{lower}, both formulae, 
    see \code{\link[stats]{step}} for details.}
  
  \item{steps}{the maximum number of steps to be considered 
    (default is 1000), see \code{\link[stats]{step}}.}
  
  \item{test}{character keyword specifying whether to compute the large
    sample Chi-squared statistic (with asymptotic Chi-squared distribution)
    or the finite sample F statistic (with approximate F distribution), see
    \code{\link[lmtest]{waldtest}}.}
  
  \item{trace}{if positive, information is printed during the running of
    \code{step}, see \code{\link[stats]{step}}.}
  
  \item{use.fitted.param}{logical scalar controlling whether fitted
    values of \code{param} (and \code{aniso} are used as initial values
    when variogram parameters are fitted afresh for \code{add}ing and 
    \code{drop}ping terms from the model (default \code{TRUE}).}
  
  \item{vcov}{a function for estimating the covariance matrix of the
    regression coefficients, see \code{\link[lmtest]{waldtest}}.}
  
  \item{verbose}{positive integer controlling logging of diagnostic
    messages to the console during model fitting, see \code{\link{georob}}
    (default \code{0}).}
    
  \item{warn}{logical scalar controlling whether warnings should be
  suppressed.}
    
  
  \item{\ldots}{additional arguments passed to methods (see in particular
    \code{waldtest.default}).}
  
}
\details{

  For a non-robust fit the  function \code{deviance} returns the residual deviance 
  \deqn{(\mbox{\boldmath $Y -X \widehat{\beta}$})^{\mathrm{T}}
    (\widehat{\tau}^2 \mbox{\boldmath$I$\unboldmath} + 
    \mbox{\boldmath$\Gamma$\unboldmath}_{\widehat{\theta}})^{-1}
    (\mbox{\boldmath $Y -X \widehat{\beta}$})
    }{(Y - X hat\beta)^T(\tau^2 I + \Gamma_hat\theta)^{-1}(Y - X hat\beta)
  }
  (see \code{\link{georob-package}} for an explanation of the notation). 
  For a robust fit the deviance is not defined. The function then computes with a warning
  the deviance of an equivalent Gaussian model with heteroscedastic nugget 
  \eqn{\tau^2/\mbox{\boldmath $w$\unboldmath}}{\tau/w} where \eqn{\mbox{\boldmath $w$\unboldmath}}{w} are 
  the \dQuote{robustness weights} \code{rweights}, see \code{\link{georobObject}}.

  \code{logLik} returns the the maximized (restricted) loglikelihood.  For
  a robust fit, the loglikelihood is not defined.  The function then
  returns \code{NA} if \code{REML==TRUE} or computes the likelihood of an
  equivalent Gaussian model with heteroscedastic nugget (see above).
  
  The methods \code{extractAIC}, \code{add1}, \code{drop1} and \code{step}
  are used for stepwise model building.  If \code{fixed==TRUE} (default)
  then the variogram parameters are kept fixed at the values fitted for
  fixed-effects model in \code{object}.  For \code{fixed==FALSE} the
  variogram parameters are fitted for each model tested by \code{add1} and
  \code{drop1}.  Then either the variogram parameters in
  \code{object$initial.objects} (\code{use.fitted.param==FALSE} or the
  fitted parameters of \code{object} (use.fitted.param==TRUE) are used
  as initial values.
  
  In addition, the functions of the \R package \pkg{multcomp} can be used
  to test general linear hypotheses about the fixed effects of the model.
  
}


\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{georobObject}} for a description of the class \code{georob};
  \code{\link{georobMethods}} for further methods for the class \code{georob}.
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

deviance(r.logzn.reml)
logLik(r.logzn.reml)

waldtest(r.logzn.reml, .~. + ffreq)
waldtest(r.logzn.reml, .~. + ffreq, fixed = FALSE)

step(r.logzn.reml, ~ sqrt(dist) + ffreq + soil)

## robust REML fit 
r.logzn.rob <- update(r.logzn.reml, tuning.psi = 1)

deviance(r.logzn.rob)
logLik(r.logzn.rob)
logLik(r.logzn.rob, REML=TRUE)

step(r.logzn.rob, ~ sqrt(dist) + ffreq + soil, fixed=FALSE, trace=2)
}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
