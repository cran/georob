% 2015-08-28 A. Papritz
% R CMD Rdconv -t html -o bla.html control.georob.Rd ; open bla.html; R CMD Rd2pdf --force control.georob.Rd; 
\encoding{macintosh}
\name{control.georob}
\alias{control.georob}
\alias{param.transf}
\alias{fwd.transf}
\alias{dfwd.transf}
\alias{bwd.transf}
\alias{control.rq}
\alias{control.nleqslv}
\alias{control.optim}
\alias{control.nlminb}

\title{Tuning Parameters for georob}

\description{
  This page documents parameters used to control \code{\link{georob}}.  It
  describes the arguments of the functions \code{control.georob},
  \code{param.transf}, \code{fwd.transf}, \code{dfwd.transf},
  \code{bwd.transf}, \code{control.rq}, \code{control.nleqslv},
  \code{control.nlminb} and \code{control.optim}, which all serve to
  control the behaviour of \code{\link{georob}}.}

\usage{
control.georob(ml.method = c( "REML", "ML" ), reparam = TRUE,
    maximizer = c( "nlminb", "optim" ), initial.param = TRUE,
    initial.fixef = c("lmrob", "rq", "lm"), bhat = NULL, 
    min.rweight = 0.25,
    param.tf = param.transf(), fwd.tf = fwd.transf(), 
    deriv.fwd.tf = dfwd.transf(), bwd.tf = bwd.transf(), 
    safe.param = 1.e12, psi.func = c("logistic", "t.dist", "huber"),
    tuning.psi.nr = 1000, 
    irwls.initial = TRUE, irwls.maxiter = 50, 
    irwls.ftol = 1.e-5, force.gradient = FALSE,
    min.condnum = 1.e-12, zero.dist = sqrt(.Machine[["double.eps"]]),
    error.family.estimation = c("gaussian", "long.tailed"),
    error.family.cov.effects = c("gaussian", "long.tailed"),
    error.family.cov.residuals = c("long.tailed", "gaussian"),
    cov.bhat = FALSE, full.cov.bhat = FALSE, cov.betahat = TRUE, 
    cov.bhat.betahat = FALSE, 
    cov.delta.bhat = TRUE, full.cov.delta.bhat = TRUE,
    cov.delta.bhat.betahat = TRUE,
    cov.ehat = TRUE, full.cov.ehat = FALSE,
    cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
    aux.cov.pred.target = FALSE, 
    hessian = TRUE,
    rq = control.rq(), lmrob = lmrob.control(),
    nleqslv = control.nleqslv(), 
    optim = control.optim(), nlminb = control.nlminb(), 
    pmm = control.pmm(), ...)
    
param.transf(variance = "log", snugget = "log", nugget = "log", scale = "log", 
    alpha = c( 
      RMaskey = "log", RMdewijsian = "logit2", RMfbm = "logit2", RMgencauchy = "logit2", 
      RMgenfbm = "logit2", RMlgd = "identity", RMqexp = "logit1", RMstable = "logit2"
    ), 
    beta = c( RMdagum = "logit1", RMgencauchy = "log", RMlgd = "log" ), 
    delta = "logit1", gamma = c( RMcauchy = "log", RMdagum = "logit1" ), 
    kappa = "logit3", lambda = "log", mu = "log", nu = "log",
    f1 = "log", f2  ="log", omega = "identity", phi = "identity", zeta = "identity")
  
fwd.transf(...)

dfwd.transf(...)

bwd.transf(...)

control.rq(tau = 0.5, rq.method = "br", rq.alpha = 0.1, ci = FALSE, iid = TRUE, 
    interp = TRUE, tcrit = TRUE, rq.beta = 0.99995, eps = 1e-06, 
    Mm.factor = 0.8, max.bad.fixup = 3, ...)
               
control.nleqslv(method = c("Broyden", "Newton"), 
    global = c( "dbldog", "pwldog", "qline", "gline", "none" ),
    xscalm = c( "fixed", "auto" ), control = list(ftol = 1e-04), ...)
    
control.optim(method = c("BFGS", "Nelder-Mead", "CG", 
        "L-BFGS-B", "SANN", "Brent"), lower = -Inf, upper = Inf, 
    control = list(reltol = 1e-05), ...)
    
control.nlminb(control = list( rel.tol = 1.e-5 ), lower = -Inf, 
    upper = Inf, ...)

}

\arguments{
  \item{ml.method}{character keyword defining whether non-robust maximum
  likelihood (\code{ML}) or restricted maximum likelihood (\code{REML}
  default) estimates will be computed (ignored if 
  \code{tuning.psi <= tuning.psi.nr}).}

  \item{reparam}{logical.  If \code{TRUE} (default) the reparametrized
  variance parameters \eqn{\sigma_Z^2}, \eqn{\eta} and \eqn{\xi} are
  estimated by Gaussian (RE)ML, otherwise the original parameters
  \eqn{\tau^2}, \eqn{\sigma_{\mathrm{n}}^2}{\sigma_n^2} and \eqn{\sigma^2}
  (cf.  subsection \emph{Estimating variance parameters by Gaussian
  (RE)ML}, section \emph{Details} of \code{\link{georob}}).}

  \item{maximizer}{character keyword defining the Gaussian (restricted)
  loglikelihood is maximized by \code{\link[stats]{nlminb}} (default) or
  \code{\link[stats]{optim}}.}

  \item{initial.param}{logical, controlling whether initial values of
  variogram parameters are computed for solving the estimating equations of
  the variogram and anisotropy parameters.  If \code{initial.param = TRUE}
  (default) robust initial values of parameters are computed by discarding
  outlying observations based on the \dQuote{robustness weights} of the
  initial fit of the regression model by \code{\link[robustbase]{lmrob}}
  and fitting the spatial linear model by Gaussian REML to the pruned data
  set.  For \code{initial.param = FALSE} no initial parameter values are
  computed and the estimating equations are solved with the initial values
  passed by \code{param} and \code{aniso} to \code{georob} (see
  \emph{Details} of \code{\link{georob}}.}
  
  \item{initial.fixef}{character keyword defining whether the function
  \code{\link[robustbase]{lmrob}} or \code{\link[quantreg]{rq}} is used to
  compute robust initial estimates of the regression parameters
  \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} (default \code{"lmrob"}).
  If the fixed effects model matrix has not full columns rank, then
  \code{\link[stats]{lm}} is used to compute initial values of the
  regression coefficients.}
  
  \item{bhat}{initial values for the spatial random effects
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB}, with
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}=\mbox{\boldmath$0$\unboldmath}}{hatB=0}
  if \code{bhat} is equal to \code{NULL} (default).}
  
  \item{min.rweight}{positive numeric.  \dQuote{Robustness weight} of
  the initial \code{\link[robustbase]{lmrob}} fit that observations must
  exceed to be used for computing robust initial estimates of variogram
  parameters by setting \code{initial.param = TRUE} (see
  \code{\link{georob}}; default \code{0.25}).}
  
  \item{param.tf}{a function such as \code{param.transf}, which returns a
  named vector of character strings that define the transformations to be
  applied to the variogram parameters for model fitting, see
  \emph{Details}.}
  
  \item{fwd.tf}{a function such as \code{fwd.transf}, which returns a named
  list of invertible functions to be used to transform variogram
  parameters, see \emph{Details}.}

  \item{deriv.fwd.tf}{a function such as \code{dfwd.transf}, which
  returns a named list of functions corresponding to the first derivatives
  of \code{fwd.tf}, see \emph{Details}.}

  \item{bwd.tf}{a function such as \code{bwd.transf}, which returns the
  named list of inverse functions corresponding to \code{fwd.tf}, see
  \emph{Details}.}

  \item{safe.param}{maximum acceptable value for any variogram parameter.
  If trial parameter values generated by \code{\link[stats]{nlminb}}
  \code{\link[stats]{optim}} or \code{\link[nleqslv]{nleqslv}} exceed
  \code{safe.param} then an error is signalled to force \code{optim} or
  \code{nleqslv} to update the trial values (default \code{1.e12}).}

  \item{psi.func}{character keyword defining what \eqn{\psi_c}-function should be
  used for robust model fitting.  Possible values are \code{"logistic"} (a
  scaled and shifted logistic cdf, default), \code{"t.dist"} (re-descending
  \eqn{\psi_c}-function associated with Student \eqn{t}-distribution with
  \eqn{c} degrees of freedom) and \code{"huber"} (Huber's
  \eqn{\psi_c}-function).}

  \item{tuning.psi.nr}{positive numeric.  If \code{tuning.psi} is less than
  \code{tuning.psi.nr} then the model is fitted robustly by solving the
  robustified estimating equations, and for \code{tuning.psi} equal to or
  larger than \code{tuning.psi.nr} the Gaussian (restricted) loglikelihood is
  maximized (default \code{1000}).}
  
  \item{irwls.initial}{logical.  If \code{TRUE} (default) the estimating
  equations of \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
  \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} are always solved by
  IRWLS from the initial estimates of
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} and
  \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta}.  If
  \code{FALSE} then IRWLS starts from respective estimates computed for the
  variogram parameter estimates of the previous iteration of \code{nleqslv}
  or \code{optim}.}
  
  \item{irwls.maxiter}{positive integer equal to the maximum number of
  IRWLS iterations to solve the estimating equations of
  \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
  \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} (default \code{50}).}
   
  \item{irwls.ftol}{numeric convergence criterion for IRWLS. Convergence is
  assumed if the objective function changes in one IRWLS iteration does not
  exceed \code{ftol}.}
  
  \item{force.gradient}{logical controlling whether the estimating
  equations or the gradient of the Gaussian restricted loglikelihood are
  evaluated even if all variogram parameters are fixed (default
  \code{FALSE}).}

  \item{min.condnum}{positive numeric.  Minimum acceptable ratio of smallest to
  largest singular value of the model matrix
  \eqn{\mbox{\boldmath$X$\unboldmath}}{X} (default \code{1.e-12}).}
  
  \item{zero.dist}{positive numeric equal to the maximum distance, separating two
  sampling locations that are still considered as being coincident.}
  
  \item{error.family.estimation}{character keyword, defining the
  probability distribution for \eqn{\varepsilon} (default:
  \code{"gaussian"}) that is used to approximate the covariance of
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB}, see
  \emph{Details}.}
  
  \item{error.family.cov.effects}{character keyword, defining the
  probability distribution for \eqn{\varepsilon} (default:
  \code{"gaussian"}) that is used to approximate the covariances of
  \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hatB},
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} and
  \eqn{\mbox{\boldmath$B$\unboldmath}-\widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB},
  see \emph{Details}.}
  
  \item{error.family.cov.residuals}{character keyword, defining the
  probability distribution for \eqn{\varepsilon} (default:
  \code{"long.tailed"}) that is used to approximate the covariances of
  \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}=\mbox{\boldmath$Y$\unboldmath}
  - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}} -
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{hat\epsilon=Y-X hat\beta -
  hatB} and \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}} +
  \widehat{\mbox{\boldmath$B$\unboldmath}} =\mbox{\boldmath$Y$\unboldmath}
  - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\epsilon+ hatB=Y-X
  hat\beta}, see \emph{Details}.}
  
  \item{cov.bhat}{logical controlling whether the covariances of
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} are returned by
  \code{georob} (default \code{FALSE}).}
  
  \item{full.cov.bhat}{logical controlling whether the full covariance
  matrix (\code{TRUE}) or only the variance vector of
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} is returned (default
  \code{FALSE}).}
  
  \item{cov.betahat}{logical controlling whether the covariance matrix of
  \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta} is returned
  (default \code{TRUE}).}
  
  \item{cov.bhat.betahat}{logical controlling whether the covariance matrix
  of \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} and
  \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta} is returned
  (default \code{FALSE}).}
  
  \item{cov.delta.bhat}{logical controlling whether the covariances of
  \eqn{\mbox{\boldmath$B$\unboldmath}-
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB} are returned (default \code{TRUE}).}
  
  \item{full.cov.delta.bhat}{logical controlling whether the full covariance
  matrix (\code{TRUE}) or only the variance vector of
  \eqn{\mbox{\boldmath$B$\unboldmath}-
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB} is returned (default \code{TRUE}).}
  
  \item{cov.delta.bhat.betahat}{logical controlling whether the covariance
  matrix of \eqn{\mbox{\boldmath$B$\unboldmath}-
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB} and
  \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta} is returned
  (default \code{TRUE}).}
  
  \item{cov.ehat}{logical controlling whether the covariances of
  \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}=\mbox{\boldmath$Y$\unboldmath}
  - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}} -
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{hat\epsilon=Y-X hat\beta -
  hatB} are returned (default \code{TRUE}).}
  
  \item{full.cov.ehat}{logical controlling whether the full covariance
  matrix (\code{TRUE}) or only the variance vector of
  \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}=\mbox{\boldmath$Y$\unboldmath}
  - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}} -
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{hat\epsilon=Y-X hat\beta -
  hatB} is returned (default \code{FALSE}).}
  
  \item{cov.ehat.p.bhat}{logical controlling whether the covariances of
  \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}} +
  \widehat{\mbox{\boldmath$B$\unboldmath}} =\mbox{\boldmath$Y$\unboldmath}
  - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\epsilon+ hatB=Y-X
  hat\beta} are returned (default \code{FALSE}).}
  
  \item{full.cov.ehat.p.bhat}{logical controlling whether the full
  covariance matrix (\code{TRUE}) or only the variance vector
  of \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}
  + \widehat{\mbox{\boldmath$B$\unboldmath}}
  =\mbox{\boldmath$Y$\unboldmath} - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\epsilon+ hatB=Y-X
  hat\beta} is returned (default \code{FALSE}).}
  
  \item{aux.cov.pred.target}{logical controlling whether a covariance term
  required for the back-transformation of kriging predictions of
  log-transformed data is returned (default \code{FALSE}).}
      
  \item{hessian}{logical scalar controlling whether for Gaussian (RE)ML the
  Hessian should be computed at the MLEs.}
  
  \item{rq}{a list of arguments passed to \code{rq} or a function such as
  \code{control.rq} that generates such a list (see
  \code{\link[quantreg]{rq}} for allowed arguments).}
  
  \item{lmrob}{a list of arguments passed to the \code{control} argument of
  \code{lmrob} or a function such as
  \code{lmrob.control} that generates such a list (see
  \code{\link[robustbase]{lmrob.control}} for allowed arguments).}
  
  \item{nleqslv}{a list of arguments passed to
  \code{nleqslv} or a function such as
  \code{control.nleqslv} that generates such a list (see
  \code{\link[nleqslv]{nleqslv}} for allowed arguments).}
  
  \item{nlminb}{a list of arguments passed to \code{nlminb}
  or a function such as \code{control.nlminb} that generates such a list
  (see \code{\link[stats]{nlminb}} for allowed arguments).}
  
  \item{optim}{a list of arguments passed to \code{optim} or a function
  such as \code{control.optim} that generates such a list (see
  \code{\link[stats]{optim}} for allowed arguments).}
  
  \item{pmm}{a list of arguments, passed  e.g. to \code{\link{pmm}} or a
  function such as \code{control.pmm} that generates such a list
  (see \code{\link{control.pmm}} for allowed arguments).}
  
  \item{...}{for \code{fwd.transf}, \code{dfwd.transf} and
  \code{bwd.transf} a named vectors of functions, extending the definition
  of transformations for variogram parameters (see \emph{Details}).}
  
  \item{variance, snugget, nugget, scale, alpha, beta, delta, gamma,
  kappa, lambda, mu, nu}{character strings with names of transformation
  functions of the variogram parameters.}
  
  \item{f1, f2, omega, phi, zeta}{character
  strings with names of transformation functions of the variogram
  parameters.}

  \item{tau, rq.method, rq.alpha, ci, iid, interp, tcrit}{arguments passed
  as \code{...} to \code{\link[quantreg]{rq}}.}
  
  \item{rq.beta, eps, Mm.factor, max.bad.fixup}{arguments passed as
  \code{...} to \code{\link[quantreg]{rq}}.}
  
  \item{method, global, xscalm, control, lower, upper, reltol,
  rel.tol}{arguments passed to related arguments of
  \code{\link[nleqslv]{nleqslv}}, \code{\link[stats]{nlminb}} and
  \code{\link[stats]{optim}}, respectively.}
      
}
%   \item{}

\details{

  \subsection{Parameter transformations}{

    The arguments \code{param.tf}, \code{fwd.tf}, \code{deriv.fwd.tf},
    \code{bwd.tf} define the transformations of the variogram parameters for
    RE(ML) estimation.  Implemented are currently \code{"log"},
    \code{"logit1"}, \code{"logit2"}, \code{"logit3"} (various variants of
    logit-transformation, see code of function \code{fwd.transf}) and \code{"identity"} (= no)
    transformations.  These are the possible values that the many arguments
    of the function \code{param.transf} accept (as quoted character strings)
    and these are the names of the list components returned by
    \code{fwd.transf}, \code{dfwd.transf} and \code{bwd.transf}.  Additional
    transformations can be implemented by:
  
    \enumerate{
    
      \item Extending the function definitions by arguments like
      
      \code{fwd.tf = fwd.transf(c(my.fun = function(x) your transformation))},\cr
      \code{deriv.fwd.tf = dfwd.transf(c(my.fun = function(x) your derivative))},\cr
      \code{bwd.tf = bwd.transf(c(my.fun = function(x) your back-transformation))},
      
      \item Assigning to a given argument of \code{param.transf} the name of
      the new function, e.g.\cr \code{variance = "my.fun"}.  }
    
    Note the values given for the arguments of \code{param.transf} must match
    the names of the functions returned by \code{fwd.transf},
    \code{dfwd.transf} and \code{bwd.transf}. 
    
  }
  
  \subsection{Approximation of covariances of fixed and random efffects and
  residuals}{
    
    The robustified estimating equations of robust REML depend on the
    covariances of \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB}.
    These covariances (and the covariances of
    \eqn{\mbox{\boldmath$B$\unboldmath}-\widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB},
    \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta},
    \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}}{hat\epsilon},
    \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}} +
    \widehat{\mbox{\boldmath$B$\unboldmath}}}{hat\epsilon+ hatB}) are
    approximated by expressions that in turn depend on the variances of
    \eqn{\varepsilon}, \eqn{\psi(\varepsilon/\tau)} and the expectation of
    \eqn{\psi'(\varepsilon/\tau) (= \partial / \partial \varepsilon
    \psi(\varepsilon/\tau))}.  The arguments
    \code{error.family.estimation}, \code{error.family.cov.effects} and
    \code{error.family.cov.residuals} control what parametric distribution
    for \eqn{\varepsilon} is used to compute the latter quantities.
    Possible options are: \code{"gaussian"} or \code{"long.tailed"}.  In the
    latter case the pdf of \eqn{\varepsilon} is assumed to be proportional
    to \eqn{1/\tau exp(-\rho(\varepsilon/\tau))} where \eqn{\psi(x)=\rho'(x)}.
  
  }
  
}

\examples{
\dontrun{
data( meuse )

r.logzn.rob <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "RMexp",
    param = c( variance = 0.15, nugget = 0.05, scale = 200 ),
    tuning.psi = 1, control = control.georob(cov.bhat = TRUE, 
    cov.ehat.p.bhat = TRUE, initial.fixef = "rq"), verbose = 2)
  
qqnorm(rstandard(r.logzn.rob, level = 0)); abline(0, 1)
qqnorm(ranef(r.logzn.rob, standard = TRUE)); abline(0, 1)
}
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{georobObject}} for a description of the class \code{georob};
  \code{\link{plot.georob}} for display of RE(ML) variogram estimates;
  \code{\link{predict.georob}} for computing robust kriging predictions; and finally
  \code{\link{georobMethods}} for further methods for the class \code{georob}.
}

% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
% 
% }

\keyword{models}
\keyword{spatial}
\keyword{robust}