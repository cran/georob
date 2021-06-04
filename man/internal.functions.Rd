% 2021-06-04 A. Papritz
% R CMD Rdconv -t html -o bla.html internal.functions.Rd ; open bla.html; R CMD Rd2pdf --force internal.functions.Rd

\encoding{utf8}
\name{internal.functions}
\alias{check.newdata}
\alias{covariances.fixed.random.effects}
\alias{crpsnorm}
% \alias{crpspd.resp.rob}
\alias{estimate.zhat}
\alias{estimating.equations.theta}
\alias{estimating.equations.z}
\alias{f.aux.Qstar}
\alias{f.aux.RSS}
\alias{f.aux.Valphaxi}
\alias{f.aux.add1.drop1}
\alias{f.aux.crpsnorm}
\alias{f.aux.eeq}
\alias{f.aux.gcr}
\alias{f.aux.gradient.nll}
\alias{f.aux.gradient.npll}
\alias{f.aux.print.gradient}
\alias{f.aux.tf.param.fwd}
\alias{f.call.set_allfitxxx_to_false}
\alias{f.call.set_allxxx_to_fitted_values}
\alias{f.call.set_onefitxxx_to_value}
\alias{f.call.set_onexxx_to_value}
\alias{f.call.set_x_to_value}
\alias{f.call.set_x_to_value_in_fun}
\alias{f.diag}
\alias{f.psi.function}
\alias{f.reparam.bkw}
\alias{f.reparam.fwd}
\alias{f.robust.uk}
\alias{f.stop.cluster}
\alias{georob.fit}
\alias{getCall.georob}
\alias{gradient.negative.loglikelihood}
\alias{K}
\alias{lmrob.control}
\alias{likelihood.calculations}
\alias{negative.log-likelihood}
\alias{partial.derivatives.variogram}
% \alias{ppd.resp.rob}
% \alias{qpd.resp.rob}
\alias{safe_pchisq}
\alias{simple.kriging.weights}
% \alias{var.pd.resp.rob}
\alias{update.zhat}

%
\title{Internal Functions of Package \pkg{georob}}
%
\description{
The unexported internal functions
  \itemize{
    \item{\code{check.newdata}}
    \item{\code{covariances.fixed.random.effects}}
    \item{\code{crpsnorm}}
    % \item{\code{crpspd.resp.rob}}
    \item{\code{estimate.zhat}}
    \item{\code{estimating.equations.theta}}
    \item{\code{estimating.equations.z}}
    \item{\code{f.aux.Qstar}}
    \item{\code{f.aux.RSS}}
    \item{\code{f.aux.Valphaxi}}
    \item{\code{f.aux.add1.drop1}}
    \item{\code{f.aux.crpsnorm}}
    \item{\code{f.aux.eeq}}
    \item{\code{f.aux.gcr}}
    \item{\code{f.aux.gradient.nll}}
    \item{\code{f.aux.gradient.npll}}
    \item{\code{f.aux.print.gradient}}
    \item{\code{f.aux.tf.param.fwd}}
    \item{\code{f.call.set_allfitxxx_to_false}}
    \item{\code{f.call.set_allxxx_to_fitted_values}}
    \item{\code{f.call.set_onefitxxx_to_value}}
    \item{\code{f.call.set_onexxx_to_value}}
    \item{\code{f.call.set_x_to_value}}
    \item{\code{f.call.set_x_to_value_in_fun}}
    \item{\code{f.diag}}
    \item{\code{f.psi.function}}
    \item{\code{f.reparam.bkw}}
    \item{\code{f.reparam.fwd}}
    \item{\code{f.robust.uk}}
    \item{\code{f.stop.cluster}}
    \item{\code{georob.fit}}
    \item{\code{getCall.georob}}
    \item{\code{gradient.negative.loglikelihood}}
    \item{\code{K}}
    \item{\code{lmrob.control}}
    \item{\code{likelihood.calculations}}
    \item{\code{negative.log-likelihood}}
    \item{\code{partial.derivatives.variogram}}
    % \item{\code{ppd.resp.rob}}
    % \item{\code{qpd.resp.rob}}
    \item{\code{safe_pchisq}}
    \item{\code{simple.kriging.weights}}
    % \item{\code{var.pd.resp.rob}}
    \item{\code{update.zhat}}
}
are not intended for direct use.  However, as any unexported function, they
can be accessed by typing \code{georob:::\var{function-name}}.

In addition, the function \code{K}, which is re-exported by the package
\pkg{georob}, is required for computing block Kriging predictions by the
function \code{f.point.block.cov} of the package \pkg{constrainedKriging}.

Furthermore, the (re-exported) function
\code{\link[robustbase]{lmrob.control}} allows to pass tuning parameters to
the function \code{\link[robustbase]{lmrob}} of the package \pkg{robustbase}.}

\author{
  Andreas Papritz \email{andreas.papritz@env.ethz.ch}\cr
  with contributions by Cornelia Schwierz.
}

\seealso{
  \code{\link{georob}} for (robust) fitting of spatial linear models;

  \code{\link{georobObject}} for a description of the class \code{georob};

  \code{\link{profilelogLik}} for computing profiles of Gaussian likelihoods;

  \code{\link{plot.georob}} for display of (RE)ML variogram estimates;

  \code{\link{control.georob}} for controlling the behaviour of \code{georob};

  \code{\link{georobModelBuilding}} for stepwise building models of class \code{georob};

  \code{\link{cv.georob}} for assessing the goodness of a fit by \code{georob};

  \code{\link{predict.georob}} for computing robust Kriging predictions;

%  \code{\link{validate.predictions}} for validating Kriging predictions;

  \code{\link{lgnpp}} for unbiased back-transformation of Kriging prediction
  of log-transformed data; and finally

  \code{\link{georobMethods}} for further methods for the class \code{georob},

  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
  for robust estimation and modelling of sample variograms.
}

\keyword{internal}

