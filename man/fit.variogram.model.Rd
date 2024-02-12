% 2024-02-9 A. Papritz
% R CMD Rdconv -t html -o bla.html fit.variogram.model.Rd ; open bla.html; R CMD Rd2pdf --force fit.variogram.model.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  fit.variogram.model.Rd

\encoding{utf8}
\name{fit.variogram.model}
\alias{control.fit.variogram.model}
\alias{fit.variogram.model}
\alias{print.fitted.variogram}
\alias{summary.fitted.variogram}
\alias{print.summary.fitted.variogram}
\alias{lines.fitted.variogram}

\title{Fitting Model Functions to Sample Variograms}

\description{
The function \code{fit.variogram.model} fits a variogram model to a sample
variogram by (weighted) non-linear least squares.  The function
\code{control.fit.variogram.model} generates a list with steering parameters which
control \code{fit.variogram.model}.  There are \code{print}, \code{summary}
and \code{lines} methods for summarizing and displaying fitted variogram
models.}

\usage{
fit.variogram.model(sv,
    variogram.model = c("RMexp", "RMaskey", "RMbessel", "RMcauchy",
        "RMcircular", "RMcubic", "RMdagum", "RMdampedcos", "RMdewijsian",
        "RMfbm", "RMgauss", "RMgencauchy", "RMgenfbm", "RMgengneiting",
        "RMgneiting", "RMlgd", "RMmatern", "RMpenta", "RMqexp",
        "RMspheric", "RMstable", "RMwave", "RMwhittle"),
    param, fit.param = default.fit.param()[names(param)],
	  aniso = default.aniso(), fit.aniso = default.fit.aniso(),
    variogram.object = NULL,
    max.lag = max(sv[["lag.dist"]]), min.npairs = 30,
    weighting.method = c("cressie", "equal", "npairs"),
    control = control.fit.variogram.model(),
    verbose = 0)

control.fit.variogram.model(maximizer = c("nlminb", "optim"),
    param.tf = param.transf(), fwd.tf = fwd.transf(),
    deriv.fwd.tf = dfwd.transf(), bwd.tf = bwd.transf(),
    hessian = TRUE, optim = control.optim(), nlminb = control.nlminb())

\method{print}{fitted.variogram}(x, digits = max(3, getOption("digits") - 3), ...)

\method{summary}{fitted.variogram}(object, correlation = FALSE, signif = 0.95, ...)

\method{lines}{fitted.variogram}(x, what = c("variogram", "covariance", "correlation"),
    from = 1.e-6, to, n = 501, xy.angle = 90, xz.angle = 90,
    col = 1:length(xy.angle), pch = 1:length(xz.angle), lty = "solid", ...)

}

\arguments{
  \item{sv}{an object of class \code{sample.variogram}, see
  \code{\link{sample.variogram}}.}

  \item{variogram.model}{a character keyword defining the variogram model
  to be fitted.  Currently, most basic variogram models provided formerly
  by the now archived package \pkg{RandomFields} can be fitted (see
  \emph{Details} of \code{\link{georob}} and \code{\link{gencorr}}).}

  \item{param}{a named numeric vector with initial values of the variogram
  parameters.  The following parameter names are allowed (see
  \emph{Details} of \code{\link{georob}} and \code{\link{georobPackage}} for
  information about the parametrization of variogram models):

  \itemize{

    \item{\code{variance}: variance (sill \eqn{\sigma^2}) of the
    auto-correlated component of the Gaussian random field
    \eqn{B(\boldsymbol{s})}{B(s)}.}

    \item{\code{snugget}: variance
    (spatial nugget \eqn{\sigma^2_{\mathrm{n}}}{\sigma^2_n})
    of the seemingly spatially uncorrelated component of
    \eqn{B(\boldsymbol{s})}{B(s)}
    (micro-scale spatial variation; default value\cr  \code{snugget = 0}).}

    \item{\code{nugget}: variance (nugget \eqn{\tau^2}) of the
    independent errors
    \eqn{\varepsilon(\boldsymbol{s})}{\epsilon(s)}.}

    \item{\code{scale}: range parameter (\eqn{\alpha}) of the variogram.}

    \item{names of additional variogram parameters such as the smoothness
    parameter \eqn{\nu} of the Whittle-\enc{Matérn}{Matern} model (see
    \code{\link{gencorr}} and \code{\link{param.names}}).} }

  }

  \item{fit.param}{a named logical vector (or a function such as
  \code{\link{default.fit.param}} that creates this vector) with the same
  names as used for \code{param}, defining which parameters are adjusted
  (\code{TRUE}) and which are kept fixed at their initial values
  (\code{FALSE}) when fitting the model.}

  \item{aniso}{a named numeric vector with initial values (or a function such as
  \code{\link{default.aniso}} that creates this vector) for fitting
  geometrically anisotropic variogram models.  The names of \code{aniso}
  are matched against the following names (see \emph{Details} and
  \code{\link{georobPackage}} for information about the parametrization of
  variogram models):

  \itemize{

      \item{\code{f1}: ratio \eqn{f_1} of lengths of second and first
        semi-principal axes of an ellipsoidal surface with constant
        semi-variance in \eqn{\mathrm{I}\!\mathrm{R}^3}{R^3} (default \code{f1 = 1}).}

      \item{\code{f2}: ratio \eqn{f_2} of lengths of third and first
        semi-principal axes of the semi-variance ellipsoid (default \code{f2 = 1}).}

      \item{\code{omega}: azimuth in degrees of the first semi-principal axis
        of the semi-variance ellipsoid (default \code{omega = 90}).}

      \item{\code{phi}: 90 degrees minus altitude of the first semi-principal axis
        of the semi-variance ellipsoid (default \code{phi = 90}).}

      \item{\code{zeta}: angle in degrees between the second semi-principal
        axis and the direction of the line defined by the intersection
        between the \eqn{x}-\eqn{y}-plane and the plane orthogonal to the
        first semi-principal axis of the semi-variance ellipsoid through the
        origin (default \code{zeta = 0}).}
    }
  }

  \item{fit.aniso}{a named logical vector (or a function such as
  \code{\link{default.fit.aniso}} that creates this vector) with the same
  names as used for \code{aniso}, defining which parameters are adjusted
  (\code{TRUE}) and which are kept fixed at their initial values
  (\code{FALSE}) when fitting the model.}

  \item{variogram.object}{an optional list that defines a possibly nested
  variogram model.  Each component is itself a list with the following
  components:

  \itemize{

      \item{\code{variogram.model}}: a character keyword defining the variogram
      model, see respective argument above.

      \item{\code{param}}: a named numeric vector with initial values of the
      variogram parameters, see respective argument above.

      \item{\code{fit.param}}: a named logical vector defining which parameters are
      adjusted, see respective argument above.

      \item{\code{aniso}}: a named numeric vector with initial values for fitting
      geometrically anisotropic variogram models, see respective argument
      above.

      \item{\code{fit.param}}: a named logical vector defining which anisotropy
      parameters are adjusted, see respective argument above.

  }

  Note that the arguments \code{variogram.model}, \code{param},
  \code{fit.param}, \code{aniso} and\cr \code{fit.aniso} are ignored when
  \code{variogram.object} is passed to \code{fit.variogram.model}. }

  \item{max.lag}{a positive numeric defining the maximum lag distance to be
  used for fitting or plotting variogram models (default all lag classes).}

  \item{min.npairs}{a positive integer defining the minimum number of data
  pairs required so that a lag class is used for fitting a variogram
  model (default \code{30}).}

  \item{weighting.method}{a character keyword defining the weights for
  non-linear least squares.  Possible values are:

    \itemize{

      \item \code{"equal"}: no weighting ,

      \item \code{"npairs"}: weighting by number of data pairs in a lag class,

      \item \code{"cressie"}: \dQuote{Cressie's weights} (default, see
      \cite{Cressie, 1993, sec.  2.6.2}).

    }
  }

  \item{verbose}{a positive integer controlling logging of diagnostic
  messages to the console during model fitting.}

  \item{control}{a list with the components \code{maximizer},
  \code{param.tf}, \code{fwd.tf}, \code{bwd.tf}, \code{hessian},
  \code{optim} and \code{nlminb}  or a function such as
  \code{control.fit.variogram.model} that generates such a list.
  See \code{\link{control.georob}} for information on
  \code{maximizer}, \code{param.tf}, \code{fwd.tf},
  \code{bwd.tf}, \code{hessian}, \code{optim} and \code{nlminb}.}

  \item{maximizer}{a character keyword defining the optimizer for nonlinear
  least squares.  Possible values are \code{\link[stats]{nlminb}} (default)
  or \code{\link[stats]{optim}}.}

  \item{hessian}{a logical scalar controlling whether the Hessian should be
  computed at the nonlinear least squares estimates.}

  \item{param.tf}{a function such as \code{\link{param.transf}}, which
  returns a named vector of character strings that define the
  transformations to be applied to the variogram parameters for model
  fitting, see \code{\link{control.georob}}.}

  \item{fwd.tf}{a function such as \code{\link{fwd.transf}}, which returns
  a named list of invertible functions to be used to transform variogram
  parameters, see \code{\link{control.georob}}.}

  \item{deriv.fwd.tf}{a function such as \code{dfwd.transf}, which
  returns a named list of functions corresponding to the first derivatives
  of \code{fwd.tf}, see \code{\link{control.georob}}.}

  \item{bwd.tf}{a function such as \code{\link{bwd.transf}}, which returns
  the named list of inverse functions corresponding to \code{fwd.tf}, see
  see \code{\link{control.georob}}.}

  \item{nlminb}{a list of arguments passed to \code{nlminb} or a function
  such as \code{\link{control.nlminb}} that generates such a list (see
  \code{\link[stats]{nlminb}} for allowed arguments).}

  \item{optim}{a list of arguments passed to \code{optim} or a function
  such as \code{\link{control.optim}} that generates such a list (see
  \code{\link[stats]{optim}} for allowed arguments).}

  \item{object, x}{an object of class \code{fitted.variogram}.}

  \item{digits}{a positive integer indicating the number of decimal digits
  to print.}

  \item{correlation}{a logical scalar controlling whether the correlation matrix of
  the fitted variogram parameters is computed (default \code{FALSE}).}

  \item{signif}{a numeric with the confidence level for computing
  confidence intervals for variogram parameters (default \code{0.95}).}

  \item{what}{a character keyword with the quantity that should be
  displayed (default \code{"variogram"}).}

  \item{from}{a numeric with the minimal lag distance used in plotting
  variogram models.}

  \item{to}{a numeric with the maximum lag distance used in plotting
  variogram models (default: largest lag distance of current plot).}

  \item{n}{a positive integer specifying the number of equally spaced lag
  distances for which semi-variances are evaluated in plotting variogram
  models (default \code{501}).}

  \item{xy.angle}{a numeric vector with azimuth angles (in degrees,
  clockwise positive from north) in \eqn{x}-\eqn{y}-plane for which
  semi-variances should be plotted.}

  \item{xz.angle}{a numeric vector with angles in \eqn{x}-\eqn{z}-plane (in
  degrees, clockwise positive from zenith to south) for which semi-variances
  should be plotted.}

  \item{col}{a vector with colours of curves to distinguish curves relating
  to different azimuth angles in \eqn{x}-\eqn{y}-plane.}

  \item{pch}{a vector with the plotting symbols added to lines to
  distinguish curves relating to different angles in
  \eqn{x}-\eqn{z}-plane.}

  \item{lty}{a vector with the line types for plotting variogram models.}

  \item{\dots}{additional arguments passed to methods.}

}
\details{

  The parametrization of geometrically anisotropic variograms is
  described in detail in \code{\link{georobPackage}}, and the section
  \emph{Details} of \code{\link{georob}} describes how the parameter
  estimates are constrained to permissible ranges.  The same
  mechanisms are used in \code{fit.variogram.model}.

  The method \code{summary} computes confidence intervals of the estimated
  variogram and anisotropy parameters from the Hessian matrix of the residual
  sums of squares, based on the asymptotic normal distribution of least
  squares estimates.  Note that the Hessian matrix with respect to the
  \emph{transformed} variogram and anisotropy parameters is used for this.
  Hence the inverse Hessian matrix is the covariance matrix of the
  transformed parameters, confidence intervals are first computed for the
  transformed parameters and the limits of these intervals are transformed
  back to the original scale of the parameters.  Optionally, \code{summary}
  reports the correlation matrix of the \emph{transformed} parameters, also
  computed from the Hessian matrix.

}


\value{

  The function \code{fit.variogram.model} generates an object of class
  \code{fitted.variogram} which is a list with the following components:

  \item{sse}{the value of the object function (weighted residual sum of
  squares) evaluated at the solution.}

  \item{variogram.object}{the estimated parameters of a possibly nested
  variograms model.  This is a list that contains for each variogram model
  structure the following components:

  \itemize{

      \item{\code{variogram.model}}: the name of the fitted parametric
      variogram model.

      \item{\code{param}}: a named numeric vector with the (estimated)
      variogram parameters.

      \item{\code{fit.param}}: a named logical vector with the flags
      defining what variogram parameters were estimated.

      \item{\code{isotropic}: logical indicating whether an isotropic
      variogram was fitted.}

      \item{\code{aniso}: a named numeric vector with the (estimated)
      anisotropy parameters.}

      \item{\code{fit.aniso}}: a named logical vector with the flags
      defining what anisotropy parameters were estimated.

      \item{\code{sincos}: a list with \code{sin} and \code{cos} of the
      angles \eqn{\omega}, \eqn{\phi} and \eqn{\zeta} that define the
      orientation of the anisotropy ellipsoid.}

      \item{\code{rotmat}: the matrix \eqn{(\boldsymbol{C}_1,
      \boldsymbol{C}_2, \boldsymbol{C}_3)}{(C_1, C_2, C_3)} (see
      \code{\link{georobPackage}}).}

      \item{\code{sclmat}: a vector with the elements 1, \eqn{1/f_1},
      \eqn{1/f_2} (see \code{\link{georobPackage}}).  }

    }

  }


  \item{param.tf}{a character vector listing the transformations of the
  variogram parameters used for model fitting.}

  \item{fwd.tf}{a list of functions for variogram parameter transformations.}

  \item{bwd.tf}{a list of functions for \emph{inverse} variogram parameter
  transformations.}

  \item{converged}{a logical scalar indicating whether numerical
  maximization by \code{\link[stats]{nlminb}} or \code{\link[stats]{optim}}
  converged.}

  \item{convergence.code}{a diagnostic integer issued by
  \code{\link[stats]{nlminb}} or \code{\link[stats]{optim}} (component
  \code{convergence}) about convergence.}

  \item{iter}{a named integer vector of length two with the number of
  function and gradient evaluations by \code{\link[stats]{nlminb}} or
  \code{\link[stats]{optim}}.}

  \item{call}{the matched call.}

  \item{residuals}{a numeric vector with the residuals, that is the sample
  semi-variance minus the fitted values.}

  \item{fitted}{a numeric vector with the modelled semi-variances.}

  \item{weights}{a numeric vector with the weights used for fitting.}

  \item{hessian.tfpa}{a symmetric matrix with the Hessian at the solution
  with respect to the transformed variogram and anisotropy parameters
  (missing if \code{hessian} is false).  This Hessian is used by
  \code{summary.fitted.variogram} to compute confidence intervals
  for the estimated parameters.}

  \item{hessian.ntfpa}{a symmetric matrix with the Hessian at the solution
  with respect to the non-transformed variogram and anisotropy parameters
  (missing if \code{hessian} is false).
  }

  The function \code{control.fit.variogram.model} returns a list with
  parameters to steer \cr \code{fit.variogram.model}, see arguments of
  the function above for its components.

  The method \code{print.fitted.variogram} invisibly returns the fitted
  variogram model unchanged.

  The method \code{summary.fitted.variogram} returns an object of class
  \code{summary.fitted.variogram} which is a list containing a subset of
  the components of the fitted variogram object (\code{call},
  \code{residuals}, \code{weights}, \code{converged},
  \code{convergence.code}, \code{iter}, \code{sse},
  \code{variogram.object}), the matrix \code{param.aniso} with the
  estimated values of the variogram parameters along with the bounds of the
  confidence intervals and optionally the correlation matrix
  \code{cor.tf.param} of the estimated transformed parameters.  There is a
  \code{print} method for objects of class \code{summary.fitted.variogram}
  which returns invisibly the object unchanged.

  The method \code{lines.fitted.variogram} is called for its side effects
  and returns the value \code{NULL} invisibly.


}

\references{

  Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}, Wiley, New
  York, \doi{10.1002/9781119115151}.

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

  \code{\link{georobMethods}} for further methods for the class \code{georob};

  \code{\link{predict.georob}} for computing robust Kriging predictions;
%  \code{\link{validate.predictions}} for validating Kriging predictions;

  \code{\link{lgnpp}} for unbiased back-transformation of Kriging prediction
  of log-transformed data;

  \code{\link{georobSimulation}} for simulating realizations of a Gaussian process
  from model fitted by \code{georob}.
%  ; and finally
%  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
%  for robust estimation and modelling of sample variograms.
}

\examples{
data(wolfcamp)

## fitting an isotropic IRF(0) model
r.sv.iso <- sample.variogram(pressure~1, data = wolfcamp,
    locations = ~x + y, lag.dist.def = seq(0, 200, by = 15))
plot(r.sv.iso, type = "l")

if(interactive()){
  ## example is run only in interactive session because cpu times exceeds 5 s
  r.irf0.iso <- fit.variogram.model(r.sv.iso, variogram.model = "RMfbm",
      param = c(variance = 100, nugget = 1000, scale = 1., alpha = 1.),
      fit.param = default.fit.param(scale = FALSE, alpha = TRUE))
  summary(r.irf0.iso, correlation = TRUE)
  lines(r.irf0.iso, line.col = "red")
}

## fitting an anisotropic IRF(0) model
r.sv.aniso <- sample.variogram(pressure~1, data = wolfcamp,
    locations = ~x + y, lag.dist.def = seq(0, 200, by = 15),
    xy.angle.def = c(0., 22.5, 67.5, 112.5, 157.5, 180.))
plot(r.sv.aniso, type = "l")

if(interactive()){
  ## example is run only in interactive session because cpu times exceeds 5 s
  r.irf0.aniso <- fit.variogram.model(r.sv.aniso, variogram.model = "RMfbm",
      param = c(variance = 100, nugget = 1000, scale = 1., alpha = 1.5),
      fit.param = default.fit.param(scale = FALSE, alpha = TRUE),
      aniso = default.aniso(f1 = 0.4, omega = 135.),
      fit.aniso = default.fit.aniso(f1 = TRUE, omega = TRUE),
      control = control.fit.variogram.model(
        maximizer = "optim",
        optim = control.optim(
          method = "BFGS", hessian = TRUE, control = list(maxit = 5000)
        )
      ))
  summary(r.irf0.aniso, correlation = TRUE)

  lines(r.irf0.aniso, xy.angle = seq(0, 135, by = 45))
}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
