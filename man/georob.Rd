% 2024-03-09 A. Papritz
% R CMD Rdconv -t html -o bla.html georob.Rd ; open bla.html; R CMD Rd2pdf --force georob.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  georob.Rd

\encoding{utf8}
\name{georob}
\alias{georob}
%
\title{Robust Fitting of Spatial Linear Models}
%
\description{
  The function \code{georob} fits a linear model with spatially correlated
  errors to geostatistical data that are possibly contaminated by
  independent outliers.  The regression coefficients and the parameters of
  the variogram model are estimated by robust or Gaussian restricted
  maximum likelihood (REML) or by Gaussian maximum likelihood (ML).}
\usage{
georob(formula, data, subset, weights, na.action, model = TRUE,
    x = FALSE, y = FALSE, contrasts = NULL, offset, locations,
    variogram.model = c("RMexp", "RMaskey", "RMbessel", "RMcauchy",
        "RMcircular", "RMcubic", "RMdagum", "RMdampedcos", "RMdewijsian",
        "RMfbm", "RMgauss", "RMgencauchy", "RMgenfbm", "RMgengneiting",
        "RMgneiting", "RMlgd", "RMmatern", "RMpenta", "RMqexp",
        "RMspheric", "RMstable", "RMwave", "RMwhittle"),
    param, fit.param = default.fit.param()[names(param)],
	  aniso = default.aniso(), fit.aniso = default.fit.aniso(),
    variogram.object = NULL,
    tuning.psi = 2, control = control.georob(),
    verbose = 0, ...)
}
%
\arguments{
  \item{formula}{a symbolic description of the regression model for the
  external drift to be fit (mandatory argument).  See
  \code{\link[stats]{lm}} and \code{\link[stats]{formula}} for more
  details.}

  \item{data}{an optional data frame, a
  \code{\link[sp]{SpatialPointsDataFrame}}, list or environment (or another
  object coercible by \code{\link{as.data.frame}} to a data frame)
  containing the variables in the model and the coordinates where the data
  was recorded.  If not found in \code{data}, the variables are taken from
  \code{environment(formula)}, typically the environment from which
  \code{georob} is called.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{weights}{an optional vector of weights to be used in the fitting
    process, currently ignored.  }

  \item{na.action}{a function which indicates what should happen when the
    data contain \code{NA}s.  The default is set by the \code{na.action}
    argument of \code{\link{options}}, and is \code{\link{na.fail}} if that is
    unset.  The \dQuote{factory-fresh} default is \code{\link{na.omit}}.
    Another possible value is \code{NULL}, no action.  Value
    \code{\link{na.exclude}} can be useful.}

  \item{model, x, y}{logical scalars.  If \code{TRUE} the corresponding
  components of the fit (the model frame, the model matrix, the response)
  are returned.  The model frame is augmented by the coordinates.}

  \item{contrasts}{an optional list.  See the \code{contrasts.arg} of
    \code{\link{model.matrix.default}}.}

  \item{offset}{this optional argument can be used to specify an \emph{a
  priori} known component to be included in the linear predictor during
  fitting.  An \code{\link{offset}} term can be included in the formula
  instead or as well, and if both are specified their sum is used.}

  \item{locations}{a one-sided formula defining the variables
  that are used as coordinates of the locations were the data was
  recorded (mandatory argument).}

  \item{variogram.model}{a character keyword defining the variogram model
  to be fitted.  Currently, most basic variogram models provided formerly
  by the now archived package \pkg{RandomFields} can be fitted (see
  \emph{Details} and \code{\link{gencorr}}).}

  \item{param}{a named numeric vector with initial values of the
  variogram parameters (mandatory argument).  The names of \code{param}
  are matched against the following names (see \emph{Details} and
  \code{\link{georobPackage}} for information about the parametrization of
  variogram models):

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
      \code{\link{gencorr}} and \code{\link{param.names}}).}

    }
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

      \item{\code{phi}: 90 degrees minus latitude of the first semi-principal axis
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

      \item{\code{variogram.model}}: a mandatory character keyword defining
      the variogram model, see respective argument above.

      \item{\code{param}}: a mandatory named numeric vector with initial
      values of the variogram parameters, see respective argument above.

      \item{\code{fit.param}}: an optional named logical vector defining
      which parameters are adjusted, see respective argument above.

      \item{\code{aniso}}: an optional named numeric vector with initial
      values for fitting geometrically anisotropic variogram models, see
      respective argument above.

      \item{\code{fit.param}}: an optional named logical vector defining
      which anisotropy parameters are adjusted, see respective argument
      above.

  }

  Note that the arguments \code{variogram.model}, \code{param},
  \code{fit.param}, \code{aniso} and\cr \code{fit.aniso} are ignored when
  \code{variogram.object} is passed to \code{georob}. }


  \item{tuning.psi}{positive numeric.  The tuning constant \eqn{c} of the
    \eqn{\psi_c}-function of the robust REML algorithm.}

  \item{control}{a list specifying parameters that control the behaviour of
  \code{georob}.  Use the function \code{\link{control.georob}} and see its
  help page for the components of \code{control}.}

  \item{verbose}{positive integer controlling logging of diagnostic
  messages to the console during model fitting.  \code{verbose = 0} largely
  suppresses such messages and \code{verbose = 4} asks for most verbose
  output (see \code{control} arguments of \code{\link[nleqslv]{nleqslv}},
  \code{\link[stats]{nlminb}} and \code{\link[stats]{optim}} and
  \code{\link{control.georob}} for information how to fine tuning
  diagnostic output generated by \code{nleqslv}, \code{nlminb} and
  \code{optim}).}

  \item{\dots}{further arguments passed to function (e.g. \code{object.}
  used internally for updating \code{georob} objects). }

}

\details{\code{georob} fits a spatial linear model by robust (\cite{\enc{Künsch}{Kuensch} et al., 2011}, \cite{\enc{Künsch}{Kuensch}
  et al., in preparation}) or Gaussian RE(ML) (\cite{Harville, 1977}).
  \code{\link{georobPackage}} describes the employed model and briefly
  sketches the robust REML estimation and the robust external drift Kriging
  method.  Here, we describe further details of \code{georob}.

  \subsection{Implemented variograms}{

    Currently, most basic variogram models provided formerly by the now
    archived package \pkg{RandomFields} can be fitted by \code{georob} (see
    argument \code{variogram.model} and \code{\link{gencorr}} for a list of
    implemented models).  Some of these models have in addition to
    \code{variance}, \code{snugget}, \code{nugget} and \code{scale} further
    parameters.  Initial values of these parameters (\code{param}) and
    fitting flags (\code{fit.param}) must be passed to \code{georob} by the
    same names as used for the models \code{RM...} in
    \code{\link{gencorr}}.  Use the function \code{\link{param.names}} to
    list additional parameters of a given variogram.model.

    The arguments \code{fit.param} and \code{fit.aniso} are used to control
    what variogram and anisotropy parameters are estimated and which are
    kept at the constant initial values.  The functions\cr
    \code{\link{default.fit.param}} and \code{\link{default.fit.aniso}} set
    reasonable default values for these arguments.  Note, as an aside, that
    the function \code{\link{default.aniso}} sets (default) values of the
    anisotropy parameters for an isotropic variogram.

    }

  \subsection{Estimating parameters of power function variogram}{

    The intrinsic variogram model \code{RMfbm} is over-parametrized when
    both the \code{variance} (plus possibly \code{snugget}) and the
    \code{scale} are estimated.  Therefore, to estimate the parameters of
    this model, \code{scale} must be kept fixed at an arbitrary value by
    using \code{fit.param["scale"] = FALSE}.}


  \subsection{Estimating parameters of geometrically anisotropic variograms}{

    The subsection \strong{Model} of \code{\link{georobPackage}} describes
    how such models are parametrized and gives definitions the various
    elements of \code{aniso}.  Some additional remarks might be helpful:

    \itemize{

      \item The first semi-principal axis points into the direction with
      the farthest reaching auto-correlation, which is described by the range
      parameter \code{scale} (\eqn{\alpha}).

      \item The ranges in the direction of the second and third
      semi-principal axes are given by \eqn{f_1\alpha} and \eqn{f_2
      \alpha}, with \eqn{0 < f_2 \leq f_1 \leq 1}{0 < f_2 <= f_1 <= 1}.

      \item The default values for \code{aniso} (\eqn{f_1=1}, \eqn{f_2=1})
        define an isotropic variogram model.

      \item Valid ranges for the angles characterizing the orientation of
      the semi-variance ellipsoid are (in degrees): \eqn{\omega} [0, 180],
      \eqn{\phi} [0, 180], \eqn{\zeta} [-90, 90].  } }


  \subsection{Estimating variance of micro-scale variation}{

    Simultaneous estimation of the variance of the micro-scale variation
    (\code{snugget}, \eqn{\sigma_\mathrm{n}^2}{sigma^2_n}), appears seemingly
    as spatially uncorrelated with a given sampling design, and of the
    variance (\code{nugget}, \eqn{\tau^2}) of the independent errors
    requires that for some locations
    \eqn{\boldsymbol{s}_i}{s_i} replicated observations are
    available.  Locations less or equal than \code{zero.dist} apart are
    thereby considered as being coincident (see
    \code{\link{control.georob}}).

  }


  \subsection{Constraining estimates of variogram parameters}{

    Parameters of variogram models can vary only within certain bounds (see
    \code{\link{param.bounds}} and \code{\link{gencorr}}
    for allowed ranges).  \code{georob} uses three mechanisms to constrain
    parameter estimates to permissible ranges:

    \enumerate{

        \item \emph{Parameter transformations}: By default, all variance
        (\code{variance}, \code{snugget}, \code{nugget}), the range
        \code{scale}, the anisotropy parameters \code{f1} and \code{f2} and
        many of the additional parameters are log-transformed before
        solving the estimating equations or maximizing the restricted
        log-likelihood and this warrants that the estimates are always
        positive (see \code{\link{control.georob}} for detailed explanations
        how to control parameter transformations).

        \item \emph{Checking permissible ranges}: The additional parameters
        of the variogram models such as the smoothness parameter \eqn{\nu}
        of the Whittle-\enc{Matérn}{Matern} model are forced to stay in the permissible
        ranges by signalling an error to \code{nleqslv}, \code{nlminb} or
        \code{optim} if the current trial values are invalid.  These
        functions then graciously update the trial values of the parameters
        and carry their task on.  However, it is clear that such a
        procedure likely gets stuck at a point on the boundary of the
        parameter space and is therefore just a workaround for avoiding
        runtime errors due to invalid parameter values.

        \item \emph{Exploiting the functionality of \code{nlminb} and
        \code{optim}}: If a spatial model is fitted non-robustly, then the
        arguments \code{lower}, \code{upper} (and \code{method} of
        \code{\link[stats]{optim}}) can be used to constrain the parameters
        (see \code{\link{control.optim}} how to pass them to \code{optim}).
        For \code{optim} one has to use the arguments \code{method =
        "L-BFGS-B"}, \code{lower = \var{l}}, \code{upper = \var{u}}, where
        \var{l} and \var{u} are numeric vectors with the lower and upper
        bounds of the \emph{transformed} parameters in the order as they
        appear in\cr \code{c(variance, snugget, nugget, scale,
        \ldots)[fit.param], aniso[fit.aniso])},\cr where \code{\ldots} are
        additional parameters of isotropic variogram models (use \cr
        \code{param.names(variogram.model)} to display the names and the
        order of the additional parameters for \code{variogram.model}).
        For \code{nlminb} one has to use the arguments \code{lower =
        \var{l}}, \code{upper = \var{u}}, where \var{l} and \var{u} are
        numeric vectors as described above.

    }

  }

  \subsection{Computing robust initial estimates of parameters for robust
  REML}{

    To solve the robustified estimating equations for
    \eqn{\boldsymbol{B}}{B} and
    \eqn{\boldsymbol{\beta}}{\beta} the following initial
    estimates are used:

    \itemize{

      \item \eqn{ \widehat{\boldsymbol{B}}=
      \boldsymbol{0},}{hatB=0,} if this turns out to be
      infeasible, initial values can be passed to \code{georob} by the
      argument \code{bhat} of \code{\link{control.georob}}.

      \item \eqn{\widehat{\boldsymbol{\beta}}}{hat\beta} is
      either estimated robustly by the function
      \code{\link[robustbase]{lmrob}}, \code{\link[quantreg]{rq}} or
      non-robustly by \code{\link[stats]{lm}} (see argument
      \code{initial.fixef} of \code{\link{control.georob}}).

    }

    Finding the roots of the robustified estimating equations of the
    variogram and anisotropy parameters is more sensitive to a good choice
    of initial values than maximizing the Gaussian (restricted)
    log-likelihood with respect to the same parameters.  If the initial
    values for \code{param} and \code{aniso} are not sufficiently close to
    the roots of the system of nonlinear equations, then
    \code{\link[nleqslv]{nleqslv}} may fail to find them.
    Setting \code{initial.param = TRUE} (see \code{\link{control.georob}})
    allows one to find initial values that are
    often sufficiently close to the roots so that
    \code{\link[nleqslv]{nleqslv}} converges.  This is achieved by:

    \enumerate{

      \item Initial values of the regression parameters are computed by
        \code{\link[robustbase]{lmrob}} irrespective of the choice for
        \code{initial.fixef} (see \code{\link{control.georob}}).

      \item Observations with \dQuote{robustness weights} of the
      \code{lmrob} fit, satisfying\cr
      \eqn{\psi_c(\widehat{\varepsilon}_i/\widehat{\tau})/(\widehat{\varepsilon}_i/\widehat{\tau})
      \leq \mbox{\code{min.rweight}}}{\psi_c(r_i/s)/(r_i/s)
      <=min.rweight}, are discarded (see
      \code{\link{control.georob}}).

      \item The model is fit to the pruned data set by Gaussian REML using
        \code{\link[stats]{nlminb}} or \code{\link[stats]{optim}}.

      \item The resulting estimates of the variogram parameters
        (\code{param}, \code{aniso}) are used as initial estimates for the
        subsequent robust fit of the model by \code{\link[nleqslv]{nleqslv}}.

    }

    Note that for step 3 above, initial values of \code{param} and
    \code{aniso} must be provided to \code{georob}.

  }

  \subsection{Estimating variance parameters by Gaussian (RE)ML}{

    Unlike robust REML, where robustified estimating equations are solved
    for the variance parameters \code{nugget} (\eqn{\tau^2}),
    \code{variance} (\eqn{\sigma^2}), and possibly \code{snugget}
    (\eqn{\sigma_{\mathrm{n}}^2}{\sigma_n^2}), for Gaussian (RE)ML the
    variances can be re-parametrized to

    \itemize{
      \item the signal variance
      \eqn{\sigma_B^2 = \sigma^2 + \sigma_{\mathrm{n}}^2}{\sigma_B^2 = \sigma^2+\sigma_n^2},

      \item the inverse relative nugget
      \eqn{\eta = \sigma_B^2 / \tau^2} and

      \item the relative auto-correlated signal variance
      \eqn{\xi = \sigma^2/\sigma_B^2}{\xi = \sigma^2 / \sigma_B^2}.
    }

    \code{georob} maximizes then a (restricted) \emph{profile
    log-likelihood} that depends only on \eqn{\eta}, \eqn{\xi},
    \eqn{\alpha}, \ldots, and \eqn{\sigma_B^2} is estimated by an explicit
    expression that depends on these parameters (e.g. \cite{Diggle and
    Ribeiro, 2006, p.  113}).  This is usually more efficient than
    maximizing the (restricted) log-likelihood with respect to the original
    variance parameters \eqn{\tau^2},
    \eqn{\sigma_{\mathrm{n}}^2}{\sigma_n^2} and \eqn{\sigma^2}.
    \code{georob} chooses the parametrization automatically, but the user
    can control it by the argument \code{reparam} of the function
    \code{\link{control.georob}}.

  }
}

\value{

An object of class \code{georob} representing a robust (or Gaussian) (RE)ML
fit of a spatial linear model.  See
\code{\link{georobObject}} for the components of the fit.

}

\references{
  Diggle, P. J. and Ribeiro, P. J. R. (2006) Model-based Geostatistics.
  Springer, New York, \doi{10.1007/978-0-387-48536-2}.

  Harville, D. A. (1977) Maximum likelihood approaches to variance
  component estimation and to related problems, \emph{Journal of the
  American Statistical Association}, \bold{72}, 320--340,
  \doi{10.1080/01621459.1977.10480998}.

  \enc{Künsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A. (in
  preparation) Robust Geostatistics.

  \enc{Künsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A.
  (2011) Robust estimation of the external drift and the variogram of
  spatial data.  Proceedings of the ISI 58th World Statistics Congress of
  the International Statistical Institute.  \doi{10.3929/ethz-a-009900710}
  }

\author{
  Andreas Papritz \email{papritz@retired.ethz.ch}\cr
  with contributions by Cornelia Schwierz.
}

\seealso{
  \code{\link{georobPackage}} for a description of the model and a brief summary of the algorithms;

%  \code{\link{georob}} for (robust) fitting of spatial linear models;
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
  from model fitted by \code{georob}; and finally

  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
  for robust estimation and modelling of sample variograms.
}

\examples{
################
## meuse data ##
################
data(meuse)

## Gaussian REML fit
r.logzn.reml <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "RMexp",
    param = c(variance = 0.15, nugget = 0.05, scale = 200),
    tuning.psi = 1000)
summary(r.logzn.reml, correlation = TRUE)
plot(r.logzn.reml, lag.dist.def = seq(0, 2000, by = 100))

## robust REML fit
if(interactive()){
  ## example is run only in interactive session because cpu times exceeds 5 s
  r.logzn.rob <- update(r.logzn.reml, tuning.psi = 1)

  summary(r.logzn.rob, correlation = TRUE)
  lines(r.logzn.rob, col = "red")
}


###################
## wolfcamp data ##
###################
data(wolfcamp)

## fitting isotropic IRF(0) model

r.irf0.iso <- georob(pressure ~ 1, data = wolfcamp, locations = ~ x + y,
    variogram.model = "RMfbm",
    param = c(variance = 10, nugget = 1500, scale = 1, alpha = 1.5),
    fit.param = default.fit.param(scale = FALSE, alpha = TRUE),
    tuning.psi = 1000)

summary(r.irf0.iso)
plot(r.irf0.iso, lag.dist.def = seq(0, 200, by = 7.5))

## fitting anisotropic IRF(0) model
if(interactive()){
  ## example is run only in interactive session because cpu times exceeds 5 s
  r.irf0.aniso <- georob(pressure ~ 1, data = wolfcamp, locations = ~ x + y,
      variogram.model = "RMfbm",
      param = c(variance = 5.9, nugget = 1450, scale = 1, alpha = 1),
      fit.param = default.fit.param(scale = FALSE, alpha = TRUE),
      aniso = default.aniso(f1 = 0.51, omega = 148.),
      fit.aniso = default.fit.aniso(f1 = TRUE, omega = TRUE),
      tuning.psi = 1000)
  summary(r.irf0.aniso)
  plot(r.irf0.aniso, lag.dist.def = seq(0, 200, by = 7.5),
      xy.angle.def = c(0, 22.5, 67.5, 112.5, 157.5, 180.),
      add = TRUE, col = 2:5)

  pchisq(2*(r.irf0.aniso[["loglik"]] - r.irf0.iso[["loglik"]]), 2, lower = FALSE)
}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
