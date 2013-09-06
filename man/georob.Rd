% 2013-09-06 A. Papritz
% R CMD Rdconv -t html -o bla.html georob.Rd ; open bla.html; R CMD Rd2pdf --force georob.Rd; 

\encoding{macintosh}
\name{georob}
\alias{georob}
%
\title{Robust Fitting of Spatial Linear Models}
%
\description{
  The function \code{georob} fits a linear model with spatially correlated
  errors to geostatistical data that are possibly contaminated by
  independent outliers.  The regression coefficients and the parameters of
  the variogram model are estimated by robust or Gaussian Restricted
  Maximum Likelihood (REML).  }
\usage{
georob(formula, data, subset, weights, na.action, model = TRUE, 
    x = FALSE, y = FALSE, contrasts = NULL, offset, locations, 
    variogram.model = c("exponential", "bessel", "cauchy", "cauchytbm", 
        "circular", "cubic", "dagum", "dampedcosine", "DeWijsian", "fractalB", 
        "gauss", "genB", "gencauchy", "gengneiting", "gneiting", "lgd1",
        "matern", "penta", "power", "qexponential", "spherical", "stable",
        "wave", "whittle"), 
    param, 
    fit.param = c(variance = TRUE, snugget = FALSE, nugget = TRUE, 
        scale = TRUE, a = FALSE, alpha = FALSE, beta = FALSE, delta = FALSE, 
        gamma = FALSE, lambda = FALSE, n = FALSE, nu = FALSE)[names(param)], 
    aniso = c(f1 = 1, f2 = 1, omega = 90, phi = 90, zeta = 0), 
    fit.aniso = c(f1 = FALSE, f2 = FALSE, omega = FALSE, 
        phi = FALSE, zeta = FALSE), 
    tuning.psi = 2, initial.param = c("exclude", "minimize", "no"), 
%     root.finding = c("nleqslv", "bbsolve"),
    control = georob.control(...),
    verbose = 0, ...)
}
%
\arguments{
  \item{formula}{a symbolic description of the regression model to be fit.  See
    \code{\link[stats]{lm}} and \code{\link[stats]{formula}} for more details.}
    
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
    
  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components
    of the fit (the model frame, the model matrix, the response) are
    returned. The model frame is augmented by the coordinates.}
    
  \item{contrasts}{an optional list.  See the \code{contrasts.arg} of
    \code{\link{model.matrix.default}}.}
    
  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.  An \code{\link{offset}} term can be included in the
    formula instead or as well, and if both are specified their sum is used.}
    
  \item{locations}{a one-sided formula defining the variables that are used
    as coordinates of the locations were the data was recorded.}
    
  \item{variogram.model}{a character keyword defining the variogram model
    to be fitted.  Currently, most basic variogram models
    provided by the package \pkg{RandomFields} can be fitted (see \emph{Details} 
    and \code{\link[RandomFields]{Variogram}}).}
  
  \item{param}{a named numeric vector with initial values of the variogram
    parameters. The following parameter names are allowed 
    (see \emph{Details} and \code{\link{georobIntro}} for information about the 
    parametrization of variogram models): 
    
    \itemize{
      \item{\code{variance}: variance (sill \eqn{\sigma^2}) of the
      auto-correlated component of the Gaussian random field 
      \eqn{B(\mbox{\boldmath$s$\unboldmath})}{B(s)}.}
      
      \item{\code{snugget}: variance 
      (spatial nugget \eqn{\sigma^2_{\mathrm{n}}}{\sigma^2_n})
      of the seemingly spatially uncorrelated component of 
      \eqn{B(\mbox{\boldmath$s$\unboldmath})}{B(s)} 
      (micro-scale spatial variation; default value\cr  \code{snugget = 0}).}
      
      \item{\code{nugget}: variance (nugget \eqn{\tau^2}) of the
      independent errors
      \eqn{\varepsilon(\mbox{\boldmath$s$\unboldmath})}{\epsilon(s)}.}
     
      \item{\code{scale}: range parameter (\eqn{\alpha}) of the variogram.}
     
      \item{names of additional variogram parameters such as the smoothness
      parameter \eqn{\nu} of the Whittle-MatŽrn model as passed by the
      argument \code{param} to \code{\link[RandomFields]{Variogram}}, see
      \code{\link{param.names}}.}
       
    }
  }
  \item{fit.param}{a named logical vector with the same names as used for 
    \code{param}, defining which parameters are adjusted (\code{TRUE}) and which 
    are kept fixed at their initial values (\code{FALSE}) when fitting the model.}
    
  \item{aniso}{a named numeric vector with initial values for fitting
    geometrically anisotropic variogram models. The following parameter names are allowed 
    (see \emph{Details} and \code{\link{georobIntro}} for information about the 
    parametrization of variogram models): 
    
    \itemize{
    
      \item{\code{f1}: ratio \eqn{f_1} of lengths of second and first
        semi-principal axes of an ellipsoidal surface with constant
        semivariance in \eqn{\mathrm{I}\!\mathrm{R}^3}{R^3} (default \code{f1 = 1}).}
      
      \item{\code{f2}: ratio \eqn{f_2} of lengths of third and first
        semi-principal axes of the semivariance ellipsoid (default \code{f2 = 1}).}
      
      \item{\code{omega}: azimuth in degrees of the first semi-principal axis 
        of the semivariance ellipsoid (default \code{omega = 90}).}
      
      \item{\code{phi}: 90 degrees minus altitude of the first semi-principal axis 
        of the semivariance ellipsoid (default \code{phi = 90}).}
     
      \item{\code{zeta}: angle in degrees between the second semi-principal
        axis and the direction of the line defined by the intersection
        between the \eqn{x}-\eqn{y}-plane and the plane orthogonal to the
        first semi-principal axis of the semivariance ellipsoid through the
        origin (default \code{zeta = 0}).} 
    } 
  }

  \item{fit.aniso}{a named logical vector with the same names as used for 
    \code{aniso}, defining which parameters are adjusted (\code{TRUE}) and which 
    are kept fixed at their initial values (\code{FALSE}) when fitting the model.}
    
  \item{tuning.psi}{positive numeric.  The tuning constant \eqn{c} of the
    \eqn{\psi_c}-function of the robust REML algorithm.}
    
  \item{initial.param}{character, controlling whether initial values of
  parameters are computed for solving the estimating equations of the
  variogram and anisotropy parameters.  
  
  If \code{initial.param = "minimize"} (defaullt) robust initial values are
  computed by minimizing the sum of the squared robustified estimating
  equations using \code{\link[stats]{optim}} (see \emph{Details}).
  If \code{initial.param = "exclude"} robust initial values of parameters are
  computed by discarding outlying observations based on the
  \dQuote{robustness weights} of the initial fit of the regression model by
  \code{\link[robustbase]{lmrob}} and fitting the spatial linear model by
  Gaussian REML to the pruned data set (see \emph{Details}).  
  For \code{initial.param = "no"} no initial parameter values are computed
  and the estimating equations are solved with the initial values passed by
  \code{param} and \code{aniso} to \code{georob}.}
  
%   \item{root.finding}{character string defining what solver is used 
%   for solving the system of nonlinear estimating equations.  Possible
%   values are:
%   
%     \itemize{
%       \item{\code{nleqslv} (default), which invokes
%       \code{\link[nleqslv]{nleqslv}} to solve the estimating equations of the
%       variogram and anisotropy parameters in combination with iterated
%       re-weighted least squares for solving the estimating equations of 
%       \eqn{ \widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} and 
%       \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta}.}
%       
%       \item{\code{bbsolve} which invokes \code{\link[BB]{BBsolve}} for solving
%       the estimating equations of \eqn{ \widehat{\mbox{\boldmath$\xi$\unboldmath}}}{hatxi} 
%       and of the variogram and anisotropy parameters directly.
%       }
%     }
%   
%   }
    
  \item{control}{a list specifying parameters that control the behaviour of
  \code{georob}.  Use the function \code{\link{georob.control}} and see its
  help page for the components of \code{control}.}
  
  \item{verbose}{positive integer controlling logging of diagnostic
  messages to the console during model fitting.  \code{verbose = 0} largely
  suppresses such messages and \code{verbose = 4} asks for most verbose
  output (see \code{control} arguments of \code{\link[nleqslv]{nleqslv}}
  and \code{\link[stats]{optim}} and \code{\link{georob.control}} for
  information how to fine tuning diagnostic output generated by
  \code{nleqslv} and \code{optim}).}
    
  \item{\dots}{can be used to specify list components of \code{control}
  directly.}
  
}

\details{\code{georob} fits a spatial linear model by robust or Gaussian REML 
  (K?nsch et al., 2011, K?nsch et al., in preparation). \code{\link{georobIntro}} 
  describes the employed model and briefly sketches the robust REML estimation 
  and the robust external-drift kriging method. Here, we describe further details 
  of \code{georob}.
  
  \subsection{Implemented variogram models}{
  
    Currently, all except one (model \code{hyperbolic}) basic variogram
    models of the package \pkg{RandomFields}, valid in
    \eqn{\mathrm{I}\!\mathrm{R}^d}{R^d}, \eqn{d > 1}, can be fitted by
    \code{georob} (see argument \code{variogram.models} for a list of
    implemented models).  Some of these models have in addition to
    \code{variance}, \code{snugget}, \code{nugget} and \code{scale} further
    parameters.  Initial values of these parameters (\code{param}) and
    fitting flags (\code{fit.param}) must be passed to \code{georob} by the
    \emph{same names as used for the argument} \code{param} of
    \code{\link[RandomFields]{Variogram}}.  Use the function
    \code{\link{param.names}} to list additional parameters of a given
    variogram.model.
    
    }
    
  \subsection{Estimation of variance of micro-scale variation}{
  
    Simultaneous estimation of the variance of the micro-scale variation
    (\code{snugget}, \eqn{\sigma_\mathrm{n}^2}{sigma^2_n}), which appears
    as seemingly uncorrelated with a given sampling design, and of the
    variance (\code{nugget}, \eqn{\tau^2}) of the independent errors
    requires that for some locations
    \eqn{\mbox{\boldmath$s$\unboldmath}_i}{s_i} replicated observations are
    available.  Locations less or equal than \code{zero.dist} apart are
    thereby considered as being coincient (see
    \code{\link{georob.control}}).
    
  }
    
  \subsection{Fitting intrinsic variogram models}{ 
  
    The intrinsic variogram model \code{fractalB} is overparametrized when
    both the \code{variance} (plus possibly \code{snugget}) and the
    \code{scale} are fitted.  Therefore, to estimate the parameters of this
    model \code{scale} must be kept fixed at an arbitray value by using
    \code{fit.param["scale"] = FALSE}.}
    
  \subsection{Fitting geometrically anisotropic variogram models}{
  
    The subsection \strong{Model} of \code{\link{georobIntro}} describes
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
      the semivariance ellipsoid are (in degrees): \eqn{\omega} [0, 180],
      \eqn{\phi} [0, 180], \eqn{\zeta} [-90, 90].  } }
  
  \subsection{Constraining estimates of variogram parameters}{
  
    Parameters of variogram models can vary only within certain bounds (see
    \code{\link{param.bounds}} and \code{\link[RandomFields]{Variogram}}
    for allowed ranges).  \code{georob} uses three mechanisms to constrain
    parameter estimates to permissible ranges:
    
    \enumerate{ 
      
      \item \emph{Parameter transformations}: By default, all variance
        (\code{variance}, \code{snugget}, \code{nugget}), the range
        \code{scale} and the anisotropy parameters \code{f1} and
        \code{f2} are log-transformed before solving the estimating
        equations or maximizing the restricted loglikelihood and this
        warrants that the estimates are always positive (see
        \code{\link{georob.control}} for controlling parameter
        transformations).
      
        \item \emph{Checking permissible ranges}: The additional parameters
        of the variogram models such as the smoothness parameter \eqn{\nu}
        of the Whittle-Mat?rn model are forced to stay in the permissible
        ranges by signalling an error to \code{nleqslv} or \code{optim} if
        the current trial values are invalid. These functions then
        graciously update the trial values of the parameters and carry
        their task on.  However, it is clear that such a procedure likely
        gets stuck at a point on the boundary of the parameter space and is
        therefore just a workaround for avoiding runtime errors due to
        invalid parameter values.
        
        \item \emph{Exploiting the functionality of \code{optim}}: If a
        spatial model is fitted non-robustly, then the arguments
        \code{lower}, \code{upper} and \code{method} of \code{\link[stats]{optim}}
        can be used to constrain the parameters (see
        \code{\link{optim.control}} how to pass them to \code{optim}).  To
        achieve this one has to use the arguments \code{method = "L-BFGS-B"},
        \code{lower = \var{l}}, \code{upper = \var{u}}, where \var{l} and
        \var{u} are numeric vectors with the lower and upper bounds of the
        \emph{transformed} parameters in the order as they appear in\cr
        \code{c( c(variance, snugget, nugget, scale, \ldots)[fit.param], aniso[fit.aniso])},\cr
        where \code{\ldots} are additional
        parameters of isotropic variogram models (use \cr
        \code{param.names(variogram.model)} to display the
        names and the order of the additional parameters for
        \code{variogram.model}).
      
    } 
      
  }
  
  \subsection{Computing robust initial estimates of parameters for robust
  REML}{
  
    To solve the robustified estimating equations for
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} the following initial
    estimates are used:
    
    \itemize{
    
      \item \eqn{ \widehat{\mbox{\boldmath$B$\unboldmath}}=
      \mbox{\boldmath$0$\unboldmath},}{hatB=0,} if this turns out to be
      unfeasible, initial values can be passed to \code{georob} by the
      argument \code{bhat} of \code{\link{georob.control}}.
      
      \item \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta} is
      either estimated robustly by the function \code{\link[robustbase]{lmrob}} or
      \code{\link[quantreg]{rq}} (see argument \code{initial.method} of
      \code{\link{georob.control}}).
      
    }
  
    Finding the roots of the robustified estimating equations of the
    variogram and anisotropy parameters is more sensitive to a good choice
    of initial values than maximizing the Gaussian restricted loglikelihood
    with respect to the same parameters.  Two options are implemented to
    get good initial values that are often sufficiently close to the roots
    so that \code{\link[nleqslv]{nleqslv}} converges:
    
    Setting \code{initial.param = "minimize"} invokes
    \code{\link[stats]{optim}} to minimize the \emph{sum of squared
    estimating equations}.  The required accuracy of the initial estimates
    are best controlled by the argument \code{abstol} of
    \code{\link[stats]{optim}}, e.g. by using the argument
    \code{control = georob.control(optim = optim.control(optim.control = list(abstol = 1.e-6)))}.
    
    Setting \code{initial.param = "exclude"} has the following effects:
    
    \enumerate{
    
      \item Initial values of the regression parameters are computed by
        \code{\link[robustbase]{lmrob}} irrespective of the choice for 
        \code{initial.method} (see \code{\link{georob.control}}).
        
      \item Observations with \dQuote{robustness weights} of the
      \code{lmrob} fit, satisfying\cr
      \eqn{\psi_c(\widehat{\varepsilon}_i/\widehat{\tau})/(\widehat{\varepsilon}_i/\widehat{\tau})
      \leq \mbox{\code{min.rweight}}}{\psi_c(r_i/s)/(r_i/s)
      <=min.rweight}, are discarded (see
      \code{\link{georob.control}}).
      
      \item The model is fit to the pruned data set by Gaussian REML using 
        \code{\link[stats]{optim}}.
      
      \item The resulting estimates of the variogram parameters
        (\code{param}, \code{aniso}) are used as initial estimates for the
        subsequent robust fit of the model by \code{\link[nleqslv]{nleqslv}}.
        
    }
    
  }
    
}
\value{

An object of class \code{georob} representing a robust (or Gaussian) REML
fit of a spatial linear model.  See
\code{\link{georobObject}} for the components of the fit.

}
\references{
  \enc{K?nsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A. (in
  preparation) Robust Geostatistics.

  \enc{K?nsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A. (2011) Robust
  estimation of the external drift and the variogram of spatial data.
  Proceedings of the ISI 58th World Statistics Congress of the International
  Statistical Institute.
  \url{http://isi2011.congressplanner.eu/pdfs/650268.pdf}
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}\cr
  \url{http://www.step.ethz.ch/people/scientific-staff/andreas-papritz}\cr
  with contributions  by Cornelia Schwierz.
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georobObject}} for a description of the class \code{georob};
  \code{\link{plot.georob}} for display of REML variogram estimates;
  \code{\link{georob.control}} for controlling the behaviour of \code{georob};
  \code{\link{cv.georob}} for assessing the goodness of a fit by \code{georob}; 
  \code{\link{predict.georob}} for computing robust kriging predictions; and finally
  \code{\link{georobModelBuilding}} for stepwise building models of class \code{georob};
  \code{\link{georobMethods}} for further methods for the class \code{georob}.
}

\examples{
\dontrun{
################
## meuse data ##
################
data( meuse )

## Gaussian REML fit
r.logzn.reml <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "exponential",
    param = c( variance = 0.15, nugget = 0.05, scale = 200 ),
    tuning.psi = 1000)
summary(r.logzn.reml, correlation = TRUE)

## robust REML fit 
r.logzn.rob <- update(r.logzn.reml, tuning.psi = 1)
    
summary(r.logzn.rob, correlation = TRUE)

plot(r.logzn.reml, lag.class.def = seq( 0, 2000, by = 100 ))
lines(r.logzn.rob, col = "red")


###################
## wolfcamp data ##
###################
data(wolfcamp, package = "geoR")
d.wolfcamp <- data.frame(x = wolfcamp[[1]][,1], y = wolfcamp[[1]][,2],
    pressure = wolfcamp[[2]])

## fitting isotropic IRF(0) model
  
r.irf0.iso <- georob(pressure ~ 1, data = d.wolfcamp, locations = ~ x + y, 
    variogram.model = "fractalB",
    param = c( variance = 10, nugget = 1500, scale = 1, alpha = 1.5 ),
    fit.param = c( variance = TRUE, nugget = TRUE, scale = FALSE, alpha = TRUE),
    tuning.psi = 1000)
  
summary(r.irf0.iso)

## fitting isotropic IRF(0) model
  
r.irf0.aniso <- georob(pressure ~ 1, data = d.wolfcamp, locations = ~ x + y, 
    variogram.model = "fractalB",
    param = c( variance = 5.9, nugget = 1450, scale = 1, alpha = 1 ),
    fit.param = c( variance = TRUE, nugget = TRUE, scale = FALSE, alpha = TRUE),
    aniso = c( f1 = 0.51, f2 = 1, omega = 148, phi = 90, zeta = 0 ),
    fit.aniso = c( f1 = TRUE, f2 = FALSE, omega = TRUE, phi = FALSE, zeta = FALSE ), 
    tuning.psi = 1000)
summary(r.irf0.aniso)

plot(r.irf0.iso, lag.class.def = seq(0, 200, by = 7.5))
plot(r.irf0.aniso, lag.class.def = seq(0, 200, by = 7.5), 
    xy.angle.def = c(0, 22.5, 67.5, 112.5, 157.5, 180.), 
    add = TRUE, col = 2:5)
    
pchisq( 2*(r.irf0.aniso[["loglik"]] - r.irf0.iso[["loglik"]]), 2, lower = FALSE )}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
