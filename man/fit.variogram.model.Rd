% 2013-06-12 A. Papritz
% R CMD Rdconv -t html -o bla.html fit.variogram.model.Rd ; open bla.html; R CMD Rd2pdf --force fit.variogram.model.Rd; 

\encoding{macintosh}
\name{fit.variogram.model}
\alias{fit.variogram.model}
\alias{print.fitted.variogram}
\alias{summary.fitted.variogram}
\alias{print.summary.fitted.variogram}
\alias{lines.fitted.variogram}

\title{Fitting Model Functions to Sample Variograms}

\description{
The function \code{fit.variogram.model} fits a variogram model to a sample
variogram by weighted non-linear least squares.  There are \code{print},
\code{summary} and
\code{lines} methods for summarizing and displaying fitted variogram
models.}

\usage{
fit.variogram.model(sv, 
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
    max.lag = max(sv[["lag.dist"]]), min.npairs = 30, 
    weighting.method = c("cressie", "equal", "npairs"), hessian = TRUE, 
    verbose = 0, ...)
    
\method{print}{fitted.variogram}(x, digits = max(3, getOption("digits") - 3), ...)

\method{summary}{fitted.variogram}(object, correlation = FALSE, signif = 0.95, ...)
    
\method{lines}{fitted.variogram}(x, from = 1.e-6, to, n = 501, xy.angle = 90, xz.angle = 90, 
    col = 1:length(xy.angle), pch = 1:length(xz.angle), lty = "solid", ...)
    
}

\arguments{
  \item{sv}{an object of class \code{sample.variogram}, see
  \code{\link{sample.variogram}}.}

  \item{variogram.model}{a character keyword defining the variogram model to
  be fitted.  Currently, most basic variogram models provided by the
  package \pkg{RandomFields} can be fitted (see \emph{Details} of
  \code{\link{georob}} and \code{\link[RandomFields]{Variogram}}).}
  
  \item{param}{a named numeric vector with initial values of the variogram
  parameters.  The following parameter names are allowed (see
  \emph{Details} of \code{\link{georob}} and \code{\link{georobIntro}} for
  information about the parametrization of variogram models):
    
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
     
       \item{names of additional variogram parameters such as the
       smoothness parameter \eqn{\nu} of the Whittle-Mat?rn model as passed
       by the argument \code{param} to
       \code{\link[RandomFields]{Variogram}}, see \code{\link{param.names}}.} } }
  
  \item{fit.param}{a named logical vector with the same names as used for 
    \code{param}, defining which parameters are adjusted (\code{TRUE}) and which 
    are kept fixed at their initial values (\code{FALSE}) when fitting the model.}
    
  \item{aniso}{a named numeric vector with initial values for fitting
    geometrically anisotropic variogram models. The following parameter names are allowed 
    (see \emph{Details} of \code{\link{georob}} and \code{\link{georobIntro}} for information about the 
    parametrization of variogram models): 
    
    \itemize{
    
      \item{\code{f1}: ratio \eqn{f_1} of lengths of second and first
        second semi-principal axes of an ellipsoidal surface with constant
        semivariance in \eqn{\mathrm{I}\!\mathrm{R}^3}{R^3} (default \code{f1 = 1}).}
      
      \item{\code{f2}: ratio \eqn{f_2} of lengths of third and first
        semi-principal axes of the semivariance ellipsoid (default \code{f2 = 1}).}
      
      \item{\code{omega}: azimuth in degrees of first semi-principal axis 
        of the semivariance ellipsoid (default \code{omega = 90}).}
      
      \item{\code{phi}: 90 degrees minus altitude of first semi-principal axis 
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
      
      \item \code{"cressie"}: \dQuote{Cressie's weights} (default, see Cressie, 1993,
      sec.  2.6.2).
      
    }
  }

  \item{hessian}{logical controlling whether the hessian is computed by
  \code{\link[stats]{optim}}.}

  \item{verbose}{positive integer controlling logging of diagnostic
  messages to the console during model fitting.}

  \item{object, x}{an object of class \code{fitted.variogram}.}
  
  \item{digits}{positive integer indicating the number of decimal digits to print.}

  \item{correlation}{logical controlling whether the correlation matrix of
  the fitted variogram parameters is computed (default \code{FALSE}).}
  
  \item{signif}{confidence level for computing confidence intervals for
  variogram parameters (default \code{0.95}).}
  
  \item{from}{numeric, minimal lag distance used in plotting variogram models.}

  \item{to}{numeric, maximum lag distance used in plotting variogram
  models (default: largest lag distance of current plot).}

  \item{n}{positive integer specifying the number of equally spaced lag
  distances for which semivariances are evaluated in plotting variogram
  models (default \code{501}).}

  \item{xy.angle}{numeric (vector) with azimuth angles (in degrees,
  clockwise positive from north) in \eqn{x}-\eqn{y}-plane for which
  semivariances should be plotted.}

  \item{xz.angle}{numeric (vector) with angles in \eqn{x}-\eqn{z}-plane (in
  degrees, clockwise positive from zenith to south) for which semivariances
  should be plotted.}

  \item{col}{color of curves to distinguish curves relating to
  different azimuth angles in \eqn{x}-\eqn{y}-plane.}

  \item{pch}{type of plotting symbols added to lines to distinguish
  curves relating to different angles in \eqn{x}-\eqn{z}-plane.}

  \item{lty}{line type for plotting variogram models.}

  \item{\dots}{additional arguments passed to \code{\link[stats]{optim}} or
  to methods.}

}
\details{The parametrization of geometrically anisotropic variograms is
described in detail in \code{\link{georobIntro}}, and the section
\emph{Details} of \code{\link{georob}} describes how the parameter
estimates are constrained to permissible ranges.  The same
mechanisms are used in \code{fit.variogram.model}.

}


\value{

  The function \code{fit.variogram.model} generates an object of class
  \code{fitted.variogram} which is a list with the following components:

  \item{sse}{the value of the object function (weighted residual sum of
  squares) evaluated at the solution.}

  \item{variogram.model}{the name of the fitted parametric variogram
  model.}
   
  \item{param}{a named vector with the (estimated) variogram parameters of
  the fitted model.}

  \item{aniso}{a list with the following components: 
   
    \itemize{

      \item{\code{isotropic}: logical indicating whether an isotropic
       variogram was fitted.}
       
      \item{\code{aniso}: a named numeric vector with the (estimated)
      anisotropy parameters of the fitted model.}
       
      \item{\code{sincos}: a list with \code{sin} and \code{cos} of the
      angles \eqn{\omega}, \eqn{\phi} and \eqn{\zeta} that define the
      orientation of the anisotropy ellipsoid.}

      \item{\code{rotmat}: the matrix
      \eqn{(\mbox{\boldmath$C$\unboldmath}_1,
      \mbox{\boldmath$C$\unboldmath}_2,
      \mbox{\boldmath$C$\unboldmath}_3)}{(C_1, C_2, C_3)} (see
      \code{\link{georobIntro}}).}
      
      \item{\code{sclmat}: a vector with the elements 1, \eqn{1/f_1},
      \eqn{1/f_2} (see \code{\link{georobIntro}}). }

    }
  }
  
  \item{param.tf}{a character vector listing the transformations of the
  variogram parameters used for model fitting.}
    
  \item{fwd.tf}{a list of functions for variogram parameter transformations.}
  
  \item{bwd.tf}{a list of functions for \emph{inverse} variogram parameter
  transformations.}
    
  \item{converged}{logical indicating whether numerical maximization by
  \code{optim} converged.}

    \item{convergence.code}{a diagnostic integer issued by
    \code{\link{optim}} (component \code{convergence}) about convergence.}

  \item{call}{the matched call.}  
  
  
  \item{residuals}{a numeric vector with the residuals, that is the sample
  semivariance minus the fitted values.}
  
  \item{fitted}{a numeric vector with the modelled semivariances.}

  \item{weights}{a numeric vector with the weights used for fitting.}
  
  \item{hessian}{a symmetric matrix giving an estimate of the Hessian at
  the solution (missing if \code{hessian} is false).}
    
}

\references{

  Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}.  New York:
  John Wiley & Sons.

}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}.
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{sample.variogram}} for computing sample variograms.
}

\examples{
data(wolfcamp, package = "geoR")

## fitting an isotropic IRF(0) model
r.sv.iso <- sample.variogram(wolfcamp[["data"]], locations = wolfcamp[[1]], 
    lag.class.def = seq(0, 200, by = 15))

r.irf0.iso <- fit.variogram.model(r.sv.iso, variogram.model = "fractalB",
    param = c(variance = 100, nugget = 1000, scale = 1., alpha = 1.),
    fit.param = c( variance = TRUE, nugget = TRUE, scale = FALSE, alpha = TRUE),
    method = "Nelder-Mead", hessian = FALSE, control = list(maxit = 5000))  
summary(r.irf0.iso, correlation = TRUE)

\dontrun{
plot( r.sv.iso, type = "l")
lines( r.irf0.iso, line.col = "red")}

## fitting an anisotropic IRF(0) model
r.sv.aniso <- sample.variogram(wolfcamp[["data"]],
    locations = wolfcamp[[1]], lag.class.def = seq(0, 200, by = 15),
    xy.angle.def = c(0., 22.5, 67.5, 112.5, 157.5, 180.))
\dontrun{
plot(r.sv.aniso, type = "l")}

r.irf0.aniso <- fit.variogram.model(r.sv.aniso, variogram.model = "fractalB",
    param = c(variance = 100, nugget = 1000, scale = 1., alpha = 1.5),
    fit.param = c(variance = TRUE, nugget = TRUE, scale = FALSE, alpha = TRUE),
    aniso = c(f1 = 0.4, f2 = 1., omega = 135, phi = 90., zeta = 0.),
    fit.aniso = c(f1 = TRUE, f2 = FALSE, omega = TRUE, phi = FALSE, zeta = FALSE),
    method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 5000))
summary(r.irf0.aniso, correlation = TRUE)

\dontrun{
lines(r.irf0.aniso, xy.angle = seq( 0, 135, by = 45))}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
