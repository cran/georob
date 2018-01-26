% 2016-02-25 A. Papritz
% R CMD Rdconv -t html -o bla.html plot.georob.Rd ; open bla.html; R CMD Rd2pdf --force plot.georob.Rd; 

\encoding{macintosh}
\name{plot.georob}
\alias{plot.georob}
\alias{lines.georob}

\title{Plot Methods for Class \code{georob}}

\description{The \code{plot} and \code{lines} methods for class
\code{georob} plot the variogram model, estimated by (robust) restricted
maximum likelihood.
\code{plot.georob} computes and plots in addition the
sample variogram of the (robust) regression residuals and can be used to
generate residual diagnostics plots (Tukey-Anscombe plot, normal QQ plots
of residuals and random effects).

}


\usage{
\method{plot}{georob}(x, what = c( "variogram", "covariance", "correlation", 
    "ta", "sl", "qq.res", "qq.ranef" ), add = FALSE, lag.dist.def, 
    xy.angle.def = c(0, 180), xz.angle.def = c(0, 180), max.lag = Inf, 
    estimator = c("mad", "qn", "ch", "matheron"), mean.angle = TRUE, 
    level = what != "ta", smooth = what == "ta" || what == "sl", 
    id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
    label.pos = c(4,2), col, pch, xlab, ylab, main, lty = "solid", ...)

\method{lines}{georob}(x, what = c("variogram", "covariance", "correlation"),
    from = 1.e-6, to, n = 501, xy.angle = 90, xz.angle = 90, 
    col = 1:length(xy.angle), pch = 1:length(xz.angle), lty = "solid", ...)
}

\arguments{
  \item{x}{an object of class \code{georob}, see \code{\link{georobObject}}.}
  
  \item{what}{character keyword for the quantity that should be displayed.
  Possible values are:
  
    \itemize{
    
      \item \code{"variogram"}: the variogram
    
      \item \code{"covariance"}: the covariance function
    
      \item \code{"correlation"}: the correlation function
      
      \item \code{"scale-location"}: square root of absolute regression
      residuals plotted against fitted values (Scale-Location plot)

      \item \code{"ta"}: regression residuals plotted against fitted values
      (Tukey-Anscombe plot)

      \item \code{"qq.res"}: normal Q-Q plot of standardized errors
      \eqn{\hat{\varepsilon}}{hat\epsilon}
      
      \item \code{"qq.ranef"}: normal Q-Q plot of standardized random
      effects \eqn{\hat{B}}{hatB}
    
    }
  
  }
  
  \item{add}{logical controlling whether a new plot should be
  generated (\code{FALSE}, default) or whether the information should be
  added to the current plot (\code{TRUE}).}

  \item{lag.dist.def}{an optional numeric scalar defining a constant bin
  width for grouping the lag distances or an optional numeric vector with
  the upper bounds of a set of contiguous bins for computing the sample
  variogram of the regression residuals, see
  \code{\link{sample.variogram}}.  If missing then the sample variogram is
  not computed.}
  
  \item{xy.angle.def}{an numeric vector defining angular classes in the
  horizontal plane for computing directional variograms.
  \code{xy.angle.def} must contain an ascending sequence of azimuth angles
  in degrees from north (positive clockwise to south), see
  \code{\link{sample.variogram}}.  Omnidirectional variograms are computed
  with the default \code{c(0,180)}.}
  
  \item{xz.angle.def}{an numeric vector defining angular classes in the
  \eqn{x}-\eqn{z}-plane for computing directional variograms.
  \code{xz.angle.def} must contain an ascending sequence of angles in
  degrees from zenith (positive clockwise to nadir), see
  \code{\link{sample.variogram}}.  Omnidirectional variograms are computed
  with the default \code{c(0,180)}.}
  
  \item{max.lag}{positive numeric defining the largest lag distance for
  which semi-variances should be computed (default no restriction).}
  
  \item{estimator}{character keyword defining the estimator for computing
  the sample variogram.  Possible values are:
  
    \itemize{
    
      \item \code{"qn"}: Genton's robust
      \code{\link[robustbase]{Qn}}-estimator (default, Genton, 1998),
      
      \item \code{"mad"}: Dowd's robust MAD-estimator (Dowd, 1984),
      
      \item \code{"matheron"}: non-robust method-of-moments estimator,
      
      \item \code{"ch"}: robust Cressie-Hawkins estimator (Cressie and
      Hawkins, 1980).
    
    }
  }
  
  \item{mean.angle}{logical controlling whether the mean lag vector (per
  combination of lag distance and angular class) is computed from the mean
  angles of all the lag vectors falling into a given class (\code{TRUE},
  default) or from the mid-angles of the respective angular classes
  (\code{FALSE}).}
  
  \item{level}{an integer giving the level for extracting the residuals
  from \code{object} for \code{what = "ta"} or \code{what = "qq.res"}.
  \code{level = 0} (default for \code{what == "ta"}) extracts the
  regression residuals \eqn{\widehat{B}(\mbox{\boldmath$s$\unboldmath}) +
  \widehat{\varepsilon}(\mbox{\boldmath$s$\unboldmath})}{hatB(s) +
  hat\epsilon(s)} and \code{level = 1} (default for \code{what ==
  "qq.res"}) only the estimated errors
  \eqn{\widehat{\varepsilon}(\mbox{\boldmath$s$\unboldmath})}{hat\epsilon(s)}.}
  
  \item{smooth}{logical controlling whether a
  \code{\link[stats]{loess.smooth}} is added to the Tukey-Anscombe plot
  (default \code{TRUE}.}
  
  \item{id.n}{number of points to be labelled in each plot, starting
  with the most extreme (see \code{\link[robustbase]{plot.lmrob}}).}

  \item{labels.id}{vector of labels, from which the labels for extreme
  points will be chosen (see \code{\link[robustbase]{plot.lmrob}}).
  \code{NULL} uses observation numbers.}

  \item{cex.id}{magnification of point labels (see
  \code{\link[robustbase]{plot.lmrob}}).}

  \item{label.pos}{positioning of labels, for the left half and right half
  of the graph respectively (see \code{\link[robustbase]{plot.lmrob}}).}
  
  \item{from}{numeric, minimal lag distance for plotting variogram models.}

  \item{to}{numeric, maximum lag distance for plotting variogram
  models (default: largest lag distance of current plot).}

  \item{n}{positive integer specifying the number of equally spaced lag
  distances for which semi-variances are evaluated in plotting variogram
  models (default \code{501}).}

  \item{xy.angle}{numeric (vector) with azimuth angles (in degrees,
  clockwise positive from north) in \eqn{x}-\eqn{y}-plane for which
  semi-variances should be plotted.}

  \item{xz.angle}{numeric (vector) with angles in \eqn{x}-\eqn{z}-plane (in
  degrees, clockwise positive from zenith to south) for which semi-variances
  should be plotted.}

  \item{col}{optional color of points and curves to distinguish items
  relating to different azimuth angles in \eqn{x}-\eqn{y}-plane.}

  \item{pch}{optional symbol for points and curves to distinguish items
  relating to different azimuth angles in \eqn{x}-\eqn{z}-plane.}

  \item{lty}{line type for plotting variogram models.}

  \item{xlab, ylab, main}{test annotation, see
  \code{\link[graphics]{plot}}.}
  
  \item{...}{additional arguments passed to
  \code{\link{plot.sample.variogram}}, \code{\link[stats]{loess.smooth}}
  and graphical methods.}
  
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}.
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  
  \code{\link{georobObject}} for a description of the class \code{georob};
  
  \code{\link{profilelogLik}} for computing profiles of Gaussian likelihoods;
  
%  \code{\link{plot.georob}} for display of RE(ML) variogram estimates;
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
\dontrun{
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

## robust REML fit 
r.logzn.rob <- update(r.logzn.reml, tuning.psi = 1)
    
summary(r.logzn.rob, correlation = TRUE)

plot(r.logzn.reml, lag.dist.def = seq(0, 2000, by = 100))
lines(r.logzn.rob, col = "red")}
}


\keyword{models}
\keyword{spatial}
\keyword{robust}
