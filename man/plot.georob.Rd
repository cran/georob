% 2012-12-14 A. Papritz
% R CMD Rdconv -t html -o bla.html plot.georob.Rd ; open bla.html; R CMD Rd2pdf --force plot.georob.Rd; 

\encoding{macintosh}
\name{plot.georob}
\alias{plot.georob}
\alias{lines.georob}

\title{Plot Methods for class \code{georob}}

\description{The \code{plot} and \code{lines} methods for class
\code{georob} plot the variogram model, estimated by (robust) restricted
maximum likelihood.
\code{plot.georob} computes and plots in addition the
sample variogram of the (robust) regression residuals.

}


\usage{
\method{plot}{georob}(x, type, plot.sv = TRUE, add = FALSE, lag.class.def, 
    xy.angle.def = c(0, 180), xz.angle.def = c(0, 180), 
    max.lag = Inf, estimator = c("mad", "qn", "ch", "matheron"), 
    mean.angle = TRUE, col, pch, lty, ...)

\method{lines}{georob}(x, from = 1.e-6, to, n = 501, xy.angle = 90, xz.angle = 90, 
    col = 1:length(xy.angle), pch = 1:length(xz.angle), lty = "solid", ...)
}

\arguments{
  \item{x}{an object of class \code{georob}, see \code{\link{georobObject}}.}
  
  \item{type}{the type of plot for display of the sample variogram,
  see \code{\link[graphics]{plot}}.}
  
  \item{plot.sv}{logical controlling if the sample variogram of the
  regression residuals, \eqn{\widehat{B}(\mbox{\boldmath$s$\unboldmath}) +
  \widehat{\varepsilon}(\mbox{\boldmath$s$\unboldmath})}{hatB(s) +
  hat\epsilon{s}} should be added to the plot (default \code{TRUE}).}
  
  \item{add}{logical controlling whether a new plot should be
  generated (\code{FALSE}, default) or whether the information should be
  added to the current plot (\code{TRUE}).}

  \item{lag.class.def}{a numeric scalar defining a constant bin width for
  grouping the lag distances or a numeric vector with the upper bounds of a
  set of contiguous bins.}
  
  \item{xy.angle.def}{an numeric vector defining angular classes
  in the horizontal plane for computing directional variograms.
  \code{xy.angle.def} must contain an ascending sequence of azimuth angles
  in degrees from north (positive clockwise to south), see \emph{Details}.
  Omnidirectional variograms are computed with the default 
  \code{c(0,180)}.}
  
  \item{xz.angle.def}{an numeric vector defining angular classes
  in the \eqn{x}-\eqn{z}-plane for computing directional variograms.
  \code{xz.angle.def} must contain an ascending sequence of angles in
  degrees from zenith (positive clockwise to nadir), see
  \emph{Details}.  Omnidirectional variograms are computed with the default
  \code{c(0,180)}.}
  
  \item{max.lag}{positive numeric defining the largest lag distance for
  which semivariances should be computed (default no restriction).}
  
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
  
  \item{from}{numeric, minimal lag distance for plotting variogram models.}

  \item{to}{numeric, maximum lag distance for plotting variogram
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

  \item{\ldots}{additional arguments passed to
  \code{\link{plot.sample.variogram}} and methods.} 
  
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}.
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{georobObject}} for a description of the class \code{georob};
  \code{\link{sample.variogram}} for computing sample variograms.
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
lines(r.logzn.rob, col = "red")}
}


\keyword{models}
\keyword{spatial}
\keyword{robust}
