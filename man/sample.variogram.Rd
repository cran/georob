% 2021-06-04 A. Papritz
% R CMD Rdconv -t html -o bla.html sample.variogram.Rd ; open bla.html; R CMD Rd2pdf --force sample.variogram.Rd;

\encoding{utf8}
\name{sample.variogram}
\alias{sample.variogram}
\alias{sample.variogram.default}
\alias{sample.variogram.formula}
\alias{sample.variogram.georob}
\alias{summary.sample.variogram}
\alias{print.summary.sample.variogram}
\alias{plot.sample.variogram}

\title{Computing (Robust) Sample Variograms of Spatial Data}

\description{The function \code{sample.variogram} computes the
sample (empirical) variogram of a spatial variable by the method-of-moment
and three robust estimators.  Both omnidirectional and direction-dependent
variograms can be computed, the latter for observation locations in a
three-dimensional domain.  There are \code{summary} and \code{plot}
methods for summarizing and displaying sample variograms.}


\usage{

sample.variogram(object, ...)

\method{sample.variogram}{default}(object, locations, lag.dist.def,
    xy.angle.def = c(0, 180), xz.angle.def = c(0, 180), max.lag = Inf,
    estimator = c("qn", "mad", "matheron", "ch"), mean.angle = TRUE, ...)

\method{sample.variogram}{formula}(object, data, subset, na.action,
    locations, lag.dist.def, xy.angle.def = c(0, 180),
    xz.angle.def = c(0, 180), max.lag = Inf,
    estimator = c("qn", "mad", "matheron", "ch"), mean.angle = TRUE, ...)

\method{sample.variogram}{georob}(object, lag.dist.def,
    xy.angle.def = c(0, 180), xz.angle.def = c(0, 180), max.lag = Inf,
    estimator = c("qn", "mad", "matheron", "ch"), mean.angle = TRUE, ...)

\method{summary}{sample.variogram}(object, ...)

\method{plot}{sample.variogram}(x, type = "p", add = FALSE,
    xlim = c(0, max(x[["lag.dist"]])),
    ylim = c(0, 1.1 * max(x[["gamma"]])), col, pch, lty, cex = 0.8,
    xlab = "lag distance", ylab = "semivariance",
    annotate.npairs = FALSE, npairs.pos = 3, npairs.cex = 0.7,
    legend = nlevels(x[["xy.angle"]]) > 1 || nlevels(x[["xz.angle"]]) > 1,
    legend.pos = "topleft", ...)

}

\arguments{
  \item{object}{a numeric vector with the values of the response for which
  the sample variogram should be computed
  (\code{sample.variogram.default}), a formula, specifying in its left part
  the response variable (right part of formula is ignored,
  \code{sample.variogram.formula}), an object of class \code{georob}\cr
  (\code{sample.variogram.georob}) or an object of class\cr
  \code{sample.variogram} (\code{summary.sample.variogram}).}

  \item{locations}{a numeric matrix with the coordinates of the locations
  where the response was observed (\code{sample.variogram.default}) or a
  one-sided formula specifying the coordinates
  (\code{sample.variogram.formula}).  The matrix may have an arbitrary
  number of columns for an omnidirectional variogram, but at most 3 columns
  if a directional variogram is computed.}

  \item{data}{an optional data frame, list or environment (or another
  object coercible by\cr \code{\link{as.data.frame}} to a data frame)
  containing the response variable and the coordinates where the data
  was recorded.  If not found in \code{data}, the variables are taken from
  \code{environment(formula)}, typically the environment from which
  \code{sample.variogram} is called.}

  \item{subset}{an optional vector specifying a subset of observations
  to be used for estimating the variogram.}

  \item{na.action}{a function which indicates what should happen when the
  data contain \code{NA}s.  The default is set by the \code{na.action}
  argument of \code{\link{options}}, and is \code{\link{na.fail}} if that is
  unset.  The \dQuote{factory-fresh} default is \code{\link{na.omit}}.
  Another possible value is \code{NULL}, no action.  Value
  \code{\link{na.exclude}} can be useful.}

  \item{lag.dist.def}{a numeric scalar defining a constant bin
  width for grouping the lag distances or a numeric vector with the bounds
  of a set of contiguous bins (upper bounds of bins except for the first
  element of \code{lag.dist.def} which is the lower bound of the first
  bin). This argument is mandatory.}

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
  which semi variances should be computed (default no restriction).}

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

  \item{x}{an object of class \code{sample.variogram}.}

  \item{type, xlim, ylim, xlab, ylab}{see respective arguments of
  \code{\link[graphics]{plot.default}}.}

  \item{add}{logical controlling whether a new plot should be
  generated (\code{FALSE}, default) or whether the information should be
  added to the current plot (\code{TRUE}).}

  \item{col}{the color of plotting symbols for distinguishing semi variances
  for angular classes in the \eqn{x}-\eqn{y}-plane.}

  \item{pch}{the type of plotting symbols for distinguishing semi variances
  for angular classes in the \eqn{x}-\eqn{z}-plane.}

  \item{lty}{the line type.}

  \item{cex}{character expansion factor for plotting symbols.}

  \item{annotate.npairs}{logical controlling whether the plotting symbols
  should be annotated by the number of data pairs per lag class.}

  \item{npairs.pos}{integer defining the position where text annotation
  about number of pairs should be plotted, see
  \code{\link[graphics]{text}}.}

  \item{npairs.cex}{numeric defining the character expansion for text
  annotation about number of pairs.}

  \item{legend}{logical controlling whether a
  \code{\link[graphics]{legend}} should be plotted.}

  \item{legend.pos}{a character keyword defining where to place the
  legend, see \code{\link[graphics]{legend}} for possible values.}

  \item{...}{additional arguments passed to
  \code{\link[graphics]{plot.formula}}.}

}

\details{
  The angular classes in the \eqn{x}-\eqn{y}- and \eqn{x}-\eqn{z}-plane are
  defined by vectors of ascending angles on the half circle.  The \eqn{i}th
  angular class is defined by the vector elements, say \var{l} and \var{u},
  with indices \eqn{i} and \eqn{i+1}.  A lag vector belongs to the
  \eqn{i}th angular class if its azimuth (or angle from zenith), say
  \eqn{\varphi}{\phi}, satisfies \eqn{ l < \varphi \leq u}{l < \phi <= u}.
  If the first and the last element of \code{xy.angle.def} or
  \code{xz.angle.def} are equal to \code{0} and \code{180} degrees,
  respectively, then the first and the last angular class are
  \dQuote{joined}, i.e., if there are \eqn{K} angles, there will be only
  \eqn{K-2} angular classes and the first class is defined by the interval
  ( \var{xy.angle.def}[\var{K}-1]-180, \var{xy.angle.def}[2] ] and the last
  class by ( \var{xy.angle.def}[\var{K}-2], \var{xy.angle.def}[\var{K}-1]].
}

\value{An object of class \code{sample.variogram}, which is a data frame
with the following components:

\tabular{ll}{
  \code{lag.dist} \tab the mean lag distance of the lag class, \cr
  \code{xy.angle} \tab the angular class in the \eqn{x}-\eqn{y}-plane, \cr
  \code{xz.angle} \tab the angular class in the \eqn{x}-\eqn{z}-plane, \cr
  \code{gamma} \tab the estimated semi-variance of the lag class, \cr
  \code{npairs} \tab the number of data pairs in the lag class, \cr
  \code{lag.x} \tab the \eqn{x}-component of the mean lag vector of the lag class, \cr
  \code{lag.x} \tab the \eqn{y}-component of the mean lag vector of the lag class, \cr
  \code{lag.z} \tab the \eqn{z}-component of the mean lag vector of the lag class.\cr
}

}


\references{
  Cressie, N. and Hawkins, D. M. (1980) Robust Estimation of the Variogram: I.
  \emph{Mathematical Geology}, \bold{12}, 115--125.

  Dowd, P. A. (1984) The variogram and Kriging: Robust and resistant
  estimators.  In \emph{Geostatistics for Natural Resources
  Characterization},  Verly, G., David, M., Journel, A. and Marechal, A. (Eds.)
  Dordrecht: D. Reidel Publishing Company, Part I, 1, 91--106.

  Genton, M. (1998) Highly Robust Variogram Estimation.  \emph{Mathematical
  Geology}, \bold{30}, 213--220.

}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}.
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms;

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

\dontrun{
plot(r.sv.iso, type = "l")}

## fitting an anisotropic IRF(0) model
r.sv.aniso <- sample.variogram(pressure~1, data = wolfcamp,
    locations = ~x + y, lag.dist.def = seq(0, 200, by = 15),
    xy.angle.def = c(0., 22.5, 67.5, 112.5, 157.5, 180.))
\dontrun{
plot(r.sv.aniso, type = "l", add = TRUE, col = 2:5)}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
