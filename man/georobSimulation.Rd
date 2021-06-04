% 2021-06-04 A. Papritz
% R CMD Rdconv -t html -o bla.html georobSimulation.Rd ; open bla.html; R CMD Rd2pdf --force georobSimulation.Rd; open georobSimulation.pdf

\encoding{utf8}
\name{georobSimulation}
\alias{georobSimulation}
\alias{condsim}
\alias{control.condsim}

\title{Simulating Realizations of Gaussian Processes from Object of Class \code{georob}}

\description{This page documents the function \code{condsim} that
simulates (un)conditional realizations of Gaussian processes from the
parameters of a spatial linear model estimated by the function
\code{\link{georob}}.}

\usage{
  condsim(object, newdata, nsim, seed,
    type =  c("response", "signal"), locations, trend.coef = NULL,
    variogram.model = NULL, param = NULL, aniso = NULL, variogram.object = NULL,
    control = control.condsim(), verbose = 0)

  control.condsim(use.grid = FALSE, grid.refinement = 2.,
    condsim = TRUE, include.data.sites = FALSE, means = FALSE,
    trend.covariates = FALSE, covariances = FALSE,
    ncores = detectCores(), pcmp = control.pcmp())
}

\arguments{

  \item{object}{an object of class \code{georob} (mandatory argument), see
    \code{\link{georobObject}}.}

  \item{newdata}{a mandatory data frame,
  \code{\link[sp]{SpatialPointsDataFrame}},
  \code{\link[sp]{SpatialPixelsDataFrame}},
  \code{\link[sp]{SpatialGridDataFrame}},
  \code{\link[sp]{SpatialPoints}},
  \code{\link[sp]{SpatialPixels}} or
  \code{\link[sp]{SpatialGrid}} object,
  with the coordinates of points for which simulations are computed
  and in which to look for variables required for computing fitted values
  or Kriging predictions,  see \code{\link{predict.georob}}.}

  \item{nsim}{number of (condititional) realizations to compute (mandatory
  argument).}

  \item{seed}{integer seed to initialize random number generation,
  see \code{\link[base]{set.seed}} (mandatory argument).}

  \item{type}{character keyword defining what target quantity should be
  simulated.  Possible values are

    \itemize{

      \item \code{"signal"}: the \dQuote{signal}
      \eqn{Z(\mbox{\boldmath$s$\unboldmath}) =
      \mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath})^\mathrm{T}
      \mbox{\boldmath$\beta$\unboldmath} +
      B(\mbox{\boldmath$s$\unboldmath})}{Z(s) = x(s)^T \beta + B(s)} of
      the process,

      \item \code{"response"}: the observations
      \eqn{Y(\mbox{\boldmath$s$\unboldmath}) =
      Z(\mbox{\boldmath$s$\unboldmath}) +
      \varepsilon(\mbox{\boldmath$s$\unboldmath}),}{Y(s)=Z(s) + \epsilon(s),}  (default),

    }

  see \code{\link{georobIntro}} for details on the model specification.}

  \item{locations}{an optional one-sided formula specifying what variables
  of \code{newdata} are the coordinates of the points for which simulated
  values are computed (default:
  \code{object[["locations.objects"]][["locations"]]}).}

  \item{trend.coef}{an optional numeric vector with the coefficients of the
  trend model to be used for computing the (conditional) mean function of
  the random process see \emph{Details}.}

  \item{variogram.model}{an optional character keyword defining the
  variogram model to be used for the simulations, see \code{\link{georob}}
  and \emph{Details}.}

  \item{param}{an optional named numeric vector with values of the
  variogram parameters used for the simulations, see \code{\link{georob}}
  and \emph{Details}.}

  \item{aniso}{an optional named numeric vector with values of anisotropy
  parameters of a variogram used for the simulations, see
  \code{\link{georob}} and \emph{Details}.}

  \item{variogram.object}{an optional list that defines a possibly nested
  variogram model used for the simulations, see \code{\link{georob}} and
  \emph{Details}.}

  \item{control}{a list with the components \code{use.grid},
  \code{grid.refinement}, \code{condsim},\cr \code{include.data.sites},
  \code{means}, \code{trend.covariates}, \code{covariances},
  \code{ncores}, and \code{pcmp} or a function such as
  \code{control.condsim} that generates such a list, see
  arguments of \code{control.condsim} for details.}

  \item{verbose}{positive integer controlling logging of diagnostic
  messages to the console.\cr \code{verbose = 0} (default) suppresses
  such messages.}

  \item{use.grid}{logical scalar (default \code{FALSE}) to control whether
  (conditional) realizations are computed for a rectangular grid instead of
  the coordinates of points contained in \code{newdata}, see \emph{Details}.}

  \item{grid.refinement}{numeric scalar that defines a factor by which the
  minimum differences of the coordinates between any pair of points in
  \code{newdata} are divided to setup the simulation grid, should be > 1
  (default 2), see \emph{Details}.}

  \item{condsim}{logical scalar (default \code{TRUE}) to control whether
  conditional (\code{TRUE}) or unconditional simulations (\code{FALSE}) are
  computed.
%   In both cases, the parameters of the mean and auto-correlation
%   functions of the underlying Gaussian processes are assumed to be
%   indentical to the parameters of the spatial linear model estimated by
%   \code{\link{georob}} (or the are taken from the optional arguments
%   \code{trend.coef}, \code{variogram.model}, \code{param}, \code{aniso},
%   \code{variogram.object}).  Hence, the uncertainty in the estimated mean
%   and auto-correlation functions are not taken into account,
%   i.~e.\``plug-in'' simulationd are computed.
  }

  \item{include.data.sites}{logical scalar, to control whether
  (conditionally) simulated values are returned also for the points
  of the original data set used to estimate the model parameters.}

  \item{means}{logical scalar, to control whether
  the (un)conditional means are included in the output.}

  \item{trend.covariates}{logical scalar, to control whether
  the covariates required for the trend model are included in the output.}

  \item{covariances}{logical scalar, to control whether the covariances
  between the points of the original data set used to estimate the model
  parameters (\code{attr gcvmat.d.d}, \code{\link{compress}}ed matrix) and
  the covariances between the simulation and the original data points
  (\code{attr gcvmat.s.d}, matrix) are returned as attributes of the
  output.  Note that these covariances are only returned if \cr
  \code{use.grid == TRUE & condsim == TRUE}.}

  \item{ncores}{positive integer controlling how many cores are used for
  parallelized computations, defaults to all cores.}

  \item{pcmp}{a list of arguments, passed  e.g. to \code{\link{pmm}} or a
  function such as \code{control.pcmp} that generates such a list
  (see \code{\link{control.pcmp}} for allowed arguments).}

}


\details{\code{condsim} (conditionally) simulates from a
Gaussian process that has a linear mean function with parameters
\eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} and an
auto-correlation structure characterized by a parametric variogram
model and variogram parameters \eqn{\tau^2}{\tau^2} and
\eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta} (see
\code{\link{georobIntro}} for the employed parametrization of the
spatial linear model).  The parameters of the mean and
auto-correlation function are either taken from the the spatial linear
model estimated by \code{\link{georob}} and passed by the argument
\code{object} to \code{condsim} or from the optional arguments
\code{trend.coef} (\eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta})
and \code{variogram.model}, \code{param}, \code{aniso} or
\code{variogram.object} (\eqn{\tau^2}{\tau^2},
\eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta}).  Note that in the
former case the uncertainty in the estimated mean and auto-correlation
parameters is not taken into account.

Simulated values are computed for the points in \code{newdata} by the
function \code{\link[RandomFields]{RFsimulate}} of the package
\pkg{RandomFields}.  Both unconditional and conditional simulations
can be computed.  In the latter cases, the simulated values are always
conditioned to the response data used to fit the spatial linear model
by \code{\link{georob}} and contained in \code{object}.

\subsection{Unconditional simulation}{

  Unconditional realizations are either computed for the exact locations of
  the points in \code{newdata} (\code{use.grid == FALSE}), irrespective of
  the fact whether these are arranged on a regular grid, or for the
  (approximate) locations of the points in \code{newdata} matched to a
  rectangular simulation grid (\code{use.grid == TRUE}).  The latter
  approach may be substantially faster for large problems because the
  simulation algorithm implemented in
  \code{\link[RandomFields]{RFsimulate}} for grids is faster than for
  arbitrary geometries of the simulation points.


	For \code{use.grid == TRUE}, a rectangular grid is constructed from the
	coordinates of the points in \code{newdata} and \code{object}.  The
	spacing of the grid is equal to the minimum distance between any pair of
	points in \code{newdata}, divided by \code{grid.refinement}.  The data
	related to the points in \code{newdata} (covariates for the trend model)
	and of the data in \code{object} (response values, covariates) are then
	assigned to the nodes of the grid that are closest to the respective
	points.  If the same grid node is assigned to several points in
	\code{newdata} (or in \code{object}) then the data of the respective
	points are averaged.  If the same node is assigned to a point in
	\code{object} and \code{newdata} then the point in \code{object} is kept
	and the concerned point in \code{newdata} is omitted.

}

\subsection{Conditional simulation}{

  Simulations are conditioned to data either by exploiting the respective
  built-in functionality of \code{\link[RandomFields]{RFsimulate}}
  (\code{use.grid == FALSE}) or by the Kriging method (\code{use.grid ==
  TRUE}, see \enc{Chilès}{Chiles} and Delfiner, 1999, sec.  7.3).  The
  latter approach may again be faster for large problems because it
  exploits the larger speed of unconditional simulations for rectangular
  grids.

}


\subsection{Parallelized computations}{

  \code{condsim} uses the packages \pkg{parallel}, \pkg{snow} and
  \pkg{snowfall} for parallelized computation of simulations.  If there
  are \eqn{m} realizations to simulate, the task is split into
  \code{ceiling(m/ncores)} sub-tasks that are then distributed to
  \code{ncores} CPUs.  Evidently, \code{ncores = 1} suppresses parallel
  execution.  By default, the function uses all available CPUs as returned
  by \code{\link[parallel]{detectCores}}.

}



}


\value{
  The output generated by \code{condsim} is an object of a ``similar''
  class as \code{newdata} (data frame,\cr
  \code{\link[sp]{SpatialPointsDataFrame}},
  \code{\link[sp]{SpatialPixelsDataFrame}},
  \code{\link[sp]{SpatialGridDataFrame}}, \cr
  \code{\link[sp]{SpatialPolygonsDataFrame}}).

  The data frame or the
  \code{data} slot of the \code{Spatial...DataFrame} objects
  have the following components:

  \itemize{

    \item the coordinates of the prediction points (only present if
    \code{newdata} is a data frame).

    \item \code{expct}: optionally the (un)conditional means.

    \item optionally the covariates required for the trend model.

    \item \code{sim.1}, \code{sim.2}, ...: the (un)conditionally
    simulated realizations.

  }

}



\references{

  \enc{Chilès}{Chiles}, J.-P. and Delfiner, P. (1999) \emph{Geostatistics:
  Modeling Spatial Uncertainty}, John Wiley & Sons, New York.

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

%  \code{\link{georobSimulation}} for simulating realizations of a Gaussian process
%  from model fitted by \code{georob}; and finally
  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
  for robust estimation and modelling of sample variograms.
}

\examples{
\dontrun{

data(meuse)
data(meuse.grid)

## convert to SpatialGridDataFrame
meuse.grid.sgdf <- meuse.grid
coordinates(meuse.grid.sgdf) <- ~ x + y
gridded(meuse.grid.sgdf) <- TRUE
fullgrid(meuse.grid.sgdf) <- TRUE

## Gaussian REML fit
r.logzn.reml <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "RMexp",
    param = c(variance = 0.15, nugget = 0.05, scale = 200),
    tuning.psi = 1000)

## Conditional simulations
r.sim <- condsim(r.logzn.reml, newdata = meuse.grid.sgdf, nsim = 100, seed = 1,
    control = control.condsim(ncores = 1))
## using multiple cores
## r.sim <- condsim(r.logzn.reml, newdata = meuse.grid.sgdf, nsim = 100, seed = 1)
str(r.sim, max=2)

## Display
spplot(r.sim, zcol = "sim.1", at = seq(3.5, 8.5, by = 0.5))
spplot(r.sim, zcol = "sim.2", at = seq(3.5, 8.5, by = 0.5))

library(lattice)
levelplot(sim.1 ~ x + y, as.data.frame(r.sim), aspect = "iso", at = seq(3.5, 8.5, by = 0.5),
  panel = function(x, y, z, subscripts, data.points, ... ){
    panel.levelplot( x, y, z, subscripts, ...)
    panel.xyplot(data.points$x, data.points$y, col = 1)
  }, data.points = meuse[, c("x", "y")]
)
levelplot(sim.2 ~ x + y, as.data.frame(r.sim), aspect = "iso", at = seq(3.5, 8.5, by = 0.5),
  panel = function(x, y, z, subscripts, data.points, ... ){
    panel.levelplot( x, y, z, subscripts, ...)
    panel.xyplot(data.points$x, data.points$y, col = 1)
  }, data.points = meuse[, c("x", "y")]
)

}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
