% 2024-02-09 A. Papritz
% R CMD Rdconv -t html -o bla.html georobSimulation.Rd ; open bla.html;
% R CMD Rd2pdf --force georobSimulation.Rd; open georobSimulation.pdf
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  georobSimulation.Rd

\encoding{utf8}
\name{georobSimulation}
\alias{georobSimulation}
\alias{condsim}
\alias{control.condsim}

\title{Simulating Realizations of Gaussian Processes}

\description{This page documents the function \code{condsim} that
simulates (un)conditional realizations of Gaussian processes from the
parameters of a spatial linear model estimated by the function
\code{\link{georob}}.}

\usage{
  condsim(object, newdata, nsim, seed, type =  c("response", "signal"),
    locations, trend.coef = NULL,
    variogram.model = NULL, param = NULL, aniso = NULL, variogram.object = NULL,
    control = control.condsim(), verbose = 0)

  control.condsim(use.grid = FALSE, grid.refinement = 2.,
    condsim = TRUE, ce.method = c( "standard", "approximate" ),
    ce.grid.expansion = 1., include.data.sites = FALSE,
    means = FALSE, trend.covariates = FALSE, covariances = FALSE,
    ncores = 1, mmax = 10000, pcmp = control.pcmp())
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

  \item{nsim}{a positive interger with the number of (condititional)
  realizations to compute (mandatory argument).}

  \item{seed}{an integer seed to initialize random number generation,
  see \code{\link[base]{set.seed}} (mandatory argument).}

  \item{type}{a character keyword defining what target quantity should be
  simulated.  Possible values are

    \itemize{

      \item \code{"signal"}: the \dQuote{signal}
      \eqn{Z(\boldsymbol{s}) =
      \boldsymbol{x}(\boldsymbol{s})^\mathrm{T}
      \boldsymbol{\beta} +
      B(\boldsymbol{s})}{Z(s) = x(s)^T \beta + B(s)} of
      the process,

      \item \code{"response"}: the observations
      \eqn{Y(\boldsymbol{s}) =
      Z(\boldsymbol{s}) +
      \varepsilon(\boldsymbol{s}),}{Y(s)=Z(s) + \epsilon(s)}  (default),

    }

  see \code{\link{georobPackage}} for details on the model specification.}

  \item{locations}{an optional one-sided formula specifying what variables
  of \code{newdata} are the coordinates of the points for which simulations
  are computed (default\cr
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
  \code{grid.refinement}, \code{condsim},
  \code{ce.method}, \code{ce.grid.expansion},
  \code{include.data.sites},
  \code{means}, \code{trend.covariates},\cr  \code{covariances},
  \code{ncores}, \code{mmax} and \code{pcmp} or a function such as
  \code{control.condsim} that generates such a list, see
  arguments of \code{control.condsim} for details.}

  \item{verbose}{a positive integer controlling logging of diagnostic
  messages to the console.\cr \code{verbose = 0} (default) suppresses
  such messages.}

  \item{use.grid}{a logical scalar (default \code{FALSE}) to control whether
  (conditional) realizations are computed for a rectangular grid instead of
  the coordinates of points contained in \code{newdata}, see \emph{Details}.}

  \item{grid.refinement}{a numeric that defines a factor by which the
  minimum differences of the coordinates between any pair of points in
  \code{newdata} are divided to setup the simulation grid, should be > 1
  (default 2), see \emph{Details}.}

  \item{condsim}{a logical scalar to control whether conditional
  (\code{TRUE} default) or unconditional simulations (\code{FALSE}) are
  computed.  }

  \item{ce.method}{a character keyword to select the method to simulate
  realizations by the circulant embedding algorithm, see \emph{Details}.}

  \item{ce.grid.expansion}{a numeric with the factor by which the
  dimensions of the simulation grid is expanded in the circulant embedding
  algorithm. Should be \eqn{\ge 1} (default 1).}

  \item{include.data.sites}{a logical scalar, to control whether
  (conditionally) simulated values are computed also for the points of the
  original data set used to estimate the model parameters and contained in
  \code{object}.}

  \item{means}{a logical scalar, to control whether
  the (un)conditional means are included in the output.}

  \item{trend.covariates}{a logical scalar, to control whether
  the covariates required for the trend model are included in the output.}

  \item{covariances}{a logical scalar, to control whether the covariances
  between the points of the original data set used to estimate the model
  parameters (\code{attr gcvmat.d.d}, \cr \code{\link{compress}}ed matrix) and
  the covariances between the simulation and the original data points
  (\code{attr gcvmat.s.d}, matrix) are returned as attributes of the
  output.  Note that these covariances are only returned if
  \code{condsim = TRUE}.}

  \item{ncores}{a positive integer controlling how many cores are
  used for parallelized computations, defaults to 1.}

  \item{mmax}{a positive integer equal to the maximum number (default
  \code{10000}) of prediction items, computed in sub-tasks executed in
  parallel, see section \emph{Details} of \code{\link{predict.georob}}.}

  \item{pcmp}{a list of arguments, passed  e.g. to \code{\link{pmm}} or a
  function such as \code{control.pcmp} that generates such a list
  (see \code{\link{control.pcmp}} for allowed arguments).}

}


\details{

 \subsection{General}{

   \code{condsim} (conditionally) simulates from a Gaussian process that
   has a linear mean function with parameters
   \eqn{\boldsymbol{\beta}}{\beta} and an auto-correlation structure
   characterized by a parametric variogram model and variogram parameters
   \eqn{\tau^2}{\tau^2} and \eqn{\boldsymbol{\theta}}{\theta} (see
   \code{\link{georobPackage}} for the employed parametrization of the
   spatial linear model).  The parameters of the mean and auto-correlation
   function are either taken from the spatial linear model estimated by
   \code{\link{georob}} and passed by the argument
   \code{object} to \code{condsim} or from the optional arguments
   \code{trend.coef} (\eqn{\boldsymbol{\beta}}{\beta})
   and \code{variogram.model}, \code{param}, \code{aniso} or \cr
   \code{variogram.object} (\eqn{\tau^2}{\tau^2},
   \eqn{\boldsymbol{\theta}}{\theta}).

   Simulated values are computed for the points in \code{newdata} and
   optionally also for the data points in \code{object} if
   \code{include.data.sites = TRUE}.  Both unconditional and conditional
   simulations can be computed.  In the latter cases, the simulated values
   are always conditioned to the response data used to fit the spatial
   linear model by \code{\link{georob}} and contained in \code{object}.

}

\subsection{Unconditional simulation}{

  Unconditional realizations are either computed for the exact locations of
  the points in \code{newdata} (\code{use.grid = FALSE}), irrespective of
  the fact whether these are arranged on a regular grid or not.
  Simulations are then generated by the Cholesky matrix decomposition
  method (e.g. \cite{\enc{Chilès}{Chiles} and Delfiner, 1999, sec.
  7.2.2}).

  For \code{use.grid = TRUE} the points in \code{newdata} are matched to a
  rectangular simulation grid and the simulations are generated for all
  nodes of this grid by the circulant embedding method (\cite{Davis and
  Bryant, 2013}; \cite{Dietrich and Newsam, 1993}; \cite{Wood and Chan,
  1994}).  For large problems this approach may be substantially faster and
  less memory demanding than the Cholesky matrix decomposition method.

	For circulant embedding, first a rectangular simulation grid is
	constructed from the coordinates of the points in \code{newdata} and
	\code{object}.  The spacings of the simulation grid is equal to the
	minimum coordinate differences between any pair of points in
	\code{newdata}, divided by \code{grid.refinement}.  The spatial extent of
	the simulation grid is chosen such that it covers the bounding boxes of
	all points in \code{newdata} and \code{object}.  The points in
	\code{newdata} and \code{object} are then matched to the closest nodes of
	the simulation grid.  If the same node is assigned to a point in
	\code{object} and \code{newdata} then the point in \code{object} is kept
	and the concerned point in \code{newdata} is omitted.

  The rectangular simulation grid is then expanded to the larger circulant
  embedding grid, and the eigenvalues of the so-called \emph{base matrix}
  (= first row of the covariance matrix of the nodes of the circulant
  embedding grid with block circulant structure, see \cite{Davies and Bryant,
  2013)} are computed by fast discrete Fourier transform (FFT).  It may
  happen that some of the eigenvalues of the base matrix are negative.  The
  standard circulant embedding algorithm then fails.

  Two approaches are implemented in \code{condsim} to handle this
  situation:

  \itemize{

    \item First, one may use the \emph{approximate circulant embedding}
    method by choosing \code{ce.method = "approximate"}.  This sets the
    negative eigenvalues of the base matrix to zero and scales the
    eigenvalues, see \cite{Chan and Wood (1994, sec.  4}, choice \eqn{\rho =
    \rho_2}).

    \item Second, one may attempt to avoid the problem of negative
    eigenvalues by increasing the size of the simulation (and circulant
    embedding) grids.  This can be achieved by choosing a value \eqn{> 1}
    for the argument \code{ce.grid.expansion}, see respective parts in
    \cite{Dietrich and Newsam (1993, sec.  4)} and \cite{Wood and Chan
    (1994, sec.  3)}.

  }

  Note that the dimension of the simulation and embedding grids are chosen
  such that the number of nodes is a highly composite integer.  This allows
  efficient FFT.

}

\subsection{Conditional simulation}{

  For both the Cholesky matrix decomposition and the circulant embedding
  approach, simulations are conditioned to data by the Kriging method,
  see \cite{\enc{Chilès}{Chiles} and Delfiner, 1999, sec. 7.3}.

}


\subsection{Parallelized computations}{

  \code{condsim} uses the packages \pkg{parallel} and \pkg{snowfall} for
  parallelized computations.  Three tasks can be executed in parallel:

  \itemize{

    \item Computation of (generalized correlations), see
    \code{\link{control.pcmp}} how to do this.

    \item Computation of Kriging predictions required for conditional
    simulations, see section \emph{Details} of
    \code{\link{predict.georob}}.

    \item Fast Fourier transform of realizations of standard normal
    deviates generated for the nodes of the base matrix (see
    \cite{Davies and Bryant, 2013, steps 3--5 of algorithm)}.  If there are
    \code{nsim} realizations to simulate, the task is split into
    \code{ceiling(nsim / ncores)} sub-tasks that are then distributed to
    \code{ncores} CPUs.  Evidently, \code{ncores = 1} (default) suppresses
    parallel execution.

  }
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

  The function \code{control.condsim} returns a list with parameters to
  steer \code{condsim}, see arguments above.

}



\references{

  \enc{Chilès}{Chiles}, J.-P., Delfiner, P. (1999) \emph{Geostatistics
  Modeling Spatial Uncertainty}, Wiley, New York,
  \doi{10.1002/9780470316993}.

  Davies, T. M., Bryant, D. (2013) On circulant embedding for gaussian
  random fields in R, \emph{Journal of Statistical Software}, \bold{55},
  1--21, \doi{10.18637/jss.v055.i09}.

  Dietrich, C. R., Newsam, G. N. (1993) A fast and exact method for
  multidimensional gaussian stochastic simulations, \emph{Water Resources
  Research}, \bold{9}, 2861--2869, \doi{10.1029/93WR01070}.

  Wood, A. T. A., Chan, G. (1994) Simulation of stationary gaussian
  processes in \eqn{[0,1]^d}, \emph{Journal of Computational and Graphcal
  Statistics}, \bold{3}, 409--432, \doi{10.2307/1390903}.

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

%  \code{\link{georobSimulation}} for simulating realizations of a Gaussian process
%  from model fitted by \code{georob}; and finally
  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
  for robust estimation and modelling of sample variograms.
}

\examples{
data(meuse)
data(meuse.grid)

## convert to SpatialGridDataFrame
meuse.grid.sgdf <- meuse.grid
coordinates(meuse.grid.sgdf) <- ~ x + y
gridded(meuse.grid.sgdf) <- TRUE
fullgrid(meuse.grid.sgdf) <- TRUE

## Gaussian REML fit
r.logzn.reml <- georob(log(zinc) ~ sqrt(dist), data = meuse, 
    locations = ~ x + y, variogram.model = "RMexp",
    param = c(variance = 0.15, nugget = 0.05, scale = 200),
    tuning.psi = 1000)

## Unconditional simulations using circulant embedding on rectangular
## simulation grid
r.sim.1 <- condsim(r.logzn.reml, newdata = meuse.grid.sgdf, nsim = 2, seed = 1,
    control = control.condsim(use.grid = TRUE, condsim = FALSE))
spplot(r.sim.1, zcol = "sim.1", at = seq(3.5, 8.5, by = 0.5))

## Conditional simulations using circulant embedding
if(interactive()){
  ## example is run only in interactive session because cpu times exceeds 5 s
  r.sim.2 <- condsim(r.logzn.reml, newdata = meuse.grid.sgdf, nsim = 2, seed = 1,
      control = control.condsim(use.grid = FALSE, condsim = TRUE))
  spplot(r.sim.2, zcol = "sim.2", at = seq(3.5, 8.5, by = 0.5))
}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
