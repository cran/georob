% 2024-02-09 A. Papritz
% R CMD Rdconv -t html -o bla.html predict.georob.Rd ; open bla.html; R CMD Rd2pdf --force predict.georob.Rd;

\encoding{utf8}
\name{predict.georob}
\alias{predict.georob}
\alias{control.predict.georob}

\title{Predict Method for Robustly Fitted Spatial Linear Models}

\description{Robust and customary external drift Kriging prediction
based on a spatial linear models fitted by \code{georob}. The
\code{predict} method for the class \code{georob} computes fitted values, point
and block Kriging predictions as
well as model terms for display by \code{\link[stats]{termplot}}.
}

\usage{
\method{predict}{georob}(object, newdata, type =  c("signal", "response", "trend", "terms"),
    terms = NULL, se.fit = TRUE, signif = 0.95, locations,
    variogram.model = NULL, param = NULL, aniso = NULL, variogram.object = NULL,
    control = control.predict.georob(), verbose = 0, ...)

control.predict.georob(full.covmat = FALSE, extended.output = FALSE,
    mmax = 10000, ncores = pcmp[["max.ncores"]], pwidth = NULL, pheight = NULL,
    napp = 1, pcmp = control.pcmp())

}


\arguments{
  \item{object}{an object of class \code{"georob"} (mandatory argument),
  see \code{\link{georobObject}}.}

  \item{newdata}{an optional data frame,
  \code{\link[sp]{SpatialPointsDataFrame}},
  \code{\link[sp]{SpatialPixelsDataFrame}},
  \code{\link[sp]{SpatialGridDataFrame}},
  \code{\link[sp]{SpatialPolygonsDataFrame}} or an (optional) object of class
  \code{\link[sp]{SpatialPoints}},
  \code{\link[sp]{SpatialPixels}} or
  \code{\link[sp]{SpatialGrid}},
  in which to look for variables
  with which to compute fitted values or Kriging predictions, see \emph{Details}.\cr If
  \code{newdata} is a \code{\link[sp]{SpatialPolygonsDataFrame}} then
  block Kriging predictions are computed, otherwise point Kriging
  predictions.}

  \item{type}{a character keyword defining what target quantity should be
  predicted (computed).  Possible values are

    \itemize{

      \item \code{"signal"}: the \dQuote{signal}
      \eqn{Z(\boldsymbol{s}) =
      \boldsymbol{x}(\boldsymbol{s})^\mathrm{T}
      \boldsymbol{\beta} +
      B(\boldsymbol{s})}{Y(s) = x(s)^T \beta + B(s)} of
      the process (default),

      \item \code{"response"}: the observations
      \eqn{Y(\boldsymbol{s}) =
      Z(\boldsymbol{s}) +
      \varepsilon(\boldsymbol{s}),}{Y(s)=Z(s) + \epsilon(s),}

      \item \code{"trend"}: the external drift
      \eqn{\boldsymbol{x}(\boldsymbol{s})^\mathrm{T}
      \boldsymbol{\beta},}{x(s)^T \beta,}

      \item \code{"terms"}: the model terms.
    }
  }

  \item{terms}{If \code{type = "terms"}, which terms (default is all terms).}

  \item{se.fit}{a logical scalar, only used if \code{type} is equal to
  \code{"terms"}, see \code{\link[stats]{predict.lm}}.}


  \item{signif}{a positive numeric equal to the tolerance or confidence level
  for computing respective intervals.  If \code{NULL} no intervals are
  computed.}

  \item{locations}{an optional one-sided formula specifying what variables
  of \code{newdata} are the coordinates of the prediction points\cr
  (default: \code{object[["locations.objects"]][["locations"]]}).}

  \item{variogram.model}{an optional character keyword defining the
  variogram model to be used for Kriging, see \code{\link{georob}} and
  \emph{Details}.}

  \item{param}{an optional named numeric vector with values of the
  variogram parameters used for Kriging, see \code{\link{georob}} and
  \emph{Details}.}

  \item{aniso}{an optional named numeric vector with values of anisotropy
  parameters of a variogram used for Kriging, see \code{\link{georob}} and
  \emph{Details}.}

  \item{variogram.object}{an optional list that defines a possibly nested
  variogram model used for Kriging, see \code{\link{georob}} and
  \emph{Details}.}

  \item{control}{a list with the components \code{full.covmat},
  \code{extended.output}, \code{mmax}, \code{ncores}, \code{pwidth},
  \code{pheight}, \code{napp} and \code{pcmp} or a function such as
  \code{control.predict.georob} that generates such a list.}

  \item{full.covmat}{a logical scalar controlling whether the full
  covariance matrix of the prediction errors is returned (\code{TRUE}) or
  only the vector with its diagonal elements (\code{FALSE}, default), see
  \emph{Value} for an explanation of the effect of \code{full.covmat}.}

  \item{extended.output}{a logical scalar controlling whether the covariance
  matrices of the Kriging predictions and of the data should be computed, see
  \emph{Details} (default \code{FALSE}).}

  \item{mmax}{a positive integer equal to the maximum number (default
  \code{10000}) of prediction items, computed in a sub-task in parallelized
  computations, see \emph{Details}.}

  \item{ncores}{a positive integer controlling how many cores are used for
  parallelized computations, see \emph{Details}.}

  \item{pwidth, pheight, napp}{numeric scalars, used to tune numeric
  integration of semi-variances for block Kriging, see
  \code{\link[constrainedKriging]{preCKrige}}.}

  \item{pcmp}{a list of arguments passed to \code{\link{pmm}} and other
  functions that carry out parallelized computations or a function such as
  \code{\link{control.pcmp}} that generates such a list (see
  \code{\link{control.pcmp}} for allowed arguments).}

  \item{verbose}{a positive integer controlling logging of diagnostic
  messages to the console.\cr  \code{verbose = 0} (default) largely suppresses
  such messages.}

  \item{...}{arguments passed to \code{control.predict.georob}.}

}

\details{
  If \code{newdata} is an object of class \code{SpatialPoints},
  \code{SpatialPixels} or \code{SpatialGrid} then the drift model may only
  use the coordinates as covariates (universal Kriging).  Furthermore the
  names used for the coordinates in \code{newdata} must be the same as in
  \code{data} when creating \code{object} (argument \code{locations} of
  \code{predict.georob} should not be used).  Note that the result returned
  by \code{predict.georob} is then an object of class
  \code{SpatialPointsDataFrame}, \code{SpatialPixelsDataFrame} or
  \code{SpatialGridDataFrame}.

  The \code{predict} method for class \code{georob} uses the packages
  \pkg{parallel} and \pkg{snowfall} for parallelized
  computation of Kriging predictions.  If there are \eqn{m} items to
  predict, the task is split into \code{ceiling(m/mmax)} sub-tasks that are
  then distributed to \code{ncores} CPUs.  Evidently, \code{ncores = 1}
  suppresses parallel execution.  By default, the function uses all
  available CPUs as returned by \code{\link[parallel]{detectCores}}.  \cr
  Note that if \code{full.covmat} is \code{TRUE} \code{mmax} must exceed
  \eqn{m} (and parallel execution is not possible).

  The argument \code{extended.output = TRUE} is used to compute all
  quantities required for (approximately) unbiased back-transformation of
  Kriging predictions of log-transformed data to the original scale of the
  measurements by \code{\link{lgnpp}}.  In more detail, the following items
  are computed:

  \itemize{

    \item \code{trend}: the fitted values,
    \eqn{\boldsymbol{x}(\boldsymbol{s})\mathrm{^T}\widehat{\boldsymbol{\beta}}}{x(s)^T hat\beta},

    \item \code{var.pred}: the variances of the Kriging predictions,
    \eqn{\mathrm{Var}_{\hat{\theta}}[\widehat{Y}(\boldsymbol{s})]}{Var[hatY(s)]} or
    \eqn{\mathrm{Var}_{\hat{\theta}}[\widehat{Z}(\boldsymbol{s})]}{Var[hatS(s)]},

    \item \code{cov.pred.target}: the covariances between the predictions and the
    prediction targets,\cr
    \eqn{\mathrm{Cov}_{\hat{\theta}}[\widehat{Y}(\boldsymbol{s}),Y(\boldsymbol{s})]}{Cov[hatY(s),Y(s)]} or
    \eqn{\mathrm{Cov}_{\hat{\theta}}[\widehat{Z}(\boldsymbol{s}),Z(\boldsymbol{s})]}{Cov[hatS(s),Z(s)]},


    \item \code{var.target}: the variances of the prediction targets
    \eqn{\mathrm{Var}_{\hat{\theta}}[Y(\boldsymbol{s})]}{Var[Y(s)]} or
    \eqn{\mathrm{Var}_{\hat{\theta}}[Z(\boldsymbol{s})]}{Var[Z(s)]}.

  }

  Note that the component \code{var.pred} is also  present if
  \code{type} is equal to \code{"trend"}, irrespective of the choice for \code{extended.output}.
  This component contains then the variances of the fitted values.

}

\value{

  The method \code{predict.georob} returns, depending on its arguments, the
  following objects:

  If \code{type} is equal to \code{"terms"} then a vector, a matrix, or a
  list with prediction results along with bounds and standard errors, see
  \code{\link[stats]{predict.lm}}.  Otherwise, the structure and contents
  of the output generated by \code{predict.georob} are determined by the
  class of \code{newdata} and the logical flags \code{full.covmat} and
  \code{extended.output}:


  If \code{full.covmat} is \code{FALSE} then the result is an object of a "similar"
  class as \code{newdata} (data frame,
  \code{\link[sp]{SpatialPointsDataFrame}},
  \code{\link[sp]{SpatialPixelsDataFrame}}
  \code{\link[sp]{SpatialGridDataFrame}}, \cr
  \code{\link[sp]{SpatialPolygonsDataFrame}}).

  The data frame or the
  \code{data} slot of the \code{Spatial...DataFrame} objects
  have the following components:

  \itemize{

    \item the coordinates of the prediction points (only present if
    \code{newdata} is a data frame).

    \item \code{pred}: the Kriging predictions (or fitted values).

    \item \code{se}: the root mean squared prediction errors (Kriging
    standard errors).

    \item \code{lower}, \code{upper}: the limits of tolerance/confidence
    intervals,

    \item \code{trend}, \code{var.pred}, \code{cov.pred.target},
    \code{var.target}: only present if \code{extended.output} is \code{TRUE},
    see \emph{Details}.

  }

  If \code{full.covmat} is \code{TRUE} then \code{predict.georob} returns a list
  with the following components:

  \itemize{

    \item \code{pred}: a data frame or a \code{Spatial...DataFrame} object
    as described above for\cr \code{full.covmat = FALSE}.

    \item \code{mse.pred}: the full covariance matrix of the prediction errors,
    \eqn{Y(\boldsymbol{s})-\widehat{Y}(\boldsymbol{s})}{Y(s)-hatY(s)} or
    \eqn{Z(\boldsymbol{s})-\widehat{Z}(\boldsymbol{s})}{s(s)-hatS(s)}
    see \emph{Details}.

    \item \code{var.pred}: the full covariance matrix of the
    Kriging predictions, see \emph{Details}.

    \item \code{cov.pred.target}: the full covariance matrix of the
    predictions and the prediction targets, see \emph{Details}.

    \item \code{var.target}: the full covariance matrix of the
    prediction targets, see \emph{Details}.

  }

  The function \code{control.predict.georob} returns a list with control
  parameters to steer \code{predict.georob}, see arguments of the
  function above for its components.

}

\references{
  Nussbaum, M., Papritz, A., Baltensweiler, A. and Walthert, L. (2014)
  Estimating soil organic carbon stocks of {Swiss} forest soils by robust
  external-drift kriging.  \emph{Geoscientific Model Development},
  \bold{7}, 1197--1210.  \doi{10.5194/gmd-7-1197-2014}.

  \enc{Künsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A. (2011) Robust
  estimation of the external drift and the variogram of spatial data.
  Proceedings of the ISI 58th World Statistics Congress of the International
  Statistical Institute.
  \doi{10.3929/ethz-a-009900710}
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

%  \code{\link{predict.georob}} for computing robust Kriging predictions;
%  \code{\link{validate.predictions}} for validating Kriging predictions;
  \code{\link{lgnpp}} for unbiased back-transformation of Kriging prediction
  of log-transformed data;

  \code{\link{georobSimulation}} for simulating realizations of a Gaussian process
  from model fitted by \code{georob}; and finally

  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
  for robust estimation and modelling of sample variograms.
}

\examples{
data(meuse)

data(meuse.grid)
coordinates(meuse.grid) <- ~x+y
meuse.grid.pixdf <- meuse.grid
gridded(meuse.grid.pixdf) <- TRUE

data(meuse.blocks, package = "constrainedKriging")

r.logzn.rob <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "RMexp", param = c(variance = 0.15, nugget = 0.05, scale = 200),
    tuning.psi = 1., control = control.georob(cov.bhat = TRUE, full.cov.bhat = TRUE))

## point predictions of log(Zn)
r.pred.points.1 <- predict(r.logzn.rob, newdata = meuse.grid.pixdf,
    control = control.predict.georob(extended.output = TRUE))
str(r.pred.points.1, max = 3)

## back-transformation of point predictions
r.backtf.pred.points <- lgnpp(r.pred.points.1)
str(r.backtf.pred.points, max = 3)

spplot(r.backtf.pred.points, zcol = "lgn.pred", main = "Zn content")

## predicting mean Zn content for whole area
if(interactive()){
  ## example is run only in interactive session because cpu times exceeds 5 s
  ## recompute point predictions with argument full.covmat = TRUE
  r.pred.points.2 <- predict(r.logzn.rob, newdata = meuse.grid.pixdf,
      control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE))
  str(r.pred.points.2, max = 3)
  r.block <- lgnpp(r.pred.points.2, is.block = TRUE, all.pred = r.backtf.pred.points@data)
  r.block
}

## block predictions of log(Zn)
if(interactive()){
  ## example is run only in interactive session because cpu times exceeds 5 s
  r.pred.block <- predict(r.logzn.rob, newdata = meuse.blocks,
      control = control.predict.georob(extended.output = TRUE,
          pwidth = 75, pheight = 75, mmax = 50))
  r.backtf.pred.block <- lgnpp(r.pred.block, newdata = meuse.grid)

  spplot(r.backtf.pred.block, zcol = "lgn.pred", main = "block means Zn content")
}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
