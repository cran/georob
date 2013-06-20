% 2013-06-12 A. Papritz
% R CMD Rdconv -t html -o bla.html predict.georob.Rd ; open bla.html; R CMD Rd2pdf --force predict.georob.Rd; 

\encoding{macintosh}
\name{predict.georob}
\alias{predict.georob}
\title{Predict Method for Robustly Fitted Spatial Linear Models}

\description{Robust and customary external drift kriging prediction 
based on a spatial linear models fitted by \code{georob}. The
\code{predict} method for the class \code{georob} computes fitted values, point
and block kriging predictions as
well as model terms for display by \code{\link[stats]{termplot}}.
}

\usage{
\method{predict}{georob}(object, newdata, type =  c("signal", "response", "trend", "terms"), 
    terms = NULL, se.fit = TRUE, signif = 0.95, mmax = 10000, locations, 
    full.covmat = FALSE, pwidth = NULL, pheight = NULL, napp = 1, 
    extended.output = FALSE, ncores = detectCores(), verbose = 0, ...)
}

\arguments{
  \item{object}{an object of class \code{"georob"}, see \code{\link{georobObject}}.}
  
  \item{newdata}{an optional data frame, 
  \code{\link[sp]{SpatialPointsDataFrame}}, 
  \code{\link[sp]{SpatialPixelsDataFrame}}, 
  \code{\link[sp]{SpatialGridDataFrame}} or 
  \code{\link[sp]{SpatialPolygonsDataFrame}} in which to look for variables
  with which to compute fitted values or kriging predictions.\cr If
  \code{newdata} is a \code{\link[sp]{SpatialPolygonsDataFrame}} then
  block kriging predictions are computed, otherwise point kriging
  predictions.}
  
  \item{type}{character keyword defining what target quantity should be
  predicted (computed).  Possible values are
  
    \itemize{
    
      \item \code{"signal"}: the \dQuote{signal}
      \eqn{S(\mbox{\boldmath$s$\unboldmath}) =
      \mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath})^\mathrm{T}
      \mbox{\boldmath$\beta$\unboldmath} +
      B(\mbox{\boldmath$s$\unboldmath})}{Y(s) = x(s)^T \beta + B(s)} of
      the process (default),
      
      \item \code{"response"}: the observations
      \eqn{Y(\mbox{\boldmath$s$\unboldmath}) =
      S(\mbox{\boldmath$s$\unboldmath}) +
      \varepsilon(\mbox{\boldmath$s$\unboldmath}),}{Y(s)=S(s) + \epsilon(s),}
      
      \item \code{"trend"}: the external drift
      \eqn{\mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath})^\mathrm{T}
      \mbox{\boldmath$\beta$\unboldmath},}{x(s)^T \beta,}
      
      \item \code{"terms"}: the model terms.
    }
  }

  \item{terms}{If \code{type = "terms"}, which terms (default is all terms).}
  
  \item{se.fit}{logical, only used if \code{type} is equal to
  \code{"terms"}, see \code{\link[stats]{predict.lm}}.}
  
  
  \item{signif}{positive numeric equal to the tolerance or confidence level
  for computing respective intervals.}
  
  \item{mmax}{integer equal to the maximum number (default \code{10000}) of
  prediction items, computed in a sub-task, see \emph{Details}.}
  
  \item{locations}{an optional one-sided formula specifying what variables
  of \code{newdata} are the coordinates of the prediction points\cr
  (default: \code{object[["locatons.objects"]]$locations}).}
  
  \item{full.covmat}{logical controlling whether the full covariance matrix
  of the prediction errors is returned (\code{TRUE}) or only the vector
  with its diagonal elements (\code{FALSE}, default), see \emph{Value} for
  an explanation of the effect of \code{full.covmat}.}
  
  \item{pwidth, pheight, napp}{numeric scalars, used to tune numeric
  integration of semivariances for block kriging, see
  \code{\link[constrainedKriging]{preCKrige}}.}
  
  \item{extended.output}{logical controlling whether the covariance
  matrices of the kriging predictions and of the data should be computed, see
  \emph{Details} (default \code{FALSE}).}
  
  \item{ncores}{positive integer controlling how many cores are used for
  parallelized computations, see \emph{Details}.}

  \item{verbose}{positive integer controlling logging of diagnostic
  messages to the console.\cr  \code{verbose = 0} (default) largely suppresses
  such messages.}
  
  \item{\ldots}{additional arguments passed to method, currently not used.}
  
}

\details{
  The \code{predict} method for class \code{georob} uses the package
  \pkg{parallel} for parallelized computation of kriging predictions.  If
  there are \eqn{m} items to predict, the task is split into
  \code{ceiling(m/mmax)} sub-tasks that are then distributed to
  \code{ncores} CPUs.  Evidently, \code{ncores = 1} suppresses parallel
  execution.  By default, the function uses all available CPUs as returned
  by \code{\link[parallel]{detectCores}}.  \cr Note that if
  \code{full.covmat} is \code{TRUE} \code{mmax} must exceed \eqn{m} (and parallel
  execution is not possible).
  
  The argument \code{extended.output = TRUE} is used to compute all
  quantities required for (approximately) unbiased back-transformation of
  kriging predictions of log-transformed data to the original scale of the
  measurements by \code{\link{lgnpp}}.  In more detail, the following items
  are computed:
  
  \itemize{    
  
    \item \code{trend}: the fitted values,
    \eqn{\mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath})\mathrm{^T}\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{x(s)^T hat\beta},
    
    \item \code{var.pred}: the variances of the kriging predictions, 
    \eqn{\mathrm{Var}_{\hat{\theta}}[\widehat{Y}(\mbox{\boldmath$s$\unboldmath})]}{Var[hatY(s)]} or 
    \eqn{\mathrm{Var}_{\hat{\theta}}[\widehat{S}(\mbox{\boldmath$s$\unboldmath})]}{Var[hatS(s)]},

    \item \code{aux.cov.pred.target}: the covariances between the predictions and the
    prediction targets,\cr
    \eqn{\mathrm{Cov}_{\hat{\theta}}[\widehat{Y}(\mbox{\boldmath$s$\unboldmath}),Y(\mbox{\boldmath$s$\unboldmath})]}{Cov[hatY(s),Y(s)]} or 
    \eqn{\mathrm{Cov}_{\hat{\theta}}[\widehat{S}(\mbox{\boldmath$s$\unboldmath}),S(\mbox{\boldmath$s$\unboldmath})]}{Cov[hatS(s),S(s)]},
    
    
    \item \code{var.target}: the variances of the prediction targets
    \eqn{\mathrm{Var}_{\hat{\theta}}[Y(\mbox{\boldmath$s$\unboldmath})]}{Var[Y(s)]} or 
    \eqn{\mathrm{Var}_{\hat{\theta}}[S(\mbox{\boldmath$s$\unboldmath})]}{Var[S(s)]}.
    
  }
  
  Note that the component \code{var.pred} is also  present if 
  \code{type} is equal to \code{"trend"}, irrespective of the choice for \code{extended.output}.  
  This component contains then the variances of the fitted values.
  
}

\value{
  If \code{type} is equal to \code{"terms"} then a vector, a matrix, or a
  list with prediction results along with bounds and standard errors, see
  \code{\link[stats]{predict.lm}}.  Otherwise, the structure and contents
  of the output generated by \code{predict.georob} are determined by the
  class of \code{newdata} and the logical flags \code{full.covmat} and
  \code{extended.output}:
  
  
  If \code{full.covmat} is \code{FALSE} then the result is an object of the same
  class as \code{newdata} (data frame,
  \code{\link[sp]{SpatialPointsDataFrame}},
  \code{\link[sp]{SpatialPixelsDataFrame}}
  \code{\link[sp]{SpatialGridDataFrame}}, \cr
  \code{\link[sp]{SpatialPolygonsDataFrame}}).  The data frame or the
  \code{data} slot of the \code{Spatial...DataFrame} objects
  have the following components:
  
  \itemize{
  
    \item the coordinates of the prediction points (only present if
    \code{newdata} is a data frame).
    
    \item \code{pred}: the kriging predictions (or fitted values).
    
    \item \code{se}: the root mean squared prediction errors (kriging
    standard errors).
    
    \item \code{lower}, \code{upper}: the limits of tolerance/confidence
    intervals,
    
    \item \code{trend}, \code{var.pred}, \code{aux.cov.pred.target},
    \code{var.target}: only present if \code{extend.output} is \code{TRUE},
    see \emph{Details}.
    
  } 
  
  If \code{full.covmat} is \code{TRUE} then \code{predict.georob} returns a list
  with the following components:
  
  \itemize{
  
    \item \code{pred}: a data frame or a \code{Spatial...DataFrame} object
    as described above for\cr \code{full.covmat = FALSE}.
    
    \item \code{mse.pred}: the full covariance matrix of the prediction errors, 
    \eqn{Y(\mbox{\boldmath$s$\unboldmath})-\widehat{Y}(\mbox{\boldmath$s$\unboldmath})}{Y(s)-hatY(s)} or 
    \eqn{S(\mbox{\boldmath$s$\unboldmath})-\widehat{S}(\mbox{\boldmath$s$\unboldmath})}{s(s)-hatS(s)} 
    see \emph{Details}.
    
    \item \code{var.pred}: the full covariance matrix of the
    kriging predictions, see \emph{Details}.
    
    \item \code{aux.cov.pred.target}: the full covariance matrix of the
    predictions and the prediction targets, see \emph{Details}.
    
    \item \code{var.target}: the full covariance matrix of the
    prediction targets, see \emph{Details}.
    
  }
  
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}}

\references{
  Nussbaum, M., Papritz, A., Baltensweiler, A. and Walthert, L. (2012)
  \emph{Organic Carbon  Stocks of Swiss Forest Soils},
  Institute of Terrestrial Ecosystems, ETH Z?rich and 
  Swiss Federal Institute for Forest, Snow and Landscape Research
  (WSL), pp. 51. 
  \url{http://e-collection.library.ethz.ch/eserv/eth:6027/eth-6027-01.pdf}

  \enc{K?nsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A. (2011) Robust
  estimation of the external drift and the variogram of spatial data.
  Proceedings of the ISI 58th World Statistics Congress of the International
  Statistical Institute.
  \url{http://isi2011.congressplanner.eu/pdfs/650268.pdf}
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{georobObject}} for a description of the class \code{georob}.
}

\examples{
\dontrun{
data(meuse )

data(meuse.grid)
meuse.grid.pixdf <- SpatialPixelsDataFrame(points = meuse.grid[, 1:2],
    data = meuse.grid[, -(1:2)])
  
library(constrainedKriging)
data(meuse.blocks)

r.logzn.rob <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "exponential", param = c( variance = 0.15, nugget = 0.05, scale = 200 ),
    tuning.psi = 1., control = georob.control(cov.bhat = TRUE, full.cov.bhat = TRUE,
        cov.bhat.betahat = TRUE, aux.cov.pred.target = TRUE))
        
## point predictions of log(Zn)
r.pred.points <- predict(r.logzn.rob, newdata = meuse.grid.pixdf, extended.output = TRUE,
    full.covmat = TRUE)
str(r.pred.points$pred@data)

## back-transformation of point predictions
r.backtf.pred.points <- lgnpp(r.pred.points)
str(r.backtf.pred.points$pred@data)

spplot(r.backtf.pred.points[["pred"]], zcol = "lgn.pred", main = "Zn content")

## predicting mean Zn content for whole area
r.block <- lgnpp(r.pred.points, is.block = TRUE, all.pred = r.backtf.pred.points[["pred"]])
r.block

## block predictions of log(Zn)
r.pred.block <- predict(r.logzn.rob, newdata = meuse.blocks, extended.output = TRUE,
    pwidth = 75, pheight = 75)
r.backtf.pred.block <- lgnpp(r.pred.block, newdata = meuse.grid)

spplot(r.backtf.pred.block, zcol = "lgn.pred", main = "block means Zn content")}

}

\keyword{models}
\keyword{spatial}
\keyword{robust}
