% 2012-12-14 A. Papritz
% R CMD Rdconv -t html -o bla.html validate.predictions.Rd ; open bla.html; R CMD Rd2pdf --force validate.predictions.Rd; 

\encoding{macintosh}
\name{validate.predictions}
\alias{validate.predictions}
\alias{plot.cv.georob}
\alias{print.cv.georob}
\alias{rstudent.cv.georob}
\alias{summary.cv.georob}
\alias{print.summary.cv.georob}

\title{Summary Statistics of (Cross-)Validation Prediction Errors}

\description{
  Functions to compute and plot summary statistics of prediction errors to
  (cross-)validate fitted spatial linear models by the criteria proposed by
  Gneiting et al.  (2007) for assessing probabilistic forecasts.}

\usage{
validate.predictions(data, pred, se.pred, 
    statistic = c("pit", "mc", "bs", "st"), ncutoff = NULL)
    
\method{plot}{cv.georob}(x, type = c("sc", "lgn.sc", "ta", "qq", "pit", "mc", "bs"), 
    ncutoff = NULL, add = FALSE, col, pch, lty, main, xlab, ylab, ...)
 
\method{print}{cv.georob}(x, digits = max(3, getOption("digits") - 3), ...)   
 
\method{rstudent}{cv.georob}(model, ...)   
 
\method{summary}{cv.georob}(object, se = FALSE, ...)   
}

\arguments{
  \item{data}{a numeric vector with observations about a response.}

  \item{pred}{a numeric vector with predictions for the response.}

  \item{se.pred}{a numeric vector with prediction standard errors.}

  \item{statistic}{character keyword defining what statistic of the
  prediction errors should be computed.  
  Possible values are (see \emph{Details}):
  
    \itemize{
    
      \item \code{"pit"}: probability integral transform (default),
      
      \item \code{"mc"}: average predictive distribution (marginal
      calibration),
      
      \item \code{"bs"}: Brier score,
      
      \item \code{"st"}: mean and dispersion statistics of (standardized)
      prediction errors.
    
    }
  }
  
  \item{ncutoff}{positive integer (\eqn{N}) giving the number of quantiles,
  for which CDFs are evaluated (\code{type = "mc"}), or the number of
  thresholds for which the Brier score is computed (\code{type = "bs"}),
  see \emph{Details} (default: \code{min(500, length(data))}).}
  
  \item{x, model, object}{objects of class \code{cv.georob}.}
  
  \item{digits}{positive integer indicating the number of decimal digits to print.}

  \item{type}{character keyword defining what type of plot is created by the
  \code{plot.cv.georob}.  Possible values are:
  
    \itemize{
    
      \item \code{"sc"}: a scatterplot of the (possibly log-transformed) response
      vs.  the respective predictions (default).
      
      \item \code{"lgn.sc"}: a scatterplot of the untransformed response
      against back-\cr transformed predictions of the log-transformed response.
      
      \item \code{"ta"}: Tukey-Anscombe plot (plot of standardized prediction
      errors vs. predictions).
      
      \item \code{"qq"}: normal QQ plot of standardized prediction errors.
      
      \item \code{"pit"}: histogram of probability integral transform, see
      \emph{Details}.
      
      \item \code{"mc"}: a marginal calibration plot, see \emph{Details},
      
      \item \code{"bs"}: a plot of Brier score vs. threshold, see
      \emph{Details}.
    
    }
  
  }
  
  \item{se}{logical controlling if the standard errors of the averaged
  continuous ranked probability score and of the mean and dispersion
  statistics of the prediction errors (see \emph{Details}) are computed
  from the respective values computed for the \eqn{K} cross-validation
  subsets (default: \code{FALSE}).}
  
  \item{add}{logical controlling wether the current high-level plot is
  added to an existing graphics without cleaning the frame before (default:
  \code{FALSE}).}
  
  \item{main, xlab, ylab}{title and axes labels of plot.}
  
  \item{col, pch, lty}{color, symbol and line type.}
  
  \item{...}{additional arguments passed to the methods.}
  
}

\details{

  \code{validate.predictions} computes the items required to evaluate (and
  plot) the diagnostic criteria proposed by Gneiting et al.  (2007) for
  assessing the \emph{calibration} and the \emph{sharpness} of
  probabilistic predictions.  To this aim, \code{validate.predictions} uses
  the assumption that the prediction errors
  \eqn{Y(\mbox{\boldmath$s$\unboldmath})-\widehat{Y}(\mbox{\boldmath$s$\unboldmath})}{Y(s)-hatY(s)}
  follow normal distributions with zero mean and standard deviations equal
  to the kriging standard errors.  This assumption is used to compute
  
  \itemize{
  
    \item the \emph{probability integral transform} (PIT),
    
    \deqn{\mathrm{PIT}_i = F_i(y_i),}{PIT_i = F_i(y_i),}
    
    where \eqn{F_i(y_i)} denotes the predictive CDF evaluted at \eqn{y_i},
    the value of the \eqn{i}th (cross-)validation datum,
    
    \item the \emph{average predictive CDF} 
    
    \deqn{\bar{F}_n(y)=1/n \sum_{i=1}^n F_i(y),}{barF_n(y)=1/n \sum_{i=1}^n F_i(y),} 
    
    where \eqn{n} is the number of (cross-)validation observations and the
    \eqn{F_i} are evaluated at \eqn{N} quantiles equal to the set of
    distinct \eqn{y_i} (or a subset of size \eqn{N} of them),
    
    \item the \emph{Brier Score}
    
    \deqn{\mathrm{BS}(y) = 1/n \sum_{i=1}^n \left(F_i(y) - I(y_i \leq y) \right)^2,}{BS(y) = 1/n \sum_{i=1}^n (F_i(y) - I(y_i \leq y) )^2,}
    
    where \eqn{I(x)} is the indicator function for the event \eqn{x}, and
    the Brier score is again evaluated at the unique values of the (cross-)validation
    observations (or a subset of size \eqn{N} of them),
    
    \item the \emph{averaged continuous ranked probability score}, CRPS, a
    strictly proper scoring criterion to rank predictions, which is related
    to the Brier score by

    \deqn{\mathrm{CRPS} = \int_{-\infty}^\infty \mathrm{BS}(y) \,dy.}{CRPS = \int_{-\infty}^\infty BS(y) dy.}

  }
  
  Gneiting et al.  (2007) proposed the following plots to validate
  probabilistic predictions:
  
  \itemize{
  
    \item A histogram of the PIT values.  For ideal predictions, with
    observed coverages of prediction intervals matching  nominal
    coverages, the PIT values have an uniform distribution.
  
   \item Plots of \eqn{\bar{F}_n(y)}{barF_n(y)} and of the empirical CDF of
   the data, say \eqn{\widehat{G}_n(y)}{hat{G}_n(y)}, and of their
   difference, \eqn{\bar{F}_n(y)-\widehat{G}_n(y)}{barF_n(y)-hat{G}_n(y)}
   vs \eqn{y}.  The forecasts are said to be \emph{marginally calibrated}
   if \eqn{\bar{F}_n(y)}{barF_n(y)} and \eqn{\widehat{G}_n(y)}{hat{G}_n(y)}
   match.
   
   \item A plot of \eqn{\mathrm{BS}(y)}{BS(y)} vs.  \eqn{y}.  Probabilistic
   predictions are said to be \emph{sharp} if the area under this curve,
   which equals CRPS, is minimized.
    
  }
  
  The \code{plot} method for class \code{cv.georob} allows to create
  these plots, along with scatterplots of observations and predictions, 
  Tukey-Anscombe and normal QQ plots of the standardized prediction
  errors.
  
  \code{summary.cv.georob} computes the mean and dispersion statistics
  of the (standardized) prediction errors (by a call to
  \code{validate.prediction} with argument \code{statistic = "st"}, see
  \emph{Value}) and the averaged continuous ranked probability score
  (\code{crps}).  If present in the \code{cv.georob} object, the error
  statistics are also computed for the errors of the unbiasedly
  back-transformed predictions of a log-transformed response.  If \code{se}
  is \code{TRUE} then these statistics are evaluated separately for the
  \eqn{K} cross-validation subsets and the standard errors of the means of
  these statistics are returned as well.
  
  The \code{print} method for class \code{cv.georob} returns the mean
  and dispersion statistics of the (standardized) prediction errors.
  
  The method \code{rstudent} returns for class \code{cv.georob} the
  standardized prediction errors.
  
}


\value{
  Depending on the argument \code{statistic}, the function
  \code{validate.predictions} returns
  
  \itemize{
  
    \item a numeric vector of PIT values if \code{statistic} is equal to \code{"pit"},
    
    \item a named numeric vector with summary statistics of the
    (standardized) prediction errors if \code{statistic} is equal to \code{"st"}.  The
    following statistics are computed:

    \tabular{rl}{%
      \code{me}     \tab mean prediction error \cr
      \code{mede}   \tab median prediction error \cr
      \code{rmse}   \tab root mean squared prediction error \cr
      \code{made}   \tab median absolute prediction error \cr
      \code{qne}    \tab Qn dispersion measure of prediction errors 
      (see \code{\link[robustbase]{Qn}}) \cr
      \code{msse}   \tab mean squared standardized prediction error \cr
      \code{medsse} \tab median squared standardized prediction error \cr
    }
    
    \item a data frame if \code{statistic} is equal to \code{"mc"} or
    \code{"bs"} with the components (see \emph{Details}):
    
    \tabular{rl}{ 
    \code{z} \tab the sorted unique (cross-)validation
      observations (or a subset of size 
    \code{ncutoff} of them)\cr
      \code{ghat} \tab the empirical CDF of the (cross-)validation
      observations \eqn{\widehat{G}_n(y)}{hat{G}_n(y)}\cr 
    \code{fbar} \tab  the average predictive distribution \eqn{\bar{F}_n(y)}{barF_n(y)}\cr
    \code{bs} \tab the Brier score \eqn{\mathrm{BS}(y)}{BS(y)}\cr }
  }
  
  The function \code{rstudent.cv.georob} returns a numeric vector with
  the standardized cross-validation prediction errors.
  
}


\references{
 Gneiting, T., Balabdaoui, F. and Raftery, A. E.(2007) Probabilistic
 forecasts, calibration and sharpness. \emph{Journal of the Royal Statistical
 Society Series B} \bold{69}, 243--268.

}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}}

\seealso{ 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{cv.georob}} for assessing the goodness of a fit by \code{georob}.
}

\examples{
\dontrun{
data( meuse )

r.logzn <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "exponential",
    param = c( variance = 0.15, nugget = 0.05, scale = 200 ),
    tuning.psi = 1)

r.logzn.cv.1 <- cv(r.logzn, seed = 1, lgn = TRUE )
r.logzn.cv.2 <- cv(r.logzn, formula = .~. + ffreq, seed = 1, lgn = TRUE )

summary(r.logzn.cv.1, se = TRUE)
summary(r.logzn.cv.2, se = TRUE)

op <- par(mfrow = c(2,2))
plot(r.logzn.cv.1, type = "lgn.sc")
plot(r.logzn.cv.2, type = "lgn.sc", add = TRUE, col = "red")
abline(0, 1, lty= "dotted")
plot(r.logzn.cv.1, type = "ta")
plot(r.logzn.cv.2, type = "ta", add = TRUE, col = "red")
abline(h=0, lty= "dotted")
plot(r.logzn.cv.2, type = "mc", add = TRUE, col = "red")
plot(r.logzn.cv.1, type = "bs")
plot(r.logzn.cv.2, type = "bs", add = TRUE, col = "red")
legend("topright", lty = 1, col = c( "black", "red"), bty = "n",
    legend = c("log(Zn) ~ sqrt(dist)", "log(Zn) ~ sqrt(dist) + ffreq"))
par(op)}
}

\keyword{models}
\keyword{spatial}
