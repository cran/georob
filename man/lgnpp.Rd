% 2024-01-26 A. Papritz
% R CMD Rdconv -t html -o bla.html lgnpp.Rd ; open bla.html; R CMD Rd2pdf --force lgnpp.Rd;

\encoding{utf8}
\name{lgnpp}
\alias{lgnpp}

\title{Unbiased Back-Transformations for Log-normal Kriging}

\description{
  The function \code{lgnpp} back-transforms point or block Kriging
  predictions of a log-transformed response variable computed by
  \code{\link{predict.georob}}.  Alternatively, the function averages
  log-normal point Kriging predictions for a block and approximates the mean
  squared prediction error of the block mean.}



\usage{
lgnpp(object, newdata, locations, is.block = FALSE, all.pred = NULL,
    extended.output = FALSE)
}

\arguments{
  \item{object}{an object with Kriging predictions of a log-transformed
  response variable as obtained by
  \code{\link{predict}(\var{georob-object, ...})}.}

  \item{newdata}{an optional object as passed as argument \code{newdata} to
  \code{\link{predict.georob}}, see \emph{Details}.}

  \item{locations}{an optional one-sided formula specifying what variables
  of \code{newdata} are the coordinates of the prediction points, see
  \code{\link{predict.georob}}.}

  \item{is.block}{an optional logical scalar (default \code{FALSE}) specifying
  whether point predictions contained in \code{object} are considered to
  belong to a single block and should be averaged after
  back-transformation.  Ignored if \code{object} contains block Kriging
  predictions, see \emph{Details}.}

  \item{all.pred}{an optional positive integer or an object as obtained by
  \code{lgnpp(predict(\var{georob-object, ...}))}, see \emph{Details}.}

  \item{extended.output}{a logical scalar controlling whether the
  covariance matrix of the errors of the back-transformed point predictions
  is added as an attribute to the result, see \emph{Details}.}

}
\details{
 The function \code{lgnpp} performs three tasks:

  \subsection{1. Back-transformation of point Kriging predictions of a
  log-transformed response}{

    The usual, marginally unbiased back-transformation for log-normal point
    Kriging is used:

    \deqn{\widehat{U}(\boldsymbol{s}) = \exp( \widehat{Z}(\boldsymbol{s}) +
      1/2 (  \mathrm{Var}_{\hat{\theta}}[ Z(\boldsymbol{s})]
      - \mathrm{Var}_{\hat{\theta}}[\widehat{Z}(\boldsymbol{s})])),}{
      hatU(s)=exp(hatS(s) + 1/2(hat\sigma^2_n + hat\sigma^2 - Var[hatS(s)])),}

    \deqn{\mathrm{Cov}_{\hat{\theta}}[
      U(\boldsymbol{s}_i) - \widehat{U}(\boldsymbol{s}_i),
      U(\boldsymbol{s}_j) - \widehat{U}(\boldsymbol{s}_j)
      ] = \mu_{\hat{\theta}}(\boldsymbol{s}_i) \mu_{\hat{\theta}}(\boldsymbol{s}_j)
    }{
      Cov[ U(s_i)-hatU(s_i), U(s_j)-hatU(s_j)] = \mu(s_i) \mu(s_j)
    }
    \deqn{
      \times \{
        \exp(\mathrm{Cov}_{\hat{\theta}}[Z(\boldsymbol{s}_i),Z(\boldsymbol{s}_j)])
        -2\exp(\mathrm{Cov}_{\hat{\theta}}[\widehat{Z}(\boldsymbol{s}_i),Z(\boldsymbol{s}_j)])
        +\exp(\mathrm{Cov}_{\hat{\theta}}[\widehat{Z}(\boldsymbol{s}_i),\widehat{Z}(\boldsymbol{s}_j)])
      \},
    }{
      *\{
      exp(Cov[Y(s_i), Y(s_j)]) - 2 exp(Cov[hatY(s_i), Y(s_j)]) + exp(Cov[hatY(s_i), hatY(s_j)])
      \},
    }

    where \eqn{\widehat{Z}}{hatY} and \eqn{\widehat{U}}{hatU} denote the
    log- and back-transformed predictions of the signal,
    and

    \deqn{\mu_{\hat{\theta}}(\boldsymbol{s}) \approx
      \exp(\boldsymbol{x}(\boldsymbol{s})\mathrm{^T}\widehat{\boldsymbol{\beta}}
      + 1/2 \mathrm{Var}_{\hat{\theta}}[Z(\boldsymbol{s})]).
    }{\mu(s) = exp( x(s)^T hat\beta + 1/2 Var[Z(s)]).}

    The expressions for the required covariance terms can be found in the
    Appendices of \cite{Nussbaum et al.  (2014)}.  Instead of the signal
    \eqn{Z(\boldsymbol{s})}{Z(s)}, predictions of the
    log-transformed response \eqn{Y(\boldsymbol{s})}{Z(s)}
    or the estimated trend
    \eqn{\boldsymbol{x}(\boldsymbol{s})^\mathrm{T}\widehat{\boldsymbol{\beta}}}{x(s)^T
    hat\beta} of the log-transformed data can be back-transformed (see
    \code{\link{georobPackage}}).  The
    above transformations are used if \code{object}
    contains point Kriging predictions (see \code{predict.georob},
    \emph{Value}) and if \code{is.block = FALSE} and \code{all.pred} is
    missing.

  }

  \subsection{2. Back-transformation of block Kriging predictions of a
  log-transformed response}{

  Block Kriging predictions of a log-transformed response variable are
  back-transformed by the approximately unbiased transformation proposed
  by \cite{Cressie (2006, Appendix C)}

  \deqn{\widehat{U}(A) = \exp( \widehat{Z}(A) + 1/2 \{
    \mathrm{Var}_{\hat{\theta}}[Z(\boldsymbol{s})] + \widehat{\boldsymbol{\beta}}\mathrm{^T}
    \boldsymbol{M}(A) \widehat{\boldsymbol{\beta}} -
    \mathrm{Var}_{\hat{\theta}}[\widehat{Z}(A)]
    \}),
  }{
    hatU(A) = exp( hatS(A) + 1/2 \{Var[Z(s)] + hat\beta^T M(A) hat\beta - Var[hatS(A)]\}),
  }

  \deqn{\mathrm{E}_{\hat{\theta}}[\{U(A) - \widehat{U}(A))^2] = \mu_{\hat{\theta}}(A)^2 \{
    \exp(\mathrm{Var}_{\hat{\theta}}[Z(A)]) - 2 \exp(\mathrm{Cov}_{\hat{\theta}}[\widehat{Z}(A),Z(A)]) + \exp(\mathrm{Var}_{\hat{\theta}}[\widehat{Z}(A)])
    \}
  }{
    E[ \{ U(A)- hatU(A) \}^2 ] = \mu(A)^2 \{
    \exp(Var[Z(A)]) - 2 \exp(Cov[Z(A),hatS(A)]) + \exp(Var[hatS(A)]) \}
  }

  where \eqn{\widehat{Z}(A)} and \eqn{\widehat{U}(A)} are the log- and
  back-transformed predictions of the block mean \eqn{U(A)}, respectively,
  \eqn{\boldsymbol{M}(A)}{M(A)} is the spatial
  covariance matrix of the covariates

  \deqn{ \boldsymbol{M}(A) = 1/|A| \int_A
  ( \boldsymbol{x}(\boldsymbol{s}) - \boldsymbol{x}(A) )
  ( \boldsymbol{x}(\boldsymbol{s}) - \boldsymbol{x}(A) )\mathrm{^T} \,d\boldsymbol{s}
  }{M(A) = 1/|A| int_B (x(s)-x(A) (x(s)-x(A))^T) ds}

  within the block \eqn{A} where

  \deqn{ \boldsymbol{x}(A) = 1/|A| \int_A \boldsymbol{x}(\boldsymbol{s}) \,d\boldsymbol{s}
  }{x(A) = 1/|A| int_B x(s) ds}

  and

  \deqn{ \mu_{\hat{\theta}}(A) \approx \exp(\boldsymbol{x}(A)\mathrm{^T}
    \widehat{\boldsymbol{\beta}} + 1/2 \mathrm{Var}_{\hat{\theta}}[Z(A)]).
  }{
    \mu(A) = exp( x(A)^T hat\beta + 1/2 Var[Z(A)]).
  }

  This back-transformation is based on the assumption that both the point data
  \eqn{U(\boldsymbol{s})}{U(s)} and the block means
  \eqn{U(A)} follow  log-normal laws, which strictly cannot hold.  But
  for small blocks the assumption works well as the bias and the loss of
  efficiency caused by this assumption are small (\cite{Cressie, 2006};
  \cite{Hofer et al., 2013}).

  The above formulae are used by \code{lgnpp} if \code{object} contains
  block Kriging predictions in the form of a
  \code{\link{SpatialPolygonsDataFrame}}.  To approximate
  \eqn{\boldsymbol{M}(A)}{M(A)}{}, one needs the covariates
  on a fine grid for the whole study domain in which the blocks lie.  The
  covariates are passed \code{lgnpp} as argument \code{newdata}, where
  \code{newdata} can be any spatial data frame accepted by
  \code{predict.georob}.  For evaluating
  \eqn{\boldsymbol{M}(A)}{M(A)} the geometry of the blocks
  is taken from the \code{polygons} slot of the\cr
  \code{SpatialPolygonsDataFrame} passed as \code{object} to \code{lgnpp}.

  }

  \subsection{3. Back-transformation and averaging of point Kriging predictions
  of a log-transformed response}{

  \code{lgnpp} allows as a further option to back-transform and
  \emph{average} point Kriging predictions passed as \code{object} to the
  function.  One then assumes that the predictions in \code{object} refer
  to points that lie in \emph{a single} block.  Hence, one uses the
  approximation

  \deqn{\widehat{U}(A) \approx \frac{1}{K} \sum_{s_i \in A} \widehat{U}(\boldsymbol{s}_i)
  }{hatU(A) = 1/K sum_{s_i in A} hatU(s_i)}

  to predict the block mean \eqn{U(A)}, where \eqn{K} is the number of
  points in \eqn{A}. The mean squared prediction error can be approximated by

  \deqn{\mathrm{E}_{\hat{\theta}}[\{U(A) - \widehat{U}(A)\}^2] \approx \frac{1}{K^2}
    \sum_{s_i \in A} \sum_{s_j \in A}
    \mathrm{Cov}_{\hat{\theta}}[
    U(\boldsymbol{s}_i) - \widehat{U}(\boldsymbol{s}_i),
    U(\boldsymbol{s}_j) - \widehat{U}(\boldsymbol{s}_j)
    ].
    }{E[\{U(A) - hatU(A)\}^2] = 1/K^2 sum_{s_i in A} sum_{s_j in A}
    Cov[ U(s_i)-hatU(s_i), U(s_j)-hatU(s_j)].
  }

 In most instances, the evaluation of the above double sum is not feasible
 because a large number of points is used to discretize the block \eqn{A}.
 \code{lgnpp} then uses the following approximations to compute the mean
 squared error (see also Appendix E of \cite{Nussbaum et al., 2014}):

  \itemize{

    \item Point prediction results are passed as \code{object} to \code{lgnpp}
    only for a \emph{random sample of points in \eqn{A}} (of size \eqn{k}),
    for which the evaluation of the above double sum is feasible.

    \item The prediction results for the \emph{complete set of points}
    within the block are passed as argument \code{all.pred} to
    \code{lgnpp}.  These results are used to compute \eqn{\widehat{U}(A)}{U(A)}.

    \item The mean squared error is then approximated by

    \deqn{
      \mathrm{E}_{\hat{\theta}}[\{U(A) - \widehat{U}(A)\}^2] \approx
      \frac{1}{K^2} \sum_{s_i \in A} \mathrm{E}_{\hat{\theta}}[ \{ U(\boldsymbol{s}_i) - \widehat{U}(\boldsymbol{s}_i)\}^2]
    }{E[\{U(A) - hatU(A)\}^2] = 1/K^2 sum_{s_i in A} E[ \{ U(s_i) - hatU(s_i)\}^2] }

    \deqn{+ \frac{K-1}{K k (k-1)} \sum_{s_i \in \mathrm{sample}}\sum_{s_j \in \mathrm{sample}, s_j \neq s_i}
    \mathrm{Cov}_{\hat{\theta}}[
    U(\boldsymbol{s}_i) - \widehat{U}(\boldsymbol{s}_i),
    U(\boldsymbol{s}_j) - \widehat{U}(\boldsymbol{s}_j)
    ].
    }{ + (K-1)/(K k (k-1)) sum_{s_i in sample} sum_{s_j in sample, s_j != s_i}
      Cov[ U(s_i)-hatU(s_i), U(s_j)-hatU(s_j)]
    }

    The first term of the RHS (and \eqn{\widehat{U}(A)}{hatU(A)}) can be
    computed from the point Kriging results contained in \code{all.pred},
    and the double sum is evaluated from the full covariance matrices of
    the predictions and the respective targets, passed to \code{lgnpp} as
    \code{object} (one has to use the arguments
    \code{control=control.predict.georob(full.covmat=TRUE)} for
    \code{\link{predict.georob}} when computing the point Kriging
    predictions stored in \code{object}).

    \item If the prediction results are not available for the complete set
    of points in \eqn{A} then \code{all.pred} may be equal to \eqn{K}.  The
    block mean is then approximated by

    \deqn{\widehat{U}(A) \approx \frac{1}{k} \sum_{s_i \in \mathrm{sample}}
    \widehat{U}(\boldsymbol{s}_i) }{hatU(A) = sum_{s_i in
    A} hatU(s_i)}

    and the first term of the RHS of the expression for the mean squared
    error by

      \deqn{ \frac{1}{kK} \sum_{s_i \in \mathrm{sample}} \mathrm{E}_{\hat{\theta}}[ \{
      U(\boldsymbol{s}_i) -
      \widehat{U}(\boldsymbol{s}_i)\}^2].
    }{ 1(k K) sum_{s_i in sample} E[ \{ U(s_i) - hatU(s_i)\}^2]. }

    \item By drawing samples repeatedly and passing the related Kriging
    results as \code{object} to \code{lgnpp}, one can reduce the error of
    the approximation of the mean squared error.

  }
  }
}

\value{
  If \code{is.block} is \code{FALSE} and \code{all.pred} is equal to
  \code{NULL} an updated object of the same class as \code{object} (see
  section \emph{Value} of \code{\link{predict.georob}}).  The data frame
  with the point or block Kriging predictions is complemented by
  \code{lgnpp} with the following new components:

  \itemize{

    \item \code{lgn.pred}: the back-transformed Kriging predictions of a
    log-transformed response.

    \item \code{lgn.se}: the standard  errors of the
    back-transformed predictions.

    \item \code{lgn.lower}, \code{lgn.upper}: the bounds of the
    back-transformed prediction intervals.

  }

  If \code{is.block} is \code{TRUE} or \code{all.pred} not equal to
  \code{NULL} a named numeric vector with two elements:

  \itemize{

    \item \code{mean}: the back-transformed block Kriging estimate, see
    \emph{Details}.

    \item \code{se}: the (approximated) block Kriging standard error, see
    \emph{Details}.
  }

  If \code{extended.output} is \code{TRUE} then the vector is supplemented
  with the attribute \code{mse.lgn.pred} that contains the full covariance
  matrix of the back-transformed point prediction errors.}

\references{
  Cressie, N. (2006) Block Kriging for Lognormal Spatial Processes.
  \emph{Mathematical Geology}, \bold{38}, 413--443,
  \doi{10.1007/s11004-005-9022-8}.

  Hofer, C., Borer, F., Bono, R., Kayser, A. and Papritz, A. 2013.
  Predicting topsoil heavy metal content of parcels of land: An empirical
  validation of customary and constrained lognormal block Kriging and
  conditional simulations.  \emph{Geoderma}, \bold{193--194}, 200--212,
  \doi{10.1016/j.geoderma.2012.08.034}.

  Nussbaum, M., Papritz, A., Baltensweiler, A. and Walthert, L. (2014)
  Estimating soil organic carbon stocks of {Swiss} forest soils by robust
  external-drift kriging.  \emph{Geoscientific Model Development},
  \bold{7}, 1197--1210.  \doi{10.5194/gmd-7-1197-2014}.

}

\author{
   Andreas Papritz \email{papritz@retired.ethz.ch}.
}

\seealso{
  \code{\link{georobPackage}} for a description of the model and a brief summary of the algorithms;

  \code{\link{georob}} for (robust) fitting of spatial linear models;

  \code{\link{predict.georob}} for computing robust Kriging predictions.
}

\examples{
\donttest{
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
r.pred.points <- predict(r.logzn.rob, newdata = meuse.grid.pixdf,
    control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE))
str(r.pred.points$pred@data)

## back-transformation of point predictions
r.backtf.pred.points <- lgnpp(r.pred.points)
str(r.backtf.pred.points$pred@data)

spplot(r.backtf.pred.points[["pred"]], zcol = "lgn.pred", main = "Zn content")

## predicting mean Zn content for whole area
r.block <- lgnpp(r.pred.points, is.block = TRUE, all.pred = r.backtf.pred.points[["pred"]])
r.block

## block predictions of log(Zn)
r.pred.block <- predict(r.logzn.rob, newdata = meuse.blocks,
    control = control.predict.georob(extended.output = TRUE,
        pwidth = 75, pheight = 75))
r.backtf.pred.block <- lgnpp(r.pred.block, newdata = meuse.grid)

spplot(r.backtf.pred.block, zcol = "lgn.pred", main = "block means Zn content")}

}

\keyword{models}
\keyword{spatial}
\keyword{robust}
