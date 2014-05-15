% 2012-12-13 A. Papritz
% R CMD Rdconv -t html -o bla.html lgnpp.Rd ; open bla.html; R CMD Rd2pdf --force lgnpp.Rd; 

\encoding{macintosh}
\name{lgnpp}
\alias{lgnpp}

\title{Unbiased Back-Transformations for Lognormal Kriging}

\description{
  The function \code{lgnpp} back-transforms point or block kriging
  predictions of a log-transformed response variable computed by
  \code{\link{predict.georob}}.  Alternatively, the function averages
  lognormal point kriging predictions for a block and approximates the mean
  squared prediction error of the block mean.}



\usage{
lgnpp(object, newdata, locations, is.block = FALSE, all.pred = NULL, 
    extended.output = FALSE)
}

\arguments{
  \item{object}{an object with kriging predictions of a log-transformed
  response variable as obtained by
  \code{\link{predict}(\var{georob-object, ...})}.}
  
  \item{newdata}{an optional object as passed as argument \code{newdata} to
  \code{\link{predict.georob}}, see \emph{Details}.}
  
  \item{locations}{an optional one-sided formula specifying what variables
  of \code{newdata} are the coordinates of the prediction points, see
  \code{\link{predict.georob}}.}
  
  \item{is.block}{an optional logical (default \code{FALSE}) specifying
  whether point predictions contained in \code{object} are considered to
  belong to a single block and should be averaged after
  back-transformation.  Ignored if \code{object} contains block kriging
  predictions, see \emph{Details}.}
  
  \item{all.pred}{an optional positive integer or an object as obtained by
  \code{lgnpp(predict(\var{georob-object, ...}))}, see \emph{Details}.}
  
  \item{extended.output}{logical controlling whether the covariance matrix
  of the errors of the back-transformed point predictions is added as an
  attribute to the result, see \emph{Details}.}
  
}
\details{
 The function \code{lgnpp} performs three tasks:
 
  \subsection{1. Back-transformation of point kriging predictions of a
  log-transformed response}{
 
    The usual, marginally unbiased back-transformation for lognormal point
    kriging is used:
 
    \deqn{\widehat{Z}(\mbox{\boldmath$s$\unboldmath}) = \exp( \widehat{S}(\mbox{\boldmath$s$\unboldmath}) + 
      1/2 (  \mathrm{Var}_{\hat{\theta}}[ S(\mbox{\boldmath$s$\unboldmath})]
      - \mathrm{Var}_{\hat{\theta}}[\widehat{S}(\mbox{\boldmath$s$\unboldmath})])),}{
      hatZ(s)=exp(hatS(s) + 1/2(hat\sigma^2_n + hat\sigma^2 - Var[hatS(s)])),}
      
    \deqn{\mathrm{Cov}_{\hat{\theta}}[
      Z(\mbox{\boldmath$s$\unboldmath}_i) - \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_i),
      Z(\mbox{\boldmath$s$\unboldmath}_j) - \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_j) 
      ] = \mu_{\hat{\theta}}(\mbox{\boldmath$s$\unboldmath}_i) \mu_{\hat{\theta}}(\mbox{\boldmath$s$\unboldmath}_j)
    }{
      Cov[ Z(s_i)-hatZ(s_i), Z(s_j)-hatZ(s_j)] = \mu(s_i) \mu(s_j)
    }
    \deqn{
      \times \{
        \exp(\mathrm{Cov}_{\hat{\theta}}[S(\mbox{\boldmath$s$\unboldmath}_i),S(\mbox{\boldmath$s$\unboldmath}_j)])
        -2\exp(\mathrm{Cov}_{\hat{\theta}}[\widehat{S}(\mbox{\boldmath$s$\unboldmath}_i),S(\mbox{\boldmath$s$\unboldmath}_j)])
        +\exp(\mathrm{Cov}_{\hat{\theta}}[\widehat{S}(\mbox{\boldmath$s$\unboldmath}_i),\widehat{S}(\mbox{\boldmath$s$\unboldmath}_j)])
      \},
    }{
      *\{
      exp(Cov[Y(s_i), Y(s_j)]) - 2 exp(Cov[hatY(s_i), Y(s_j)]) + exp(Cov[hatY(s_i), hatY(s_j)])
      \},  
    }
    
    where \eqn{\widehat{S}}{hatY} and \eqn{\widehat{Z}}{hatZ} denote the
    log- and back-transformed predictions of the signal,
    and
    
    \deqn{\mu_{\hat{\theta}}(\mbox{\boldmath$s$\unboldmath}) \approx 
      \exp(\mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath})\mathrm{^T}\widehat{\beta} 
      + 1/2 \mathrm{Var}_{\hat{\theta}}[S(\mbox{\boldmath$s$\unboldmath})]).    
    }{\mu(s) = exp( x(s)^T hat\beta + 1/2 Var[S(s)]).}
    
    The expressions for the required covariance terms can be found in the
    Appendices of Nussbaum et al.  (2012).  Instead of the signal
    \eqn{S(\mbox{\boldmath$s$\unboldmath})}{S(s)}, predictions of the
    log-transformed response \eqn{Y(\mbox{\boldmath$s$\unboldmath})}{S(s)}
    or the estimated trend
    \eqn{\mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath})^\mathrm{T}\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{x(s)^T
    hat\beta} of the log-transformed data can be back-transformed (see
    \code{\link{georobIntro}}).  The
    above transformations are used if \code{object}
    contains point kriging predictions (see \code{predict.georob},
    \emph{Value}) and if \code{is.block = FALSE} and \code{all.pred} is
    missing.
    
  }
 
  \subsection{2. Back-transformation of block kriging predictions of a
  log-transformed response}{
  
  Block kriging predictions of a log-transformed response variable are
  back-transformed by the approximately unbiased transformation proposed
  by Cressie (2006)
  
  \deqn{\widehat{Z}(B) = \exp( \widehat{S}(B) + 1/2 \{
    \mathrm{Var}_{\hat{\theta}}[S(\mbox{\boldmath$s$\unboldmath})] + \widehat{\mbox{\boldmath$\beta$\unboldmath}}\mathrm{^T}
    \mbox{\boldmath$M$\unboldmath}(B) \widehat{\mbox{\boldmath$\beta$\unboldmath}} -
    \mathrm{Var}_{\hat{\theta}}[\widehat{S}(B)]  
    \}),
  }{
    hatZ(B) = exp( hatS(B) + 1/2 \{Var[S(s)] + hat\beta^T M(B) hat\beta - Var[hatS(B)]\}),
  }
  
  \deqn{\mathrm{E}_{\hat{\theta}}[\{Z(B) - \widehat{Z}(B))^2] = \mu_{\hat{\theta}}(B)^2 \{
    \exp(\mathrm{Var}_{\hat{\theta}}[S(B)]) - 2 \exp(\mathrm{Cov}_{\hat{\theta}}[\widehat{S}(B),S(B)]) + \exp(\mathrm{Var}_{\hat{\theta}}[\widehat{S}(B)])
    \}
  }{
    E[ \{ Z(B)- hatZ(B) \}^2 ] = \mu(B)^2 \{
    \exp(Var[S(B)]) - 2 \exp(Cov[S(B),hatS(B)]) + \exp(Var[hatS(B)]) \}
  }
  
  where \eqn{\widehat{S}(B)} and \eqn{\widehat{Z}(B)} are the log- and
  back-transformed predictions of the block mean \eqn{Z(B)}, respectively,
  \eqn{\mbox{\boldmath$M$\unboldmath}(B)}{M(B)} is the spatial
  covariance matrix of the covariates
  
  \deqn{ \mbox{\boldmath$M$\unboldmath}(B) = 1/|B| \int_B 
  ( \mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}) - \mbox{\boldmath$x$\unboldmath}(B) )
  ( \mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}) - \mbox{\boldmath$x$\unboldmath}(B) )\mathrm{^T} \,d\mbox{\boldmath$s$\unboldmath}
  }{M(B) = 1/|B| int_B (x(s)-x(B) (x(s)-x(B))^T) ds}
  
  with 
  \deqn{ \mbox{\boldmath$x$\unboldmath}(B) = 1/|B| \int_B \mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}) \,d\mbox{\boldmath$s$\unboldmath}  
  }{x(B) = 1/|B| int_B x(s) ds}
  
  and 
  
  \deqn{ \mu_{\hat{\theta}}(B) \approx \exp(\mbox{\boldmath$x$\unboldmath}(B)\mathrm{^T}
    \widehat{\mbox{\boldmath$\beta$\unboldmath}} + 1/2 \mathrm{Var}_{\hat{\theta}}[S(B)]). 
  }{
    \mu(B) = exp( x(B)^T hat\beta + 1/2 Var[S(B)]).
  }
  
  These equations are based on the assumption that both the point data
  \eqn{Z(\mbox{\boldmath$s$\unboldmath})}{Z(s)} and the block means
  \eqn{Z(B)} follow  lognormal laws, which strictly cannot hold.  But
  for small blocks the assumption works well as the bias and the loss of
  efficiency caused by this assumption are small (Cressie, 2006;
  Hofer et al., 2013). 
  
  The above formulae are used by \code{lgnpp} if \code{object} contains
  block kriging predictions in the form of a
  \code{\link{SpatialPolygonsDataFrame}}.  To approximate
  \eqn{\mbox{\boldmath$M$\unboldmath}(B)}{M(B)}, one needs the covariates
  on a fine grid within the block \eqn{B}.  The covariates are passed to
  \code{lgnpp} as argument \code{newdata}, where \code{newdata} can be any
  spatial data frame accepted by \code{predict.georob}.  For evaluating
  \eqn{\mbox{\boldmath$M$\unboldmath}(B)}{M(B)} the geometry of the blocks
  is taken from the \code{polygons} slot of the
  \code{SpatialPolygonsDataFrame} passed as \code{object} to \code{lgnpp}.
  
  }
 
  \subsection{3. Backtransformation and averaging of point kriging predictions
  of a log-transformed response}{
  
  \code{lgnpp} allows as a further option to back-transform and
  \emph{average} point kriging predictions passed as \code{object} to the
  function.  One then assumes that the predictions refer to points
  that lie in \emph{a single} block.  Hence, one uses the approximation
  
  \deqn{\widehat{Z}(B) \approx \frac{1}{K} \sum_{s_i \in B} \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_i)
  }{hatZ(B) = 1/K sum_{s_i in B} hatZ(s_i)}
  
  to predict the block mean \eqn{Z(B)}, where \eqn{K} is the number of
  points in \eqn{B}. The mean squared error can be approximated by
  
  \deqn{\mathrm{E}_{\hat{\theta}}[\{Z(B) - \widehat{Z}(B)\}^2] \approx \frac{1}{K^2}
    \sum_{s_i \in B} \sum_{s_j \in B} 
    \mathrm{Cov}_{\hat{\theta}}[
    Z(\mbox{\boldmath$s$\unboldmath}_i) - \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_i),
    Z(\mbox{\boldmath$s$\unboldmath}_j) - \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_j) 
    ].
    }{E[\{Z(B) - hatZ(B)\}^2] = 1/K^2 sum_{s_i in B} sum_{s_j in B} 
    Cov[ Z(s_i)-hatZ(s_i), Z(s_j)-hatZ(s_j)].  
  }
    
 In most instances, the evaluation of the above double sum is not feasible
 because a large number of points is used to discretize the block \eqn{B}.
 lgnpp then uses the following approximation for the mean squared error
 (see also Appendix E of Nussbaum et al., 2012):
  
  \itemize{
  
    \item Prediction results are passed as \code{object} to \code{lgnpp}
    only for a \emph{random sample of points in \eqn{B}} (of size \eqn{k}),
    for which the evaluation of the above double sum is feasible.
    
    \item The prediction results for the \emph{complete set of points}
    within the block are passed as argument \code{all.pred} to
    \code{lgnpp}.  These results are used to compute \eqn{\widehat{Z}(B)}{Z(B)}.
    
    \item The mean squared error is then approximated by  
    
    \deqn{
      \mathrm{E}_{\hat{\theta}}[\{Z(B) - \widehat{Z}(B)\}^2] \approx 
      \frac{1}{K^2} \sum_{s_i \in B} \mathrm{E}_{\hat{\theta}}[ \{ Z(\mbox{\boldmath$s$\unboldmath}_i) - \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_i)\}^2]
    }{E[\{Z(B) - hatZ(B)\}^2] = 1/K^2 sum_{s_i in B} E[ \{ Z(s_i) - hatZ(s_i)\}^2] }
    
    \deqn{+ \frac{K-1}{K k (k-1)} \sum_{s_i \in \mathrm{sample}}\sum_{s_j \in \mathrm{sample}, s_j \neq s_i}
%     \mathrm{Cov}_{\hat{\theta}}[
    Z(\mbox{\boldmath$s$\unboldmath}_i) - \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_i),
    Z(\mbox{\boldmath$s$\unboldmath}_j) - \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_j) 
    ].
    }{ + (K-1)/(K k (k-1)) sum_{s_i in sample} sum_{s_j in sample, s_j != s_i} 
      Cov[ Z(s_i)-hatZ(s_i), Z(s_j)-hatZ(s_j)]
    }
    
    The first term of the RHS can be computed from the point kriging
    results contained in \code{all.pred}, and the double sum is evaluated
    from the full covariance matrices of the predictions and the respective
    targets, passed to \code{lgnpp} as \code{object}.
    
    \item If the prediction results are not available for the complete set
    of points in \eqn{B} then \code{all.pred} may be equal to \eqn{K}.  The
    block mean is then approximated by
    
    \deqn{\widehat{Z}(B) \approx \frac{1}{k} \sum_{s_i \in \mathrm{sample}}
    \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_i) }{hatZ(B) = sum_{s_i in
    B} hatZ(s_i)}
    
    and the first term of the RHS of the expression for the mean squared
    error by
    
      \deqn{ \frac{1}{kK} \sum_{s_i \in \mathrm{sample}} \mathrm{E}_{\hat{\theta}}[ \{
      Z(\mbox{\boldmath$s$\unboldmath}_i) -
      \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_i)\}^2].
    }{ 1(k K) sum_{s_i in sample} E[ \{ Z(s_i) - hatZ(s_i)\}^2]. }
    
    \item By drawing samples repeatedly and passing the related kriging
    results as \code{object} to \code{lgnpp}, one can reduce the error of
    the approximation of the mean squared error.
    
  }
  }
}
 
\value{
  If \code{is.block} is\code{FALSE} and \code{all.pred} is equal to
  \code{NULL} an updated object of the same class as \code{object} (see
  section \emph{Value} of \code{\link{predict.georob}}).  The data frame
  with the point or block kriging predictions is complemented by
  \code{lgnpp} with the following new components:
  
  \itemize{
  
    \item \code{lgn.pred}: the back-transformed kriging predictions of a
    log-transformed response.
    
    \item \code{lgn.se}: the standard  errors of the
    back-transformed predictions.
    
    \item \code{lgn.lower}, \code{lgn.upper}: the bounds of the
    back-transformed prediction intervals.
  
  }
  
  If \code{is.block} is \code{TRUE} or \code{all.pred} not equal to
  \code{NULL} a named numeric vector with two elements:
  
  \itemize{ 
    
    \item \code{mean}: the back-transformed block kriging estimate, see
    \emph{Details}.
      
    \item \code{mse}: the (approximated) block kriging variance, see
    \emph{Details}.
  }
  
  If \code{extended.output} is \code{TRUE} then the vector is supplemented
  with the attribute \code{mse.lgn.pred} that contains the full covariance
  matrix of the back-transformed point prediction errors.}

\references{
  Cressie, N. (2006) Block Kriging for Lognormal Spatial Processes.
  \emph{Mathematical Geology}, \bold{38}, 413--443.
  
  Hofer, C., Borer, F., Bono, R., Kayser, A. and Papritz, A. 2013. Predicting
  topsoil heavy metal content of parcels of land: An empirical validation
  of customary and constrained lognormal block kriging and conditional
  simulations. \emph{Geoderma}, \bold{193--194}, 200--212.

  Nussbaum, M., Papritz, A., Baltensweiler, A. and Walthert, L. (2012)
  \emph{Organic Carbon  Stocks of Swiss Forest Soils},
  Institute of Terrestrial Ecosystems, ETH Zurich and 
  Swiss Federal Institute for Forest, Snow and Landscape Research
  (WSL), pp. 51. 
  \url{http://e-collection.library.ethz.ch/eserv/eth:6027/eth-6027-01.pdf}
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}.
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{predict.georob}} for computing robust kriging predictions.
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
