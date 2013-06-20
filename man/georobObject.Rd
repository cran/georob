% 2012-12-13 A. Papritz
% R CMD Rdconv -t html -o bla.html georobObject.Rd ; open bla.html; R CMD Rd2pdf --force georobObject.Rd; 

\encoding{macintosh}
\name{georobObject}
\alias{georobObject}
\title{Fitted georob Object}
\description{
  An object of class \code{georob} as returned by \code{\link{georob}} and
  representing a (robustly) fitted spatial linear model.  Objects of this
  class have methods for cross-validation (see
  \code{\link{cv.georob}}), for computing (robust) kriging predictions
  (see \code{\link{predict.georob}}), for plotting (see
  \code{\link{plot.georob}}) and for common generic functions (see
  \code{\link{georobMethods}}).  }

\value{
  The following components are included in a \code{georob} object:
   
  \item{loglik}{the maximized restricted Gaussian loglikelihood of a
  non-robust REML fit or \code{NA} for a robust fit if 
  \code{tuning.psi} is less than \code{1000}.}
   
  \item{variogram.model}{the name of the fitted parametric variogram
  model.}
   
  \item{param}{a named numeric vector with the (estimated) variogram
  parameters.}
   
  \item{aniso}{a list with the following components: 
   
    \itemize{

      \item{\code{isotropic}: logical indicating whether an isotropic
       variogram was fitted.}
       
      \item{\code{aniso}: a named numeric vector with the (estimated)
      anisotropy parameters.}
       
      \item{\code{sincos}: a list with \code{sin} and \code{cos} of the
      angles \eqn{\omega}, \eqn{\phi} and \eqn{\zeta} that define the
      orientation of the anisotropy ellipsoid.}

      \item{\code{rotmat}: the matrix
      \eqn{(\mbox{\boldmath$C$\unboldmath}_1,
      \mbox{\boldmath$C$\unboldmath}_2,
      \mbox{\boldmath$C$\unboldmath}_3)}{(C_1, C_2, C_3)} (see
      \code{\link{georobIntro}}).}
      
      \item{\code{sclmat}: a vector with the elements 1, \eqn{1/f_1},
      \eqn{1/f_2} (see \code{\link{georobIntro}}). }

    }
  }
  
  \item{gradient}{a named numeric vector with the estimating equations
  (robust REML) evaluated at the solution or the gradient of the maximized
  restricted loglikelihood (Gaussian REML).}

  \item{psi.func}{the name of the parametric \code{psi.function} used for the
  robust model fit (see \code{\link{georob.control}}).}

  \item{tuning.psi}{the value of the tuning constant \eqn{c} of the
  \eqn{\psi_c}-function.}

  \item{coefficients}{a named vector with the estimated regression coefficients.}

  \item{fitted.values}{a named vector with the fitted values of the external drift 
  \eqn{\mbox{\boldmath$X$\unboldmath}\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{X hat\beta}.}

  \item{bhat}{a named vector with the predicted spatial random effects
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} at the data locations.}

  \item{residuals}{a named vector with the residuals 
  \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}=\mbox{\boldmath$Y$\unboldmath} - 
  \mbox{\boldmath$X$\unboldmath} \widehat{\mbox{\boldmath$\beta$\unboldmath}} -
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{hat\epsilon=Y-X hat\beta - hatB}.}
   
  \item{rweights}{a named numeric vector with the \dQuote{robustness weights}  
  \eqn{\psi(\widehat{\varepsilon}_i/\widehat{\tau}) / 
    (\widehat{\varepsilon}_i/\widehat{\tau})}{
    \psi(hat\epsilon_i/hat\tau) / (hat\epsilon_i/hat\tau)}.}
    
  \item{converged}{logical indicating whether numerical maximization of
    the restricted loglikelihood by \code{optim} or root finding by
    \code{nleqslv} converged.}

  \item{convergence.code}{a diagnostic integer issued by \code{\link{optim}}
    (component \code{convergence}) or \code{\link{nleqslv}} (component
    \code{termcd}) about convergence.}

  \item{iter}{a named integer vector of length two, indicating either

    \itemize{

      \item the number of function and gradient evaluations when maximizing
      the restricted Gaussian loglikelihood by \code{\link{optim}}, or
      
      \item the number of function and Jacobian evaluations when solving
      the robustified estimating equations by \code{\link{nleqslv}}.
      
    }
  }
 
  \item{Tmat}{the compressed design matrix for replicated observations at
  coincident locations (integer vector that contains for each observation
  the row index of the respective unique location).}
  
  \item{cov}{a list with covariance matrices (or diagonal variance
  vectors).  Covariance matrices are stored in \emph{compressed form} (see
  \code{\link{compress}}) and can be expanded to square matrices by
  \code{\link{expand}}.  What \code{cov} actually contains depends on the
  flags passed to \code{georob} for computing covariances (see
  \code{\link{georob.control}}).  Possible components are:
  
  \itemize{
  
      \item \code{cov.bhat}: the covariances of 
      \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB}.
      
      \item \code{cov.betahat}: the covariances of 
      \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta}.
      
      \item \code{cov.bhat.betahat}: the covariances of 
      \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} and
      \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta}.
      
      \item \code{cov.delta.bhat}: the covariances of 
      \eqn{\mbox{\boldmath$B$\unboldmath}- \widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB}.
      
      \item \code{cov.delta.bhat.betahat}: the covariances of 
      \eqn{\mbox{\boldmath$B$\unboldmath}- \widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB}
      and 
      \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta}.
      
      \item \code{cov.ehat}: the covariances of 
      \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}=\mbox{\boldmath$Y$\unboldmath} - 
      \mbox{\boldmath$X$\unboldmath} \widehat{\mbox{\boldmath$\beta$\unboldmath}} -
      \widehat{\mbox{\boldmath$B$\unboldmath}}}{hat\epsilon=Y-X hat\beta - hatB}.
      
      \item \code{cov.ehat.p.bhat}: the covariances of 
      \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}} + 
      \widehat{\mbox{\boldmath$B$\unboldmath}}
      =\mbox{\boldmath$Y$\unboldmath} - 
      \mbox{\boldmath$X$\unboldmath} \widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\epsilon+ hatB=Y-X hat\beta}.
      
      \item \code{aux.cov.pred.target}: a covariance term required for the
      back-trans-\cr formation of kriging predictions of log-transformed data.
      
    }
  
  }
     
  \item{param.tf}{a character vector indicating the transformations of the
  variogram parameters for model fitting.}
    
  \item{fwd.tf}{a list of functions for variogram parameter transformations.}
  
  \item{bwd.tf}{a list of inverse functions for \emph{inverse} variogram
  parameter transformations.}
    
  \item{hessian}{a symmetric matrix giving an estimate of the Hessian at
  the solution if the model was fitted non-robustly with the argument
  \code{hessian = TRUE} (see \code{\link{georob.control}}).  Missing
  otherwise.}
    
  \item{expectations}{a named numeric vector with the expectations of
  \eqn{\partial \psi_c(x)/\partial x}{d\psi_c'(x)/dx} (\code{dpsi}) and
  \eqn{\psi_c^2(x)} (\code{psi2}) with respect to a standard normal
  distribution.}
     
  \item{Valpha.objects}{a list of matrices in \emph{compressed form} with
  the components:
 
    \itemize{
    
      \item \code{Valpha.inverse}: the inverse of 
        \eqn{\mbox{\boldmath$V$\unboldmath}_\alpha}{V_\alpha}.
        
      \item \code{Valpha.ilcf}: the inverse of the lower triangular factor
        \eqn{\mbox{\boldmath$L$\unboldmath}}{L} of the Cholesky
        decomposition \eqn{\mbox{\boldmath$V$\unboldmath}_\alpha =
        \mbox{\boldmath$L$\unboldmath}
        \mbox{\boldmath$L$\unboldmath}\mathrm{^T}}{V_\alpha=L L^T}.
    
      \item \code{Valpha.ucf}: the upper triangular Cholesky factor
        \eqn{\mbox{\boldmath$L$\unboldmath}\mathrm{^T}}{L^T}.
        
      \item \code{gcr.constant}: the constant \eqn{\gamma_0} (see expression for
      \eqn{\mbox{\boldmath$V$\unboldmath}_\alpha}{Valpha} in section
      \bold{Model} of \code{\link{georobIntro}}.
    
    }
 
  }
  
  \item{locations.object}{a list  with 2 components:
 
    \itemize{
    
      \item \code{locations}: a formula indicating the coordinates of the
      measurement locations.
        
      \item \code{locations.coords}: a numeric matrix with the coordinates
      of the measurement locations.
    
    }
 
  }
  
  \item{initial.objects}{a list with 5 components:
 
    \itemize{
    
      \item \code{coefficients}: initial estimates of
      \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} computed either by
      \code{\link[robustbase]{lmrob}} or \code{\link[quantreg]{rq}}.
      
      \item \code{bhat}: initial predictions of
      \eqn{\mbox{\boldmath$B$\unboldmath}}{B}.
      
      \item \code{param}: numeric vector with initial estimates of the
      variogram parameters, either computed (\code{initial.param = TRUE})
      or as passed to \code{georob} (\code{initial.param = FALSE}).
      
      \item \code{fit.param}: logical vector indicating which variogram
      parameters were fitted.
      
      \item \code{aniso}: numeric vector with initial estimates of the
      anisotropy parameters, either either computed (\code{initial.param = TRUE})
      or as passed to \code{georob} (\code{initial.param = FALSE}).
      
      \item \code{fit.aniso}: logical vector indicating which anisotropy
      parameters were fitted.
      
      \item \code{initial.param}: logical flag indicating whether robust
      initial variogram parameter estimates were computed (see
      \code{\link{georob}}, \emph{Details}).
            
    }
 
  }
  
  \item{MD}{optionally a matrix of robust distances in the space spanned by
  \eqn{\mbox{\boldmath$X$\unboldmath}}{X} (see argument \code{compute.rd}
  of \code{\link[robustbase]{lmrob.control}} and
  \code{\link{georob.control}}).}
  
  \item{model, x, y}{if requested the model frame, the model matrix and the
  response, respectively.}  
  
  \item{\code{na.action}, \code{offset}, \code{contrasts}, \code{xlevels},
  \code{rank}, \code{df.residual}, \code{call}, \code{terms}}{further
  components of the fit as described for an object of class
  \code{\link[stats]{lm}}.}
  
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{georob.control}} for controlling the behaviour of \code{georob};
  \code{\link{plot.georob}} for display of REML variogram estimates;
  \code{\link{cv.georob}} for assessing the goodness of a fit by \code{georob}; 
  \code{\link{predict.georob}} for computing robust kriging predictions; and finally
  \code{\link{georobMethods}} for further methods for the class \code{georob}.
}


\keyword{models}
\keyword{spatial}
\keyword{robust}
