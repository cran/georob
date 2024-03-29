% 2024-03-09 A. Papritz
% R CMD Rdconv -t html -o bla.html georobObject.Rd ; open bla.html; R CMD Rd2pdf --force georobObject.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  georobObject.Rd
\encoding{utf8}
\name{georobObject}
\alias{georobObject}

\title{Fitted georob Object}

\description{
  An object of class \code{georob} as returned by \code{\link{georob}} and
  representing a (robustly) fitted spatial linear model.  Objects of this
  class have methods for model building (see
  \code{\link{georobModelBuilding}}) and cross-validation (see
  \code{\link{cv.georob}}), for computing (robust) Kriging predictions (see
  \code{\link{predict.georob}}), for plotting (see
  \code{\link{plot.georob}}) and for common generic functions (see
  \code{\link{georobMethods}}).  }

\value{
  A \code{georob} object is a list with following components:

  \item{loglik}{the maximized (restricted) Gaussian log-likelihood of a
  non-robust (RE)ML fit or \code{NA} for a robust fit if
  \code{tuning.psi} is less than \code{tuning.psi.nr}.}

  \item{variogram.object}{the estimated parameters of a possibly nested
  variograms model.  This is a list that contains for each variogram model
  structure the following components:

  \itemize{

      \item{\code{variogram.model}}: the name of the fitted parametric variogram
      model.

      \item{\code{param}}: a named numeric vector with the (estimated) variogram
      parameters.

      \item{\code{fit.param}}: a named logical vector with the flags defining what
      variogram parameters were estimated.

      \item{\code{isotropic}: logical indicating whether an isotropic
      variogram was fitted.}

      \item{\code{aniso}: a named numeric vector with the (estimated)
      anisotropy parameters.}

      \item{\code{fit.aniso}}: a named logical vector with the flags defining what
      anisotropy parameters were estimated.

      \item{\code{sincos}: a list with \code{sin} and \code{cos} of the
      angles \eqn{\omega}, \eqn{\phi} and \eqn{\zeta} that define the
      orientation of the anisotropy ellipsoid (see \code{\link{georobPackage}}).}

      \item{\code{rotmat}: the matrix
      \eqn{(\boldsymbol{C}_1,
      \boldsymbol{C}_2,
      \boldsymbol{C}_3)}{(C_1, C_2, C_3)} (see
      \code{\link{georobPackage}}).}

      \item{\code{sclmat}: a vector with the elements 1, \eqn{1/f_1},
      \eqn{1/f_2} (see \code{\link{georobPackage}}). }

    }

  }

  \item{gradient}{a named numeric vector with the estimating equations
  (robust REML) or the gradient of the maximized (restricted) log-likelihood
  (Gaussian (RE)ML) evaluated at the solution.}

  \item{tuning.psi}{the value of the tuning constant \eqn{c} of the
  \eqn{\psi_c}-function.}

  \item{coefficients}{a named vector with the estimated regression coefficients
  \eqn{\widehat{\boldsymbol{\beta}}}{hat\beta}.}

  \item{fitted.values}{a named vector with the fitted values of the
  external drift
  \eqn{\boldsymbol{X}\widehat{\boldsymbol{\beta}}}{X
  hat\beta}.}

  \item{bhat}{a named vector with the predicted spatial random effects
  \eqn{\widehat{\boldsymbol{B}}}{hatB} at the data locations.}

  \item{residuals}{a named vector with the residuals
  \eqn{\widehat{\boldsymbol{\varepsilon}}=\boldsymbol{Y} -
  \boldsymbol{X} \widehat{\boldsymbol{\beta}} -
  \widehat{\boldsymbol{B}}}{hat\epsilon=Y-X hat\beta - hatB}.}

  \item{rweights}{a named numeric vector with the \dQuote{robustness weights}
  \eqn{\psi_c(\widehat{\varepsilon}_i/\widehat{\tau}) /
    (\widehat{\varepsilon}_i/\widehat{\tau})}{
    \psi_c(hat\epsilon_i/hat\tau) / (hat\epsilon_i/hat\tau)}.}

  \item{converged}{a logical scalar indicating whether numerical maximization of
  the (restricted)\cr log-likelihood by \code{nlminb} or \code{optim} or root
  finding by \code{nleqslv} converged.}

  \item{convergence.code}{a diagnostic integer issued by
  \code{\link[stats]{nlminb}}, \code{\link[stats]{optim}} (component
  \code{convergence}) or \code{\link[nleqslv]{nleqslv}} (component
  \code{termcd}) about convergence.}

  \item{iter}{a named integer vector of length two, indicating either

    \itemize{

      \item the number of function and gradient evaluations when maximizing
      the (restricted) Gaussian log-likelihood by \code{\link[stats]{nlminb}}
      or \code{\link[stats]{optim}}, or

      \item the number of function and Jacobian evaluations when solving
      the robustified estimating equations by
      \code{\link[nleqslv]{nleqslv}}.

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
  \code{\link{control.georob}}).  Possible components are:

  \itemize{

      \item \code{cov.bhat}: the covariances of
      \eqn{\widehat{\boldsymbol{B}}}{hatB}.

      \item \code{cov.betahat}: the covariances of
      \eqn{\widehat{\boldsymbol{\beta}}}{hat\beta}.

      \item \code{cov.delta.bhat}: the covariances of
      \eqn{\boldsymbol{B}- \widehat{\boldsymbol{B}}}{B-hatB}.

      \item \code{cov.delta.bhat.betahat}: the covariances of
      \eqn{\boldsymbol{B}- \widehat{\boldsymbol{B}}}{B-hatB}
      and
      \eqn{\widehat{\boldsymbol{\beta}}}{hat\beta}.

      \item \code{cov.ehat}: the covariances of
      \eqn{\widehat{\boldsymbol{\varepsilon}}=\boldsymbol{Y} -
      \boldsymbol{X} \widehat{\boldsymbol{\beta}} -
      \widehat{\boldsymbol{B}}}{hat\epsilon=Y-X hat\beta - hatB}.

      \item \code{cov.ehat.p.bhat}: the covariances of
      \eqn{\widehat{\boldsymbol{\varepsilon}} +
      \widehat{\boldsymbol{B}}
      =\boldsymbol{Y} -
      \boldsymbol{X} \widehat{\boldsymbol{\beta}}}{hat\epsilon+ hatB=Y-X hat\beta}.

      \item \code{cov.pred.target}: a covariance term required for the
      back-trans-\cr formation of Kriging predictions of log-transformed data.

    }

  }

  \item{expectations}{a named numeric vector with the expectations of
  \eqn{\partial \psi_c(x)/\partial x}{d\psi_c'(x)/dx} and \eqn{\psi_c^2(x)}
  with respect to a standard normal distribution (\code{exp.gauss.dpsi},
  \code{var.gauss.psi}) and the expectations of
  \eqn{\varepsilon^2}{\epsilon^2} and \eqn{\psi_c^2(x)} with respect to the
  long-tailed distribution of \eqn{\varepsilon}{\epsilon}
  (\code{var.f0.eps}, \code{var.f0.psi}).}

  \item{Valphaxi.objects}{a list of matrices in \emph{compressed form} with
  (among others) the following components:

    \itemize{

      \item \code{Valpha}: a list with the (generalized) correlation
      matrices (\code{Valpha}) of the nested variogram models structures
      along with the constants (\code{gcr.constant}) added to the respective
      semi-variance matrices.

      \item \code{Valphaxi}: the (generalized) correlation matrix
      \eqn{\boldsymbol{V}_{\alpha,\xi} =
      \boldsymbol{\Gamma}_{\alpha,\xi} /
      (\sigma_{\mathrm{n}}^2+\sigma^2 )}{V_{\alpha, \xi} = \Gamma_\theta /
      (\sigma_n^2+\sigma^2)} that includes the spatial nugget effect.

      \item \code{Valphaxi.inverse}: the inverse of
      \eqn{\boldsymbol{V}_{\alpha,\xi}}{V_{\alpha, \xi}}.

      \item \code{log.det.Valphaxi}:
      \eqn{\log(\det(\boldsymbol{V}_{\alpha,\xi}))}{log(det(V_{\alpha, \xi}))}.

    }
  }

  \item{zhat.objects}{a list of matrices in (partly) \emph{compressed form} with
  the following components:

    \itemize{

      \item \code{Aalphaxi}: the matrix
      \eqn{(\boldsymbol{X}^T
      \boldsymbol{V}_{\alpha,\xi}^{-1}\boldsymbol{X})^{-1}
      \boldsymbol{X}^T\boldsymbol{V}_{\alpha,\xi}^{-1}
      }{(X^T V_{\alpha, \xi}^-1 X)^-1 X^T V_{\alpha, \xi}^-1 }.

      \item \code{Palphaxi}: the matrix
      \eqn{\boldsymbol{I}-
      \boldsymbol{X} \boldsymbol{A}_{\alpha,\xi}
      }{I - X A_{\alpha, \xi}}.

      \item \code{Valphaxi.inverse.Palphaxi}: the matrix
      \eqn{\boldsymbol{V}^{-1}_{\alpha,\xi}
      \boldsymbol{P}_{\alpha,\xi} }{V^-1_{\alpha, \xi}
      P_{\alpha, \xi}}.

    }
  }
  \item{locations.object}{a list with 3 components:

    \itemize{

      \item \code{locations}: a formula indicating the coordinates of the
      measurement locations.

      \item \code{coordinates}: a numeric matrix with the coordinates
      of the measurement locations.

      \item \code{lag.vectors}: a numeric matrix with the lag vectors
      between any distinct pairs of measurement locations.

    }

  }

  \item{initial.objects}{a list with 3 components:

    \itemize{

      \item \code{coefficients}: initial estimates of
      \eqn{\boldsymbol{\beta}}{\beta} computed either by
      \code{\link[robustbase]{lmrob}} or \code{\link[quantreg]{rq}}.

      \item \code{bhat}: initial predictions of
      \eqn{\boldsymbol{B}}{B}.

      \item \code{variogram.object}: the initial values of the parameters
      of a possibly nested variograms model.  This is a list with the same
      structure as described above for the component
      \code{variogram.object}.

    }

  }

  \item{hessian.tfpa}{a symmetric matrix with the Hessian (observed
  Fisher information) at the solution with respect to the transformed
  variogram and anisotropy parameters if the model was fitted non-robustly
  with the argument \code{hessian = TRUE} (see
  \code{\link{control.georob}}).  Missing otherwise.  This Hessian is used
  by \code{\link{summary.georob}} to compute confidence intervals for the
  estimated parameters.}

  \item{hessian.ntfpa}{a symmetric matrix with the Hessian (observed
  Fisher information) at the solution with respect to the non-transformed
  variogram and anisotropy parameters if the model was fitted non-robustly
  with the argument \code{hessian = TRUE} (see
  \code{\link{control.georob}}).  Missing otherwise.}

  \item{control}{a list with control parameters generated by
  \code{\link{control.georob}}.}

  \item{MD}{optionally a matrix of robust distances in the space spanned by
  \eqn{\boldsymbol{X}}{X} (see argument \code{compute.rd}
  of \code{\link[robustbase]{lmrob.control}} and
  \code{\link{control.georob}}).}

  \item{model, x, y}{if requested the model frame, the model matrix and the
  response, respectively.}

  \item{\code{na.action}, \code{offset}, \code{contrasts}, \code{xlevels},
  \code{rank}, \code{df.residual}, \code{call}, \code{terms}}{further
  components of the fit as described for an object of class
  \code{\link[stats]{lm}}.}

}

\author{
   Andreas Papritz \email{papritz@retired.ethz.ch}.
}

\seealso{
  \code{\link{georobPackage}} for a description of the model and a brief summary of the algorithms;

  \code{\link{georob}} for (robust) fitting of spatial linear models;

%  \code{\link{georobObject}} for a description of the class \code{georob};
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
  from model fitted by \code{georob}; and finally

  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
  for robust estimation and modelling of sample variograms.
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
