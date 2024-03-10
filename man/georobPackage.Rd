% 2024-03-09 A. Papritz
% R CMD Rd2pdf --force georobPackage.Rd;  R CMD Rdconv -t html -o bla.html georobPackage.Rd ; open bla.html;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  georobPackage.Rd
\encoding{utf8}
\name{georob-package}
\alias{georob-package}
\alias{georobPackage}
\alias{georobIntro}
\docType{package}
%
\title{The georob Package}
%
\description{This is a summary of the features and functionality of
  \pkg{georob}, a package in \R for customary and robust geostatistical
  analyses.}
%

\details{

% intro

  \pkg{georob} is a package for customary and robust analyses of
  geostatistical data.
  Such data, say \eqn{y_i=y(\boldsymbol{s}_i)}{y(s_i)}, are
  recorded at a set of locations,
  \eqn{\boldsymbol{s}_i}{s_i}, \eqn{i=1,2, \ldots, n}, in a
  domain \eqn{G \in \mathrm{I}\!\mathrm{R}^d}{G in R^d}, \eqn{d \in (1,2,3)}{d in (1,2,3)}, along
  with covariate information
  \eqn{x_j(\boldsymbol{s}_i)}{x_j(s_i)}, \eqn{j=1,2,
  \ldots, p}.

	%   model

  \subsection{Model}{
    We use the following model for the data
    \eqn{y_i=y(\boldsymbol{s}_{i})}{y_i=y(s_i)}:
    \deqn{Y(\boldsymbol{s}_i) =
      Z(\boldsymbol{s}_i) + \varepsilon =
      \boldsymbol{x}(\boldsymbol{s}_i)^\mathrm{T}
         \boldsymbol{\beta} +
      B(\boldsymbol{s}_i) +
      \varepsilon_i,}{%
      Y(s_i) = Z(s_i) + \epsilon = x(s_i)^T \beta + B(s_i) + \epsilon_i,}
    where
    \eqn{Z(\boldsymbol{s}_i)=\boldsymbol{x}(\boldsymbol{s}_i)^\mathrm{T}
         \boldsymbol{\beta} +
      B(\boldsymbol{s}_i)}{Z(s_i) = x(s_i)^T \beta + B(s_i)} is the so-called signal,
    \eqn{\boldsymbol{x}(\boldsymbol{s}_i)^\mathrm{T}
      \boldsymbol{\beta}}{x(s_i)^T\beta}
    is the external drift,
    \eqn{\{B(\boldsymbol{s})\}}{\{B(s)\}} is an unobserved stationary or
    intrinsic spatial Gaussian random field with zero mean, and
    \eqn{\varepsilon_i}{\epsilon_i} is an
    \emph{i.i.d} error from a possibly long-tailed distribution with scale parameter
    \eqn{\tau} (\eqn{\tau^2} is usually called nugget effect).
    In vector form the model is written as
    \deqn{\boldsymbol{Y} = \boldsymbol{X} \boldsymbol{\beta} + \boldsymbol{B} + \boldsymbol{\varepsilon},}{%
      Y = X \beta + B + \epsilon,}
    where \eqn{\boldsymbol{X}}{X} is the model matrix with the
    rows
    \eqn{\boldsymbol{x}(\boldsymbol{s}_i)^\mathrm{T}}{%
      x^T(s_i)}.

    The (generalized) covariance matrix of the vector of
    spatial Gaussian random effects
    \eqn{\boldsymbol{B}}{B}
    is denoted by
    \deqn{\mathrm{E}[
      \boldsymbol{B}\,
      \boldsymbol{B}^\mathrm{T}] =
      \boldsymbol{\Gamma}_\theta =
      \sigma_{\mathrm{n}}^2\boldsymbol{I} +
      \sigma^2\boldsymbol{V}_\alpha =
      \sigma_B^2 \, \boldsymbol{V}_{\alpha,\xi} =
      \sigma_B^2 \, ((1-\xi) \, \boldsymbol{I} +
      \xi\, \boldsymbol{V}_\alpha )
      ,}{%
      E[B B^T] = \Gamma_\theta = \sigma_n^2 I+\sigma^2 V_\alpha =
      \sigma_B^2 V_{\alpha,\xi} = \sigma_B^2 ((1-\xi) I - \xi V_\alpha)
      ,}
    where
    \eqn{\sigma_{\mathrm{n}}^2}{\sigma_n^2}
    is the variance of seemingly uncorrelated micro-scale variation in
    \eqn{B(\boldsymbol{s})}{B(s)}
    that cannot be resolved with the chosen sampling design,
    \eqn{\boldsymbol{I}} is the identity matrix,
    \eqn{\sigma^2} is the variance of the captured auto-correlated variation in
    \eqn{B(\boldsymbol{s})}{B(s)},
    \eqn{\sigma_B^2=\sigma_{\mathrm{n}}^2+\sigma^2}{\sigma_B^2=\sigma_n^2+\sigma^2}
    is the signal variance, and
    \eqn{\xi=\sigma^2/\sigma_B^2}{\xi=\sigma^2/\sigma_B^2}.
    To estimate both
    \eqn{\sigma_{\mathrm{n}}^2}{\sigma_n^2} and \eqn{\tau^2} (and not only their sum), one needs
    replicated measurements for some of the
    \eqn{\boldsymbol{s}_i}{s_i}.

    We define
    \eqn{\boldsymbol{V}_{\alpha}}{V_{\alpha}}
    to be the (generalized) correlation matrix with elements
    \deqn{(\boldsymbol{V}_{\alpha})_{ij} =
      \gamma_0 - \gamma(|\boldsymbol{A}\;(
      \boldsymbol{s}_i-\boldsymbol{s}_j)|),}{
      (V_{\alpha})_ij = \gamma_0 - \gamma(|A (s_i - s_j )|),}
    where the constant \eqn{\gamma_0} is chosen large enough so that
    \eqn{\boldsymbol{V}_{\alpha}}{V_{\alpha}}
    is positive definite,
    \eqn{\gamma(\cdot)}{v()} is a valid stationary or intrinsic variogram, and
    \eqn{\boldsymbol{A} =
      \boldsymbol{A}(\alpha, f_1, f_2; \omega, \phi, \zeta)}{
      A = A(\alpha, f_1, f_2; \omega, \phi, \zeta)}
    is a matrix that is used to model geometrically anisotropic auto-correlation.
    In more detail,
    \eqn{\boldsymbol{A}}{A}
    maps an arbitrary point on an ellipsoidal surface with constant semi-variance in
    \eqn{ \mathrm{I}\!\mathrm{R}^3}{R^3},
    centred on the origin, and having lengths of semi-principal axes,
    \eqn{\boldsymbol{p}_1}{p_1},
      \eqn{\boldsymbol{p}_2}{p_2},
      \eqn{\boldsymbol{p}_3}{p_3},
    equal to
    \eqn{|\boldsymbol{p}_1|=\alpha}{|p_1|=\alpha},
      \eqn{|\boldsymbol{p}_2|=f_1\,\alpha}{|p_2|=f_1 \alpha} and
      \eqn{|\boldsymbol{p}_3|=f_2\,\alpha}{|p_3|=f_2 \alpha},
      \eqn{0 < f_2 \leq f_1 \leq 1}{0 < f_2 <= f_1 <= 1},
    respectively, onto the surface of the unit ball centred on the origin.

    The orientation of the ellipsoid is defined by the three angles
    \eqn{\omega}, \eqn{\phi} and \eqn{\zeta}:
    \describe{

      \item{\eqn{\omega}}{is the azimuth of \eqn{\boldsymbol{p}_1}{p_1}
        (= angle between north and the projection
        of
        \eqn{\boldsymbol{p}_1}{p_1}
        onto the \eqn{x}-\eqn{y}-plane,
        measured from north to south positive clockwise in degrees),
      }

      \item{\eqn{\phi}}{is 90 degrees minus the latitude of
        \eqn{\boldsymbol{p}_1}{p_1}
        (= angle between the zenith and
        \eqn{\boldsymbol{p}_1}{p_1},
        measured from zenith to nadir positive clockwise in degrees), and
      }

      \item{\eqn{\zeta}}{is the angle between
        \eqn{\boldsymbol{p}_2}{p_2}
        and the direction of the line, say \eqn{y^\prime}{y'},
        defined by the intersection between the
        \eqn{x}-\eqn{y}-plane and the plane orthogonal to
        \eqn{\boldsymbol{p}_1}{p_1} running through the origin
        (\eqn{\zeta} is measured from \eqn{y^\prime}{y'} positive counter-clockwise in degrees).
      }
    }

    The transformation matrix is given by
    \deqn{\boldsymbol{A}=
      \left(\begin{array}{ccc}
        1/\alpha & 0 & 0\\
        0 & 1/(f_1\,\alpha) & 0\\
        0 & 0 & 1/(f_2\,\alpha) \\
      \end{array}\right)
      (
        \boldsymbol{C}_1,
        \boldsymbol{C}_2,
        \boldsymbol{C}_3,
      )
    }{A=diag(1/\alpha, 1/(f_1\,\alpha),1/(f_2\,\alpha)) (C_1, C_2, C_3)}
    where
    \deqn{\boldsymbol{C}_1^\mathrm{T} = ( \sin\omega \sin\phi, -\cos\omega \cos\zeta - \sin\omega \cos\phi \sin\zeta, \cos\omega \sin\zeta - \sin\omega \cos\phi \cos\zeta )
    }{C_1^T=( sin\omega sin\phi, -cos\omega cos\zeta - sin\omega cos\phi sin\zeta, cos\omega sin\zeta - sin\omega cos\phi cos\zeta )
    }
    \deqn{\boldsymbol{C}_2^\mathrm{T} = ( \cos\omega \sin\phi, \sin\omega \cos\zeta - \cos\omega \cos\phi \sin\zeta, -\sin\omega \sin\zeta - \cos\omega \cos\phi\cos\zeta )
    }{C_2^T=( cos\omega sin\phi, sin\omega cos\zeta - cos\omega cos\phi sin\zeta, -sin\omega sin\zeta - cos\omega cos\phicos\zeta )
    }
    \deqn{\boldsymbol{C}_3^\mathrm{T} = (\cos\phi, \sin\phi \sin\zeta, \sin\phi \cos\zeta )
    }{C_3^T=(cos\phi, sin\phi sin\zeta, sin\phi cos\zeta )
    }
    To model geometrically anisotropic variograms in
    \eqn{ \mathrm{I}\!\mathrm{R}^2}{R^2}
    one has to set \eqn{\phi=90} and \eqn{f_2 = 1},
    and for \eqn{f_1 = f_2 = 1}
    one obtains the model for isotropic auto-correlation
    with range parameter \eqn{\alpha}.
    Note that for isotropic auto-correlation the software processes data for
    which \eqn{d} may exceed 3.

    Two remarks are in order:

    \enumerate{

      \item Clearly, the (generalized) covariance matrix of the observations
      \eqn{\boldsymbol{Y}}{Y} is given by
      \deqn{\mathrm{Cov}[\boldsymbol{Y},\boldsymbol{Y}^\mathrm{T}]
        = \tau^2 \boldsymbol{I} + \boldsymbol{\Gamma}_\theta.
      }{
        Cov[Y, Y^T] = \tau^2 I + \Gamma_\theta.
      }

      \item Depending on the context, the term \dQuote{variogram
      parameters} denotes sometimes all parameters of a geometrically
      anisotropic variogram model, but in places only the parameters of an
      isotropic variogram model, i.e. \eqn{\sigma^2, \ldots, \alpha,
      \ldots} and \eqn{f_1, \ldots, \zeta} are denoted by the term
      \dQuote{anisotropy parameters}.  In the sequel
      \eqn{\boldsymbol{\theta}}{\theta} is used to denote all
      variogram and anisotropy parameters except the nugget effect
      \eqn{\tau^2}.

    }


  }

%   estimation

  \subsection{Estimation}{
    The unobserved spatial random effects
    \eqn{\boldsymbol{B}}{B} at the data locations
    \eqn{\boldsymbol{s}_i}{s_i}
    and the model parameters
    \eqn{\boldsymbol{\beta}}{\beta}, \eqn{\tau^2} and
    \eqn{\boldsymbol{\theta}^\mathrm{T} =
      (\sigma^2, \sigma_{\mathrm{n}}^2, \alpha, \ldots,  f_{1}, f_{2},
        \omega, \phi, \zeta)
    }{\theta^T =
      (\sigma^2, \sigma^2_n, \tau^2, \alpha, f_1, f_2,
        \omega, \phi, \zeta, \ldots)
    }
    are unknown and are estimated in \pkg{georob} either by Gaussian
    (\cite{Harville, 1977}) or robust (\cite{\enc{Künsch}{Kuensch} et al., 2011})
    restricted maximum likelihood (REML) or
    Gaussian maximum likelihood (ML). Here \var{\ldots}
    denote further parameters of the variogram such as the smoothness parameter
    of the Whittle-\enc{Matérn}{Matern} model.

    In brief, the robust REML method is based on the insight that for
    given \eqn{\boldsymbol{\theta}}{\theta} and \eqn{\tau^2} the
    Kriging predictions (= BLUP) of
    \eqn{\boldsymbol{B}}{B} and the generalized least
    squares (GLS = ML) estimates of
    \eqn{\boldsymbol{\beta}}{\beta} can be obtained
    simultaneously by maximizing
    \deqn{
      - \sum_i
        \left(
          \frac{
            y_i -
            \boldsymbol{x}(\boldsymbol{s}_i)^\mathrm{T}
            \boldsymbol{\beta} -
            B(\boldsymbol{s}_i)
          }{\tau}
        \right)^2 -
        \boldsymbol{B}^{\mathrm{T}}
        \boldsymbol{\Gamma}^{-1}_\theta
        \boldsymbol{B}
    }{
      - \sum_i ( y_i - x^T(s_i)\beta - B(s_i) )^2 / \tau^2
      - B^T \Gamma_\theta^-1 B
    }
    with respect to
    \eqn{\boldsymbol{B}}{B} and
    \eqn{\boldsymbol{\beta}}{\beta}, e.g.
    \cite{Harville (1977)}.
    %
    Hence, the BLUP of \eqn{\boldsymbol{B}}{B},
    ML estimates of \eqn{\boldsymbol{\beta}}{\beta},
    \eqn{\boldsymbol{\theta}}{\theta} and \eqn{\tau^2}
    are obtained by maximizing
    \deqn{
      - \log(\det(
        \tau^2 \boldsymbol{I} +
        \boldsymbol{\Gamma}_\theta
      )) -
      \sum_i
        \left(
          \frac{
            y_i -
            \boldsymbol{x}(\boldsymbol{s}_i)^\mathrm{T}
            \boldsymbol{\beta} -
            B(\boldsymbol{s}_i)
          }{\tau}
        \right)^2 -
        \boldsymbol{B}^{\mathrm{T}}
        \boldsymbol{\Gamma}^{-1}_\theta
        \boldsymbol{B}
    }{
      - log(det(  \tau^2 I + \Gamma_\theta) )
      - \sum_i ( y_i - x^T(s_i)\beta - B(s_i) )^2 / \tau^2
      - B^T \Gamma_\theta^-1 B
    }
    jointly with respect to
    \eqn{\boldsymbol{B}}{B},
    \eqn{\boldsymbol{\beta}}{\beta},
    \eqn{\boldsymbol{\theta}}{\theta} and  \eqn{\tau^2}
    or by solving the respective estimating equations.

    The estimating equations can then by robustified by
    \itemize{
      \item replacing the standardized errors, say
        \eqn{\varepsilon_i/\tau = (
          y_i - \boldsymbol{x}(\boldsymbol{s}_i)^\mathrm{T}
          \boldsymbol{\beta} -
          B(\boldsymbol{s}_i)
        ) / \tau}{\epsilon_i/\tau = (y_i - x^T(s_i)\beta - B(s_i)) / \tau},
        by a bounded or re-descending \eqn{\psi}-function,
        \eqn{\psi_c(\varepsilon_i/\tau)}{\psi_c(\epsilon_i/\tau)},
        of them (e.g. \cite{Maronna et al, 2006, chap. 2}) and by
      \item introducing suitable bias correction terms for Fisher consistency at
        the Gaussian model,
    }
    see \cite{\enc{Künsch}{Kuensch} et al. (2011)} for details.
    The robustified estimating equations
    are solved numerically by a combination of iterated re-weighted least squares
    (IRWLS) to estimate \eqn{\boldsymbol{B}}{B} and
    \eqn{\boldsymbol{\beta}}{\beta} for given
    \eqn{\boldsymbol{\theta}}{\theta} and \eqn{\tau^2}
    and nonlinear root finding by the function
    \code{\link[nleqslv]{nleqslv}} of the \R package \pkg{nleqslv}
    to get \eqn{\boldsymbol{\theta}}{\theta} and  \eqn{\tau^2}.
    The robustness of the procedure is controlled by the tuning parameter \eqn{c}
    of the \eqn{\psi_c}-function. For \eqn{c \ge 1000}{c>=1000} the algorithm computes
    Gaussian (RE)ML estimates and customary plug-in Kriging predictions.
    Instead of solving the Gaussian (RE)ML estimating equations, our software then
    maximizes the Gaussian (restricted) log-likelihood using  \code{\link[stats]{nlminb}} or
    \code{\link[stats]{optim}}.

    \pkg{georob} uses variogram models that were provided formerly by the
    now archived \R package \pkg{RandomFields} and are now implemented in
    the function \code{\link{gencorr}} of \pkg{georob}.
    Currently, estimation of the parameters of the following models is
    implemented:

    \code{"RMaskey"}, \code{"RMbessel"}, \code{"RMcauchy"},
    \code{"RMcircular"}, \code{"RMcubic"}, \code{"RMdagum"}, \cr
    \code{"RMdampedcos"}, \code{"RMdewijsian"}, \code{"RMexp"} (default),
    \code{"RMfbm"}, \code{"RMgauss"}, \cr \code{"RMgencauchy"},
    \code{"RMgenfbm"}, \code{"RMgengneiting"}, \code{"RMgneiting"},
    \code{"RMlgd"}, \cr \code{"RMmatern"}, \code{"RMpenta"}, \code{"RMqexp"},
    \code{"RMspheric"}, \code{"RMstable"}, \code{"RMwave"}, \cr
    \code{"RMwhittle"}.

    For most variogram parameters, closed-form expressions of \eqn{\partial
    \gamma/ \partial \theta_i}{d\gamma/d\theta_i} and \eqn{\partial \gamma/
    \partial \tau^2}{d\gamma/d\tau^2} are used in the computations.
    However, for the parameter \eqn{\nu} of the models \code{"RMbessel"},
    \code{"RMmatern"} and \code{"RMwhittle"} \eqn{\partial \gamma/ \partial
    \nu}{d\gamma/d\nu} is evaluated numerically by the function
    \code{\link[stats]{numericDeriv}}, and this results in an increase in
    computing time when \eqn{\nu} is estimated.  }

  %     prediction

  \subsection{Prediction}{
    Customary and robust plug-in external drift point Kriging predictions
    can be computed for an non-sampled location
    \eqn{\boldsymbol{s}_0}{s_0}
    from the covariates
    \eqn{\boldsymbol{x}(\boldsymbol{s}_0)}{
      x(s_0)},
    the estimated parameters
    \eqn{\widehat{\boldsymbol{\beta}}}{hat\beta},
    \eqn{\widehat{\boldsymbol{\theta}}}{hat\theta}
    and the predicted random effects
    \eqn{\widehat{\boldsymbol{B}}}{hatB}
    by
    \deqn{
      \widehat{Y}(\boldsymbol{s}_0) = \widehat{Z}(\boldsymbol{s}_0) =
      \boldsymbol{x}(\boldsymbol{s}_0)^\mathrm{T}
      \widehat{\boldsymbol{\beta}} +
      \boldsymbol{\gamma}^\mathrm{T}_{\widehat{\theta}}(\boldsymbol{s}_0)
      \boldsymbol{\Gamma}^{-1}_{\widehat{\theta}}
      \widehat{\boldsymbol{B}},
    }{
      hatY(s_0) = hatZ(s_0) = x^T(s_0) hat\beta + \gamma^T_hat\theta(s_0) \Gamma_hat\theta^-1 hatB,
    }

    where
    \eqn{\boldsymbol{\Gamma}_{\widehat{\theta}}}{\Gamma_hat\theta}
    is the estimated (generalized) covariance matrix of
    \eqn{\boldsymbol{B}}{B} and
    \eqn{\boldsymbol{\gamma}_{\widehat{\theta}}(\boldsymbol{s}_0)}{
    \gamma^T_hat\theta(s_0)}
    is the vector with the estimated (generalized) covariances between
    \eqn{\boldsymbol{B}}{B} and
    \eqn{B(\boldsymbol{s}_0)}{B(s_0)}.
    Kriging variances can be computed as well, based on approximated covariances of
    \eqn{\widehat{\boldsymbol{B}}}{hatB} and
    \eqn{\widehat{\boldsymbol{\beta}}}{hat\beta}
    (see \cite{\enc{Künsch}{Kuensch} et al., 2011}, and appendices of
    \cite{Nussbaum et al., 2014}, for details).

    The package \pkg{georob} provides in addition software for computing
    customary and robust external drift \emph{block} Kriging predictions.
    The required integrals of the generalized covariance function are
    computed by functions of the \R package \pkg{constrainedKriging}.  }

	% functionality

  \subsection{Functionality}{
    For the time being, the functionality of \pkg{georob} is limited
    to geostatistical analyses of \emph{single} response variables.
    No software is currently available for customary and robust
    multivariate geostatistical analyses.
    \pkg{georob} offers functions for:
    \enumerate{

      \item Robustly fitting a spatial linear model to data that are
      possibly contaminated by independent errors from a long-tailed
      distribution by robust REML (see functions \code{\link{georob}} ---
      which also fits such models efficiently by Gaussian (RE)ML ---
      \code{\link{profilelogLik}} and \code{\link{control.georob}}).

      \item Extracting estimated model components (see
      \code{\link{residuals.georob}}, \code{\link{rstandard.georob}},
      % \code{\link{rstudent.georob}},
      \code{\link{ranef.georob}}).

      \item Robustly estimating sample variograms and for fitting variogram
      model functions to them (see \code{\link{sample.variogram}} and
      \code{\link{fit.variogram.model}}).

      \item Model building by forward and backward selection of covariates
      for the external drift (see \code{\link{waldtest.georob}},
      \code{\link{step.georob}}, \code{\link{add1.georob}},
      \code{\link{drop1.georob}}, \code{\link{extractAIC.georob}}, \cr
      \code{\link{logLik.georob}}, \code{\link{deviance.georob}}).  For a
      robust fit, the log-likelihood is not defined.  The function then
      computes the (restricted) log-likelihood of an equivalent Gaussian
      model with heteroscedastic nugget (see \code{\link{deviance.georob}}
      for details).

      \item Assessing the goodness-of-fit and predictive power of the model
      by \var{K}-fold cross-validation (see \code{\link{cv.georob}} and
      \code{\link{validate.predictions}}).

      \item Computing customary and robust external drift point and block
      Kriging predictions (see \code{\link{predict.georob}},
      \code{\link{control.predict.georob}}).

      \item Unbiased back-transformation of both point and block Kriging
      predictions of log-transformed data to the original scale of the
      measurements (see \code{\link{lgnpp}}).

      \item Computing unconditional and conditional Gaussian simulations
      from a fitted spatial linear model (see \code{\link{condsim}}).

    }
  }
}

\references{
  Harville, D. A. (1977) Maximum likelihood approaches to variance
  component estimation and to related problems, \emph{Journal of the
  American Statistical Association}, \bold{72}, 320--340,
  \doi{10.1080/01621459.1977.10480998}.

  \enc{Künsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A.
  (in preparation) Robust Geostatistics.

  \enc{Künsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A.
  (2011) Robust estimation of the external drift and the variogram of
  spatial data.  Proceedings of the ISI 58th World Statistics Congress of
  the International Statistical Institute.  \doi{10.3929/ethz-a-009900710}

  Maronna, R. A., Martin, R. D. and Yohai, V. J. (2006) \emph{Robust
  Statistics Theory and Methods}, Wiley, Hoboken, \doi{10.1002/0470010940}.

  Nussbaum, M., Papritz, A., Baltensweiler, A. and Walthert, L. (2014)
  Estimating soil organic carbon stocks of {Swiss} forest soils by robust
  external-drift kriging.  \emph{Geoscientific Model Development},
  \bold{7}, 1197--1210.  \doi{10.5194/gmd-7-1197-2014}.

}

\author{
   Andreas Papritz \email{papritz@retired.ethz.ch}.
}

\seealso{
%  \code{\link{georobPackage}} for a description of the model and a brief summary of the algorithms;
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
  from model fitted by \code{georob}; and finally

  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
  for robust estimation and modelling of sample variograms.
}

%
\keyword{models}
\keyword{spatial}
\keyword{robust}
\keyword{package}
