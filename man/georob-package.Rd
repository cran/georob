% 2016-02-26 A. Papritz
% R CMD Rd2pdf --force georob-package.Rd;  R CMD Rdconv -t html -o bla.html georob-package.Rd ; open bla.html; 
% Rd2latex("georob-package.Rd", "../vignettes/georob-package.tex", outputEncoding="macintosh")

\encoding{macintosh}
\name{georob-package}
\alias{georob-package}
\alias{georobIntro}
\docType{package}
%
\title{The georob Package}
%
\description{This is a summary of the features and functionality of
  \pkg{georob}, a package in \R for robust geostatistical analyses.}
%
\details{

% intro
 
  \pkg{georob} is a package for robust analyses of geostatistical data.
  Such data, say \eqn{y_i=y(\mbox{\boldmath$s$\unboldmath}_i)}{y(s_i)}, are
  recorded at a set of locations,
  \eqn{\mbox{\boldmath$s$\unboldmath}_i}{s_i}, \eqn{i=1,2, \ldots, n}, in a
  domain \eqn{G \in \mathrm{I}\!\mathrm{R}^d}{G in R^d}, \eqn{d \in (1,2,3)}{d in (1,2,3)}, along
  with covariate information
  \eqn{x_j(\mbox{\boldmath$s$\unboldmath}_i)}{x_j(s_i)}, \eqn{j=1,2,
  \ldots, p}.
  
%   model
    
  \subsection{Model}{
    We use the following model for the data 
    \eqn{y_i=y(\mbox{\boldmath$s$\unboldmath}_{i})}{y_i=y(s_i)}:
    \deqn{Y(\mbox{\boldmath$s$\unboldmath}_i) = 
      Z(\mbox{\boldmath$s$\unboldmath}_i) + \varepsilon =
      \mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_i)^\mathrm{T}
         \mbox{\boldmath$\beta$\unboldmath} +
      B(\mbox{\boldmath$s$\unboldmath}_i) +
      \varepsilon_i,}{%
      Y(s_i) = Z(s_i) + \epsilon = x(s_i)^T \beta + B(s_i) + \epsilon_i,} 
    where 
    \eqn{Z(\mbox{\boldmath$s$\unboldmath}_i)=\mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_i)^\mathrm{T}
         \mbox{\boldmath$\beta$\unboldmath} +
      B(\mbox{\boldmath$s$\unboldmath}_i)}{Z(s_i) = x(s_i)^T \beta + B(s_i)} is the so-called signal, 
    \eqn{\mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_i)^\mathrm{T}
      \mbox{\boldmath$\beta$\unboldmath}}{x(s_i)^T\beta} 
    is the external drift, 
    \eqn{\{B(\mbox{\boldmath$s$\unboldmath})\}}{\{B(s)\}} is an unobserved stationary or
    intrinsic spatial Gaussian random field with zero mean, and
    \eqn{\varepsilon_i}{\epsilon_i} is an
    \emph{i.i.d} error from a possibly long-tailed distribution with scale parameter
    \eqn{\tau} (\eqn{\tau^2} is usually called nugget effect).
    In vector form the model is written as
    \deqn{\mbox{\boldmath $Y = X \beta + B + \varepsilon$\unboldmath},}{%
      Y = X \beta + B + \epsilon,} 
    where \eqn{\mbox{\boldmath $X$\unboldmath}}{X} is the model matrix with the
    rows 
    \eqn{\mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_i)^\mathrm{T}}{%
      x^T(s_i)}.
      
    The (generalized) covariance matrix of the vector of  
    spatial Gaussian random effects 
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} 
    is denoted by 
    \deqn{\mathrm{E}[
      \mbox{\boldmath$B$\unboldmath}\,
      \mbox{\boldmath$B$\unboldmath}^\mathrm{T}] =
      \mbox{\boldmath$\Gamma$\unboldmath}_\theta =
      \sigma_{\mathrm{n}}^2\mbox{\boldmath$I$\unboldmath} +
      \sigma^2\mbox{\boldmath$V$\unboldmath}_\alpha = 
      \sigma^2_Z \, \mbox{\boldmath$V$\unboldmath}_{\alpha,\xi} = 
      \sigma^2_Z \, (\xi \, \mbox{\boldmath$I$\unboldmath} + 
      (1-\xi)\, \mbox{\boldmath$V$\unboldmath}_\alpha )
      ,}{%
      E[B B^T] = \Gamma_\theta = \sigma_n^2 I+\sigma^2 V_\alpha =
      \sigma_Z^2 V_{\alpha,\xi} = \sigma_Z^2 (\xi I - (1-\xi) V_\alpha)
      ,}
    where 
    \eqn{\sigma_{\mathrm{n}}^2}{\sigma_n^2} 
    is the variance of seemingly uncorrelated micro-scale variation in 
    \eqn{B(\mbox{\boldmath$s$\unboldmath})}{B(s)} 
    that cannot be resolved with the chosen sampling design,  
    \eqn{\sigma^2} is the variance of the captured auto-correlated variation in 
    \eqn{B(\mbox{\boldmath$s$\unboldmath})}{B(s)},  
    \eqn{\sigma_Z^2=\sigma_{\mathrm{n}}^2+\sigma^2}{\sigma_Z^2=\sigma_n^2+\sigma^2} 
    is the signal variance, and 
    \eqn{\xi=\sigma_{\mathrm{n}}^2/\sigma_Z^2}{\xi=\sigma_n^2/\sigma^_Z^2}. 
    To estimate both 
    \eqn{\sigma_{\mathrm{n}}^2}{\sigma_n^2} and \eqn{\tau^2} (and not only their sum), one needs
    replicated measurements for some of the 
    \eqn{\mbox{\boldmath$s$\unboldmath}_i}{s_i}.
    
    We define 
    \eqn{\mbox{\boldmath$V$\unboldmath}_{\alpha}}{V_{\alpha}} 
    to be the matrix with elements
    \deqn{(\mbox{\boldmath$V$\unboldmath}_{\alpha})_{ij} =
      \gamma_0 - \gamma(|\mbox{\boldmath$A$\unboldmath}\;(
      \mbox{\boldmath$s$\unboldmath}_i-\mbox{\boldmath$s$\unboldmath}_j)|),}{
      (V_{\alpha})_ij = \gamma_0 - \gamma(|A (s_i - s_j )|),}
    where the constant \eqn{\gamma_0} is chosen large enough so that 
    \eqn{\mbox{\boldmath$V$\unboldmath}_{\alpha}}{V_{\alpha}}
    is positive definite,
    \eqn{\gamma(\cdot)}{v()} is a valid stationary or intrinsic variogram, and
    \eqn{\mbox{\boldmath$A$\unboldmath} = 
      \mbox{\boldmath$A$\unboldmath}(\alpha, f_1, f_2; \omega, \phi, \zeta)}{
      A = A(\alpha, f_1, f_2; \omega, \phi, \zeta)} 
    is a matrix that is used to model geometrically anisotropic auto-correlation.
    In more detail, 
    \eqn{\mbox{\boldmath$A$\unboldmath}}{A} 
    maps an arbitrary point on an ellipsoidal surface with constant semi-variance in 
    \eqn{ \mathrm{I}\!\mathrm{R}^3}{R^3},
    centred on the origin, and having lengths of semi-principal axes, 
    \eqn{\mbox{\boldmath$p$\unboldmath}_1}{p_1}, 
      \eqn{\mbox{\boldmath$p$\unboldmath}_2}{p_2}, 
      \eqn{\mbox{\boldmath$p$\unboldmath}_3}{p_3}, 
    equal to 
    \eqn{|\mbox{\boldmath$p$\unboldmath}_1|=\alpha}{|p_1|=\alpha}, 
      \eqn{|\mbox{\boldmath$p$\unboldmath}_2|=f_1\,\alpha}{|p_2|=f_1 \alpha} and 
      \eqn{|\mbox{\boldmath$p$\unboldmath}_3|=f_2\,\alpha}{|p_3|=f_2 \alpha},
      \eqn{0 < f_2 \leq f_1 \leq 1}{0 < f_2 <= f_1 <= 1},
    respectively, onto the surface of the unit ball centred on the origin.
    
    The orientation of the ellipsoid is defined by the three angles
    \eqn{\omega}, \eqn{\phi} and \eqn{\zeta}: 
    \describe{
    
      \item{\eqn{\omega}}{is the azimuth of \eqn{\mbox{\boldmath$p$\unboldmath}_1}{p_1}
        (= angle between north and the projection  
        of
        \eqn{\mbox{\boldmath$p$\unboldmath}_1}{p_1}
        onto the \eqn{x}-\eqn{y}-plane, 
        measured from north to south positive clockwise in degrees),
      }
      
      \item{\eqn{\phi}}{is 90 degrees minus the altitude of 
        \eqn{\mbox{\boldmath$p$\unboldmath}_1}{p_1}
        (= angle between the zenith and
        \eqn{\mbox{\boldmath$p$\unboldmath}_1}{p_1},
        measured from zenith to nadir positive clockwise in degrees), and
      }
      
      \item{\eqn{\zeta}}{is the angle between 
        \eqn{\mbox{\boldmath$p$\unboldmath}_2}{p_2}
        and the direction of the line, say \eqn{y^\prime}{y'}, 
        defined by the intersection between the
        \eqn{x}-\eqn{y}-plane and the plane orthogonal to 
        \eqn{\mbox{\boldmath$p$\unboldmath}_1}{p_1} running through the origin
        (\eqn{\zeta} is measured from \eqn{y^\prime}{y'} positive counter-clockwise in degrees).
      }
    }
    
    The transformation matrix is given by
    \deqn{\mbox{\boldmath$A$\unboldmath}=
      \left(\begin{array}{ccc}
        1/\alpha & 0 & 0\\
        0 & 1/(f_1\,\alpha) & 0\\
        0 & 0 & 1/(f_2\,\alpha) \\
      \end{array}\right) 
      (
        \mbox{\boldmath$C$\unboldmath}_1,
        \mbox{\boldmath$C$\unboldmath}_2,
        \mbox{\boldmath$C$\unboldmath}_3,
      ) 
    }{A=diag(1/\alpha, 1/(f_1\,\alpha),1/(f_2\,\alpha)) (C_1, C_2, C_3)}
    where
    \deqn{\mbox{\boldmath$C$\unboldmath}_1^\mathrm{T} = (\sin\phi \sin\omega, -\cos\zeta \cos\omega + \sin\zeta \cos\phi \sin\omega, -\cos\omega \sin\zeta - \cos\zeta \cos\phi \sin\omega)
    }{C_1^T=(sin\phi sin\omega, -cos\zeta cos\omega + sin\zeta cos\phi sin\omega, -cos\omega sin\zeta - cos\zeta cos\phi sin\omega)
    }
    \deqn{\mbox{\boldmath$C$\unboldmath}_2^\mathrm{T} = (-\sin\phi \cos\omega, \cos\omega \sin\zeta \cos\phi + \cos\zeta \sin\omega, -\cos\zeta \cos\omega \cos\phi + \sin\zeta \sin\omega)
    }{C_2^T=(-sin\phi cos\omega, cos\omega sin\zeta cos\phi + cos\zeta sin\omega, -cos\zeta cos\omega cos\phi + sin\zeta sin\omega)
    }
    \deqn{\mbox{\boldmath$C$\unboldmath}_3^\mathrm{T} = (\cos\phi, -\sin\phi \sin\zeta, \cos\zeta \sin\phi)
    }{C_3^T=(cos\phi, -sin\phi sin\zeta, cos\zeta sin\phi)
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
      \eqn{\mbox{\boldmath $Y$\unboldmath}}{Y} is given by
      \deqn{\mathrm{Cov}[\mbox{\boldmath $Y$\unboldmath},\mbox{\boldmath $Y$\unboldmath}^\mathrm{T}]
        = \tau^2 \mbox{\boldmath$I$\unboldmath} + \mbox{\boldmath$\Gamma$\unboldmath}_\theta.      
      }{
        Cov[Y, Y^T] = \tau^2 I + \Gamma_\theta.
      }
      
      \item Depending on the context, the term \dQuote{variogram
      parameters} denotes sometimes all parameters of a geometrically
      anisotropic variogram model, but in places only the parameters of an
      isotropic variogram model, i.e. \eqn{\sigma^2, \ldots, \alpha,
      \ldots} and \eqn{f_1, \ldots, \zeta} are denoted by the term
      \dQuote{anisotropy parameters}.  In the sequel
      \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta} is used to denote all
      variogram and anisotropy parameters except the nugget effect
      \eqn{\tau^2}.
      
    }
   
    
  }
  
%   estimation
  
  \subsection{Estimation}{
    The unobserved spatial random effects 
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} at the data locations 
    \eqn{\mbox{\boldmath$s$\unboldmath}_i}{s_i}
    and the model parameters 
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta}, \eqn{\tau^2} and
    \eqn{\mbox{\boldmath$\theta$\unboldmath}^\mathrm{T} = 
      (\sigma^2, \sigma_{\mathrm{n}}^2, \alpha, \ldots,  f_{1}, f_{2}, 
        \omega, \phi, \zeta)
    }{\theta^T = 
      (\sigma^2, \sigma^2_n, \tau^2, \alpha, f_1, f_2, 
        \omega, \phi, \zeta, \ldots)
    } 
    are unknown and are estimated in \pkg{georob} either by Gaussian or 
    robust restricted maximum likelihood (REML) or 
    Gaussian maximum likelihood (ML). Here \var{\ldots} 
    denote further parameters of the variogram such as the smoothness parameter
    of the Whittle-\enc{Matrn}{Matern} model.
    
    In brief, the robust REML method is based on the insight that for
    given \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta} and \eqn{\tau^2} the
    kriging predictions (= BLUP) of
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and the generalized least
    squares (GLS = ML) estimates of
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} can be obtained
    simultaneously by maximizing 
    \deqn{
      - \sum_i 
        \left(
          \frac{
            y_i - 
            \mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_i)^\mathrm{T}
            \mbox{\boldmath$\beta$\unboldmath} -
            B(\mbox{\boldmath$s$\unboldmath}_i)
          }{\tau}
        \right)^2 - 
        \mbox{\boldmath$B$\unboldmath}^{\mathrm{T}} 
        \mbox{\boldmath$\Gamma$\unboldmath}^{-1}_\theta
        \mbox{\boldmath$B$\unboldmath}
    }{
      - \sum_i ( y_i - x^T(s_i)\beta - B(s_i) )^2 / \tau^2
      - B^T \Gamma_\theta^-1 B  
    }
    with respect to
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} e.g.
    \cite{Harville (1977)}.
    %
    Hence, the BLUP of \eqn{\mbox{\boldmath$B$\unboldmath}}{B}, 
    ML estimates of \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta}, 
    \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta} and \eqn{\tau^2}
    are obtained by maximizing 
    \deqn{
      - \log(\det(
        \tau^2 \mbox{\boldmath$I$\unboldmath} + 
        \mbox{\boldmath$\Gamma$\unboldmath}_\theta        
      )) -
      \sum_i 
        \left(
          \frac{
            y_i - 
            \mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_i)^\mathrm{T}
            \mbox{\boldmath$\beta$\unboldmath} -
            B(\mbox{\boldmath$s$\unboldmath}_i)
          }{\tau}
        \right)^2 - 
        \mbox{\boldmath$B$\unboldmath}^{\mathrm{T}} 
        \mbox{\boldmath$\Gamma$\unboldmath}^{-1}_\theta
        \mbox{\boldmath$B$\unboldmath}
    }{
      - log(det(  \tau^2 I + \Gamma_\theta) )
      - \sum_i ( y_i - x^T(s_i)\beta - B(s_i) )^2 / \tau^2
      - B^T \Gamma_\theta^-1 B  
    }
    jointly with respect to 
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B},
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta}, 
    \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta} and  \eqn{\tau^2}
    or by solving the respective estimating equations.
    
    The estimating equations can then by robustified by 
    \itemize{
      \item replacing the standardized residuals, say
        \eqn{\varepsilon_i/\tau}{\epsilon_i/\tau},
        by a bounded or re-descending \eqn{\psi}-function, 
        \eqn{\psi_c(\varepsilon_i/\tau)}{\psi_c(\epsilon_i/\tau)},
        of them (e.g. Marona et al, 2006, chap. 2) and by
      \item introducing suitable bias correction terms for Fisher consistency at 
        the Gaussian model,
    }
    see \cite{\enc{Knsch}{Kuensch} et al. (2011)} for details. 
    The robustified estimating equations
    are solved numerically by a combination of iterated re-weighted least squares
    (IRWLS) to estimate \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} for given 
    \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta} and \eqn{\tau^2}
    and nonlinear root finding by the function 
    \code{\link[nleqslv]{nleqslv}} of the \R package \pkg{nleqslv}
    to get \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta} and  \eqn{\tau^2}.
    The robustness of the procedure is controlled by the tuning parameter \eqn{c}
    of the \eqn{\psi_c}-function. For \eqn{c \ge 1000}{c>=1000} the algorithm computes
    Gaussian (RE)ML estimates and customary plug-in kriging predictions. 
    Instead of solving the Gaussian (RE)ML estimating equations, our software then
    maximizes the Gaussian (restricted) log-likelihood using  \code{\link[stats]{nlminb}} or
    \code{\link[stats]{optim}}.
        
    \pkg{georob} uses variogram models implemented in the \R package
    \pkg{RandomFields} (see \code{\link[RandomFields]{RMmodel}}).
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
    
    For most variogram parameters, closed-form
    expressions of \eqn{\partial \gamma/ \partial
    \theta_i}{d\gamma/d\theta_i} are used in the computations.  However,
    for the parameter \eqn{\nu} of the models \code{"RMbessel"},
    \code{"RMmatern"} and \code{"RMwhittle"} \eqn{\partial \gamma/ \partial
    \nu}{d\gamma/d\nu} is evaluated numerically by the function
    \code{\link[stats]{numericDeriv}}, and this results in an increase in
    computing time when \eqn{\nu} is estimated.  }
  
  %     prediction
  
  \subsection{Prediction}{ 
    Robust plug-in external drift point kriging predictions 
    can be computed for an non-sampled location 
    \eqn{\mbox{\boldmath$s$\unboldmath}_0}{s_0}
    from the covariates
    \eqn{\mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_0)}{
      x(s_0)},
    the estimated parameters 
    \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta},
    \eqn{\widehat{\mbox{\boldmath$\theta$\unboldmath}}}{hat\theta}
    and the predicted random effects
    \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB}
    by
    \deqn{
      \widehat{Y}(\mbox{\boldmath$s$\unboldmath}_0) = \widehat{Z}(\mbox{\boldmath$s$\unboldmath}_0) = 
      \mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_0)^\mathrm{T}
      \widehat{\mbox{\boldmath$\beta$\unboldmath}} + 
      \mbox{\boldmath$\gamma$\unboldmath}^\mathrm{T}_{\widehat{\theta}}(\mbox{\boldmath$s$\unboldmath}_0)
      \mbox{\boldmath$\Gamma$\unboldmath}^{-1}_{\widehat{\theta}}
      \widehat{\mbox{\boldmath$B$\unboldmath}},
    }{
      hatY(s_0) = hatZ(s_0) = x^T(s_0) hat\beta + \gamma^T_hat\theta(s_0) \Gamma_hat\theta^-1 hatB,
    }
    
    where
    \eqn{\mbox{\boldmath$\Gamma$\unboldmath}_{\widehat{\theta}}}{\Gamma_hat\theta}
    is the estimated (generalized) covariance matrix of
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
    \eqn{\mbox{\boldmath$\gamma$\unboldmath}_{\widehat{\theta}}(\mbox{\boldmath$s$\unboldmath}_0)}{
    \gamma^T_hat\theta(s_0)}
    is the vector with the estimated (generalized) covariances between 
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
    \eqn{B(\mbox{\boldmath$s$\unboldmath}_0)}{B(s_0)}.
    Kriging variances can be computed as well, based on approximated covariances of 
    \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} and 
    \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta}
    (see \cite{Knsch et al., 2011}, and appendices of
    \cite{Nussbaum et al., 2012}, for details).
    
    The package \pkg{georob} provides in addition software for computing
    robust external drift \emph{block} kriging predictions.  The required
    integrals of the generalized covariance function are computed by
    functions of the \R package \pkg{constrainedKriging}.
  }

% functionality
  
  \subsection{Functionality}{
    For the time being, the functionality of \pkg{georob} is limited 
    to robust geostatistical analyses of \emph{single} response variables. 
    No software is currently available for robust multivariate geostatistical 
    analyses.
    \pkg{georob} offers functions for: 
    \enumerate{
      \item Robustly fitting a spatial linear model to data that are 
      possibly contaminated by independent errors from a long-tailed
      distribution by robust REML (see functions \code{\link{georob}} --- which also
      fits such  models efficiently by Gaussian (RE)ML --- \code{\link{profilelogLik}} and 
      \code{\link{control.georob}}).
      \item Extracting estimated model components (see \code{\link{residuals.georob}}, 
      \code{\link{rstandard.georob}}, \cr 
      % \code{\link{rstudent.georob}}, 
      \code{\link{ranef.georob}}).
      \item Robustly estimating sample variograms and for fitting 
      variogram model functions to them (see 
      \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}).
      \item Model building by forward and backward selection of 
      covariates for the external drift 
      (see \code{\link{waldtest.georob}}, \code{\link{step.georob}}, 
      \code{\link{add1.georob}}, \code{\link{drop1.georob}}, 
      \code{\link{extractAIC.georob}}, \cr \code{\link{logLik.georob}}, 
      \code{\link{deviance.georob}}). 
      For  a robust fit, the log-likelihood is not defined.  The function then
      computes the (restricted) log-likelihood of an equivalent Gaussian model with
      heteroscedastic nugget (see \code{\link{deviance.georob}} for details).
      \item Assessing the goodness-of-fit and predictive power of the model by 
      \var{K}-fold cross-validation (see \code{\link{cv.georob}} and 
      \code{\link{validate.predictions}}).
      \item Computing robust external drift point and block kriging 
      predictions (see \code{\link{predict.georob}}, \code{\link{control.predict.georob}}).
      \item Unbiased back-transformation of both point and block
      kriging predictions of log-transformed data to the original scale of
      the measurements (see \code{\link{lgnpp}}).
    }
  }
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}\cr
  \url{http://www.step.ethz.ch/people/scientific-staff/andreas-papritz}
}

\references{
  Nussbaum, M., Papritz, A., Baltensweiler, A. and Walthert, L. (2012)
  \emph{Organic Carbon  Stocks of Swiss Forest Soils},
  Institute of Terrestrial Ecosystems, ETH Zurich and 
  Swiss Federal Institute for Forest, Snow and Landscape Research
  (WSL), pp. 51. 
  \url{http://e-collection.library.ethz.ch/eserv/eth:6027/eth-6027-01.pdf}

  \enc{Knsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A. (in
  preparation) Robust Geostatistics.

  \enc{Knsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A. (2011) Robust
  estimation of the external drift and the variogram of spatial data.
  Proceedings of the ISI 58th World Statistics Congress of the International
  Statistical Institute.
  \url{http://e-collection.library.ethz.ch/eserv/eth:7080/eth-7080-01.pdf}
  
  Maronna, R. A., Martin, R. D. and Yohai, V. J. (2006) Robust Statistics Theory and
  Methods, John Wiley \& Sons.
}

\seealso{
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{georobObject}} for a description of the class \code{georob};
  \code{\link{plot.georob}} for display of RE(ML) variogram estimates;
  \code{\link{control.georob}} for controlling the behaviour of \code{georob};
  \code{\link{cv.georob}} for assessing the goodness of a fit by \code{georob}; 
  \code{\link{predict.georob}} for computing robust kriging predictions; and finally
  \code{\link{georobModelBuilding}} for stepwise building models of class \code{georob};
  \code{\link{georobMethods}} for further methods for the class \code{georob},
  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}} 
  for robust estimation and modelling of sample variograms.
}

%
\keyword{models}
\keyword{spatial}
\keyword{robust}
\keyword{package}
