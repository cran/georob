% 2012-12-12 A. Papritz
% R CMD Rdconv -t html -o bla.html georob-package.Rd ; open bla.html; R CMD Rd2pdf --force georob-package.Rd; 

\encoding{macintosh}
\name{georob-package}
\alias{georob-package}
\alias{georobIntro}
\docType{package}
%
\title{The georob Package}
%
\description{This is a summary of the features of \pkg{georob}, 
  a package in \R for robust geostatistical analyses.}
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
    \deqn{Y(\mbox{\boldmath$s$\unboldmath}) = 
      S(\mbox{\boldmath$s$\unboldmath}) + 
      \varepsilon(\mbox{\boldmath$s$\unboldmath}) =
      \mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath})^\mathrm{T}
         \mbox{\boldmath$\beta$\unboldmath} +
      B(\mbox{\boldmath$s$\unboldmath}) +
      \varepsilon(\mbox{\boldmath$s$\unboldmath}),}{%
      Y(s) = S(s) + \epsilon(s) = x(s)^T \beta + B(s) + \epsilon(s),} 
    where 
    \eqn{S(\mbox{\boldmath$s$\unboldmath})=\mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath})^\mathrm{T}
         \mbox{\boldmath$\beta$\unboldmath} +
      B(\mbox{\boldmath$s$\unboldmath})}{S(s) = x(s)^T \beta + B(s)} is the so-called signal, 
    \eqn{\mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath})^\mathrm{T}
      \mbox{\boldmath$\beta$\unboldmath}}{x^T\beta} 
    is the external drift, 
    \eqn{B(\mbox{\boldmath$s$\unboldmath})}{B(s)} is a stationary or
    intrinsic spatial Gaussian random field with zero mean, and
    \eqn{\varepsilon(\mbox{\boldmath$s$\unboldmath})}{\epsilon(s)} are
    \emph{i.i.d} errors from a possibly long-tailed distribution with scale parameter
    \eqn{\tau} (\eqn{\tau^2} is usually called nugget effect).
    In vector form the model leads to 
    \deqn{\mbox{\boldmath $Y = X \beta + B + \varepsilon$\unboldmath},}{%
      Y = X \beta + B + \epsilon,} 
    where \eqn{\mbox{\boldmath $X$\unboldmath}}{X} is the model matrix consisting
    of the rows 
    \eqn{\mbox{\boldmath $x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_i)^\mathrm{T}}{%
      x^T(s_i)}.
      
    The (generalized) covariance matrix of the vector of unobserved 
    spatial Gaussian random effects 
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} 
    is denoted by 
    \deqn{\mathrm{E}[
      \mbox{\boldmath$B$\unboldmath}\,
      \mbox{\boldmath$B$\unboldmath}^\mathrm{T}] =
      \mbox{\boldmath$\Gamma$\unboldmath}_\theta =
      \sigma_{\mathrm{n}}^2\mbox{\boldmath$I$\unboldmath} +
      \sigma^2\mbox{\boldmath$V$\unboldmath}_\alpha,}{%
      E[B B^T] = \sigma_n^2 I+\sigma^2 V_\alpha,}
    where \eqn{\sigma_{\mathrm{n}}^2}{\sigma_n^2} is the variance of seemingly 
    uncorrelated micro-scale variation in 
    \eqn{B(\mbox{\boldmath$s$\unboldmath})}{B(s)} 
    that cannot be resolved with the chosen sampling design, and 
    \eqn{\sigma^2} is the variance of the captured auto-correlated variation in 
    \eqn{B(\mbox{\boldmath$s$\unboldmath})}{B(s)}. To estimate both 
    \eqn{\sigma_\mathrm{n}^2} and \eqn{\tau^2} (and not only their sum), one needs
    replicated measurements for some of the 
    \eqn{\mbox{\boldmath$s$\unboldmath}_i}.
    
    We define 
    \eqn{\mbox{\boldmath$V$\unboldmath}_\alpha}{V_\alpha} 
    to be the matrix with elements
    \deqn{(\mbox{\boldmath$V$\unboldmath}_\alpha)_{ij} =
      m - \gamma(|\mbox{\boldmath$A$\unboldmath}\;(
      \mbox{\boldmath$s$\unboldmath}_i-\mbox{\boldmath$s$\unboldmath}_j) |),}{
      (V_\alpha)_ij = \gamma_0 - \gamma(|A (s_i - s_j )|),}
    where the constant \eqn{\gamma_0} is chosen large enough so that 
    \eqn{\mbox{\boldmath$V$\unboldmath}_\alpha}{V_\alpha}
    is positive definite,
    \eqn{\gamma(\cdot)}{v()} is a valid stationary or intrinsic variogram, and
    \eqn{\mbox{\boldmath$A$\unboldmath} = 
      \mbox{\boldmath$A$\unboldmath}(\alpha, f_1, f_2; \omega, \phi, \zeta)}{
      A = A(\alpha, f_1, f_2; \omega, \phi, \zeta)} 
    is a matrix that is used to model geometrically anisotropic auto-correlation.
    In more detail, 
    \eqn{\mbox{\boldmath$A$\unboldmath}}{A} 
    maps an arbitrary point on an ellipsoidal surface with constant semivariance in 
    \eqn{ \mathrm{I}\!\mathrm{R}^d}{R^d},
    centred on the origin, and having lengths of semi-principal axes, 
    \eqn{\mbox{\boldmath$P$\unboldmath}_1}{P_1}, 
      \eqn{\mbox{\boldmath$P$\unboldmath}_2}{P_2}, 
      \eqn{\mbox{\boldmath$P$\unboldmath}_3}{P_3}, 
    equal to 
    \eqn{|\mbox{\boldmath$P$\unboldmath}_1|=\alpha}{|P_1|=\alpha}, 
      \eqn{|\mbox{\boldmath$P$\unboldmath}_2|=f_1\,\alpha}{|P_2|=f_1 \alpha} and 
      \eqn{|\mbox{\boldmath$P$\unboldmath}_3|=f_2\,\alpha}{|P_3|=f_2 \alpha},
      \eqn{0 < f_2 \leq f_1 \leq 1}{0 < f_2 <= f_1 <= 1},
    respectively, onto the surface of the unit ball centred on the origin.
    
    The orientation of the ellipsoid is defined by the three angles
    \eqn{\omega}, \eqn{\phi} and \eqn{\zeta}: 
    \describe{
    
      \item{\eqn{\omega}}{is the azimuth of \eqn{\mbox{\boldmath$P$\unboldmath}_1}{P_1}
        (= angle between north and the projection  
        of
        \eqn{\mbox{\boldmath$P$\unboldmath}_1}{P_1}
        onto the \eqn{x}-\eqn{y}-plane, 
        measured from north to south positive clockwise in degrees),
      }
      
      \item{\eqn{\phi}}{is 90 degrees minus the altitude of 
        \eqn{\mbox{\boldmath$P$\unboldmath}_1}{P_1}
        (= angle between the zenith and
        \eqn{\mbox{\boldmath$P$\unboldmath}_1}{P_1},
        measured from zenith to nadir positive clockwise in degrees), and
      }
      
      \item{\eqn{\zeta}}{is the angle between 
        \eqn{\mbox{\boldmath$P$\unboldmath}_2}{P_2}
        and the direction of the line, say \eqn{y^\prime}{y'}, 
        defined by the intersection between the
        \eqn{x}-\eqn{y}-plane and the plane orthogonal to 
        \eqn{\mbox{\boldmath$P$\unboldmath}_1}{P_1} running through the origin
        (\eqn{\zeta} is measured from \eqn{y^\prime}{y'} positive counterclockwise in degrees).
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
    Note that for isotropic auto-correlation \eqn{d} may exceed 3.
    
    Depending on the context, the term \dQuote{variogram parameters}
    denotes sometimes all parameters of a geometrically anisotropic
    variogram model, but in places only the parameters of an isotropic
    variogram model, i.e. \eqn{\sigma^2, \ldots, \alpha, \ldots} and
    \eqn{f_1, \ldots, \zeta} are denoted by the term \dQuote{anisotropy
    parameters}. 
    
  }
  
%   estimation
  
  \subsection{Estimation}{
    The  spatial random effects 
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B}
    and the  model parameters 
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} and
    \eqn{\mbox{\boldmath$\theta$\unboldmath}^\mathrm{T} = 
      (\sigma^2, \sigma_{\mathrm{n}}^2, \tau^2, \alpha, \newline \ldots,  f_{1}, f_{2}, 
        \omega, \phi, \zeta)
    }{\theta^T = 
      (\sigma^2, \sigma^2_n, \tau^2, \alpha, f_1, f_2, 
        \omega, \phi, \zeta, \ldots)
    } 
    are estimated by robust Restricted Maximum Likelihood (REML). Here \var{\ldots} 
    denote further parameters of the variogram such as the smoothness parameter
    of the Whittle-MatŽrn model.
    In brief, the robust REML method is based on the insight that 
    the kriging predictions of \eqn{\mbox{\boldmath$B$\unboldmath}}{B} 
    and the Gaussian Maximum Likelihood estimates of 
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} and
    \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta} can be obtained 
    simultaneously by maximizing 
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
    with respect to 
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B},
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta}, 
    \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta}.
    The respective estimating equations can then by robustified by 
    \itemize{
      \item replacing the standardized residuals, say
        \eqn{\widehat{\varepsilon}/\widehat{\tau}}{r/s},
        by a bounded function, 
        \eqn{\psi_c(\widehat{\varepsilon}/\widehat{\tau})}{\psi_c(r/s)},
        of them and by
      \item introducting suitable bias correction terms for Fisher consistency at 
        the Gaussian model,
    }
    see KŸnsch et al. (2011) for details. The robustified estimating equations
    are solved numerically by a combination of iterated re-weighted least squares
    (IRWLS) to estimate \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
    \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} for given 
    \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta}
    and nonlinear root finding by the function 
    \code{\link[nleqslv]{nleqslv}} of the \R package \pkg{nleqslv}
    to get \eqn{\mbox{\boldmath$\theta$\unboldmath}}{\theta}.
    The robustness of the procedure is controlled by the tuning parameter \eqn{c}
    of the \eqn{\psi_c}-function. For \eqn{c \ge 1000}{c>=1000} the algorithm computes
    Gaussian REML estimates and customary plug-in kriging predictions. 
    Instead of solving the Gaussian REML estimating equations, our software then
    maximizes the Gaussian restricted loglikelihood using 
    \code{\link[stats]{optim}}.
    
    \pkg{georob} uses variogram models implemented in the \R package 
    \pkg{RandomFields} (see \code{\link[RandomFields]{Variogram}}). Currently, 
    estimation of the parameters of the following models is implemented: 
    \code{"bessel"}, \code{"cauchy"}, \code{"cauchytbm"}, \code{"circular"},
    \code{"cubic"}, \code{"dagum"}, \code{"dampedcosine"}, \code{"DeWijsian"},\cr
    \code{"exponential"}, \code{"fractalB"}, \code{"gauss"}, \code{"genB"},
    \code{"gencauchy"}, \code{"gengneiting"}, \code{"gneiting"}, \code{"lgd1"},
    \code{"matern"}, \code{"penta"}, \code{"power"}, \code{"qexponential"},
    \code{"spherical"}, \code{"stable"}, \code{"wave"}, \code{"whittle"}.
    For most variogram parameters, closed-form expressions of 
    \eqn{\partial \gamma/ \partial \theta_i}{d\gamma/d\theta_i}
    are used  in the computations. 
    However, for the parameter \eqn{\nu} of the models 
    \code{"bessel"}, \code{"matern"} and \code{"whittle"}
    \eqn{\partial \gamma/ \partial \theta_i}{d\gamma/d\theta_i} 
    is evaluated numerically by the function 
    \code{\link[stats]{numericDeriv}}, and this results in an 
    increase in computing time when \eqn{\nu} is estimated.
  }
  
  %     prediction
  
  \subsection{Prediction}{ 
    Robust plug-in external drift point kriging predictions 
    can be computed for an unsampled location 
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
      \widehat{Y}{\mbox{\boldmath$s$\unboldmath}_0} = 
      \mbox{\boldmath$x$\unboldmath}(\mbox{\boldmath$s$\unboldmath}_0)^\mathrm{T}
      \widehat{\mbox{\boldmath$\beta$\unboldmath}} + 
      \mbox{\boldmath$\gamma$\unboldmath}^\mathrm{T}_{\widehat{\theta}}(\mbox{\boldmath$s$\unboldmath}_0)
      \mbox{\boldmath$\Gamma$\unboldmath}^{-1}_{\widehat{\theta}}
      \widehat{\mbox{\boldmath$B$\unboldmath}},
    }{
      hatY(s_0) = x^T(s_0) hat\beta + \gamma_^That\theta(s_0) \Gamma_hat\theta^-1 hatB,
    }
    
    where
    \eqn{\mbox{\boldmath$\Gamma$\unboldmath}_{\widehat{\theta}}}{\Gamma_hat\theta}
    is the estimated (generalized) covariance matrix of
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
    \eqn{\mbox{\boldmath$\gamma$\unboldmath}_{\widehat{\theta}}(\mbox{\boldmath$s$\unboldmath}_0)}{
    gamma_hat\theta(s_0)}
    is the vector with the estimated (generalized) covariances between 
    \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
    \eqn{B(\mbox{\boldmath$s$\unboldmath}_0)}{B(s_0)}.
    Kriging variances can be computed as well, based on approximated covariances of 
    \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} and 
    \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta}
    (see KŸnsch et al., 2011, and Appendices of
    Nussbaum et al., 2012, for details).
    
    The package \pkg{georob} provides in addition software for computing
    robust external drift \emph{block} kriging predictions.  The required
    integrals of the generalized covariance function are computed by
    functions of the \R package \pkg{constrainedKriging}.
  }

  \subsection{Main functionality}{
    For the time being, the functionality of \pkg{georob} is limited 
    to robust geostatistical analyses of \emph{single} response variables. 
    No software is currently available for robust multivariate geostatistical 
    analyses.
    \pkg{georob} offers functions for: 
    \enumerate{
      \item Robustly fitting a spatial linear model to data that are 
      possibly contaminated by independent errors from a long-tailed
      distribution by robust REML (see \code{\link{georob}}, which also
      fits such  models efficiently by Gaussian REML).
      \item Assessing the goodness-of-fit of the model by 
      \var{K}-fold cross-validation (see \code{\link{cv.georob}}).
      \item Computing robust external drift point and block kriging 
      predictions (see \code{\link{predict.georob}}).
      \item Unbiased back-transformation of both point and block
      kriging predictions of log-transformed data to the original scale of
      the measurements (see \code{\link{lgnpp}}).
      \item Robustly estimating sample variograms and for fitting 
      variogram model functions to them (see 
      \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}).
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
  Institute of Terrestrial Ecosystems, ETH ZŸrich and 
  Swiss Federal Institute for Forest, Snow and Landscape Research
  (WSL), pp. 51. 
  \url{http://e-collection.library.ethz.ch/eserv/eth:6027/eth-6027-01.pdf}

  \enc{KŸnsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A. (in
  preparation) Robust Geostatistics.

  \enc{KŸnsch}{Kuensch}, H. R., Papritz, A., Schwierz, C. and Stahel, W. A. (2011) Robust
  estimation of the external drift and the variogram of spatial data.
  Proceedings of the ISI 58th World Statistics Congress of the International
  Statistical Institute.
  \url{http://isi2011.congressplanner.eu/pdfs/650268.pdf}
}

\seealso{
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{georobObject}} for a description of the class \code{georob};
  \code{\link{plot.georob}} for display of REML variogram estimates;
  \code{\link{georob.control}} for controlling the behaviour of \code{georob};
  \code{\link{cv.georob}} for assessing the goodness of a fit by \code{georob}; 
  \code{\link{predict.georob}} for computing robust kriging predictions; and finally
  \code{\link{georobMethods}} for further methods for the class \code{georob},
  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}} 
  for robust estimation and modelling of sample variograms.
}

%
\keyword{spatial}
\keyword{robust}
\keyword{package}
