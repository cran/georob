% 2024-02-09 A. Papritz
% R CMD Rdconv -t html -o bla.html gencorr.Rd ; open bla.html; R CMD Rd2pdf --force gencorr.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  gencorr.Rd

\encoding{utf8}
\name{gencorr}
\alias{gencorr}
%
\title{Variogram Models}
%
\description{
  The function \code{gencorr} computes intrinsic or
  stationary isotropic generalized correlations (= negative semi-variances
  computed with sill (variance) parameter equal to 1) for a set of basic
  variogram models formerly made available by the function \code{RFfctn} of
  the now archived \R package \pkg{RandomFields}.}

\usage{
gencorr(x, variogram.model, param)
}
%
\arguments{
  \item{x}{a numeric vector with scaled lag distances, i.e. lag distances
  divided by the range parameter \code{param["scale"]}.}

  \item{variogram.model}{a character keyword defining the variogram model.
  Currently, the following models are implemented:\cr
    \code{"RMaskey"}, \code{"RMbessel"}, \code{"RMcauchy"},
    \code{"RMcircular"}, \code{"RMcubic"}, \code{"RMdagum"}, \cr
    \code{"RMdampedcos"}, \code{"RMdewijsian"}, \code{"RMexp"} (default),
    \code{"RMfbm"}, \code{"RMgauss"}, \cr \code{"RMgencauchy"},
    \code{"RMgenfbm"}, \code{"RMgengneiting"}, \code{"RMgneiting"},
    \code{"RMlgd"}, \cr \code{"RMmatern"}, \code{"RMpenta"}, \code{"RMqexp"},
    \code{"RMspheric"}, \code{"RMstable"}, \code{"RMwave"}, \cr
    \code{"RMwhittle"}, see \emph{Details}.}

  \item{param}{a named numeric vector with values of the additional
  parameters of the variogram models such as the smoothness parameter of
  the Whittle-Matérn model, see \code{\link{param.names} for the names of
  these parameters}.  Note that some variogram models (e.g.
  \code{"RMcircular"}) do not have any additional parameters.}

}

\details{

  The name and parametrization of the variogram models originate from the
  function \code{RFfctn} of \pkg{RandomFields}.  The equations and further
  information are taken (with minor modifications) from the help pages of
  the respective functions of the archived \R package \pkg{RandomFields},
  version 3.3.14 (\cite{Schlather et al., 2022}).  Note that the variance
  (sill, \code{param["variance"]}) and the range parameters
  (\code{param["scale"]}) are assumed to be equal to 1 in the following
  formulae, and \eqn{x} is the lag distance.  The variogram functions are
  stationary and are valid for any number of dimensions if not mentioned
  otherwise.

  The following  intrinsic or stationary isotropic variogram
  functions \eqn{\gamma(x)} are implemented in \code{gencorr}:

  \itemize{

    \item \code{RMaskey}
    \deqn{
      \gamma(x)= 1 - (1-x)^\alpha 1_{[0,1]}(x)
    }
    \eqn{1_{[0,1]}(x)} is the indicator function equal to 1 for \eqn{x \in
    [0,1]} and 0 otherwise.  This variogram function is
    valid for dimension \eqn{d} if \eqn{\alpha \ge (d+1)/2}.  For
    \eqn{\alpha=1} we get the well-known triangle (or tent) model, which is
    only valid on the real line.

    \item \code{RMbessel}
    \deqn{
      \gamma(x) = 1 - 2^\nu \Gamma(\nu+1) x^{-\nu} J_\nu(x)
    }
    where \eqn{\nu \ge \frac{d-2}2}{\nu \ge (d-2)/2}, \eqn{\Gamma} denotes
    the gamma function and \eqn{J_\nu} is a Bessel function of first kind.
    This models a hole effect (see \cite{\enc{Chilès}{Chiles} and Delfiner, 1999,
    p. 92}).
    An important case is \eqn{\nu=-0.5} which gives the variogram
    function
    \deqn{\gamma(x)= 1 - \cos(x)}{\gamma(x)= 1 - cos(x)}
    and which is only valid for \eqn{d=1} (this equals \code{RMdampedcos}
    for \eqn{\lambda = 0}).
    A second important case is \eqn{\nu=0.5} with variogram function
        \deqn{
      \gamma(x) = \left(1 - \frac{\sin(x)}{x}\right) 1_{x>0}
    }{\gamma(x) = (1 - sin(x)/x)  1_{x>0}}
    which is valid for \eqn{d \le 3}. This coincides with \code{RMwave}.

    \item \code{RMcauchy}
    \deqn{\gamma(x) = 1 - (1 + x^2)^{-\gamma}}

    where \eqn{\gamma > 0}.  The parameter \eqn{\gamma} determines the
    asymptotic power law.  The smaller \eqn{\gamma}, the longer the
    long-range dependence.  The generalized Cauchy Family
    (\code{RMgencauchy}) includes this family for the choice \eqn{\alpha =
    2} and \eqn{\beta = 2 \gamma}.

    \item \code{RMcircular}
    \deqn{
      \gamma(x) = 1 - \left(1 -\frac{2}{\pi} \left(x \sqrt{1-x^2} + \arcsin(x)\right)\right) 1_{[0,1]}(x)
    }{\gamma(x) = 1 - (1 -2/\pi (x \sqrt{(1-x^2)} + asin(x))) 1_{[0,1]}(x)}
    This variogram function is valid only for dimensions \eqn{d \le 2}.

    \item \code{RMcubic}
    \deqn{
      \gamma(x) = 1 - (1-7 x^2 + 8.75 x^3 - 3.5 x^5 + 0.75 x^7) 1_{[0,1]}(x)
    }
    The model is only valid for dimensions \eqn{d \le 3}.  It is a 2 times
    differentiable variogram function with compact support (see
    \cite{\enc{Chilès}{Chiles} and Delfiner, 1999, p. 84}).

    \item \code{RMdagum}
    \deqn{
      \gamma(x) = (1+x^{-\beta})^{-\gamma / \beta}
    }{\gamma(x) = (1+x^(-\beta))^(-\gamma / \beta)}
    The parameters \eqn{\beta} and \eqn{\gamma} can be varied in the
    intervals \eqn{(0,1]} and \eqn{(0,1)}, respectively.  Like the
    generalized Cauchy model (\code{RMgencauchy}) the Dagum family can be
    used to model separately fractal dimension and Hurst effect
    (see \cite{Berg et al., 2008}).

    \item \code{RMdampedcos}
    \deqn{
      \gamma(x) = 1 - \exp(-\lambda x) \cos(x)
    }{\gamma(x) = 1 - exp(-\lambda x) cos(x)}
    The model is valid for any dimension \eqn{d}.  However, depending on
    the dimension of the random field the following bound
    \eqn{
      \lambda \ge 1/{\tan(\pi/(2d))}
    }{\lambda \ge 1/{tan(\pi/(2d))}} has to be respected.
    This variogram function models a hole effect
    (see \cite{\enc{Chilès}{Chiles} and Delfiner, 1999, p. 92}).
    For \eqn{\lambda = 0} we obtain
    \deqn{\gamma(x)= 1 - \cos(x)}{\gamma(x)= 1 - cos(x)} which is only valid
    for \eqn{d=1} and corresponds to \code{RMbessel} for \eqn{\nu=-0.5}.

    \item \code{RMdewijsian}
    \deqn{
      \gamma(x) = \log(1 + x^{\alpha})
    }{\gamma(x)=log(1 + x^{\alpha}+1)}
    where \eqn{\alpha \in (0,2]}.  This is an intrinsic variogram function.
    Originally, the logarithmic model \eqn{\gamma(x)
    = \log(x)}{\gamma(x) = log(x)} was named after de Wijs and
    reflects a principle of
    similarity (see  \cite{\enc{Chilès}{Chiles} and Delfiner, 1999, p. 90}).  But
    note that \eqn{\gamma(x) = \log(x)}{\gamma(x) = log(x)} is not
    a valid variogram function.

    \item \code{RMexp}
    \deqn{\gamma(x) = 1 - e^{-x}}{\gamma(x) = 1 - exp(-x)}
    This model is a special case of the Whittle model
    (\code{RMwhittle}) if \eqn{\nu=0.5}
    and of the stable family (\code{RMstable})
    if \eqn{\nu = 1}. Moreover, it is the continuous-time analogue
    of the first order auto-regressive time series covariance structure.

    \item \code{RMfbm}
    \deqn{\gamma(x) = x^\alpha}
    where \eqn{\alpha \in (0,2)}{0 < \alpha < 2}.  This is an
    intrinsically stationary variogram function.  For \eqn{\alpha=1}
    we get a variogram function corresponding to a standard
    Brownian Motion.  For \eqn{\alpha \in (0,2)}{0 < \alpha < 2} the
    quantity \eqn{H = \frac{\alpha}{2}}{H=\alpha/2} is called Hurst index
    and determines the fractal dimension \eqn{D = d + 1 - H} of the corresponding
    Gaussian sample paths where \eqn{d} is the
    dimension of the random field
    (see \cite{\enc{Chilès}{Chiles} and Delfiner, 1999, p. 89}).

    \item \code{RMgauss}
    \deqn{\gamma(x) = 1 - e^{-x^2}}{\gamma(x)= 1 - exp(-x^2)}

    The Gaussian model has an infinitely differentiable variogram
    function.  This smoothness is artificial.  Furthermore, this often
    leads to singular matrices and therefore numerically instable
    procedures (see \cite{Stein, 1999, p. 29}).  The Gaussian model is included in
    the stable class (\code{RMstable}) for the choice \eqn{\alpha = 2}.

    \item \code{RMgencauchy}
    \deqn{\gamma(x) = 1 - (1 + x^\alpha)^{-\beta/\alpha}}
    where \eqn{\alpha \in (0,2]}{0 < \alpha \le 2} and \eqn{\beta >
    0}.  This model has a smoothness parameter \eqn{\alpha} and
    a parameter \eqn{\beta} which determines the asymptotic power law.
    More precisely, this model admits simulating random fields where
    fractal dimension \emph{D} of the Gaussian sample path and Hurst
    coefficient \emph{H} can be chosen independently (compare also with
    \code{RMlgd}): Here, we have \eqn{D = d + 1 - \alpha/2, \alpha \in
    (0,2]}{ D = d + 1 - \alpha/2, 0 < \alpha \le 2} and \eqn{ H = 1 -
    \beta/2, \beta > 0}.  The smaller \eqn{\beta}, the longer the
    long-range dependence.  The variogram function is very regular
    near the origin, because its Taylor expansion only contains even terms
    and reaches its sill slowly.  Note that the Cauchy Family
    (\code{RMcauchy}) is included in this family for the choice \eqn{\alpha
    = 2} and \eqn{\beta = 2 \gamma}.

    \item \code{RMgenfbm}
    \deqn{
      \gamma(x) = (1 + x^{\alpha})^{\delta/\alpha} - 1
    }
    where \eqn{\alpha \in (0,2)}{0 < \alpha < 2} and \eqn{\delta \in
    (0,1)}. This is an intrinsic variogram function.

    \item \code{RMgengneiting}
    This is a family of stationary models whose elements are specified by
    the two parameters \eqn{\kappa} and \eqn{\mu} with \eqn{\kappa} being a
    non-negative integer and \eqn{\mu \ge \frac{d}{2}}{\mu \ge d/2} with
    \eqn{d} denoting the dimension of the random field (the models can be
    used for any dimension).  Let \eqn{\beta = \mu + 2\kappa +1/2}.

    For \eqn{\kappa = 0} the model equals the Askey model (\code{RMaskey})
    and is therefore not implemented.

    For \eqn{\kappa = 1} the model is given by

    \deqn{
      \gamma(x) = 1 - \left(1+\beta x \right)(1-x)^{\beta} 1_{[0,1]}(x),
        \qquad \beta = \mu +2.5,
    }{
      \gamma(x) = 1 - (1+\beta x)(1-x)^\beta 1_{[0,1]}(x), \beta = \mu + 2.5,}

    If \eqn{\kappa = 2}
    \deqn{
      \gamma(x) = 1 - \left(1 + \beta x + \frac{\beta^{2} - 1}{3}
        x^{2} \right)(1-x)^{\beta} 1_{[0,1]}(x), \qquad
        \beta = \mu+4.5,
    }{
      \gamma(x) = 1 -(1 + \beta x + (\beta^2 - 1) x^2 / 3)
        (1 - x)^\beta 1_{[0,1]}(x), \beta = \mu + 4.5,
    }

    and for \eqn{\kappa = 3}
    \deqn{
      \gamma(x) = 1 - \left( 1 + \beta x + \frac{(2\beta^{2}-3)}{5} x^{2}+
        \frac{(\beta^2 - 4)\beta}{15} x^{3} \right)(1-x)^\beta 1_{[0,1]}(x),
        \beta = \mu+6.5,
    }{
      \gamma(x) = 1 - (1 + \beta x + (2 \beta^2 - 3)x^2 / 5 +
        (\beta^2 - 4) \beta x^3 / 15)(1-x)^\beta 1_{[0,1]}(x),
        \beta = \mu + 6.5
    }

    A special case of this family is \code{RMgneiting} (with
    \eqn{s = 1} there) for the choice \eqn{\kappa = 3, \mu =
    3/2}.

    \item \code{RMgneiting}

    \deqn{
      \gamma(x) = 1 - (1 + 8 s x + 25 s^2 x^2 + 32 s^3 x^3)(1-s x)^8
    }{\gamma(x) = 1 - (1 + 8 s x + 25 s^2 x^2 + 32 s^3 x^3)(1-s x)^8 }

    if \eqn{0 \le x \le \frac{1}{s}}{0 <= x <= 1/s} and
    \deqn{\gamma(x)= 1} otherwise.  Here,
    \eqn{s=0.301187465825}.   This variogram function is
    valid only for dimensions less than or equal to 3.  It is 6 times
    differentiable and has compact support.  This model is an alternative
    to \code{RMgauss} as its graph is hardly distinguishable from the graph
    of the Gaussian model, but possesses neither the mathematical nor the
    numerical disadvantages of the Gaussian model.  It is a special case of
    \code{RMgengneiting} for the choice \eqn{\kappa=3, \mu=1.5}.

    \item \code{RMlgd}
    \deqn{
      \gamma(x) = \frac{\beta}{\alpha + \beta} x^{\alpha} 1_{[0,1]}(x) +
        (1 - \frac{\alpha}{\alpha + \beta} x^{-\beta}) 1_{x>1}(x)
    }{
      \gamma(x) = \beta/(\alpha + \beta) x^(\alpha) 1_{[0,1]}(x) +
        (1 - \alpha / (\alpha + \beta)x^(-\beta)) 1_{x>1}(r)
    }
    where \eqn{\beta >0} and \eqn{0 < \alpha \le (3-d)/2}{0 < \alpha \le
    (3-d)/2}, with \eqn{d}{d} denoting the dimension of the random field.
    The model is only valid for dimension \eqn{d=1,2}.  This model admits
    simulating random fields where fractal dimension \eqn{D} of the
    Gaussian sample and Hurst coefficient \eqn{H} can be chosen
    independently (compare also \code{RMgencauchy}): Here, the random field
    has fractal dimension \eqn{D = d+1 - \alpha/2} and Hurst coefficient
    \eqn{H = 1-\beta/2} for \eqn{0< \beta \le 1}.

    \item \code{RMmatern}
    \deqn{
      \gamma(x) = 1 - \frac{2^{1-\nu}}{\Gamma(\nu)} (\sqrt{2\nu}x)^\nu
        K_\nu(\sqrt{2\nu}x)
    }{
      \gamma(x) = 1 - 2^{1- \nu} \Gamma(\nu)^{-1} (\sqrt{2\nu} x)^\nu
        K_\nu(\sqrt{2\nu} x)
    }
    where \eqn{\nu > 0} and \eqn{K_\nu} is the modified Bessel function of
    second kind.  This is one of 3 possible parametrizations (Whittle,
    \enc{Matérn}{Matern}, Handcock-Wallis) of the Whittle-\enc{Matérn}{Matern} model.  The
    Whittle-\enc{Matérn}{Matern} model is the model of choice if the smoothness of a
    random field is to be parametrized: the sample paths of a Gaussian
    random field with this covariance structure are \eqn{m} times
    differentiable if and only if \eqn{\nu > m} (see \cite{Gneiting and Guttorp,
    2010, p.  24}).  Furthermore, the fractal dimension \eqn{D} of the
    Gaussian sample paths is determined by \eqn{\nu}: We have \eqn{D = d +
    1 - \nu, \nu \in (0,1)}{D = d + 1 - \nu, 0 < \nu < 1} and \eqn{D = d}
    for \eqn{\nu > 1}{\nu > 1} where \eqn{d} is the dimension of the random
    field (see \cite{Stein, 1999, p. 32}).  If \eqn{\nu=0.5}{\nu=0.5} the \enc{Matérn}{Matern}
    model equals \code{RMexp}.  For \eqn{\nu} tending to \eqn{\infty} a
    rescaled Gaussian model \code{RMgauss} appears as limit for the
    Handcock-Wallis parametrization.

    \item \code{RMpenta}
    \deqn{
      \gamma(x) = 1 - \left(1 - \frac{22}{3}x^{2} + 33 x^{4} - \frac{77}{2} x^{5} +
        \frac{33}{2} x^{7} - \frac{11}{2} x^{9} + \frac{5}{6}x^{11}\right) 1_{[0,1]}(x)
    }{\gamma(x)= 1 - (1 - 22/3 x^{2} + 33 x^{4} - 77/2 x^{5} + 33/2 x^{7} - 11/2 x^{9} + 5/6 x^{11}) 1_{[0,1]}(x)}
    The model is only valid for dimensions \eqn{d \le 3}{d \le 3}.  It has
    a 4 times differentiable variogram function with compact support (cf.
    \cite{\enc{Chilès}{Chiles} and Delfiner, 1999, p. 84}).

    \item \code{RMqexp}
    \deqn{\gamma(x)= 1 - \frac{2 e^{-x} - \alpha e^{-2x}}{ 2 - \alpha }}
    where \eqn{\alpha \in [0,1]}{0 \le \alpha \le 1}.

    \item \code{RMspheric}
    \deqn{
      \gamma(x) = 1 - \left(1 - \frac{3}{2} x + \frac{1}{2} x^3\right)  1_{[0,1]}(x)
    }{\gamma(x) = 1 - (1 - 1.5 x + 0.5 x^3) 1_{[0,1]}(x)}
    This variogram model is valid only for dimensions less than or equal
    to 3 and has compact support.

    \item \code{RMstable}
    \deqn{\gamma(x) = 1 - e^{-x^\alpha}}{\gamma(x)= 1 - e^{-x^\alpha}}
    where \eqn{\alpha \in (0,2]}{0 < \alpha \le 2}.  The parameter
    \eqn{\alpha} determines the fractal dimension \eqn{D} of the Gaussian
    sample paths: \eqn{D = d + 1 - \alpha/2} where \eqn{d} is the dimension
    of the random field.  For \eqn{\alpha < 2} the Gaussian sample paths
    are not differentiable (see \cite{Gelfand et al., 2010, p.  25}).  The stable
    family includes the exponential model (\code{RMexp}) for \eqn{\alpha =
    1} and the Gaussian model ( \code{RMgauss}) for \eqn{\alpha = 2}.

    \item \code{RMwave}

    \deqn{
      \gamma(x) = \left(1 - \frac{\sin(x)}{x}\right) 1_{x>0}
    }{\gamma(x) = (1 - sin(x)/x)  1_{x>0}}

    The model is only valid for dimensions \eqn{d \le 3}.  It is a special
    case of \code{RMbessel} for \eqn{\nu = 0.5}.  This variogram models a
    hole effect (see \cite{\enc{Chilès}{Chiles} and Delfiner, 1999, p. 92}).

    \item \code{RMwhittle}

    \deqn{
      \gamma(x)=1 - \frac{2^{1- \nu}}{\Gamma(\nu)} x^{\nu}K_{\nu}(x)
    }{\gamma(x)=1 - 2^{1- \nu}  \Gamma(\nu)^{-1}x^{\nu}K_{\nu}(x)}

    where \eqn{\nu > 0} and \eqn{K_\nu} is the modified Bessel function of
    second kind.  This is one of 3 possible parametrizations (Whittle,
    \enc{Matérn}{Matern}, Handcock-Wallis) of the
    Whittle-\enc{Matérn}{Matern} model, for further details, see
    information for entry \code{RMmatern} above.

  }

}

\value{
A numeric vector with generalized correlations (= negative semi-variances
computed with variance parameter \code{param["variance"] = 1}).
}

\references{

  Berg, C., Mateau, J., Porcu, E. (2008) The dagum family of isotropic
  correlation functions, \emph{Bernoulli}, \bold{14}, 1134--1149,
  \doi{10.3150/08-BEJ139}.

  \enc{Chilès}{Chiles}, J.-P., Delfiner, P. (1999) \emph{Geostatistics
  Modeling Spatial Uncertainty}, Wiley, New York,
  \doi{10.1002/9780470316993}.

  Gneiting, T. (2002) Compactly supported correlation functions.
  \emph{Journal of Multivariate Analysis}, \bold{83}, 493--508,
  \doi{10.1006/jmva.2001.2056}.

  Gneiting, T., Schlather, M. (2004) Stochastic models which separate
  fractal dimension and Hurst effect.  \emph{SIAM review} \bold{46},
  269--282, \doi{10.1137/S0036144501394387}.

  Gneiting, T., Guttorp, P. (2010) Continuous Parameter Stochastic Process
  Theory, In Gelfand, A. E., Diggle, P. J., Fuentes, M., Guttrop, P. (Eds.)
  \emph{Handbook of Spatial Statistics}, CRC Press, Boca Raton, p. 17--28,
  \doi{10.1201/9781420072884}.

  Schlather M., Malinowski A., Oesting M., Boecker D., Strokorb K., Engelke
  S., Martini J., Ballani F., Moreva O., Auel J., Menck P.J., Gross S.,
  Ober U., Ribeiro P., Ripley B.D., Singleton R., Pfaff B., R Core Team
  (2022).  \emph{RandomFields: Simulation and Analysis of Random Fields}.
  \R package version 3.3.14,
  \url{https://cran.r-project.org/src/contrib/Archive/RandomFields/}.

  Stein, M. L. (1999) \emph{Interpolation of Spatial Data: Some Theory for
  Kriging}, Springer, New York, \doi{10.1007/978-1-4612-1494-6}.

}

\author{
  Andreas Papritz \email{papritz@retired.ethz.ch}}

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

  \code{\link{predict.georob}} for computing robust Kriging predictions;

%   \code{\link{validate.predictions}} for validating Kriging predictions;

  \code{\link{lgnpp}} for unbiased back-transformation of Kriging prediction
  of log-transformed data;

  \code{\link{georobSimulation}} for simulating realizations of a Gaussian process
  from model fitted by \code{georob}; and finally

  \code{\link{sample.variogram}} and \code{\link{fit.variogram.model}}
  for robust estimation and modelling of sample variograms.
}

\examples{
## scaled lag distances
x <- seq(0, 3, length.out = 100)

## generalized correlations stable model
y <- gencorr(x, variogram.model = "RMstable", param = c(alpha = 1.5))
plot(x, y)

## generalized correlations circular model
y <- gencorr(x, variogram.model = "RMcircular")
plot(x, y)
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
