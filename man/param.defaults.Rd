% 2024-02-01 A. Papritz
% R CMD Rdconv -t html -o bla.html param.defaults.Rd ; open bla.html; R CMD Rd2pdf --force param.defaults.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  param.defaults.Rd

\encoding{utf8}
\name{default.aniso}
\alias{default.aniso}
\alias{default.fit.param}
\alias{default.fit.aniso}

\title{Setting Default Values of Variogram Parameters}

\description{
  Auxiliary functions to set sensible default values for anisotropy
  parameters and for controlling what variogram and anisotropy parameters
  should be estimated.}

\usage{

default.aniso(f1 = 1., f2 = 1., omega = 90., phi = 90., zeta = 0.)

default.fit.param(
  variance = TRUE, snugget = FALSE, nugget = TRUE, scale = TRUE,
  alpha = FALSE, beta = FALSE, delta = FALSE, gamma = FALSE,
  kappa = FALSE, lambda = FALSE, mu = FALSE, nu = FALSE)

default.fit.aniso(f1 = FALSE, f2 = FALSE, omega = FALSE,
  phi = FALSE, zeta = FALSE)
}

\arguments{
  \item{variance}{variance (sill \eqn{\sigma^2}) of the
  auto-correlated component of the Gaussian random field
  \eqn{B(\boldsymbol{s})}{B(s)}.}

  \item{snugget}{variance (spatial nugget
  \eqn{\sigma^2_{\mathrm{n}}}{\sigma^2_n}) of the seemingly spatially
  uncorrelated component of \eqn{B(\boldsymbol{s})}{B(s)}
  (micro-scale spatial variation; default value \code{snugget = 0}).}

  \item{nugget}{variance (nugget \eqn{\tau^2}) of the independent error
  \eqn{\varepsilon(\boldsymbol{s})}{\epsilon(s)}.}

  \item{scale}{range parameter (\eqn{\alpha}) of the variogram.}

  \item{alpha, beta, delta, gamma, kappa, lambda, mu, nu}{names of
  additional variogram parameters such as the smoothness parameter
  \eqn{\nu} of the Whittle-\enc{Matérn}{Matern} model (see
  \code{\link{gencorr}} and \code{\link{param.names}}).}

  \item{f1}{positive ratio \eqn{f_1} of lengths of second and first
  semi-principal axes of an ellipsoidal surface with constant semi-variance
  in \eqn{\mathrm{I}\!\mathrm{R}^3}{R^3} (default \code{f1 = 1}), see
  subsection \emph{Model} of \code{\link{georobPackage}}.}

  \item{f2}{positive ratio \eqn{f_2} of lengths of third and first
  semi-principal axes of the semi-variance ellipsoid (default \code{f2 =
  1}), see subsection \emph{Model} of \code{\link{georobPackage}}.}

  \item{omega}{azimuth in degrees of the first semi-principal axis of the
  semi-variance ellipsoid (default \code{omega = 90}), see subsection
  \emph{Model} of \code{\link{georobPackage}}.}

  \item{phi}{90 degrees minus altitude of the first semi-principal axis of
  the semi-variance ellipsoid (default \code{phi = 90}), see subsection
  \emph{Model} of \code{\link{georobPackage}}.}

  \item{zeta}{angle in degrees between the second semi-principal axis and
  the direction of the line defined by the intersection between the
  \eqn{x}-\eqn{y}-plane and the plane orthogonal to the first
  semi-principal axis of the semi-variance ellipsoid through the origin
  (default \code{zeta = 0}), see subsection \emph{Model} of
  \code{\link{georobPackage}}.}

}

\value{
  Either a named numeric vector with initial values of anisotropy
  parameters (\code{default.aniso}) or named logical vectors, controlling
  what parameters should be estimated (\code{default.fit.param},\cr
  \code{default.fit.aniso}).}

\author{
   Andreas Papritz \email{papritz@retired.ethz.ch}.
}

\seealso{
  \code{\link{georobPackage}} for a description of the model and a brief summary of the algorithms;

  \code{\link{georob}} for (robust) fitting of spatial linear models.
}

\examples{
default.aniso(f1 = 0.5, omega = 45)
default.fit.param(scale=FALSE, alpha = TRUE)
default.fit.aniso(f1 = TRUE, omega = TRUE)
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
