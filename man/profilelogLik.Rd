% 2015-11-27 A. Papritz
% R CMD Rdconv -t html -o bla.html profilelogLik.Rd ; open bla.html; R CMD Rd2pdf --force profilelogLik.Rd; 

\encoding{macintosh}
\name{profilelogLik}
\alias{profilelogLik}

\title{Profile Likelihood}
\description{
  The function \code{profilelogLik} computes for an array of fixed
  variogram parameters the profile log-likelihood by maximizing the
  (restricted) log-likelihood with respect to the remaining variogram
  parameters, the fixed and random effects.  }

\usage{
  profilelogLik(object, values, use.fitted = TRUE, verbose = 0, 
    ncores = min(detectCores(), NROW(values)))
}

\arguments{
  \item{object}{an object of class \code{"georob"} (mandatory argument),
  see \code{\link{georobObject}}.}
  
   \item{values}{a \code{data.frame} or a \code{matrix} with the values of
   the variogram and anisotropy parameters that should be kept fixed
   (mandatory argument, see \code{\link{georob}} and
   \code{\link{georobIntro}} for information about the parametrization of
   variogram models).  The names of the columns of \code{values} must match
   the names of variogram and anisotropy parameters.}
  
  \item{use.fitted}{logical scalar controlling whether the fitted variogram
  parameters of \code{object} should be used as initial values (default
  \code{TRUE}) when maximizing the profile log-likelihood or the initial
  values used to generate \code{object}.}
  
  \item{verbose}{positive integer controlling logging of diagnostic
  messages to the console during model fitting, see \code{\link{georob}}.}
  
  \item{ncores}{positive integer controlling how many cores are used for
  parallelized computations, see \emph{Details}.}
}

\details{

  For robust REML fits \code{profilelogLik} returns (possibly with a
  warning) the log-likelihood of the Gaussian (RE)ML fit of the equivalent
  Gaussian spatial linear model with heteroscedastic nugget.  
  
  \bold{Note} that \emph{the data frame passed as} \code{data}
  \emph{argument to} \code{georob} \emph{must exist in the user workspace
  when calling \code{profilelogLik}}.
  
  \code{profilelogLik} uses the package \pkg{parallel} for
  parallelized computation of the profile likelihood.  By default, the
  function uses \code{NROW(values)} CPUs but not more than are physically
  available (as returned by \code{\link[parallel]{detectCores}}).
  
  \code{profilelogLik} uses the function \code{\link[stats]{update}} to
  re-estimated the model with partly fixed variogram parameters.
  Therefore, any argument accepted by \code{\link{georob}} except
  \code{data} can be changed when re-fitting the model.  Some of them (e.g.
  \code{verbose}) are explicit arguments of
  \code{profilelogLik}, but also the remaining ones can be passed by
  \code{\ldots} to the function.

}

\value{
  A \code{data.frame} with the columns of \code{values}, a column
  \code{loglik} (contains the maximized [restricted] log-likelihood),
  columns with the estimated variogram and fixed effect parameters, columns
  with the gradient of the (restricted) log-likelihood (or the roots of the
  estimating equations) and a column \code{converged}, indicating whether
  convergence has occurred \code{converged ==1} when fitting the respective
  model.  }

\author{
  Andreas Papritz \email{andreas.papritz@env.ethz.ch}
}

\seealso{
  \code{\link{georob}} for (robust) fitting of spatial linear models.
}

\examples{
\dontrun{

data(meuse)

r.logzn.ml <- georob(log(zinc)~sqrt(dist), data=meuse, locations=~x+y,
    variogram.model="RMexp", param=c(variance=0.15, nugget=0.05, scale=200),
    tuning.psi=1000, control=control.georob(ml.method="ML"))
  
r.prflik <- profilelogLik(r.logzn.ml, values=expand.grid(scale=seq(75, 600, by=25)))
plot(loglik~scale, r.prflik, type="l")
abline(v=r.logzn.ml$variogram.object[[1]]$param["scale"], lty="dotted")
abline(h=r.logzn.ml$loglik-0.5*qchisq(0.95, 1), lty="dotted")
  
}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
