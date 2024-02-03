% 2024-02-01 A. Papritz
% R CMD Rdconv -t html -o bla.html reexportedfun.Rd ; open bla.html;
% R CMD Rd2pdf --force reexportedfun.Rd; open reexportedfun.pdf
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  reexportedfun.Rd

\name{reexported-fun}
\alias{K}
\alias{lmrob.control}
% \alias{waldtest}

\title{Re-Exported Functions from \R package imports}

\description{The imported functions \code{\link[constrainedKriging]{K}},
\code{\link[robustbase]{lmrob.control}},
and \code{\link[lmtest]{waldtest}} are re-exported for ease of use without
attaching the respective packages.  }

\usage{
K(dist, model)

lmrob.control(setting, seed = NULL, nResample = 500, tuning.chi = NULL,
  bb = 0.5, tuning.psi = NULL, max.it = 50, groups = 5, n.group = 400,
  k.fast.s = 1, best.r.s = 2, k.max = 200, maxit.scale = 200, k.m_s = 20,
  refine.tol = 1e-7, rel.tol = 1e-7, scale.tol = 1e-10, solve.tol = 1e-7,
  zero.tol = 1e-10, trace.lev = 0, mts = 1000,
  subsampling = c("nonsingular", "simple"), compute.rd = FALSE,
  method = "MM", psi = "bisquare", numpoints = 10, cov = NULL,
  split.type = c("f", "fi", "fii"), fast.s.large.n = 2000,
  # only for outlierStats() :
  eps.outlier = function(nobs) 0.1 / nobs,
  eps.x = function(maxx) .Machine$double.eps^(.75)*maxx,
  compute.outlier.stats = method, warn.limit.reject = 0.5,
  warn.limit.meanrw = 0.5, \dots)

}

\arguments{

  \item{dist}{a numeric vector with distances.}

  \item{model}{an object of class \dQuote{\code{covmodel}} that defines an
  isotropic covariance model, see
  \code{\link[constrainedKriging]{covmodel}}.}

  \item{setting}{a string specifying alternative default values, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{seed}{\code{NULL} or an integer vector compatible with
  \code{\link{.Random.seed}}, see \code{\link[robustbase]{lmrob.control}}.
  }

  \item{nResample}{number of re-sampling candidates to be used to find the
  initial S-estimator, see \code{\link[robustbase]{lmrob.control}}.}

  \item{tuning.chi}{tuning constant vector for the S-estimator, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{bb}{expected value under the normal model of the \dQuote{chi}, see
  \code{\link[robustbase]{lmrob.control}}}

  \item{tuning.psi}{tuning constant vector for the redescending
  M-estimator, see \code{\link[robustbase]{lmrob.control}}.}

  \item{max.it}{integer specifying the maximum number of IRWLS iterations,
  see \code{\link[robustbase]{lmrob.control}}.}

  \item{groups}{(for the fast-S algorithm): Number of random subsets to use
  when the data set is large, see \code{\link[robustbase]{lmrob.control}}.}

  \item{n.group}{(for the fast-S algorithm): Size of each of the
  \code{groups} above, see \code{\link[robustbase]{lmrob.control}}.}

  \item{k.fast.s}{(for the fast-S algorithm): Number of local improvement
  steps (\dQuote{\emph{I-steps}}) for each re-sampling candidate, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{best.r.s}{(for the fast-S algorithm): Number of of best candidates
  to be iterated further, see \code{\link[robustbase]{lmrob.control}}.}

  \item{k.max}{(for the fast-S algorithm): maximal number of refinement
  steps for the \dQuote{fully} iterated best candidates, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{maxit.scale}{integer specifying the maximum number of C level
  \code{find_scale()} iterations (in fast-S and M-S algorithms), see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{k.m_s}{(for the M-S algorithm): specifies after how many
  unsuccessful refinement steps the algorithm stops, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{refine.tol}{(for the fast-S algorithm): relative convergence
  tolerance for the fully iterated best candidates, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{rel.tol}{(for the RWLS iterations of the MM algorithm): relative
  convergence tolerance for the parameter vector, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{scale.tol}{(for the scale estimation iterations of the S
  algorithm): relative convergence tolerance for the \code{scale}
  \eqn{\sigma(.)}, see \code{\link[robustbase]{lmrob.control}}, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{solve.tol}{(for the S algorithm): relative tolerance for inversionsee
  \code{\link[robustbase]{lmrob.control}}.}

  \item{zero.tol}{for checking 0-residuals in the S algorithm, non-negative
  number, see \code{\link[robustbase]{lmrob.control}}.}

  \item{trace.lev}{integer indicating if the progress of the MM-algorithm
  and the fast-S algorithms should be traced, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{mts}{maximum number of samples to try in subsampling algorithm, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{subsampling}{type of subsampling to be used, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{compute.rd}{a logical scalar indicating if robust distances (based on the
  MCD robust covariance estimator) are to be computed for the robust
  diagnostic plots, see \code{\link[robustbase]{lmrob.control}}.}

  \item{method}{string specifying the estimator-chain, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{psi}{string specifying the type \eqn{\psi}-function used., see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{numpoints}{number of points used in Gauss quadrature, see
  \code{\link[robustbase]{lmrob.control}}}

  \item{cov}{function or string with function name to be used to calculate
  covariance matrix estimate, see \code{\link[robustbase]{lmrob.control}}.}

  \item{split.type}{determines how categorical and continuous variables are
  split, see \code{\link[robustbase]{lmrob.control}}.}

  \item{fast.s.large.n}{minimum number of observations required to switch
  from ordinary \dQuote{fast S} algorithm to an efficient \dQuote{large n}
  strategy, see \code{\link[robustbase]{lmrob.control}}}

  \item{eps.outlier}{limit on the robustness weight below which an
  observation is considered to be an outlier, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{eps.x}{limit on the absolute value of the elements of the design
  matrix below which an element is considered zero, see
  \code{\link[robustbase]{lmrob.control}}.}

  \item{compute.outlier.stats}{vector of character strings, each valid to
  be used as \code{method} argument, see
  \code{\link[robustbase]{lmrob.control}}}

  \item{warn.limit.reject}{see \code{\link[robustbase]{lmrob.control}}.}

  \item{warn.limit.meanrw}{limit of the mean robustness per factor level
  below which (\eqn{\leq}{<=}) a warning is produced.  Set to \code{NULL}
  to disable warning.}

  \item{\dots}{some methods for the generic function
  \code{\link[lmtest]{waldtest}} require additional arguments, see
  respective help pages.}

}

\details{
  The function \code{\link[constrainedKriging]{K}} is required for
  computing block Kriging predictions by the function
  \code{\link[constrainedKriging]{f.point.block.cov}} of the package
  \pkg{constrainedKriging}.

  Furthermore, the function \code{\link[robustbase]{lmrob.control}} allows
  to pass tuning parameters to the function \code{\link[robustbase]{lmrob}}
  of the package \pkg{robustbase}, which is used for computing robust
  initial values of the regression coefficients.  }

\value{
  See help pages of \code{\link[constrainedKriging]{K}} and
  \code{\link[robustbase]{lmrob.control}} for the output generated by these
  functions.
}

\keyword{models}
\keyword{robust}
\keyword{regression}
