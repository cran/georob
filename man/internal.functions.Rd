% 2013-09-06 A. Papritz
% R CMD Rdconv -t html -o bla.html internal.functions.Rd ; open bla.html; R CMD Rd2pdf --force internal.functions.Rd 
\name{internal.functions}
\alias{compute.covariances}
\alias{compute.estimating.equations}
\alias{compute.semivariance}
\alias{dcorr.dparam}
\alias{estimate.xihat}
\alias{gcr}
\alias{georob.fit}
\alias{getCall.georob}
\alias{gradient.negative.restricted.loglikelihood}
\alias{K}
\alias{negative.restr.loglikelihood}
\alias{prepare.likelihood.calculations}
\alias{update.xihat}
%
\title{Internal Functions of Package \pkg{georob}}
%
\description{
The internal unexported functions 
\itemize{
  \item\code{compute.covariances},
  \item\code{compute.estimating.equations},
  \item\code{compute.semivariance},
  \item\code{dcorr.dparam},
  \item\code{estimate.xihat},
  \item\code{gcr},
  \item\code{georob.fit},
  \item\code{getCall.georob},
  \item\code{gradient.negative.restricted.loglikelihood},
  \item\code{negative.restr.loglikelihood},
  \item\code{prepare.likelihood.calculations},
  \item\code{update.xihat}
}
are not intended for direct use.  However, as any unexported function they
can be accessed by typing
\code{georob:::\var{function-name}}.  

In addition, the (exported) function \code{K} is required for computing
block kriging predictions by the (unexported) function
\code{f.point.block.cov} of the package
\pkg{constrainedKriging}.
}

\author{
  Andreas Papritz \email{andreas.papritz@env.ethz.ch}\cr
  \url{http://www.step.ethz.ch/people/scientific-staff/andreas-papritz}\cr
  with contributions by Cornelia Schwierz.  }

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

\keyword{internal}

