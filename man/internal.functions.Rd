% 2013-06-03 A. Papritz
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
The internal, unexported functions 
\code{compute.covariances},
\code{compute.estimating.equations},
\code{compute.semivariance},
\code{dcorr.dparam},
\code{estimate.xihat},
\code{gcr},
\code{georob.fit},
\code{getCall.georob},
\code{gradient.negative.restricted.loglikelihood},
\code{K},
\code{negative.restr.loglikelihood},
\code{prepare.likelihood.calculations},
\code{update.xihat}
are not intended for direct use.  However, as any unexported function they
can be accessed by typing
\code{georob:::\var{function-name}}.
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

