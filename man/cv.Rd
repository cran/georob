% 2021-06-04 A. Papritz
% R CMD Rdconv -t html -o bla.html cv.Rd ; open bla.html; R CMD Rd2pdf --force cv.Rd;

\encoding{utf8}
\name{cv}
\alias{cv}

\title{Generic Cross-validation}

\description{Generic function for cross-validating models.}

\usage{
cv(object, ...)
}

\arguments{

  \item{object}{any model object.}

  \item{\dots}{additional arguments as required by the methods.}
}

\value{
  will depend on the method function used; see the respective documentation.}

\author{
   Andreas Papritz \email{papritz@retired.ethz.ch}.
}

\seealso{
  \code{\link{georob}} for (robust) fitting of spatial linear models;

  \code{\link{cv.georob}} for assessing the goodness of a model fitted by
  \code{georob}.  }


\keyword{models}
\keyword{spatial}
