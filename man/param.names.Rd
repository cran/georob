% 2012-12-13 A. Papritz
% R CMD Rdconv -t html -o bla.html param.names.Rd ; open bla.html; R CMD Rd2pdf --force param.names.Rd; 

\encoding{macintosh}
\name{param.names}
\alias{param.names}
\alias{param.bounds}

\title{Names and Permissible Ranges of Variogram Parameters}

\description{
  Helper functions to query names and permissble ranges of variogram
  parameters.}

\usage{
param.names(model)

param.bounds(model, d, param)
}

\arguments{
  \item{model}{character keyword denoting a valid variogram,
  see \code{\link{georob}} and \code{\link{georobIntro}}.} 
  
  \item{d}{integer equal number of dimensions of the survey domain.}
  
  \item{param}{named numeric vector with variogram parameters, see
  \code{\link{georob}}.} 
}

\value{
  Either a character vector with the names of the additional variogram
  parameters such as the smoothness parameter of the Whittle-Matérn model
  (\code{param.names}) or a named list with the lower and upper bounds of
  permissble parameter ranges.  
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models.
}

\examples{
param.names("gengneiting")
param.bounds("gengneiting", d = 2, param = c(n = 2))
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
