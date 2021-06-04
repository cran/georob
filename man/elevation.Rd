% 2021-06-04 A. Papritz
% R CMD Rdconv -t html -o bla.html elevation.Rd ; open bla.html; R CMD Rd2pdf --force elevation.Rd;

%  hunspell -d en_GB -t -i mac  *.Rd

\encoding{utf8}

\name{elevation}
\alias{elevation}

\docType{data}

\title{Elevation Data}

\description{
  Surface elevation data taken from Davis (1972).
}

\usage{data(elevation)}


\format{
  A data frame with 52 observations on the following 3 variables:


  \describe{

    \item{\code{x}}{a numeric vector with the easting coordinate in
    multiplies of 50 feet.}

    \item{\code{y}}{a numeric vector with the northing coordinate in
    multiplies of 50 feet..}

    \item{\code{height}}{a numeric vector with the elevation in
    multiples of 10 feet.}

  }
}

\note{
  The data were imported from the package \pkg{geoR}.
}

\source{
  Davis, J.C. (1973) \emph{Statistics and Data Analysis in Geology}.  New
  York: Wiley.
}


\examples{
data(elevation)
summary(elevation)
}

\keyword{datasets}
