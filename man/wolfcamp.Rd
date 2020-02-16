% 2020-02-05 A. Papritz
% R CMD Rdconv -t html -o bla.html wolfcamp.Rd ; open bla.html; R CMD Rd2pdf --force wolfcamp.Rd; 

%  hunspell -d en_GB -t -i mac  *.Rd

\encoding{utf8}

\name{wolfcamp}
\alias{wolfcamp}

\docType{data}

\title{Wolfcamp Aquifer Data}

\description{
  Piezometric head measurements taken at the Wolfcamp Aquifer, Texas, USA.
  See \cite{Cressie (1993, p. 212--214)} for description of the scientific
  problem and the data.  Original data were converted to SI units:
  coordinates are given in kilometers and pressure heads in meters. 
}

\usage{data(wolfcamp)}


\format{
  A data frame with 85 observations on the following 3 variables:
  

  \describe{

    \item{\code{x}}{a numeric vector with the easting coordinate in
    kilometers.}

    \item{\code{y}}{a numeric vector with the northing coordinate in
    kilometers.}

    \item{\code{pressure}}{a numeric vector with the piezometric head in
    meters.}

  }
}

\note{
  The data were imported from the package \pkg{geoR}.
}

\source{
  Harper, W.V. and Furr, J.M. (1986) Geostatistical analysis of
  potentiometric data in the Wolfcamp Aquifer of the Palo Duro Basin,
  Texas.  \emph{Technical Report BMI/ONWI-587, Bettelle Memorial Institute,
  Columbus, OH}. 
}

\references{
  Cressie, N.A.C (1993) \emph{Statistics for Spatial Data}. New York: Wiley.
  
  Papritz, A. and Moyeed, R. (2001) Parameter uncertainty in spatial
  prediction: checking its importance by cross-validating the Wolfcamp
  and Rongelap data sets. \emph{GeoENV 2000: Geostatistical for
    Environmental Applications. Ed. P. Monestiez, D. Allard,
    R. Froidevaux}. Kluwer. 
}

\examples{
data(wolfcamp)
summary(wolfcamp)
}

\keyword{datasets}
