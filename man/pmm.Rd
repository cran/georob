% 2014-07-29 A. Papritz
% R CMD Rdconv -t html -o bla.html pmm.Rd ; open bla.html; R CMD Rd2pdf --force pmm.Rd; 

\encoding{macintosh}
\name{pmm}
\alias{pmm}
\alias{parallel.control}

\title{Parallelized Matrix Multiplication}

\description{
  This page documents the functions \code{pmm} for parallelized matrix
  multiplication and the function \code{parallel.control}, which controls
  the behaviour of \code{pmm} and other functions that execute task in
  parallel.}
  
\usage{
pmm(A, B, control = parallel.control())

parallel.control(pmm.ncores = 1, gradient.ncores = 1, max.ncores = detectCores(), 
    f = 2, sfstop = FALSE, allow.recursive = TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A, B}{matrices to be mulitplied.}

  \item{control}{a list of with the arguments \code{ncores}, \code{f}, and
  \code{sfstop} or a function such as \code{parallel.control} that generates such
  a list.} 
  
  \item{pmm.ncores}{number (integer, default 2) of cores used for parallelized
  matrix multiplication.}
  
  \item{gradient.ncores}{number (integer, default 2) of cores used for parallelized
  evaluation of estimating equations and gradient.}
  
  \item{max.ncores}{maximum number of cores (integer, default all cores of
  a machine) used for parallelized computations.}
  
  \item{f}{number (integer, default 2) of tasks assigned on non-Windows OS
  to each core in parallelized matrix multiplication.}
  
  \item{sfstop}{logical controlling whether the SNOW socket cluster is
  stopped after each parallelized matrix multiplication on windows OS
  (default \code{FALSE}).}
  
  \item{allow.recursive}{logical controlling whether nested parallelized
  computation should be allowed (default \code{TRUE}).}
  
  \item{...}{further arguments, currently not used.}
  
}


\value{
  \describe{
    \item{\code{pmm}:}{the matrix product \code{A \%*\% B},}
    
    \item{\code{parallel.control}:}{a list with components
    \code{pmm.ncores}, \code{gradient.ncores}, \code{max.ncores}, \code{f},
    \code{sfstop}, \code{allow.recursive}.} } }
   
\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}
}

\examples{
\dontrun{
A <- as.matrix(dist(rnorm(2000)))
B <- as.matrix(dist(rnorm(2000)))
system.time(C <- pmm(A, B, control = parallel.control(pmm.ncores = 1)))
system.time(C <- pmm(A, B, control = parallel.control(pmm.ncores = 4)))
}}

\keyword{models}
\keyword{spatial}
\keyword{robust}
