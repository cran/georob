% 2015-06-30 A. Papritz
% R CMD Rdconv -t html -o bla.html pmm.Rd ; open bla.html; R CMD Rd2pdf --force pmm.Rd; 

\encoding{macintosh}
\name{pmm}
\alias{pmm}
\alias{control.pcmp}

\title{Parallelized Matrix Multiplication}

\description{
  This page documents the function \code{pmm} for parallelized matrix
  multiplication and the function \code{control.pcmp}, which controls
  the behaviour of \code{pmm} and other functions that execute tasks in
  parallel.}
  
\usage{
pmm(A, B, control = control.pcmp())

control.pcmp(pmm.ncores = 1L, gcr.ncores = 1L, max.ncores = detectCores(), 
    f = 2L, sfstop = FALSE, allow.recursive = TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A, B}{matrices to be multiplied.}

  \item{control}{a list of with the arguments \code{ncores}, \code{f}, and
  \code{sfstop} or a function such as \code{control.pcmp} that generates such
  a list.} 
  
  \item{pmm.ncores}{number (integer, default 1) of cores used for parallelized
  matrix multiplication.}
  
  \item{gcr.ncores}{number (integer, default 1) of cores used for parallelized
  computation of semi-variance matrix.}
  
  \item{max.ncores}{maximum number of cores (integer, default all cores of
  a machine) used for parallelized computations.}
  
  \item{f}{number (integer, default 2) of tasks assigned to each core in
  parallelized operations.}
  
  \item{sfstop}{logical controlling whether the SNOW socket cluster is
  stopped after each parallelized matrix multiplication on windows OS
  (default \code{FALSE}).}
  
  \item{allow.recursive}{logical controlling whether nested parallelized
  computation should be allowed (default \code{TRUE}).}
  
  \item{...}{further arguments, currently not used.}
  
}

\details{
Parallelized matrix multiplication shortens computing time for large data
sets (\eqn{n>1000}).  However, spawning child processes requires itself
resources and increasing the number of cores for parallel matrix
multiplication does not always result in reduced computing time.  A
sensible default for the number of cores is likely
\code{pmm.ncores=2}.

Note, however, that very substantial reductions in computing time
results when one uses the \strong{OpenBLAS} library instead of the
reference BLAS library that ships with R, see
\url{http://www.openblas.net/} and R FAQ for your OS. With OpenBLAS no
gains are obtained by using more than one core for matrix
multiplication, and one should therefore use the default argument
\code{pmm.ncores = 1} for \code{control.pcmp()}.

\code{max.ncores} controls how many child processes are spawned in total.
This can be used to prevent that child processes spawn
themselves children which may result in a considerable number of child
processes.

}

\value{
  \describe{
    \item{\code{pmm}:}{the matrix product \code{A \%*\% B},}
    
    \item{\code{control.pcmp}:}{a list with components
    \code{pmm.ncores}, \code{gcr.ncores}, \code{max.ncores}, \code{f},
    \code{sfstop}, \code{allow.recursive}.} } }
   
\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}
}

\examples{
\dontrun{
A <- as.matrix(dist(rnorm(2000)))
B <- as.matrix(dist(rnorm(2000)))
system.time(C <- pmm(A, B, control = control.pcmp(ncores = 1)))
system.time(C <- pmm(A, B, control = control.pcmp(ncores = 4)))
}}

\keyword{models}
\keyword{spatial}
\keyword{robust}
