% 2024-03-09 A. Papritz
% R CMD Rdconv -t html -o bla.html pmm.Rd ; open bla.html; R CMD Rd2pdf --force pmm.Rd;

\encoding{utf8}
\name{pmm}
\alias{pmm}
\alias{control.pcmp}

\title{Parallelized Matrix Multiplication}

\description{
  This page documents the function \code{pmm} for parallelized matrix
  multiplication and the function\cr \code{control.pcmp}, which controls
  the behaviour of \code{pmm} and other functions that execute tasks in
  parallel.}

\usage{
pmm(A, B, control = control.pcmp())

control.pcmp(pmm.ncores = 1, gcr.ncores = 1, max.ncores = parallel::detectCores(),
    f = 1, sfstop = FALSE, allow.recursive = FALSE,
    fork = !identical(.Platform[["OS.type"]], "windows"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A, B}{two numeric matrices to be multiplied.}

  \item{control}{a list with the arguments \code{pmm.ncores},
  \code{gcr.ncores}, \code{max.ncores}, \code{f}, \code{sfstop} and
  \code{allow.recursive} or a function such as \code{control.pcmp} that
  generates such a list.}

  \item{pmm.ncores}{an integer (default 1) with the number of cores used
  for parallelized matrix multiplication.}

  \item{gcr.ncores}{an integer (default 1) with the number of cores used
  for parallelized computation of (generalized) covariances or
  semi-variances.}

  \item{max.ncores}{maximum number of cores (integer, default all cores of
  a machine) used for parallelized computations.}

  \item{f}{an integer (default 1) with the number of tasks assigned to each
  core in parallelized operations.}

  \item{sfstop}{a logical scalar controlling whether the SNOW socket
  cluster is stopped after each parallelized matrix multiplication on
  windows OS (default \code{FALSE}).}

  \item{allow.recursive}{a logical scalar controlling whether parallelized
  matrix multiplicaction and computation of generalized) covariances should
  be allowed by child processes running already in parallel (default
  \code{FALSE}).}

  \item{fork}{a logical scalar controlling whether forking should be used for
  parallel computations (default \code{TRUE} on unix and \code{FALSE} on
  windows operating systems).  Note that settting \code{fork = TRUE} on
  windows suppresses parallel computations.}

  \item{...}{further arguments, currently not used.}

}

\details{
Parallelized matrix multiplication shortens computing time for large data
sets (\eqn{n>1000}).  However, spawning child processes requires itself
resources and increasing the number of cores for parallelized matrix
multiplication and parallelized computation of covariances does not always
result in reduced computing time.  Furthermore, these operations may be
initiated by child processes, that were itself spawned by functions like
\code{\link{cv.georob}}, \code{\link{predict.georob}},
\code{\link{profilelogLik}}, \code{\link{add1.georob}},
\code{\link{drop1.georob}} and \code{\link{step.georob}}.  By default,
parallelized matrix multiplication and computation of covariances is then
suppressed to avoid that child processes itself spawn child processes.  To
allow parallelized matrix multipliation and parallelized computation of
covariances by child processes one has to use the argument
\code{allow.recursive = TRUE}.

Note that very substantial reductions in computing time results when one
uses the \strong{OpenBLAS} library instead of the reference BLAS library
that ships with R, see
\url{https://www.openblas.net/} and R FAQ for your OS. With OpenBLAS no
gains are obtained by using more than one core for matrix
multiplication, and one should therefore use the default arguments
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
    \code{sfstop},\cr \code{allow.recursive}.} } }

\author{
   Andreas Papritz \email{papritz@retired.ethz.ch}.
}

\seealso{
  \code{\link{georobPackage}} for a description of the model and a brief summary of the algorithms;

  \code{\link{georob}} for (robust) fitting of spatial linear models.
}


\examples{
if(interactive()){
  ## example is run only in interactive session because cpu times exceeds 5 s
  A <- as.matrix(dist(rnorm(2000)))
  B <- as.matrix(dist(rnorm(2000)))
  system.time(C <- A \%*\% B)
  system.time(C <- pmm(
      A, B, control = control.pcmp(pmm.ncores = 2L)))
}
}

\keyword{models}
\keyword{spatial}
\keyword{robust}
