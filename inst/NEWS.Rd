% 2024-02-09 A. Papritz
\name{georobNews}
\encoding{UTF-8}

\title{News for Package \pkg{georob}}

\section{Changes in version 0.3-18 (2024-02-09)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item minor revision of vignette
      
      \item substituted all \code{donttest} statements in examples with
      lengthy computations by interactive() clause which has the effect
      that executes lenghty computations are only carried out interactive
      sessions

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item correction of error in \code{predict.georob() when processing
      output}
      
      \item \code{plot.sample.variogram()} now returns NULL
      
      \item eliminated lonely item in \code{NEWS.Rd}

    }
  }
}


\section{Changes in version 0.3-17 (2024-02-01)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item eliminated dependency on functions \code{RFfctn()} and
      \code{RFsimulate()} of archived package \pkg{RandomFields} by new
      functions \code{gencorr()}, \code{sim.chol.decomp()},
      \code{sim.circulant.embedding()}

      \item added \code{on.exit(options(oldoptions))} and
      \code{on.exit(par(oldpar))} after call \code{options()} and
      \code{par()} in function code

      \item querying class attribute by \code{inherits()}

      \item changed default for argument \code{allow.recursive} of function
      \code{control.pcmp()}

      \item more efficient calculation of \code{lag.vectors} for
      geometrically anisotropic variograms

      \item thorough revision of function \code{condsim()} (new functions
      \code{sim.chol.decomp()} and \code{sim.circulant.embedding()} for
      simulating unconditional Gaussian realizations by Cholesky matrix
      decomposition and circulant embedding method)

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item corrected parameter bounds for models RMlgd, RMdampedcos, RMfbm
      and RMgenfbm

      \item corrected symbol for signal variance in package help page and
      vignette
      
      \item complemented value section of various help pages, diverse fixes 
      in order to comply with CRAN policies in Rd-files

    }
  }
}


\section{Changes in version 0.3-16 (2023-01-28)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item substituted dependency on \pkg{RandomFields}, version 3.3.14,
      by dependency on \pkg{RandomFieldsMLV3010} version 0.1.0 (= legacy
      version of \pkg{RandomFields}, version 3.0.10, with reduced
      functionality)

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item correction of error in \code{add1.georob}, \code{drop1.georob}
      and \code{step.georob} occuring in parallelized computations on
      windows OS

      \item forcing socket clusters for parallelized computations in
      \code{condsim} on non-windows OS

    }
  }
}


\section{Changes in version 0.3-15 (2023-01-24)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item new file \code{inst/NEWS.Rd}

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item new default (FALSE) for argument allow.recursive of
      \code{control.pcmp}

      \item sanitized help pages complying with KaTeX

      \item only NAMESPACE is imported for the formerly attached packages
      \pkg{graphics}, \pkg{parallel}, \pkg{snowfall}, \pkg{stats} and
      \pkg{utils}

      \item vignette and help pages that use suggested packages
      \pkg{gstat}, \pkg{lattice} and \pkg{multcomp} changed

    }
  }
}


\section{Changes in version 0.3-14 (2021-06-04)}{
  \subsection{BUG FIXES}{
    \itemize{

      \item correction of error in check of initial values of anisotropy
      parameters in function \code{georob.fit}

    }
  }
}


\section{Changes in version 0.3-13 (2020-03-27)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item functions \code{fit.variogram.model} and \code{georob} compute
      now Hessian of non-transformed variogram parameters in addition to
      Hessian of transformed parameters

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item correction of error in sanity checks of arguments of function
      \code{georob}

    }
  }
}


\section{Changes in version 0.3-12 (2020-02-16)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item new data sets \code{wolfcamp} and \code{elevation}

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item all code files renamed

      \item sanity checks of arguments in calls of \code{if} and
      \code{switch}

      \item DESCRIPTION: omission of \emph{suggest} dependency for orphaned
      package \pkg{geoR} in

      \item vignette: conditional use of objects provided by
      \emph{suggested} packages \pkg{gstat} and \pkg{multcomp}

    }
  }
}


\section{Changes in version 0.3-11 (2019-12-18)}{
  \subsection{BUG FIXES}{
    \itemize{

      \item NAMESPACE: correction of dependencies

      \item DESCRIPTION: correction of dependencies and licence information

      \item method \code{summary.georob} takes terms component now from
      georob object

      \item function \code{condsim}: sanity checks of specified function
      arguments and for calls of \code{if}

      \item function \code{condsim}: correction of errors when processing
      output of function \code{RFsimulate}

      \item changes in calls of \code{class}

    }
  }
}


\section{Changes in version 0.3-10 (2019-05-24)}{
  \subsection{BUG FIXES}{
    \itemize{

      \item method \code{summary.georob}: correction of error when
      printing confidence intervals of variogram parameters

      \item function \code{condsim}: correction of error that occurred when
      preparing output consisting of a single realization

    }
  }
}


\section{Changes in version 0.3-9 (2019-04-07)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item function \code{fit.variogram.model}: new control argument and
      possibility to choose between \code{optim} and \code{nlminb}
      optimizer for non-linear least squares

      \item new function \code{control.fit.variogram.model} for providing
      sensible defaults to argument \code{control} of
      \code{fit.variogram.model}

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item internal function \code{f.aux.gcr}: correction of error
      introduced in version 0.3-8 when computing semivariances for
      anisotropic variogram models

    }
  }
}


\section{Changes in version 0.3-8 (2019-03-29)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item function \code{condsim}: conditioning data is passed as
      SpatialPointsDataFrame to \code{RFsimulate}

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item internal function \code{f.aux.gcr}: changes for new version
      3.3 of package \pkg{RandomFields}

    }
  }
}


\section{Changes in version 0.3-7 (2018-08-27)}{
  \subsection{BUG FIXES}{
    \itemize{

      \item internal function \code{f.aux.sim.1}: elimination of grid.nodes
      argument

      \item methods \code{add1.georob}, \code{drop1.georob},
      \code{step.georob}: elimination of data argument

    }
  }
}


\section{Changes in version 0.3-6 (2018-01-26)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item improved memory management in all parallized computations

      \item method \code{predict.georob}: generates error message if
      names of coordinates in object and newdata are not the same

      \item method \code{predict.georob}: new optional arguments for
      specifying a variogram model different from the one contained in
      object

      \item new functions \code{condsim} and \code{control.condsim} for
      (un)conditional simulations

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item internal function \code{f.call.set_allxxx_to_fitted_values}
      variogram.model is also updated

    }
  }
}


\section{Changes in version 0.3-5 (2017-10-20)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item internal function \code{georob.fit}: changes to allow non-zero
      snugget if there are no replicated observations

      \item new methods \code{coef.georob} and \code{print.coef.georob}

    }
  }
}


\section{Changes in version 0.3-4 (2017-05-09)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item internal function \code{georob.fit}: issuinng warning for
      negative-definite hessian

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item functions \code{georob} and \code{control.rq}: small changes
      and error correction in interface to \code{rq} for computing robust
      initial values

    }
  }
}


\section{Changes in version 0.3-3 (2017-01-24)}{
  \subsection{BUG FIXES}{
    \itemize{

      \item function \code{fit.variogram.model}: correction of minor bug
      for fitting 3d anisotropic variogram models to sample variograms

    }
  }
}


\section{Changes in version 0.3-2 (2016-11-28)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item internal function \code{georob.fit}: checking ml.method and
      presence of intercept for intrinsic models

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item method \code{predict.georob}: correction of error when
      computing predictions for intrinsic variograms

      \item interal function \code{likelihood.calculations}: correction of
      error when printing trial parameter values if they exceed safe bounds

      \item (internal) functions \code{fit.variogram.model},
      \code{likelihood.calculations,georob.fit},
      \code{partial.derivatives.variogram}: correcting of error in 3d-rotation
      matrix for geometrically anisotropic variograms

    }
  }
}


\section{Changes in version 0.3-1 (2016-09-09)}{
  \subsection{BUG FIXES}{
    \itemize{

      \item correction of minor error in internal function
    }
  }
}


\section{Changes in version 0.3-0 (2016-09-08)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item functionality to use nested variogram models

      \item function \code{control.georob}: change of default value for
      argument error.family.cov.residuals

      \item method \code{plot.sample.variogram}: new argument lty

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item internal function \code{bwd.transf}: correction of error for
      logit-transformations

      \item method \code{cv.georob}: correction of error when computing
      standard errors of cv kriging predictions

    }
  }
}


\section{Changes in version 0.2-4 (2016-07-27 )}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item improvements in all parallelized computations

      \item methods \code{add1.georob}, \code{drop1.georob},
      \code{step.georob}: extended functionality for refitting variogram
      parameters or keeping them fixed, improved computational efficiency,
      diagnostics about convergence

      \item method \code{cv.georob}: allow recursive calls of
      \code{mclapply}

      \item function \code{georob}: allowing use of \code{lmrob} for computing
      initial fixed effects for rank-deficient model matrix

      \item function \code{profilelogLik} and method
      \code{waldtest.georob}: improved computationsl efficiency

      \item method \code{predict.georob}: SpatialPoints, SpatialPixels and
      SpatialGrid as newdata objects

      \item method \code{plot.cv.georob}: added line color and type for
      scatterplot smooths


    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item method \code{fixef.georob}: correcting wrong component name

      \item internal function \code{georob.fit}: corrected names of
      gradient components

    }
  }
}


\section{Changes in version 0.2-3 (2016-03-02)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item method \code{plot.georob}: changes for residual diagnostic plots

      \item method \code{plot.cv.georob}: minor changes for adding plots to
      existing diagnostics graphics

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item ...

    }
  }
}


\section{Changes in version 0.2-2 (2015-12-03)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item new way to control which variogram parameters are fitted

      \item generic and methods for computing sample variogram

      \item methods \code{cv.georob} and \code{predict.georob}: catching
      errors occurring during parallelized computations

      \item functions \code{validate.preditions},
      \code{fit.variogram.model}: check whether mandatory arguments were
      provided

      \item method \code{plot.georob}: changes for Tukey-Anscombe, normal
      QQ-plots and for plotting variograms

      \item issuing additional warnings in various computations

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item function \code{fit.variogram.model}: correcting error in
      setting default values for fit.param

    }
  }
}


\section{Changes in version 0.2-1 (2015-09-02)}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item Gaussian (RE)ML estimation of parameters of re-parametrized
      variogram

      \item new function \code{profilelogLik} for computing profile
      loglikelihood of variogram parameters

      \item extended functionality for transforming variogram parameters

      \item changes for improving computational efficiency

      \item functionality for optional computation of covariances for
      long-tailed error distributions

      \item many more improvements, see ChangeLog for period from
      2015-04-03 to 2015-09-02

    }
  }
  \subsection{BUG FIXES}{
    \itemize{

      \item correction of multiple errors, see ChangeLog for period from
      2015-04-03 to 2015-09-02


    }
  }
}
%
%
% \section{Changes in version x.y-z (YYYY-MM-DD)}{
%   \subsection{NEW FEATURES}{
%     \itemize{
%
%       \item ...
%
%     }
%   }
%   \subsection{BUG FIXES}{
%     \itemize{
%
%       \item ...
%
%     }
%   }
% }
