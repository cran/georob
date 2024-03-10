##    ##############################################################################
### f.aux.gcr

f.aux.gcr <-
function(
  lag.vectors, variogram.object, gcr.constant = NULL, symmetric = TRUE,
  irf.models = control.georob()[["irf.models"]],
  control.pcmp, verbose
)
{

  ##  Function computes the generalized correlation (matrix) for the lag
  ##  distances in lag.vectors.  The result is a generalized correlation matrix
  ##  that is positive definite.

  ##  cf. HRK's notes of 2011-06-17 on "Robust Kriging im intrinsischen
  ##  Fall"

  ##  2011-12-27 ap
  ##  2012-02-07 AP modified for geometrically anisotropic variograms
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-03-05 AP changes for version 3 of RandomFields
  ## 2014-08-18 AP changes for parallelized computations
  ## 2015-07-17 AP new name of function, Gaussian (RE)ML estimation for reparametrized variogram
  ## 2015-07-23 AP Valpha (correlation matrix without spatial nugget) no longer stored
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-09 AP changes for nested variogram models
  ## 2016-08-10 AP changes for isotropic variogram models
  ## 2017-12-23 AP improved memory management in parallel computations
  ## 2019-03-29 AP call to RFvariogram{RandomFields, version 3.1} replaced by
  ##               call to RFfctn{RandomFields, version 3.3}
  ## 2019-04-07 AP correction of error when computing semivariances for anisotropic variogram models
  ## 2023-11-17 AP substitution of function RFfctn{RandomFields} by new function gencorr
  ## 2023-11-17 AP elimination of calls to function RFoptions{RandomFields}
  ## 2023-12-20 AP added on.exit(options(old.opt)), deleted options(error = NULL)
  ## 2023-12-20 AP replacement of identical(class(...), ...) by inherits(..., ...)
  ## 2023-12-28 AP modified gcr.constant for RMfbm model

#### -- check consistency of arguments

  if( !is.null( gcr.constant ) ){
    if( !is.list( gcr.constant ) ||
      !identical( length(gcr.constant), length(variogram.object) )
    ) stop( "lengths of 'gcr.constant' and 'variogram.object' differ" )
  } else {
    gcr.constant <- as.list( rep( NA_real_, length(variogram.object) ) )
  }

#### -- compute generalized correlations

  res <- lapply(
    1L:length(variogram.object),
    function( i, x, gcr.constant, lag.vectors, control.pcmp ){

      variogram.model <- x[[i]][["variogram.model"]]
      param           <- x[[i]][["param"]]
      param           <- param[!names(param) %in% c( "variance", "snugget", "nugget")]
      aniso           <- x[[i]][c("aniso", "sclmat", "rotmat")]
      gcr.constant    <- gcr.constant[[i]]

      result <- list( error = TRUE )


# #### --- preparation

      if( NCOL( lag.vectors ) > 1L ){

        ### anisotropic model

        ### matrix for coordinate transformation

        A <- aniso[["sclmat"]] * aniso[["rotmat"]] / param["scale"]

        ### rotated and scaled lag distance

        scaled.lag.distance <- sqrt( rowSums( (lag.vectors %*% t(A))^2 ) )

      } else {

        ### isotropic model

        ### scaled lag distance

        scaled.lag.distance <- lag.vectors / param["scale"]

      }

      ### auxiliary function to compute generalized correlations in parallel

      f.aux <- function( i ){

        ### note: gencorr computes the negative semivariance
        ### for stationary and IRF models, required is negative semivariance
        result <- try(
          gencorr(
            x = scaled.lag.distance[s[i]:e[i]],
            variogram.model = variogram.model, param = param[-1L]
          ),
          silent = TRUE
        )

        if( inherits( result, "try-error" ) || any( is.na( result ) ) ){
          "gen.corr.error"
        } else {
          result
        }

      }


#### --- compute generalized correlations

      ## determine number of cores

      ncores <- control.pcmp[["gcr.ncores"]]

      ## definition of junks to be evaluated in parallel

      k <- control.pcmp[["f"]] * ncores
      n <- NROW(lag.vectors)
      dn <- floor( n / k )
      s <- ( (0L:(k-1L)) * dn ) + 1L
      e <- (1L:k) * dn
      e[k] <- n

      ## set default value for control of forking if missing (required for backward compatibility)

      if( is.null( control.pcmp[["fork"]] ) ){
        control.pcmp[["fork"]] <- !identical( .Platform[["OS.type"]], "windows" )
      }

      ## compute generalized correlations in parallel

      if( ncores > 1L && !control.pcmp[["fork"]] ){

        if( !sfIsRunning() ){
          old.opt <- options( error = f.stop.cluster )
          on.exit( options( old.opt ) )

          junk <- sfInit( parallel = TRUE, cpus = ncores )

          junk <- sfLibrary( georob, verbose = FALSE )

          junk <- sfExport(
            "s", "e", "scaled.lag.distance", "variogram.model", "param"
          )

         }

        Valpha <- sfLapply( 1L:k, f.aux )

        if( control.pcmp[["sfstop"]] ){
          junk <- sfStop()
        }

      } else {

        Valpha <- mclapply( 1L:k, f.aux, mc.cores = ncores )

      }

      not.ok <- any( sapply(
          Valpha,
          function( x ) identical( x, "gen.corr.error" ) || any( is.na( x ) )
        ))

      if( !not.ok ){

        Valpha <- unlist( Valpha )

        ##  compute additive constant for positive definiteness, this
        ##  implements a sufficient condition for positive definiteness of
        ##  Valpha (strong row sum criterion)

        if( is.na( gcr.constant ) ){
          if( variogram.model %in% irf.models ){
            if( identical( variogram.model, "RMfbm" ) ){
              ## cf Chiles & Delfiner, 1999, eqs 7.34 & 7.35, p. 511
              t.a <- param["alpha"]
              t.ah <- 0.5 * t.a
              t.n <- 0.5 + 0.5 * NCOL(lag.vectors)
              gcr.constant <- max(
                2.,
                gamma(0.5 + t.ah) * gamma(1. - t.ah) / sqrt(pi) *
                gamma(t.a + t.n) / gamma(1. + t.a) / gamma(t.n)
              )
            } else {
              gcr.constant <- 2.
            }
            gcr.constant <- -min( Valpha ) * gcr.constant

          } else {
            gcr.constant <- 1.
          }
        }

        ## compute generalized correlation

        Valpha <- gcr.constant + Valpha

        ## convert generalized correlation vector to symmetric correlation matrices

        if( symmetric ){
          Valpha <- list(
            diag = rep( gcr.constant, 0.5 * ( 1. + sqrt( 1. + 8L * length( Valpha ) ) ) ),
            tri = Valpha
          )
          attr( Valpha, "struc" ) <- "sym"
        }

#### --- collect results

        result[["error"]]        <- FALSE
        result[["gcr.constant"]] <- gcr.constant
        result[["Valpha"]]       <- Valpha

      } else {

        warning( "there were errors: call function with argument 'verbose' > 1" )
        if( verbose > 3. ) cat(
          "\n an error occurred when computing the negative semivariance matrix\n"
        )

      }

      return( result )

    }, x = variogram.object, gcr.constant = gcr.constant,
    lag.vectors = lag.vectors, control.pcmp = control.pcmp
  )

}


##   ##############################################################################
### partial.derivatives.variogram

partial.derivatives.variogram <-
  function(
    d.param, x,
    variogram.model, param, aniso, sincos, sclmat, rotmat,
    verbose
  )
{

  ##  Function to compute partial derivatives of generalized
  ##  correlation matrix with respect to scale and extra parameters

  ##  06 Apr 2011  C.Schwierz
  ##  2011-07-17 ap
  ##  2012-01-24 ap RMcauchytbm and RMlgd models added
  ##  2012-01-25 ap extra model parameter with same names as in Variogram{RandomFields}
  ##  2012-02-07 AP modified for geometrically anisotropic variograms
  ##  2013-06-12 AP substituting [["x"]] for $x in all lists
  ##  2014-05-15 AP changes for version 3 of RandomFields
  ##  2015-07-17 AP new name of function, scaling with 1-xi eliminated
  ##  2016-08-04 AP changes for nested variogram models
  ##  2016-11-14 AP correcting error in 3d rotation matrix for geometrically anisotropic variograms
  ##  2020-02-07 AP sanity checks for switch() and if()

  aniso.name <- names( aniso )
  alpha <- unname( param["scale"] )
  n = NCOL( x )

  if( n > 1L ){

    aux <- rotmat %*% t(x)

    ## scaled lag distance

    hs <- sqrt( colSums( ( sclmat * aux )^2 ) ) / alpha

  } else {

    hs <- x / alpha
  }

#### -- partial derivatives of scaled lag distance with respect to
  ## anisotropy parameters

  dhs.daniso <- switch(

    d.param[1],

    f1 = {
      colSums(
        ( c( 0., -1. / aniso["f1"]^2, 0. )[1L:n] * sclmat ) * aux^2
      )
    },

    f2 = {
      colSums(
        ( c( 0., 0., -1. / aniso["f2"]^2 )[1L:n] * sclmat ) * aux^2
      )
    },
    omega = {
      drotmat <- with(
        sincos,
        rbind(
          c(             co*sp,            -so*sp, 0. ),
          c(  so*cz - co*cp*sz,  co*cz + so*cp*sz, 0. ),
          c( -so*sz - co*cp*cz, -co*sz + so*cp*cz, 0. )
        )[ 1L:n, 1L:n, drop = FALSE ]
      )
      colSums(
        ( sclmat * drotmat %*% t(x) ) * ( sclmat * aux )
      )
    },

    phi = {
      drotmat <- with(
        sincos,
        rbind(
          c(     so*cp,     co*cp,    -sp ),
          c(  so*sp*sz,  co*sp*sz,  cp*sz ),
          c(  so*sp*cz,  co*sp*cz,  cp*cz )
        )[ 1L:n, 1L:n, drop = FALSE ]
      )
      colSums(
        ( sclmat * drotmat %*% t(x) ) * ( sclmat * aux )
      )
    },

    zeta = {
      drotmat <- with(
        sincos,
        rbind(
          c(                0.,               0.,     0. ),
          c(  co*sz - so*cp*cz, -so*sz - co*cp*cz,  sp*cz ),
          c(  co*cz + so*cp*sz, -so*cz + co*cp*sz, -sp*sz )
        )[ 1L:n, 1L:n, drop = FALSE ]
      )
      colSums(
        ( sclmat * drotmat %*% t(x) ) * ( sclmat * aux )
      )
    },

    NA
  ) / ( hs * alpha^2 )

#### -- partial derivative of scaled lag distance with respect to scale
  ##  parameter

  dhs.dscale <- -hs / alpha

#### -- compute derivative of generalized correlation matrix with
  ##  respect to scale and extra parameters

  result <- switch(
    variogram.model[1],

#### --- RMbessel
    RMbessel = {

      A <- unname( param["nu"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -( 2.^A * besselJ( hs, 1.+A ) * gamma( 1+A ) ) / hs^A
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( 2^A * besselJ( hs, 1+A ) * gamma(1 + A) ) / hs^A,
      #   0.
      # )


      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        nu = {
          myenv <- new.env()
          assign( "hs", hs, envir = myenv )
          assign( "nu", param["nu"], envir = myenv )
          as.vector(
            attr(
              numericDeriv(
                expr = quote(
                  2.^nu * gamma( nu+1. ) * besselJ( hs, nu ) / hs^nu
                ),
                theta = "nu",
                rho = myenv
              ),
              "gradient"
            )
          )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMbessel

#### --- RMcauchy
    RMcauchy = {

      A <- unname( param["gamma"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -2. * A * hs * ( 1.+hs^2 )^(-1.-A)

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        gamma = {
          -( 1. + hs^2 )^(-A) * log( 1. + hs^2 )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMcauchy

#     RMcauchytbm = {
#
#       A <- unname( param["alpha"] )
#       B <- unname( param["beta"] )
#       C <- unname( param["gamma"] )
#
#       ## derivative with of generalized covariance with respect to
#       ## scaled lag distance
#
#       dgc.dhs <- -(
#         B * hs^(-1.+A) * (1.+hs^A)^(-2.-B/A) * ( A + C - (-B + C) * hs^A )
#       ) / C
#       # dgc.dhs <- ifelse(
#       #   hs > 0.,
#       #   -(
#       #     B * hs^(-1.+A) * (1.+hs^A)^(-2.-B/A) * ( A + C - B * hs^A + C * hs^A)
#       #   ) / C,
#       #   if( A > 1. ){
#       #     0.
#       #   } else if( identical( A, 1 ) ){
#       #     -B * (1.+C) / C
#       #   } else {
#       #     -Inf
#       #   }
#       # )
#
#
#       switch(
#         d.param[1],
#         scale = dgc.dhs * dhs.dscale,
#         # scale = {
#         #   ( B * hs^A * (1.+hs^A)^(-2.-B/A) * (A + C + (-B+C) * hs^A ) ) / ( C * scale )
#         # },
#         alpha = {
#           ( B * (1.+hs^A)^(-2. - B/A) * (
#               -( A * hs^A * ( A + C + (-B+C) * hs^A ) * log(hs) ) +
#               ( 1. + hs^A) * (C + (-B+C) * hs^A ) * log( 1.+hs^A )
#             )
#           ) / (A^2 * C )
#         },
#         # alpha = {
#         #   ifelse(
#         #     hs > 0.,
#         #     ( B * (1.+hs^A)^(-2. - B/A) * (
#         #         -( A * hs^A * ( A + C + (-B+C) * hs^A ) * log(hs) ) +
#         #         ( 1. + hs^A) * (C + (-B+C) * hs^A ) * log( 1.+hs^A )
#         #       )
#         #     ) / (A^2 * C ),
#         #     0.
#         #   )
#         # },
#         beta = {
#           ( -( A * hs^A) - (C + (-B+C) * hs^A ) * log( 1.+hs^A ) ) /
#           ( A*C * (1.+hs^A)^( (A+B)/A ) )
#         },
#         gamma = {
#           ( B * hs^A ) / ( C^2 * (1.+hs^A)^( (A+B)/A) )
#         },
#         dgc.dhs * dhs.daniso
#       )
#     }, ##  end case RMcauchytbm

#### --- RMcircular
    RMcircular = {

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- ( -4. * sqrt( 1.-hs[sel]^2 ) ) / pi

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMcircular

#### --- RMcubic
    RMcubic = {

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- hs[sel] * ( -14. + 26.25*hs[sel] - 17.5*hs[sel]^3 + 5.25*hs[sel]^5 )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMcubic

#### --- RMdagum
    RMdagum = {

      A <- unname( param["beta"] )
      B <- unname( param["gamma"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -B / ( hs * ( 1.+hs^(-A) )^(B/A) * ( 1.+hs^A ) )
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( B / ( hs * ( 1.+ hs^(-A) )^(B/A) * (1. + hs^A ) ) ),
      #   -Inf
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse(
        #     hs > 0.,
        #     B / ( ( 1. + hs^(-A) )^(B/A) * (1. + hs^A) * scale ),
        #     0.
        #   )
        # },
        beta = {
          -( B * ( A * log(hs) + (1.+hs^A) * log( 1.+hs^(-A) ) ) ) /
          ( A^2 * ( 1.+hs^(-A) )^(B/A) * ( 1.+hs^A ) )
        },
        # beta = {
        #   ifelse(
        #     hs > 0.,
        #     -( B * ( A * log(hs) + (1.+hs^A) * log( 1.+hs^(-A) ) ) ) /
        #     ( A^2 * ( 1.+hs^(-A) )^(B/A) * ( 1.+hs^A ) ),
        #     0.
        #   )
        # },
        gamma = {
          log( 1. + hs^(-A) ) / ( A * (1. + hs^(-A) )^(B/A) )
        },
        # gamma = {
        #   ifelse(
        #     hs > 0.,
        #     log( 1. + hs^(-A) ) / ( A * (1. + hs^(-A) )^(B/A) ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMdagum

#### --- RMdampedcos
    RMdampedcos = {

      A <- unname( param["lambda"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -( ( A * cos(hs) + sin(hs) ) / exp( A*hs ) )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        lambda = {
          -exp( -A * hs ) * hs * cos( hs )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMdampedcos

#### --- RMdewijsian
    RMdewijsian = {


      A <- unname( param["alpha"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -A / ( hs + hs^(1.-A) )
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( A / ( hs + hs^(1.-A) ) ),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( A * hs^A )/( scale + hs^A * scale )
        # },
        alpha = {
          -( ( hs^A * log( hs ) ) / ( 1. + hs^A ) )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -( ( hs^A * log( hs ) ) / ( 1. + hs^A ) ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMdewijsian


#### --- RMexp
    RMexp = {

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -exp( -hs )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case exponential

#### --- RMfbm
    RMfbm = {

      A <- unname( param["alpha"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -A * hs^(-1.+A)
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( A * hs^(-1.+A) ),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   A * hs^A / scale
        # },
        alpha = {
          -hs^A * log( hs )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -( hs^A * log( hs ) ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMfbm

#### --- RMgauss
    RMgauss = {

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -2. * hs / exp( hs^2 )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgauss

#### --- RMgenfbm
    RMgenfbm = {

      A <- unname( param["alpha"] )
      B <- unname( param["delta"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -( A * B * hs^(-1.+A) * (1.+hs^A)^(-1.+B))
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( A * B * hs^(-1.+A) * (1.+hs^A)^(-1.+B)),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -B.
      #   } else {
      #     0.
      #   }
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( A * B * hs^A * (1.+hs^A)^(-1.+B) ) / scale
        # },
        alpha = {
          -( B * hs^A * (1.+hs^A)^(-1.+B) * log(hs) )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -( B * hs^A * (1.+hs^A)^(-1.+B) * log(hs) ),
        #     0.
        #   )
        # },
        delta = {
          -( (1. + hs^A )^B * log( 1. + hs^A ) )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgenfbm

#### --- RMgencauchy
    RMgencauchy = {

      A <- unname( param["alpha"] )
      B <- unname( param["beta"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -( ( B * hs^(-1.+A)) / (1.+hs^A)^((A+B)/A))
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( ( B * hs^(-1.+A)) / (1.+hs^A)^((A+B)/A)),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -B.
      #   } else {
      #     0.
      #   }
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( B * hs^A ) / ( (  1.+hs^A)^((A+B)/A) * scale )
        # },
        alpha = {
          B * ( 1. + hs^A )^(-(A+B)/A) * (
            -A * hs^A * log( hs ) +
            ( 1. + hs^A ) * log( 1. + hs^A )
          ) / A^2
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     B * ( 1. + hs^A )^(-(A+B)/A) * (
        #       -A * hs^A * log( hs ) +
        #       ( 1. + hs^A ) * log( 1. + hs^A )
        #     ) / A^2,
        #     0.
        #   )
        # },
        beta = {
          -( log( 1.+hs^A ) / ( A * (1.+hs^A)^(B/A) ) )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgencauchy


#### --- RMgengneiting
    RMgengneiting = {


      A <- unname( param["kappa"] )
      B <- unname( param["mu"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- if( identical( as.integer(A), 1L ) ){
        (2.5+B) * (1.-hs[sel])^(2.5+B) - (2.5+B) * (1.-hs[sel])^(1.5+B) * (1.+(2.5+B) * hs[sel])
      } else if( identical( as.integer(A), 2L ) ){
        (1. - hs[sel])^(4.5+B) * (4.5+B + 2./3. * (3.5+B) * (5.5+B) * hs[sel] ) -
        (4.5+B) * (1. - hs[sel])^(3.5+B) * (
          1. + hs[sel]*(4.5 + B + 6.416666666666666*hs[sel] + B/3. * (9.+B) * hs[sel] )
        )
      } else if( identical( as.integer(A), 3L ) ){
        (1. - hs[sel])^(6.5+B) * (6.5 + B + 0.8 * (5.275255128608411+B) * (7.724744871391589+B) * hs[sel] +
          0.2 * (4.5+B) * (6.5+B) * (8.5+B) * hs[sel]^2) -
        (6.5+B) * (1. - hs[sel])^(5.5+B) * (1. + (6.5+B) * hs[sel] + 0.4 * (5.275255128608411+B) *
          (7.724744871391589+B) * hs[sel]^2 + 0.2/3 * (4.5+B) * (6.5+B) * (8.5+B) * hs[sel]^3
        )
      } else {
        stop( "RMgengneiting model undefined for 'kappa' != 1:3" )
      }

      result <- rep( 0., length( hs ) )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        mu = {
          result[sel] <- if( identical( as.integer(A), 1L ) ){
            (1. - hs[sel])^(2.5+B) * (hs[sel] + (1. + (2.5+B) * hs[sel]) * log(1. - hs[sel]))
          } else if( identical( as.integer(A), 2L ) ){
            (1. - hs[sel])^(4.5+B) * (hs[sel] + 2./3. * (4.5+B) * hs[sel]^2 + (1. + hs[sel] * (
                  4.5 + B + 6.416666666666666*hs[sel] +  B/3. * (9.+B) * hs[sel]) ) * log(1. - hs[sel])
            )
          } else if( identical( as.integer(A), 3L ) ){
            (1. - hs[sel])^(6.5 + B)*
            (hs[sel] + (5.2 + 0.8*B)*hs[sel]^2 +
              0.2*(5.345299461620754 + B)*
              (7.654700538379246 + B)*hs[sel]^3 +
              (1. + hs[sel]*(6.5 + 1.*B +
                  0.4*(5.275255128608411 + B)*
                  (7.724744871391589 + B)*hs[sel] +
                  0.06666666666666667*(4.5 + B)*

                  (6.5 + B)*(8.5 + B)*hs[sel]^2))*
              log(1. - hs[sel]))
          }
          result
        },
        dgc.dhs * dhs.daniso
      )


    }, ##  end case Gengneiting


#### --- RMgneiting
    RMgneiting = {

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < (1. / 0.301187465825)
      dgc.dhs[sel] <- (1. - 0.301187465825*hs[sel])^7 * (
        -1.9957055705418814*hs[sel] -  4.207570523270417*hs[sel]^2 - 2.896611435848653*hs[sel]^3
      )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgneiting

#### --- RMlgd
    RMlgd = {

      A <- unname( param["alpha"] )
      B <- unname( param["beta"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -( A * B * hs^(-1.-B) ) / (A+B)
      sel <- hs <= 1.
      dgc.dhs[sel] <- -( A * B * hs[sel]^(-1.+A) ) / (A+B)

      # dgc.dhs <- ifelse(
      #   hs > 0.
      #   ifelse(
      #     hs <= 1.,
      #     -( A * B * hs^(-1.+A) ) / (A+B),
      #     -( A * B * hs^(-1.-B) ) / (A+B)
      #   ),
      #   if( identical( A, 1. ) ){
      #     -B / ( B + 1. )
      #   } else if( A < 1. ){
      #     -Inf
      #   }
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse(
        #     hs > 0.,
        #     ifelse(
        #       hs <= 1.,
        #       A * B * hs^A,
        #       A * B / hs^B
        #     ) / ( (A+B) * scale ),
        #     0.
        #   )
        # },
        alpha = {
          result <- B / ( (A+B)^2 * hs^B )
          sel <- hs <= 1.
          result[sel] <- -( B * hs[sel]^A * ( -1. + (A+B ) * log( hs[sel] ) ) ) / (A+B)^2
          result
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     ifelse(
        #       hs <= 1.,
        #       -( B * hs^A * ( -1. + (A+B ) * log( hs ) ) ) / (A+B)^2,
        #       B / ( (A+B)^2 * hs^B )
        #     ),
        #     0.
        #   )
        # },
        beta = {
          result <- -A * ( 1. + (A+B) * log( hs ) ) / ( (A+B)^2 * hs^B )
          sel <- hs <= 1.
          result[sel] <- -A * hs[sel]^A / (A+B)^2
          result
        },
        # beta = {
        #   ifelse(
        #     hs > 0.,
        #     ifelse(
        #       hs <= 1.,
        #       -A * hs^A / (A+B)^2,
        #       -A * ( 1. + (A+B) * log( hs ) ) / ( (A+B)^2 * hs^B )
        #     ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMlgd

#### --- RMmatern
    RMmatern = {

      A <- unname( param["nu"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -(
        2.^(1.5 - A/2.) * sqrt(A) * ( sqrt(A) * hs )^A * besselK( sqrt(2.*A)*hs, -1.+A )
      ) / gamma(A)
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -(
      #     ( 2.^(1.5 - A/2.) * sqrt(A) * ( sqrt(A) * hs )^A * besselK( sqrt(2.) * sqrt(A) * hs , -1.+A )
      #     ) / gamma(A)
      #   ),
      #   if( A < 0.5 ){
      #     -Inf
      #   } else if( identical( A, 0.5 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse(
        #     hs > 0.,
        #     ( 2.^(1.5 - A/2.) * ( sqrt(A) * hs )^(1.+A) *
        #       besselK( sqrt(2.*A) * hs, A-1.)
        #     ) / (scale * gamma(A) ),
        #     0.
        #   )
        # },
        nu = {
          myenv <- new.env()
          assign( "hs", hs, envir = myenv )
          assign( "nu", param["nu"], envir = myenv )
          as.vector(
            attr(
              numericDeriv(
                expr = quote(
                  2.^(1.-nu) / gamma(nu) *
                  ( sqrt( 2.*nu ) * hs )^nu * besselK( sqrt( 2.*nu ) * hs, nu )
                ),
                theta = "nu",
                rho = myenv
              ),
              "gradient"
            )
          )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMmatern

#### --- RMpenta
    RMpenta = {

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- ( 11. * (-1.+hs[sel])^5 * hs[sel] * (2.+hs[sel]) * ( 4. + hs[sel] * ( 18. + 5. * hs[sel] * (3.+hs[sel]) ) ) ) / 6.

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMpenta

#### --- RMaskey
    RMaskey = {

      A <- unname( param["alpha"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- -(A * (1.-hs[sel])^(-1.+A))

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        alpha = {
          result <- rep( 0., length( hs ) )
          sel <- hs < 1.
          result[sel] <- ( 1. - hs[sel] )^A * log( 1. - hs[sel] )
          result
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMaskey

#### --- RMqexp
    RMqexp = {

      A <- unname( param["alpha"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- 2. * (-A + exp(hs) ) / ( (-2.+A ) * exp(2.*hs) )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        alpha = {
          ( 2. * exp( -2.*hs ) * ( -1. + exp( hs ) ) ) / (-2.+A)^2
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMqexp


#### --- RMspheric
    RMspheric = {

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- -1.5 + 1.5 * hs[sel]^2

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMspheric

#### --- RMstable
    RMstable = {

      A <- unname( param["alpha"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -( ( A * hs^(-1.+A) ) / exp(hs^A) )
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( ( A * hs^(-1.+A) ) / exp(hs^A) ),
      #   if( A > 1. ){
      #     0.
      #   } else if( identical( A, 1. ) ){
      #     -1.
      #   } else {
      #     -Inf
      #   }
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( A * exp( -hs^A ) * hs^A ) / scale
        # },
        alpha = {
          -exp( -hs^A ) * hs^A * log( hs )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -exp( -hs^A ) * hs^A * log( hs ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMstable

#### --- RMwave
    RMwave = {

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- ( hs * cos(hs) - sin(hs) ) / hs^2
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   ( hs * cos(hs) - sin(hs) ) / hs^2,
      #   0.
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case wave

#### --- RMwhittle
    RMwhittle = {

      A <- unname( param["nu"] )

      ## derivative with of generalized covariance with respect to
      ## scaled lag distance

      dgc.dhs <- -( 2.^(1.-A) * hs^A * besselK( hs, -1.+A ) ) / gamma(A)
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( 2^(1.-A) * hs^A * besselK( hs, -1.+A ) ) / gamma(A),
      #   if( A < 0.5 ){
      #     -Inf
      #   } else if( identical( A, 0.5 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )

      switch(
        d.param[1],
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse(
        #     hs > 0.,
        #     (
        #       2.^(1-A) * h * hs^A * besselK( hs, -1.+A )
        #     ) / ( scale^2 * gamma(A) ),
        #     0.
        #   )
        #
        # },
        nu = {
          myenv <- new.env()
          assign( "hs", hs, envir = myenv )
          assign( "nu", param["nu"], envir = myenv )
          as.vector(
            attr(
              numericDeriv(
                expr = quote(
                  2.^(1.-nu) / gamma(nu) * hs^nu * besselK( hs, nu )
                ),
                theta = "nu",
                rho = myenv
              ),
              "gradient"
            )
          )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case whittle

    stop(
      paste(
        variogram.model,
        "model: derivatives for scale, extra variogram and anisotropy parameters undefined"
      )
    )

  ) ##  end switch cov.model

#### -- convert to matrix

  result <- list(
    diag = rep( 0., 0.5 * ( 1. + sqrt( 1. + 8. * length( result ) ) ) ),
    tri = result
  )
  attr( result, "struc" ) <- "sym"
  result <- expand( result )

  return( result )

}


################################################################################

f.aux.Valphaxi <- function(
  lag.vectors, variogram.object, xi, control.pcmp, verbose
){

  ## auxiliary function to compute generalized correlation matrix and
  ## related items (called by likelihood.calculations)

  ## 2014-07-29 A. Papritz
  ## 2015-07-23 AP Valpha (correlation matrix without spatial nugget) no longer stored
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-04 AP changes for nested variogram models
  ## 2023-12-21 AP replacement of identical(class(...), ...) by inherits(..., ...)

  result <- list(
    error = TRUE, Valpha = NULL, Valphaxi = NULL,
    Valphaxi.inverse = NULL, log.det.Valphaxi = NULL
  )

  Valpha <- f.aux.gcr(
    lag.vectors = lag.vectors,
    variogram.object = variogram.object,
    gcr.constant = NULL,
    control.pcmp = control.pcmp,
    verbose = verbose
  )

  if( any( sapply( Valpha, function(x) x[["error"]] ) ) ) return( result )

  Valphaxi <- list(
    diag = rowSums(
      sapply(
        1L:length(Valpha),
        function( i, x, xi ){
          xi[i] * x[[i]][["Valpha"]][["diag"]]
        }, x = Valpha, xi = xi
      )
    ) + ( 1. - sum(xi) ),
    tri = rowSums(
      sapply(
        1L:length(Valpha),
        function( i, x, xi ){
          xi[i] * x[[i]][["Valpha"]][["tri"]]
        }, x = Valpha, xi = xi
      )
    )
  )
  attr( Valphaxi, "struc" ) <- "sym"

  Valphaxi <- expand( Valphaxi )

  t.vchol <- try( chol( Valphaxi ), silent = TRUE )

  if( !inherits( t.vchol, "try-error" ) ){

    result[["error"]]            <- FALSE
    result[["Valpha"]]           <- Valpha
    result[["Valphaxi"]]         <- Valphaxi
    result[["Valphaxi.inverse"]] <- chol2inv( t.vchol )
    result[["log.det.Valphaxi"]] <- 2. * sum( log( diag( t.vchol) ) )

    attr( result[["Valphaxi"]], "struc" )         <- "sym"
    attr( result[["Valphaxi.inverse"]], "struc" ) <- "sym"

  }

  result
}


