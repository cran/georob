##  ###########################################################################

control.predict.georob <-
function(
  full.covmat = FALSE, extended.output = FALSE,
  mmax = 10000,  ncores = pcmp[["max.ncores"]],
  pwidth = NULL, pheight = NULL, napp = 1,
  pcmp = control.pcmp()
){

  ## auxiliary function to set meaningful default values for predict.georob

  ## 2014-07-29 A. Papritz
  ## 2016-07-20 AP changes for parallel computations
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check validity of arguments

  stopifnot(identical(length(full.covmat), 1L)     && is.logical(full.covmat))
  stopifnot(identical(length(extended.output), 1L) && is.logical(extended.output))

  stopifnot(identical(length(mmax), 1L)   && is.numeric(mmax)   && mmax >= 1)
  stopifnot(identical(length(ncores), 1L) && is.numeric(ncores) && ncores >= 1)
  stopifnot(identical(length(napp), 1L)   && is.numeric(napp)   && napp >= 1)

  stopifnot(is.null(pwidth)  || (identical(length(pwidth), 1L)  && is.numeric(pwidth)  && pwidth > 0))
  stopifnot(is.null(pheight) || (identical(length(pheight), 1L) && is.numeric(pheight) && pheight > 0))

  stopifnot(is.list(pcmp))

  list(
    full.covmat = full.covmat,
    extended.output = extended.output,
    mmax = mmax, ncores = ncores,
    pwidth = pwidth, pheight = pheight, napp = napp,
    pcmp = pcmp

  )

}


##  ###########################################################################
### predict.georob

predict.georob <-
function(
  object, newdata,
  type = c( "signal", "response", "trend", "terms" ),
  terms = NULL, se.fit = TRUE,
  signif = 0.95,
  locations,
  variogram.model = NULL, param = NULL, aniso = NULL,
  variogram.object = NULL,
  control = control.predict.georob(),
  verbose = 0, ...
)
{

  ## ToDos:

  ## - try fuer kritische Berechungen
  ## - Anpassung fuer Matrix Package

  ## Given a fitted georob object, the function computes either the trend
  ## or (robust) kriging predictions of the signal or the observations for
  ## newdata or extracts the fitted trend, the trend terms, the signal or
  ## the observations for the support locations if newdata is not given.
  ## both point or block predictions are computed if newdata is specified.

  ## 2011-10-07 A. Papritz
  ## 2012-01-03 AP modified for replicated observations and for parallel processing
  ## 2012-01-05 AP modified for compress storage of matrices
  ## 2012-02-07 AP modified for geometrically anisotropic variograms
  ## 2012-03-02 AP eliminated possibility for logging to file in parallel processing
  ## 2012-03-19 AP correction of error in parallel processing on Windows
  ## 2012-03-28 AP correction of error when processing newdata with NAs
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2012-10-18 AP changes for new definition of eta
  ## 2012-11-04 AP handling compressed cov.betahat
  ## 2012-11-30 AP use of SpatialGridDataFrame and SpatialPixelDataFrame for newdata
  ## 2013-01-19 AP correction of error in computing lag distance matrix between support
  ##               and prediction points
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-05-23 AP correct handling of missing observations
  ## 2013-05-06 AP changes for solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-02-18 AP correcting error in computing predictions for model with offset
  ## 2014-04-23 AP correcting error when computing predictions for data locations
  ## 2014-08-18 AP new argument control for passing tuning parameters to function
  ## 2014-08-18 AP changes for parallelized computations
  ## 2015-03-04 AP some changes for reducing computation effort
  ## 2015-06-23 AP modifications for robust prediction of response
  ## 2015-07-20 AP inactivation of modifications for robust prediction of response
  ##               (variables: rp.response, se.signal, scld.res, resscl)
  ## 2015-08-27 AP correcting error in processing output
  ## 2015-11-30 AP catching errors occurring during parallel computations
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-07-22 AP SpatialPoints, SpatialPixels and SpatialGrid as newdata objects
  ## 2016-07-27 AP correcting error when newdata is SpatialPoints, ...  objects
  ## 2016-08-11 AP changes for nested variogram models
  ## 2016-11-28 AP correcting error when computing predictions for intrinsic variograms
  ## 2017-10-24 AP error message if names of coordinates in object and newdata are not the same
  ## 2017-12-22 AP improved memory management in parallel computations
  ## 2018-01-22 AP optional specification of new variogram model
  ## 2019-12-13 AP correcting use of class() in if() and switch()
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ##  ##############################################################################

  RFoptions(newAniso=FALSE)

#### -- check arguments

  ## match arguments

  type <- match.arg( type )

  ## check validity of arguments

  if(!missing(newdata)) check.newdata(newdata)

  stopifnot(identical(length(se.fit), 1L) && is.logical(se.fit))

  stopifnot(is.null(signif) || (identical(length(signif), 1L) && is.numeric(signif) && signif > 0 && signif < 1))
  stopifnot(identical(length(verbose), 1L) && is.numeric(verbose) && verbose >= 0)

  stopifnot(is.null(param) || is.numeric(param))
  stopifnot(is.null(aniso) || is.numeric(aniso))

  stopifnot(is.list(control))
  stopifnot(is.null(variogram.object) || is.list(variogram.object))

  stopifnot(is.null(terms) || is.character(terms))
  stopifnot(is.null(variogram.model) || is.character(variogram.model))

#### -- setup or check contents of variogram.object

  if( !all(
      is.null(variogram.model), is.null(param), is.null(aniso),
      is.null(variogram.object)
    )
  ){

    cl <- object[["call"]]

    if( is.null( variogram.object ) ){

      ## either variogram.model, param, or aniso have been specified

      if( "variogram.object" %in% names(cl) ) stop(
        "variogram parameters were specified for 'object' as argument 'variogram.object'",
        " --- specifiy new variogram model in the same way"
      )

      if( !is.null(variogram.model) ){
        variogram.model <- match.arg(
          variogram.model,
          c( "RMexp", "RMaskey", "RMbessel", "RMcauchy",
            "RMcircular", "RMcubic", "RMdagum", "RMdampedcos", "RMdewijsian", "RMfbm",
            "RMgauss", "RMgencauchy", "RMgenfbm", "RMgengneiting", "RMgneiting", "RMlgd",
            "RMmatern", "RMpenta", "RMqexp", "RMspheric", "RMstable",
            "RMwave", "RMwhittle"
          )
        )
      } else variogram.model <- object[["variogram.object"]][[1]][["variogram.model"]]

      ## match names of param, aniso

      if( !is.null(param) ){
        if( is.numeric(param) ){
          tmp <- names( param )
          tmp <- sapply(tmp, function(x, choices){
              match.arg(x, choices)
            },
            choices = names( default.fit.param() )
          )
          names( param ) <- nme.param <- tmp
        } else stop( "'param' must be a numeric vector" )
      } else param <- object[["variogram.object"]][[1]][["param"]]

      fit.param <- default.fit.param( variance = FALSE, nugget = FALSE, scale = FALSE )

      if( !is.null(aniso) ){
        if( is.numeric(aniso) ){
          tmp <- names( aniso )
          tmp <- sapply(tmp, function(x, choices) match.arg(x, choices),
            choices = names( default.aniso() )
          )
          names( aniso ) <- tmp
        } else stop( "'aniso' must be a numeric vector" )
      } else aniso <- object[["variogram.object"]][[1]][["aniso"]]

      fit.aniso <- default.fit.aniso()

      ## create variogram.object

      variogram.object <- list(
        list(
          variogram.model = variogram.model,
          param = param, fit.param = fit.param,
          aniso = aniso, fit.aniso = fit.aniso
        )
      )

    } else {

      ## variogram.object has been passed to function, check contents

      if( !"variogram.object" %in% names(cl) ) stop(
        "variogram parameters were not specified for 'object' as argument 'variogram.object'",
        " --- specifiy new variogram model in the same way"
      )

      variogram.object <- lapply(
        variogram.object,
        function( y ){

          variogram.model <- y[["variogram.model"]]
          param      <- y[["param"]]
          aniso      <- y[["aniso"]]

          if( !is.null(variogram.model) ){
            y[["variogram.model"]] <- match.arg(
              variogram.model,
              c( "RMexp", "RMaskey", "RMbessel", "RMcauchy",
                "RMcircular", "RMcubic", "RMdagum", "RMdampedcos", "RMdewijsian", "RMfbm",
                "RMgauss", "RMgencauchy", "RMgenfbm", "RMgengneiting", "RMgneiting", "RMlgd",
                "RMmatern", "RMpenta", "RMqexp", "RMspheric", "RMstable",
                "RMwave", "RMwhittle"
              )
            )
          } else stop( "component 'variogram.model' missing in 'variogram.object'")

          ## match names of components

          nme.param <- NULL

          if( !is.null(param) ){
            if( is.numeric(param) ){
              tmp <- names( param )
              tmp <- sapply(tmp, function(x, choices){
                  match.arg(x, choices)
                },
                choices = names( default.fit.param() )
              )
              names( param ) <- nme.param <- tmp
            } else stop( "'param' must be a numeric vector" )
            y[["param"]] <- param
          } else stop( "component 'param' missing in 'variogram.object'")

          y[["fit.param"]] <- default.fit.param( variance = FALSE, nugget = FALSE, scale = FALSE )

          if( !is.null(aniso) ){
            if( is.numeric(aniso) ) {
              tmp <- names( aniso )
              tmp <- sapply(tmp, function(x, choices) match.arg(x, choices),
                choices = names( default.aniso() )
              )
              names( aniso ) <- tmp
            } else stop( "'aniso' must be a numeric vector" )
            y[["aniso"]] <- aniso
          } else y[["aniso"]] <- default.aniso()

          y[["fit.aniso"]] <- default.fit.aniso()

          y

        }
      )
    }

    ## replace variogram.object in object

    object[["variogram.object"]] <- variogram.object

    ## set all initial values to specified parameter values

    object[["call"]] <- f.call.set_allxxx_to_fitted_values( object )

    ## fix all variogram parameters

    cl <- object[["call"]]

    cl <- f.call.set_allfitxxx_to_false( cl )

    ## set hessian equal to FALSE and avoid computation of unneeded
    ## covariance matrices

    cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "hessian", FALSE )
    cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "cov.bhat", FALSE )
    cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "cov.ehat", FALSE )
    cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "cov.ehat.p.bhat", FALSE )

    ## set verbose argument

    cl <- f.call.set_x_to_value( cl, "verbose", 0 )

    ## update call in object

    object[["call"]] <- cl

    ## re-fit object

    object <- update( object )

  }

#### -- further preparations

  ## expand matrices

  Valphaxi.objects <- expand( object[["Valphaxi.objects"]] )
  zhat.objects     <- expand( object[["zhat.objects"]] )
  object[["cov"]]  <- expand( object[["cov"]]  )

  if( missing( locations ) ){
    locations <- object[["locations.objects"]][["locations"]]
  }

  ## check the consistency of the provided arguments

  if( !missing( newdata ) && class( newdata )[1] == "SpatialPolygonsDataFrame" ){

    ## check whether pwidth and pheight were provided

    if( is.null( control[["pwidth"]] ) || is.null( control[["pheight"]] ) ) stop(
      "'pwidth' and 'pheight' must be provided for block kriging"
    )

    ## map names of variogram models of RandomFields version 3 to version 2

    object[["variogram.object"]] <- lapply(
      object[["variogram.object"]],
      function(x){

        variogram.model <- x[["variogram.model"]]
        isotropic <- x[["isotropic"]]
        param <- x[["param"]]

        if( variogram.model %in% control.georob()[["irf.models"]] ) stop(
          "block kriging not yet implemented for unbounded variogram models"
        )
        if( !isotropic ) stop(
          "block kriging not yet implemented for anisotropic variograms"
        )

        variogram.model.v2 <- gsub("^RM", "", variogram.model )

        variogram.model.v2 <- switch(
          variogram.model.v2[1],
          askey = stop(
            "variogram model 'RMaskey' not implemented in package constrainedKriging"
          ),
          dagum = stop(
            "variogram model 'RMdagum' not implemented in package constrainedKriging"
          ),
          dewijsian = stop(
            "variogram model 'RMdewijsian' not implemented in package constrainedKriging"
          ),
          fbm = stop(
            "variogram model 'RMfbm' not implemented in package constrainedKriging"
          ),
          genfbm = stop(
            "variogram model 'RMgenfbm' not implemented in package constrainedKriging"
          ),
          dampedcos = "dampedcosine",
          exp = "exponential",
          lgd = "lgd1",
          qexp = "qexponential",
          spheric = "spherical",
          variogram.model.v2
        )

        if( identical( variogram.model.v2, "gengneiting" ) ) param[6L] <- sum( param[5L:6L] ) + 0.5

        x[["variogram.model.v2"]] <- variogram.model.v2
        x[["param"]] <- param

        x

      }
    )

  }

  if( control[["full.covmat"]] ){
    if( verbose > 0. ){
      cat(
        "\ncomputing full covariance matrix of prediction errors\n"
      )
      if(
        !missing( newdata ) && class( newdata )[1] == "SpatialPolygonsDataFrame"
      ) cat(
        "requires some computing time for block kriging, be patient ...\n"
      )
    }
  }

  if( identical( type, "terms" ) &&
    !( missing( newdata ) || is.null( newdata )) ) stop(
    "predicting terms for newdata not yet implemented"
  )

  if( !missing( newdata ) && class( newdata )[1] %in% c( "SpatialPoints", "SpatialPixels", "SpatialGrid" ) ){
    t.formula <- as.formula( paste( as.character( formula( object ) )[-2L], collapse = "" ) )
    tmp <- try(
      get_all_vars( t.formula, as.data.frame( coordinates( newdata ) ) ),
      silent = TRUE
    )
    if( identical( class( tmp ), "try-error" ) ) stop(
      "'newdata' is a SpatialPoints, SpatialPixels or SpatialGrid object\n but drift covariates are not functions of coordinates"
    )
  }

  if( control[["extended.output"]] &&
    any( sapply(object[["variogram.object"]], function(x) x[["variogram.model"]])
      %in% control.georob()[["irf.models"]] )
  ) stop( "extended output cannot be computed for intrinsic variogram models" )

  ## extract fixed effects terms of object

  tt <- terms( object )

  ## extract fixed effects design matrix of support data

  X <- model.matrix(
    tt,
    model.frame( object )
  )
  attr.assign <- attr( X, "assign" )
  X <- X[!duplicated( object[["Tmat"]] ), , drop = FALSE]
  attr( X, "assign" ) <- attr.assign

  n <- length( object[["bhat"]] )

  ## extract the coordinates of the support locations

  locations.coords <-
    object[["locations.objects"]][["coordinates"]][!duplicated( object[["Tmat"]] ), , drop = FALSE]

  #   ## extract residuals if robust predictions of response for newdata are computed
  #
  #   scld.res <- NULL
  #   se.signal <- NULL
  #   rp.response <- FALSE
  #
  #   if(
  #     object[["tuning.psi"]] < object[["control"]][["tuning.psi.nr"]] &&
  #     identical( type, "response" )
  #   ){
  #     if( control[["extended.output"]] ) warning(
  #       "variances of prediction targets (response) are underestimated"
  #     )
  #     if( !( missing( newdata ) || is.null( newdata ) ) ){
  #       rp.response <- FALSE
  #       resscl <- 1.
  #       warning( "scale factor for computing empirical distribution of residuals equals 1" )
  #       scld.res <- object[["residuals"]] / resscl
  #     }
  #   }

#### -- compute missing covariance matrices

  cov.betahat    <- is.null( object[["cov"]][["cov.betahat"]] )
  cov.delta.bhat   <- is.null( object[["cov"]][["cov.delta.bhat"]] ) ||
    !is.matrix( object[["cov"]][["cov.delta.bhat"]] )
  cov.delta.bhat.betahat <- is.null( object[["cov"]][["cov.delta.bhat.betahat"]] )

  cov.bhat    <- control[["extended.output"]] & (
    is.null( object[["cov"]][["cov.bhat"]] ) || !is.matrix( object[["cov"]][["cov.bhat"]] )
  )
  cov.bhat.betahat  <-  control[["extended.output"]] & is.null( object[["cov"]][["cov.bhat.betahat"]] )
  cov.p.t  <-  control[["extended.output"]] & is.null( object[["cov"]][["cov.pred.target"]] )

  if( any( c( cov.betahat, cov.delta.bhat, cov.delta.bhat.betahat,
        control[["extended.output"]] & ( cov.bhat || cov.bhat.betahat || cov.p.t )
      )
    )
  ){ ## cov

    r.cov <- covariances.fixed.random.effects(
      Valphaxi.objects = Valphaxi.objects[c("Valphaxi", "Valphaxi.inverse")],
      Aalphaxi = zhat.objects[["Aalphaxi"]],
      Palphaxi = zhat.objects[["Palphaxi"]],
      Valphaxi.inverse.Palphaxi = zhat.objects[["Valphaxi.inverse.Palphaxi"]],
      rweights = object[["rweights"]],
      XX = X, TT = object[["Tmat"]], TtT = as.vector( table( object[["Tmat"]] ) ),
      names.yy = rownames( model.frame( object ) ),
      nugget = object[["variogram.object"]][[1L]][["param"]]["nugget"],
      eta =  f.reparam.fwd( object[["variogram.object"]] )[[1L]][["param"]]["nugget"],
      expectations = object[["expectations"]], family = "gaussian",
      cov.bhat = cov.bhat, full.cov.bhat = cov.bhat,
      cov.betahat = cov.betahat,
      cov.bhat.betahat = cov.bhat.betahat,
      cov.delta.bhat = cov.delta.bhat, full.cov.delta.bhat = cov.delta.bhat,
      cov.delta.bhat.betahat = cov.delta.bhat.betahat,
      cov.ehat = FALSE, full.cov.ehat = FALSE,
      cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
      aux.cov.pred.target = cov.p.t,
      control.pcmp = control[["pcmp"]],
      verbose = verbose
    )

    if( r.cov[["error"]] ) stop(
      "an error occurred when computing the covariances of fixed and random effects",
    )

    if( is.null( object[["cov"]] ) ) object[["cov"]] <- list()

    if( cov.betahat )     object[["cov"]][["cov.betahat"]]    <- r.cov[["cov.betahat"]]
    if( cov.delta.bhat )  object[["cov"]][["cov.delta.bhat"]] <- r.cov[["cov.delta.bhat"]]
    if( cov.delta.bhat.betahat ) object[["cov"]][["cov.delta.bhat.betahat"]] <-
      r.cov[["cov.delta.bhat.betahat"]]
    if( control[["extended.output"]] && cov.bhat )   object[["cov"]][["cov.bhat"]] <- r.cov[["cov.bhat"]]
    if( control[["extended.output"]] && cov.bhat.betahat ) object[["cov"]][["cov.bhat.betahat"]] <-
      r.cov[["cov.bhat.betahat"]]
    if( control[["extended.output"]] && cov.p.t ) object[["cov"]][["cov.pred.target"]] <-
      r.cov[["cov.pred.target"]]

  } ## end cov

  ## compute lower cholesky factor of covariance matrix of delta = (b -
  ## bhat) and betahat - beta

  cov.delta.bhat.betahat.l <- try(
    t(
      chol(
        rbind(
          cbind(
            object[["cov"]][["cov.delta.bhat"]],
            object[["cov"]][["cov.delta.bhat.betahat"]]
          ),
          cbind(
            t( object[["cov"]][["cov.delta.bhat.betahat"]] ),
            object[["cov"]][["cov.betahat"]]
          )
        )
      )
    ), silent = TRUE
  )
  if( identical( class( cov.delta.bhat.betahat.l ), "try-error" ) ) stop(
    "covariance matrix of kriging errors 'b-bhat' and 'betahat' not positive definite"
  )

  ## compute covariance matrix of betahat and bhat for extended output

  cov.betahat.l <- try( t( chol( object[["cov"]][["cov.betahat"]] ) ) )
  if( identical( class( cov.betahat.l ), "try-error" ) ) stop(
    "covariance matrix of 'betahat' not positive definite"
  )

  if( control[["extended.output"]] ){

    ## compute covariance matrix of bhat and betahat

    cov.bhat.betahat <-  rbind(
      cbind(
        object[["cov"]][["cov.bhat"]],
        object[["cov"]][["cov.bhat.betahat"]]
      ),
      cbind(
        t( object[["cov"]][["cov.bhat.betahat"]] ),
        object[["cov"]][["cov.betahat"]]
      )
    )

    cov.p.t <- object[["cov"]][["cov.pred.target"]]

  } else {

    cov.bhat.betahat <- NULL
    cov.p.t <- NULL

  }

  ## extract signal variance, xi, nugget and gcr.constant

  tmp <- f.reparam.fwd( object[["variogram.object"]] )

  var.signal <- attr(tmp, "var.signal" )
  xi <- sapply( tmp, function(x) x[["param"]]["variance"] )

  nugget <- object[["variogram.object"]][[1L]][["param"]]["nugget"]

  gcr.constant <- lapply(
    object[["Valphaxi.objects"]][["Valpha"]],
    function(x) x[["gcr.constant"]]
  )

  ##########################

#### -- compute predictions

  if( missing( newdata ) || is.null( newdata ) ){

    ##############

#### --- terms for support locations

    ## no newdata object: compute terms for support locations
    ## code borrowed from predict.lm

    if( identical( type, "terms" ) ){

      beta <- coef( object )
      aa <- attr( X, "assign" )
      ll <- attr ( tt, "term.labels" )
      hasintercept <- attr( tt, "intercept") > 0L
      if (hasintercept) ll <- c( "(Intercept)", ll )
      aaa <- factor( aa, labels = ll )
      asgn <- split( order(aa), aaa )
      if( hasintercept ) {
        asgn$"(Intercept)" <- NULL
        avx <- colMeans( X )
        termsconst <- sum( avx * beta )
      }
      nterms <- length( asgn )

      if( nterms > 0L ){

        predictor <- matrix( ncol = nterms, nrow = length( object[["Tmat"]] ) )
        dimnames( predictor ) <- list(
          rownames( model.frame( object ) ),
          names(asgn)
        )
        if( se.fit ){
          ip <- matrix( ncol = nterms, nrow = length( object[["Tmat"]] ) )
          dimnames( ip ) <- list(
            rownames( model.frame( object ) ),
            names(asgn)
          )
        }

        if (hasintercept)  X <- sweep(X, 2L, avx, check.margin = FALSE)

        for( i in seq.int( 1L, nterms, length.out = nterms) ){

          ii <- asgn[[i]]

          predictor[ , i] <- X[object[["Tmat"]], ii, drop = FALSE] %*% beta[ii]

          if( se.fit ){
            t.cov.betahat.l <- t(
              chol( object[["cov"]][["cov.betahat"]][ ii, ii, drop = FALSE] )
            )
            ip[ , i] <- rowSums(
              ( X[object[["Tmat"]], ii, drop = FALSE] %*% t.cov.betahat.l )^2
            )
          }

        }

        if( !is.null( terms ) ){
          predictor <- predictor[ , terms, drop = FALSE]
          if( se.fit ) ip <- ip[ , terms, drop = FALSE]
        }

      } else {
        predictor <- ip <- matrix(0., length( object[["Tmat"]] ), 0L)
      }

      attr( predictor, "constant" ) <- if( hasintercept )  termsconst  else  0

      if( se.fit ){
        se <- sqrt(ip)
        if (type == "terms" && !is.null(terms))
        se <- se[, terms, drop = FALSE]
      }
      if( missing(newdata) && !is.null(na.act <- object[["na.action"]] ) ){
        predictor <- napredict( na.act, predictor )
        if (se.fit) se <- napredict( na.act, se )
      }
      result <- if( se.fit ){
        list(
          fit = predictor, se.fit = se,
          df = object[["df.residual"]],
          residual.scale = unname( sqrt( nugget ) )
        )
      } else {
        predictor
      }
      ## end "terms"

    } else {

      ##############

#### --- predictions for support locations

      ## no newdata object: compute predictions for support locations

      ## compute predictions

      pred <- switch(
        type,
        "response" = model.response( model.frame( object ) ),
        "signal"   = object[["fitted.values"]] + object[["bhat"]][object[["Tmat"]]],
        "trend"    = object[["fitted.values"]]
      )

#       var.pred        <- NULL
#       var.target      <- NULL
#       cov.pred.target <- NULL

      ## compute covariance matrix of signal

      if( control[["extended.output"]] ){
        V <- var.signal * Valphaxi.objects[["Valphaxi"]]
      }


      ## compute MSEP and (co-)variances of targets and predictions

      t.result <- switch(
        type,
        response = {  ## response
          mse.pred <- rep( 0., length( object[["Tmat"]] ) )
          if( control[["extended.output"]] ){
            var.pred <- var.target <- cov.pred.target <- rep(
              var.signal * Valphaxi.objects[["Valphaxi"]][1,1] + nugget,
              length( object[["Tmat"]] )
            )
          }
          if( control[["full.covmat"]] ){
            mse.pred <- diag( mse.pred )
            if( control[["extended.output"]] ){
              var.pred <- V
              diag( var.pred ) <- diag( var.pred ) + nugget
              var.pred <- var.target <- cov.pred.target <- var.pred[object[["Tmat"]], object[["Tmat"]]]
            }
          }
          c(
            list( mse.pred = mse.pred ),
            if( control[["extended.output"]] ) list(
              var.pred = var.pred, var.target = var.target, cov.pred.target = cov.pred.target
            )
          )
        },
        signal = {    ## signal
          aux <- cbind(
            cov.delta.bhat.betahat.l[1L:n,1L:n] - X %*% cov.delta.bhat.betahat.l[-(1L:n),1L:n],
            - X  %*% cov.delta.bhat.betahat.l[-(1L:n),-(1L:n)]
          )
          aux <- aux[object[["Tmat"]], , drop = FALSE]
          if( control[["full.covmat"]] ){
            mse.pred <- tcrossprod( aux )
          } else {
            mse.pred <- rowSums( aux^2 )
          }
          if( control[["extended.output"]] ){
            aux <- cov.bhat.betahat[1L:n, -(1L:n), drop = FALSE] %*% t(X)
            var.pred <- cov.bhat.betahat[1L:n, 1L:n, drop = FALSE] + aux + t(aux) +
              X %*% cov.bhat.betahat[-(1L:n), -(1L:n), drop = FALSE] %*% t(X)
            var.pred <- var.pred[object[["Tmat"]], object[["Tmat"]]]
            var.target <- V[object[["Tmat"]], object[["Tmat"]]]
            cov.pred.target <- pmm(
              (cov.p.t[1L:n,] + X %*% cov.p.t[-(1L:n),]),
              V, control = control[["pcmp"]]
            )
            cov.pred.target <- cov.pred.target[object[["Tmat"]], object[["Tmat"]]]
            if( !control[["full.covmat"]] ){
              var.pred <- diag( var.pred )
              var.target <- diag( var.target )
              cov.pred.target <- diag( cov.pred.target )
            }
          }
          c(
            list( mse.pred = mse.pred ),
            if( control[["extended.output"]] ) list(
              var.pred = var.pred, var.target = var.target, cov.pred.target = cov.pred.target
            )
          )
        },
        trend = {     ## trend
          aux <- X %*% cov.betahat.l
          aux <- aux[object[["Tmat"]], , drop = FALSE]
          if( control[["full.covmat"]] ){
            mse.pred <- matrix( NA_real_, length( object[["Tmat"]] ), length( object[["Tmat"]] ) )
            var.pred <- tcrossprod( aux )
          } else {
            mse.pred <- rep( NA_real_, length( object[["Tmat"]] ) )
            var.pred <- rowSums( aux^2 )
          }
          list( mse.pred = mse.pred, var.pred = var.pred )
        }
      )

      ## collect results

      pred.se <-  sqrt( f.diag( t.result[["mse.pred"]] ) )

      result <- data.frame(
        pred = pred,
        se = pred.se
      )
      if( !is.null( signif ) ){
        result <- cbind(
          result,
          data.frame(
            lower = pred + qnorm( 0.5 * ( 1.-signif[1L] ) ) * pred.se,
            upper = pred + qnorm( 0.5 * ( 1.+signif[1L] ) ) * pred.se
          )
        )
      }
      if( control[["extended.output"]] ){
        result <- cbind( result, trend = fitted( object ) )
      }


      if( !is.null(t.result[["var.pred"]]) ){
        result[["var.pred"]] <- f.diag( t.result[["var.pred"]] )
      }
      if( !is.null(t.result[["cov.pred.target"]]) ){
        result[["cov.pred.target"]]<- f.diag( t.result[["cov.pred.target"]] )
      }
      if( !is.null(t.result[["var.target"]]) ){
        result[["var.target"]]<- f.diag( t.result[["var.target"]] )
      }

      if( identical( type, "trend" ) ) result <- result[, c( "pred", "var.pred" )]

      result <- as.data.frame( napredict( object[["na.action"]], as.matrix( result ) ) )

      if( control[["full.covmat"]] ){

        result <- c(
          list(
            pred = result,
            mse.pred = napredict(
              object[["na.action"]], t( napredict( object[["na.action"]], mse.pred ) ) )
          ),
          if( !is.null(t.result[["var.pred"]]) ){
            var.pred = napredict(
              object[["na.action"]], t( napredict( object[["na.action"]], t(t.result[["var.pred"]]) ) )
            )
            dimnames( var.pred ) <- NULL
            list( var.pred = var.pred )
          },
          if( !is.null(t.result[["cov.pred.target"]]) ){
            cov.pred.target = napredict(
              object[["na.action"]], t( napredict( object[["na.action"]], t(t.result[["cov.pred.target"]]) ) )
            )
            dimnames( cov.pred.target ) <- NULL
            list( cov.pred.target = cov.pred.target )
          },
          if( !is.null(t.result[["var.target"]]) ){
            var.target = napredict(
              object[["na.action"]], t( napredict( object[["na.action"]], t(t.result[["var.target"]]) ) )
            )
            dimnames( var.target ) <- NULL
            list( var.target = var.target )
          }
        )
        dimnames( result[["mse.pred"]] ) <- NULL

        if( identical( type, "trend" ) ) result <- result[c( "pred", "var.pred" )]
      }


    }
    ## end no newdata object

  } else {

    ##############

#### --- predictions for new locations

    ## compute predictions for newdata

    Terms <- delete.response( tt )
    Terms.loc <- terms( locations )
    attr( Terms.loc, "intercept" ) <- 0

    ## get the model frame for newdata

    mf.newdata <- switch(
      class( newdata )[1],
      "data.frame" = model.frame(
        Terms, newdata, na.action = na.pass, xlev = object[["xlevels"]]
      ),
      "SpatialPoints" = ,
      "SpatialPixels" = ,
      "SpatialGrid" = model.frame(
        Terms, as.data.frame( coordinates( newdata ) ), na.action = na.pass,
        xlev = object[["xlevels"]]
      ),
      "SpatialPointsDataFrame" = ,
      "SpatialPixelsDataFrame" = ,
      "SpatialGridDataFrame" = ,
      "SpatialPolygonsDataFrame" = model.frame(
        Terms, slot( newdata, "data" ), na.action = na.pass,
        xlev = object[["xlevels"]]
      ),
      stop(
        "cannot construct model frame for class(newdata) ='",
        class( newdata )
      )
    )

    ## check whether variables that will be used to compute the
    ## predictions agree with those in object

    if( !is.null( cl <- attr(Terms, "dataClasses" ) ) )
    .checkMFClasses( cl, mf.newdata )

    ## get fixed effects design matrix for newdata

    pred.X <- model.matrix( Terms, mf.newdata,
      contrasts.arg = object[["contrasts"]] )

    ## deal with non-NULL offset

    offset <- rep( 0., NROW(pred.X) )
    if( !is.null( off.num <- attr( tt, "offset" ) ) ){
      for( i in off.num ) {
        offset <- offset + eval( attr( tt, "variables" )[[i + 1L]], newdata )
      }
    }
    if( !is.null( object[["call"]][["offset"]] ) ){
      offset <- offset + eval( object[["call"]][["offset"]], newdata )
    }

    ## get matrix of coordinates of newdata for point kriging

    pred.coords <- switch(
      class( newdata )[1],
      "data.frame" = model.matrix(
        Terms.loc,
        model.frame(
          Terms.loc, newdata, na.action = na.pass
        )
      ),
      "SpatialPoints" = ,
      "SpatialPixels" = ,
      "SpatialGrid" = ,
      "SpatialPointsDataFrame" = ,
      "SpatialPixelsDataFrame" = ,
      "SpatialGridDataFrame" = coordinates( newdata ),
      "SpatialPolygonsDataFrame" = NULL
    )

    if( !is.null( pred.coords ) &&
      !(
        NCOL( locations.coords ) == NCOL( pred.coords ) &&
        all( colnames( pred.coords ) ==
          colnames(object[["locations.objects"]][["coordinates"]]) )
      )
    ) stop(
      "inconsistent number and/or names of coordinates in 'object' and in 'newdata'"
    )

    ## number of items to predict

    m <- NROW( newdata )

    ## determine number of prediction parts

    n.part <- ceiling( m / control[["mmax"]] )
    rs <- ( 0L:(n.part-1L)) * control[["mmax"]] + 1L
    re <- ( 1L:(n.part  )) * control[["mmax"]]; re[n.part] <- m

    ncores <- min( n.part, control[["ncores"]] )

    ncores.available <- control[["pcmp"]][["max.ncores"]]
    if( sfIsRunning() ) sfStop()

    control.pcmp <- control[["pcmp"]]
    control.pcmp[["pmm.ncores"]] <- min(
      control.pcmp[["pmm.ncores"]],
      max( 1L, floor( (ncores.available - ncores) / ncores ) )
    )
    if( ncores > 1L && !control.pcmp[["allow.recursive"]] ) control.pcmp[["pmm.ncores"]] <- 1L

    if( control[["full.covmat"]] && n.part > 1L ) stop(
      "full covariance matrix of prediction errors cannot ",
      "be computed\n  if prediction problem is split into several parts\n",
      "-> increase 'mmax' to avoid splitting"
    )


    ## handle parallel processing

    ## auxiliary function to compute the predictions for one part

    f.aux <- function(
      i
#       rs, re, n.part,
#       type,
#       locations.coords, betahat, bhat, response,
#       #       rp.response, scld.res,
#       pred.X, pred.coords, offset, newdata, mf.newdata,
#       variogram.object, var.signal, xi, nugget, gcr.constant,
#       cov.delta.bhat.betahat.l, cov.betahat.l, cov.bhat.betahat, cov.p.t,
#       Valphaxi, Valphaxi.inverse,
#       pwidth, pheight, napp,
#       signif,
#       extended.output, full.covmat,
#       control.pcmp,
#       verbose
    ){

      ## objects
      ##
      ## rs, re, n.part, type, locations.coords,
      ## betahat, bhat, response, variogram.object,
      ## # rp.response, scld.res,
      ## pred.X, pred.coords, offset,
      ## newdata, mf.newdata, var.signal, xi, nugget,
      ## gcr.constant, cov.delta.bhat.betahat.l, cov.betahat.l,
      ## cov.bhat.betahat, cov.p.t,
      ## Valphaxi, Valphaxi.inverse,
      ## pwidth, pheight, napp, extended.output, full.covmat,
      ## signif, control.pcmp, verbose
      ##
      ## are taken from parent enviroment

      if( verbose > 0. )
      cat( "  predicting part ", i, " of ", n.part, "\n" )

      ## select the data for the current part

      pred.X <- pred.X[ rs[i]:re[i], , drop = FALSE]
      offset <- offset[ rs[i]:re[i] ]
      if( class(newdata)[1] %in% c( "SpatialPoints", "SpatialPixels", "SpatialGrid" ) ){
        newdata <- mf.newdata[ rs[i]:re[i], ]
      } else{
        newdata <- newdata[ rs[i]:re[i], ]
      }

      if( !is.null( pred.coords ) ) {
        pred.coords <- pred.coords[ rs[i]:re[i], , drop = FALSE]
      }

      ## compute the predictions for the current part

      result <- f.robust.uk(
        type = type, terms = terms,
        locations.coords = locations.coords,
        betahat = betahat,
        bhat = bhat,
        response = response,
        #         rp.response = rp.response,
        #         scld.res = scld.res,
        pred.X = pred.X, pred.coords = pred.coords, offset = offset, newdata = newdata,
        variogram.object = variogram.object, var.signal = var.signal,
        xi = xi, nugget = nugget, gcr.constant = gcr.constant,
        cov.delta.bhat.betahat.l = cov.delta.bhat.betahat.l,
        cov.betahat.l = cov.betahat.l,
        cov.bhat.betahat = cov.bhat.betahat,
        cov.p.t = cov.p.t,
        Valphaxi = Valphaxi, Valphaxi.inverse = Valphaxi.inverse,
        pwidth = pwidth, pheight = pheight, napp = napp,
        signif = signif,
        extended.output = extended.output,
        full.covmat = full.covmat,
        control.pcmp = control.pcmp,
        verbose = verbose
      )

      return( result )
    }

    ## prepare items to pass to function f.aux

    betahat           <- object[["coefficients"]]
    bhat              <- object[["bhat"]]
    response          <- model.response( model.frame( object ) )
    variogram.object  <- object[["variogram.object"]]
    Valphaxi          <- Valphaxi.objects[["Valphaxi"]]
    Valphaxi.inverse  <- Valphaxi.objects[["Valphaxi.inverse"]]
    pwidth            <- control[["pwidth"]]
    pheight           <- control[["pheight"]]
    napp              <- control[["napp"]]
    extended.output   <- control[["extended.output"]]
    full.covmat       <- control[["full.covmat"]]

  ## set default value for control of forking if missing (required for backward compatibility)

  if( is.null( control[["pcmp"]][["fork"]] ) ){
    control[["pcmp"]][["fork"]] <- !identical( .Platform[["OS.type"]], "windows" )
  }
    ## compute the predictions for all the parts

    if( ncores > 1L && !control[["pcmp"]][["fork"]] ){

      ## create a SNOW cluster on windows OS

      clstr <- makeCluster( ncores, type = "SOCK")
      save( clstr, file = "SOCKcluster.RData" )
      options( error = f.stop.cluster )

      ## export required items to workers

      junk <- clusterEvalQ( clstr, require( georob, quietly = TRUE ) )

      junk <- clusterExport(
        clstr,
        c( "rs", "re", "n.part",
          "type",
          "locations.coords",
          "betahat", "bhat", "response", "variogram.object",
          #       "rp.response", "scld.res",
          "pred.X", "pred.coords", "offset", "newdata", "mf.newdata",
          "var.signal", "xi", "nugget", "gcr.constant",
          "cov.delta.bhat.betahat.l", "cov.betahat.l", "cov.bhat.betahat", "cov.p.t",
          "Valphaxi", "Valphaxi.inverse",
          "pwidth", "pheight", "napp", "extended.output", "full.covmat",
          "signif",
          "control.pcmp",
          "verbose" ),
        envir = environment()
      )


      t.result <- try(
        parLapply(
          clstr,
          1L:n.part,
          f.aux
        )
      )

      f.stop.cluster( clstr )

      #         junk <- parLapply( clstr, 1L:length(clstr), function( i ) sfStop() )
      #         junk <- stopCluster( clstr )
      #         if( file.exists( "SOCKcluster.RData" ) ){
      #           file.remove( "SOCKcluster.RData" )
      #         }
      #         options( error = NULL )

    } else {

      ## fork child processes on non-windows OS

      t.result <- try(
        mclapply(
          1L:n.part,
          f.aux,
          mc.cores = ncores,
          mc.allow.recursive = control.pcmp[["allow.recursive"]]
        )
      )

    }

    has.error <- sapply(
      t.result, function( x ) identical( class(x), "try-error" )
    )

    if( any( has.error ) ){
      cat( "\nerror(s) occurred when computing kriging predictions in parallel:\n\n" )
      sapply( t.result[has.error], cat)
      cat( "\nuse 'ncores=1' and 'verbose = 1' to avoid parallel computations and to see where problem occurs\n\n" )
      stop()
    }

    ## delete items
    rm(
      betahat, bhat, response, variogram.object, Valphaxi,
      Valphaxi.inverse, pwidth, pheight, napp, extended.output, full.covmat
    )

    ## collect results of the various parts into a single list

    result <- t.result[[1L]]
    if( length( t.result ) > 1L ){
      for( i in 2L:length( t.result ) ) {
        result <- rbind( result, t.result[[i]] )
      }
    }

    ## end compute predictions for newdata

  }


#### -- prepare output

  ## complement kriging result with coordinate information on prediction
  ## targets

  if( missing( newdata ) || is.null( newdata ) ){

    coords <- napredict( object[["na.action"]], object[["locations.objects"]][["coordinates"]] )

    if( !identical( type, "terms" ) ){
      if( control[["full.covmat"]] ){
        result[["pred"]] <- data.frame( coords, result[["pred"]] )
      } else {
        result <- data.frame( coords, result )
      }
    }

  } else {

    t.pred <- if( control[["full.covmat"]] ){
      result[["pred"]]
    } else {
      result
    }

    if( class(newdata)[1] == "SpatialPolygonsDataFrame" ){

      t.pred <- SpatialPolygonsDataFrame(
        Sr = SpatialPolygons( newdata@polygons ),
        data = t.pred
      )

    } else {

      #       sel <- match( "se.signal", colnames( t.pred ) )
      #
      #       if( identical( type, "response" ) ) se.signal <- t.pred[, sel]
      #       t.pred <- t.pred[, -sel]

      t.pred <- data.frame( pred.coords, t.pred )

      if( class( newdata )[1] != "data.frame" ){
        coordinates( t.pred ) <- locations
        if( !( class( newdata )[1] %in% c( "SpatialPoints", "SpatialPointsDataFrame" ) ) ){
          gridded( t.pred ) <- TRUE
          if( !( class( newdata )[1] %in% c( "SpatialPixels", "SpatialPixelsDataFrame" ) ) ){
            fullgrid( t.pred ) <- TRUE
          }
        }
      }

    }

    if( control[["full.covmat"]] ){
      result[["pred"]] <- t.pred
    } else {
      result <- t.pred
    }

  }

  ## set attributes required for back-transformation by lgnpp

  if( !identical( type, "terms" ) ){
    if( control[["full.covmat"]] ){
      if( is.data.frame( result[["pred"]] ) ){
        attr( result[["pred"]], "variogram.object" )      <- object[["variogram.object"]]
        attr( result[["pred"]], "psi.func" )              <- object[["control"]][["psi.func"]]
        attr( result[["pred"]], "tuning.psi" )            <- object[["tuning.psi"]]
        attr( result[["pred"]], "type" )                  <- type
        #         attr( result[["pred"]], "scaled.residuals" )      <- scld.res
        #         attr( result[["pred"]], "se.signal" )             <- se.signal
      } else {
        attr( result[["pred"]]@data, "variogram.object" ) <- object[["variogram.object"]]
        attr( result[["pred"]]@data, "psi.func" )         <- object[["control"]][["psi.func"]]
        attr( result[["pred"]]@data, "tuning.psi" )       <- object[["tuning.psi"]]
        attr( result[["pred"]]@data, "type" )             <- type
        #         attr( result[["pred"]]@data, "scaled.residuals" ) <- scld.res
        #         attr( result[["pred"]]@data, "se.signal" )        <- se.signal
        if( class( result[["pred"]] )[1] == "SpatialPolygonsDataFrame" ){
          attr( result[["pred"]]@data, "coefficients" )   <- object[["coefficients"]]
          attr( result[["pred"]]@data, "terms" )          <- object[["terms"]]
          attr( result[["pred"]]@data, "locations" )      <- object[["locations.objects"]][["locations"]]
        }
      }
    } else {
      if( is.data.frame( result ) ){
        attr( result, "variogram.object" )      <- object[["variogram.object"]]
        attr( result, "psi.func" )              <- object[["control"]][["psi.func"]]
        attr( result, "tuning.psi" )            <- object[["tuning.psi"]]
        attr( result, "type" )                  <- type
        #         attr( result, "scaled.residuals" )      <- scld.res
        #         attr( result, "se.signal" )             <- se.signal
      } else {
        attr( result@data, "variogram.object" ) <- object[["variogram.object"]]
        attr( result@data, "psi.func" )         <- object[["control"]][["psi.func"]]
        attr( result@data, "tuning.psi" )       <- object[["tuning.psi"]]
        attr( result@data, "type" )             <- type
        #         attr( result@data, "scaled.residuals" ) <- scld.res
        #         attr( result@data, "se.signal" )        <- se.signal
        if( class( result )[1] == "SpatialPolygonsDataFrame" ){
          attr( result@data, "coefficients" )   <- object[["coefficients"]]
          attr( result@data, "terms" )          <- object[["terms"]]
          attr( result@data, "locations" )      <- object[["locations.objects"]][["locations"]]
        }
      }
    }
  }

  invisible( result )

}

##  ###########################################################################
### f.robust.uk

## auxiliary function for computing robust kriging predictions for a set
## of prediction targets

f.robust.uk <- function(
  type, terms,
  locations.coords, betahat, bhat, response,
  #   rp.response, scld.res,
  pred.X, pred.coords, offset, newdata,
  variogram.object, var.signal, xi, nugget, gcr.constant,
  cov.delta.bhat.betahat.l, cov.betahat.l, cov.bhat.betahat, cov.p.t,
  Valphaxi, Valphaxi.inverse,
  pwidth, pheight, napp,
  signif,
  extended.output, full.covmat,
  control.pcmp, verbose
){ ## f.robust.uk

  ## function computes robust (or Gaussian) universal point or block
  ## kriging predictions

  ## 2011-07-29 A. Papritz
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2015-06-24 AP modifications for robust prediction of response
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-11 AP changes for nested variogram models
  ## 2016-11-28 AP correcting error when computing predictions for intrinsic variograms

#### -- preparation

  n <- length( bhat )

  ## exclude prediction items with missing information

  ex <- attr( na.omit( pred.X ), "na.action" )

  if( !is.null( pred.coords ) ){
    ex <- unique( c( ex, attr( na.omit( pred.coords ), "na.action" ) ) )
  }

  if( !is.null( ex ) ) {
    ex <- ( 1L:NROW(pred.X) ) %in% sort( ex )
  } else {
    ex <- rep( FALSE, NROW(pred.X) )
  }

  if( any( !ex ) ){

#### -- compute predictions

#### --- trend surface

    t.pred <- t.trend <- drop( pred.X[!ex, , drop = FALSE ] %*% betahat ) + offset[!ex]

    if( !identical( type, "trend" ) ){

      ## compute point or block kriging predictions

      ## get covariance matrix (cov.target) of B at predictons locations and
      ## covariance matrix (gamma) between B at prediction and support
      ## locations

      if( !is.null( pred.coords ) ){

#### --- covariance objects for point kriging

        ## generalized (co-)variance (matrix) of signal at prediction points

        if( full.covmat && NROW( pred.coords[!ex, , drop = FALSE ] ) > 1L ){

          ## lag vectors for all distinct pairs

          if( !all( sapply( variogram.object, function(x) x[["isotropic"]] ) ) ){
            indices.pairs <- combn( NROW( pred.coords[!ex, , drop = FALSE ] ), 2L )
            lag.vectors <- pred.coords[!ex, , drop = FALSE ][ indices.pairs[2L,], ]
                         - pred.coords[!ex, , drop = FALSE ][ indices.pairs[1L,], ]
          } else {
            lag.vectors <- as.vector( dist( pred.coords[!ex, , drop = FALSE ] ) )
          }
          attr( lag.vectors, "ndim.coords" ) <- NCOL(pred.coords[!ex, , drop = FALSE ])

          ##  generalized covariance matrix

          Valpha <- f.aux.gcr(
            lag.vectors = lag.vectors,
            variogram.object = variogram.object,
            gcr.constant = gcr.constant,
            control.pcmp = control.pcmp,
            verbose = verbose
          )

          if( any( sapply( Valpha, function(x) x[["error"]] ) ) ) stop(
            "an error occurred when computing semivariances between prediction points"
          )

          t.var.target <- list(
            diag = var.signal * rowSums(
              sapply(
                1L:length(Valpha),
                function( i, x, xi ){
                  xi[i] * x[[i]][["Valpha"]][["diag"]]
                }, x = Valpha, xi = xi
              )
            ) + ( 1. - sum(xi) ),
            tri = var.signal * rowSums(
              sapply(
                1L:length(Valpha),
                function( i, x, xi ){
                  xi[i] * x[[i]][["Valpha"]][["tri"]]
                }, x = Valpha, xi = xi
              )
            )
          )
          attr( t.var.target, "struc" ) <- "sym"

          t.var.target <- expand( t.var.target )

        } else {

          t.var.target <- var.signal * Valphaxi[1, 1]

        }

        ## generalized covariance matrix between prediction and support
        ## points

        if( !all( sapply( variogram.object, function(x) x[["isotropic"]] ) ) ){
          indices.pairs <- expand.grid(
            1L:NROW( pred.coords[!ex, , drop = FALSE ] ),
            1L:NROW( locations.coords )
          )
          lag.vectors <- (pred.coords[!ex, , drop = FALSE ])[ indices.pairs[, 1L], ] -
            locations.coords[ indices.pairs[, 2L], ]
        } else {
          lag.vectors <- as.vector( rdist( pred.coords[!ex, , drop = FALSE ], locations.coords ) )
        }
        attr( lag.vectors, "ndim.coords" ) <- NCOL(pred.coords[!ex, , drop = FALSE ])

        ## functions of version 3 of RandomFields

        Valpha <- f.aux.gcr(
          lag.vectors = lag.vectors,
          variogram.object = variogram.object,
          gcr.constant = gcr.constant,
          symmetric = FALSE,
          control.pcmp = control.pcmp,
          verbose = verbose
        )

        if( any( sapply( Valpha, function(x) x[["error"]] ) ) ) stop(
          "an error occurred when computing semivariances between support ",
          "and prediction points"
        )

        gamma <- rowSums(
          sapply(
            1L:length(Valpha),
            function( i, x, xi ){
              xi[i] * x[[i]][["Valpha"]]
            }, x = Valpha, xi = xi
          )
        )

        ## add spatial nugget if prediction and support locations coincides

        if( sum(xi) < 1. ){
          if( NCOL(lag.vectors) > 1L ){
            sel <- rowSums(lag.vectors) == 0.
          } else {
            sel <- lag.vectors == 0.
          }
          gamma[sel] <- gamma[sel] + (1. - sum(xi) )
        }

        gamma <- var.signal * matrix( gamma, nrow = NROW( pred.coords[!ex, , drop = FALSE ] ) )


        #         print(str(gamma))
        #         stop()

      } else {

#### --- covariance objects for block kriging

        ## construct covmodel

        tmp <- lapply(
          1L:length(variogram.object),
          function(i, x, type){

            variogram.model.v2 <- x[[i]][["variogram.model.v2"]]
            param              <- x[[i]][["param"]]

            ## setup covariance model list

            t.covmodel <- covmodel(
              modelname = variogram.model.v2,
              mev = switch(
                type,
                response = 0.,
                signal = unname( if( identical(i, 1L) ) param["nugget"] else 0. )
              ),
              nugget = switch(
                type,
                response = unname( if( identical(i, 1L) ) sum( param[c("snugget", "nugget")] ) else 0. ),
                signal = unname( if( identical(i, 1L) ) param["snugget"] else 0. )
              ),
              variance = unname( param["variance"] ),
              scale = unname( param["scale"] ),
              parameter = unname(
                if( length(param) > 4L-(i-1L)*2L ){
                  param[-(1:(4L-(i-1L)*2L))]
                } else {
                  NULL
                }
              )
            )
          }, x = variogram.object, type = type
        )

        t.covmodel <- tmp[[1]]
        if( length(tmp) > 1L ){
          for( i in 2L:length(tmp) ) t.covmodel <- c( t.covmodel, tmp[[i]] )
        }
        class(t.covmodel) <- class(tmp[[1]])

         ## variances of the prediction blocks

        t.preCKrige <- preCKrige(
          newdata = newdata[!ex, , drop = FALSE ],
          model = t.covmodel,
          pwidth = pwidth, pheight = pheight, napp = napp
        )
        t.var.target <- sapply(
          t.preCKrige@covmat,
          function( x ) c( x )
        )

        if( full.covmat ){

          t.neighbours <- lapply(
            1L:length( newdata[!ex, , drop = FALSE ] ),
            function(i) integer()
          )
          t.neighbours[[1L]] <- 2L:length( newdata[!ex, , drop = FALSE ] )

          t.preCKrige.aux <- preCKrige(
            newdata = newdata[!ex, , drop = FALSE ],
            neighbours = t.neighbours,
            model = t.covmodel,
            pwidth = pwidth, pheight = pheight, napp = napp
          )

          t.se <- sqrt( t.var.target )
          t.var.target <- t.se * t( t.se *
            cov2cor( t.preCKrige.aux@covmat[[1L]] ) )

        }


        ## get rid of mev component in covariance model list

        t.covmodel <- t.preCKrige@model[
        unlist(
          lapply(
            1L:length(t.preCKrige@model),
            function( i, m ){
              m[[i]][["model"]] != "mev"
            },
            m = t.preCKrige@model
          )
        )
        ]

        ## covariances between the support points and the prediction
        ## blocks

        gamma <- t(
          sapply(
            t.preCKrige@pixconfig,
            function( x, locations, model ){
              f.point.block.cov(
                pixconfig = x,
                locations = locations,
                model = model
              )
            },
            locations = locations.coords,
            model = t.covmodel
          )
        )

      }  ## end of block krighing


      ## initialize covariance matrix of predictions and covariance matrix
      ## of predictions and observations

      t.var <- NULL

#### --- compute uk predictions

      # gammaVi <- gamma %*% Valphaxi.objects[["Valphaxi.inverse"]] / var.signal
      #       gammaVi <- pmm( gamma, Valphaxi.objects[["Valphaxi.inverse"]], control = control.pcmp ) /
      #         var.signal
      gammaVi <- pmm( gamma, Valphaxi.inverse, control = control.pcmp ) / var.signal
      t.pred <- t.pred + drop( gammaVi %*% bhat )


#### --- compute uk variance (= (co-)variances of prediction errors)

      #         aux <- cbind(
      #           gammaVi %*% cov.delta.bhat.betahat.l[1L:n, 1L:n] - pred.X[!ex, , drop = FALSE ] %*% cov.delta.bhat.betahat.l[-(1L:n), 1L:n],
      #           - pred.X[!ex, , drop = FALSE ] %*% cov.delta.bhat.betahat.l[-(1L:n), -(1L:n)]
      #         )
      aux <- cbind(
        pmm( gammaVi, cov.delta.bhat.betahat.l[1L:n, 1L:n], control = control.pcmp ) -
        pmm(
          pred.X[!ex, , drop = FALSE ], cov.delta.bhat.betahat.l[-(1L:n), 1L:n],
          control = control.pcmp
        ),
        - pred.X[!ex, , drop = FALSE ] %*% cov.delta.bhat.betahat.l[-(1L:n), -(1L:n)]
      )

      if( full.covmat ){
        #           t.mse.pred <- tcrossprod( aux ) + t.var.target - gammaVi %*% t( gamma )
        t.mse.pred <- tcrossprod( aux ) + t.var.target - pmm(
          gammaVi, t( gamma ), control = control.pcmp
        )
      } else {
        t.mse.pred <- rowSums( aux^2 ) + t.var.target - rowSums( gammaVi * gamma )
      }

      #       t.pred.se.signal <- sqrt( f.diag( t.mse.pred ) )

      if( identical( type, "response" ) && !is.null( pred.coords ) ){

        #         if( rp.response ){
        #
        #           # Gaussian mixture predictive distribution for response (robust kriging)
        #
        #           #           t.var.pd.response <- var.pd.resp.rob( t.pred, t.pred.se.signal, scld.res )
        #           t.var.pd.response <- f.diag( t.mse.pred ) +
        #             var(scld.res) * (length(scld.res) - 1L) / length(scld.res)
        #
        #         } else {

        # Gaussian predictive distribution for response (non-robust kriging)

        t.var.pd.response <- f.diag( t.mse.pred ) + unname( nugget )

        #         }

        if( full.covmat ){
          diag( t.mse.pred ) <- t.var.pd.response
          diag( t.var.target ) <- diag( t.var.target ) + unname( nugget )
        } else {
          t.mse.pred <- t.var.pd.response
          t.var.target <- t.var.target + unname( nugget )
        }

      }

      if( extended.output ){

        ## compute covariance matrix of uk predictions and
        ## covariance matrix of uk predictions and prediction targets
        ## (needed for lognormal kriging)

        aux0 <- cbind( gammaVi, pred.X[!ex, , drop = FALSE ] )
        aux1 <- pmm( aux0, cov.bhat.betahat, control = control.pcmp )
        aux2 <- pmm( gamma, t(cov.p.t), control = control.pcmp )

        if( !full.covmat ){
          t.var.pred        <- rowSums( aux0 * aux1 )
          t.cov.pred.target <- rowSums( aux0 * aux2 )
        } else {
          t.var.pred        <- pmm( aux0, t(aux1), control = control.pcmp )
          t.cov.pred.target <- pmm( aux0, t(aux2), control = control.pcmp )
        }

        #           t.var.pred <- aux0 %*% cov.bhat.betahat %*% t( aux0 )
        #           t.cov.pred.target <- aux0 %*% cov.p.t %*% t( gamma )

      }

      ## for type == "response" correct  predictions for prediction
      ## locations that coincide with data locations

      if( !is.null( pred.coords ) && identical( type, "response" ) ){

        exx <- apply(
          pred.coords[!ex, , drop = FALSE ],
          1L,
          function(x, lc){
            tmp <- colSums( abs( t(lc) - x ) ) < sqrt(.Machine$double.eps)
            if( sum( tmp ) ){
              (1L:length(tmp))[tmp][1L]
            } else NA_integer_
          },
          lc = locations.coords
        )
        exx <- unname(exx)

        sel <- !is.na( exx )

        if( length( exx[sel] ) ){

          if( extended.output ){
            tmp <- var.signal * Valphaxi
            diag(tmp) <- diag(tmp) + nugget
          }

          t.pred[sel] <- response[exx[sel]]

          if( full.covmat ){
            t.mse.pred[sel, ] <- 0.
            t.mse.pred[, sel] <- 0.
            if( extended.output ){
              t.var.pred[sel, ] <- t.cov.pred.target[sel, ]
              t.var.pred[, sel] <- t.cov.pred.target[, sel]
              t.var.pred[sel, sel] <- tmp[exx[sel], exx[sel]]
              t.cov.pred.target[sel, sel] <- tmp[exx[sel], exx[sel]]
            }
          } else {
            t.mse.pred[sel] <- 0.
            if( extended.output ){
              t.var.pred[sel] <- diag(tmp)[exx[sel]]
              t.cov.pred.target[sel] <- t.var.pred[sel]
            }
          }
        }
      }

      ## end compute kriging predictions

    } else {

#### --- compute variance of trend surface prediction

      if( full.covmat ){
        t.var.pred <- tcrossprod( pred.X[!ex, , drop = FALSE ] %*% cov.betahat.l )
        t.mse.pred <- matrix( NA_real_, nrow = NROW( t.var.pred ), ncol = NCOL( t.var.pred ) )
        t.var.target <- matrix( 0., nrow = NROW( t.var.pred ), ncol = NCOL( t.var.pred ) )
        t.cov.pred.target <- matrix( 0., nrow = NROW( t.var.pred ), ncol = NCOL( t.var.pred ) )
      } else {
        t.var.pred <- rowSums( (pred.X[!ex, , drop = FALSE ] %*% cov.betahat.l)^2 )
        t.mse.pred <- rep( NA_real_, length( t.var.pred ) )
        t.var.target <- rep( 0., length( t.var.pred ) )
        t.cov.pred.target <- rep( 0., length( t.var.pred ) )
      }

    }

  } else {

    t.pred            <- NULL
    t.trend           <- NULL
    t.mse.pred        <- NULL
    t.var.pred        <- NULL
    t.cov.pred.target <- NULL

  }

#### -- prepare output

  ## add items with missing information back

  sr <- (1L:NROW(pred.X))[!ex]

  pred <- rep( NA_real_, NROW(pred.X) )
  if( length(sr) ) pred[sr] <- t.pred

  #   pred.se.signal <- rep( NA_real_, NROW(pred.X) )
  #   if( length(sr) ) pred.se.signal[sr] <- t.pred.se.signal

  if( extended.output ){
    trend <- rep( NA_real_, NROW(pred.X) )
    if( length(sr) ) trend[sr] <- t.trend
  }

  if( full.covmat ){

    mse.pred <- matrix( rep( NA_real_, NROW(pred.X)^2 ), nrow = NROW(pred.X) )
    if( length(sr) ) mse.pred[sr, sr] <- t.mse.pred

    if( identical( type, "trend" ) || extended.output ){
      var.pred <- matrix( rep( NA_real_, NROW(pred.X)^2 ), nrow = NROW(pred.X) )
      if( length(sr) ) var.pred[sr, sr] <- t.var.pred
    }
    if( extended.output ){
      cov.pred.target <- matrix( rep( NA_real_, NROW(pred.X)^2 ), nrow = NROW(pred.X) )
      if( length(sr) ) cov.pred.target[sr, sr] <- t.cov.pred.target
      var.target <- matrix( rep( NA_real_, NROW(pred.X)^2 ), nrow = NROW(pred.X) )
      if( length(sr) ) var.target[sr, sr] <- t.var.target
    }

  } else {

    mse.pred <- rep( NA_real_, NROW(pred.X) )
    if( length(sr) ) mse.pred[sr] <- t.mse.pred
    if( identical( type, "trend" ) || extended.output ){
      var.pred <- rep( NA_real_, NROW(pred.X) )
      if( length(sr) ) var.pred[sr] <- t.var.pred
    }
    if( extended.output ){
      cov.pred.target <- rep( NA_real_, NROW(pred.X) )
      if( length(sr) ) cov.pred.target[sr] <- t.cov.pred.target
      var.target <- rep( NA_real_, NROW(pred.X) )
      if( length(sr) ) var.target[sr] <- t.var.target
    }

  }

  ## collect results

  pred.se <-  sqrt( f.diag( mse.pred ) )

  #   result <- data.frame(
  #     pred = pred,
  #     se = pred.se,
  #     se.signal = pred.se.signal
  #   )

  result <- data.frame(
    pred = pred,
    se = pred.se
  )

  if( !is.null( signif ) ){
    #     if( rp.response ){
    #       t.quantiles <- qpd.resp.rob(
    #         0.5 * ( 1. + signif[1L] * c( -1L, 1L ) ),
    #         pred, pred.se.signal, scld.res
    #       )
    #     } else {
    t.quantiles <- cbind(
      pred + qnorm( 0.5 * ( 1.-signif[1L] ) ) * pred.se,
      pred + qnorm( 0.5 * ( 1.+signif[1L] ) ) * pred.se
    )
    #     }
    colnames( t.quantiles ) <- c( "lower", "upper" )
    result <- cbind( result, t.quantiles )
  }

  if( extended.output ) result[["trend"]] <- trend

  if( identical( type, "trend" ) || extended.output ){
    result[["var.pred"]] <- f.diag( var.pred )
    if( extended.output ){
      result[["cov.pred.target"]] <- f.diag( cov.pred.target )
      result[["var.target"]]      <- f.diag( var.target )
    }
  }

  if(
    !is.null( row.names( newdata ) ) &&
    length( row.names( newdata ) ) == NROW( result )
  ) row.names( result ) <- row.names( newdata )

  if( full.covmat ){
    result <- list( pred = result, mse.pred = mse.pred )
    if( extended.output ){
      result[["var.pred"]]        <- var.pred
      result[["cov.pred.target"]] <- cov.pred.target
      result[["var.target"]]      <- var.target
    }
  }

  return( result )
  ## end robust.uk
}

##  ###########################################################################
### simple.kriging.weights

simple.kriging.weights <- function(
  pred.coords, newdata = NULL,
  object = NULL,
  support.coords, locations, variogram.object,
  type = c("response", "signal"),
  covariances = FALSE,
  control = control.predict.georob(),
  verbose = 0, ...
)
{

  ## simplified version of predict.georob that computes only the simple
  ## kriging weights

  ## 2018-01-22 A. Papritz
  ## 2019-12-13 AP correcting use of class() in if() and switch()


  ##  ###########################################################################

#### -- auxiliary function for computing generalized covariances between
  ## prediction targets and support data

  f.gamma <- function(
    type,
    support.coords, pred.coords, newdata,
    variogram.object, var.signal, xi, nugget, gcr.constant,
    pwidth, pheight, napp,
    control.pcmp, verbose
  ){

    ## simplified version of f.robust.uk that returns a list with the
    ## simple kriging weights (gammaVi) and the variances of the prediction
    ## targets (t.var.target) (both are computed as in f.robust.uk)

    ## 2018-01-22 A. Papritz

    ## exclude prediction items with missing information

    if( !is.null( pred.coords ) ){              # point kriging

      n <- NROW(pred.coords)
      ex <- attr( na.omit( pred.coords ), "na.action" )
      if(length(ex)){
        ex <- ( 1L:n ) %in% sort( ex )
      } else {
        ex <- rep( FALSE, n )
      }


    } else {

      n <- NROW(newdata)
      ex <- rep( FALSE, n )                     # block kriging

    }

    if( any( !ex ) ){

      if( !is.null( pred.coords ) ){

        ## point kriging

        ## generalized covariance matrix between prediction and support
        ## points

        if( !all( sapply( variogram.object, function(x) x[["isotropic"]] ) ) ){
          indices.pairs <- expand.grid(
            1L:NROW( pred.coords[!ex, , drop = FALSE ] ),
            1L:NROW( support.coords )
          )
          lag.vectors <- (pred.coords[!ex, , drop = FALSE ])[ indices.pairs[, 1L], ] -
          support.coords[ indices.pairs[, 2L], ]
        } else {
          lag.vectors <- as.vector( rdist( pred.coords[!ex, , drop = FALSE ], support.coords ) )
        }
        attr( lag.vectors, "ndim.coords" ) <- NCOL(pred.coords[!ex, , drop = FALSE ])

        ## functions of version 3 of RandomFields

        Valpha <-  f.aux.gcr(
          lag.vectors = lag.vectors,
          variogram.object = variogram.object,
          gcr.constant = gcr.constant,
          symmetric = FALSE,
          control.pcmp = control.pcmp,
          verbose = verbose
        )

        if( any( sapply( Valpha, function(x) x[["error"]] ) ) ) stop(
          "an error occurred when computing semivariances between support ",
          "and prediction points"
        )

        gamma <- rowSums(
          sapply(
            1L:length(Valpha),
            function( i, x, xi ){
              xi[i] * x[[i]][["Valpha"]]
            }, x = Valpha, xi = xi
          )
        )

        ## add spatial nugget if prediction and support locations coincides

        if( sum(xi) < 1. ){
          if( NCOL(lag.vectors) > 1L ){
            sel <- abs( rowSums(lag.vectors) ) < control.georob()[["zero.dist"]]
          } else {
            sel <- lag.vectors < control.georob()[["zero.dist"]]
          }
          gamma[sel] <- gamma[sel] + (1. - sum(xi) )
        }

        gamma <- var.signal * gamma

        ## add nugget for prediction points that match a support point if
        ## weights are computed for response

        if( identical( type, "response" ) ){
          if( NCOL(lag.vectors) > 1L ){
            sel <- abs( rowSums(lag.vectors) ) < control.georob()[["zero.dist"]]
          } else {
            sel <- lag.vectors < control.georob()[["zero.dist"]]
          }
          gamma[sel] <- gamma[sel] + unname( nugget )
        }

        ## convert to matrix

        gamma <- matrix( gamma, nrow = NROW( pred.coords[!ex, , drop = FALSE ] ) )

      } else {

        ## block kriging

        ## construct covmodel

        tmp <- lapply(
          1L:length(variogram.object),
          function(i, x, type){

            variogram.model.v2 <- x[[i]][["variogram.model.v2"]]
            param              <- x[[i]][["param"]]

            ## setup covariance model list

            t.covmodel <- covmodel(
              modelname = variogram.model.v2,
              mev = switch(
                type,
                response = 0.,
                signal = unname( if( identical(i, 1L) ) param["nugget"] else 0. )
              ),
              nugget = switch(
                type,
                response = unname( if( identical(i, 1L) ) sum( param[c("snugget", "nugget")] ) else 0. ),
                signal = unname( if( identical(i, 1L) ) param["snugget"] else 0. )
              ),
              variance = unname( param["variance"] ),
              scale = unname( param["scale"] ),
              parameter = unname(
                if( length(param) > 4L-(i-1L)*2L ){
                  param[-(1:(4L-(i-1L)*2L))]
                } else {
                  NULL
                }
              )
            )
          }, x = variogram.object, type = type
        )

        t.covmodel <- tmp[[1]]
        if( length(tmp) > 1L ){
          for( i in 2L:length(tmp) ) t.covmodel <- c( t.covmodel, tmp[[i]] )
        }
        class(t.covmodel) <- class(tmp[[1]])

        ## variances of the prediction blocks

        t.preCKrige <- preCKrige(
          newdata = newdata[!ex, , drop = FALSE ],
          model = t.covmodel,
          pwidth = pwidth, pheight = pheight, napp = napp
        )
        t.var.target <- sapply(
          t.preCKrige@covmat,
          function( x ) c( x )
        )

        ## get rid of mev component in covariance model list

        t.covmodel <- t.preCKrige@model[
        unlist(
          lapply(
            1L:length(t.preCKrige@model),
            function( i, m ){
              m[[i]][["model"]] != "mev"
            },
            m = t.preCKrige@model
          )
        )
        ]

        ## covariances between the support points and the prediction
        ## blocks

        gamma <- t(
          sapply(
            t.preCKrige@pixconfig,
            function( x, locations, model ){
              f.point.block.cov(
                pixconfig = x,
                locations = locations,
                model = model
              )
            },
            locations = support.coords,
            model = t.covmodel
          )
        )

      }  ## end of block krighing


    } else {

      gamma <- NULL

    }

    ## add items with missing information back

    sr <- (1L:n)[!ex]

    result <- matrix( rep( NA_real_, n * NCOL(gamma) ), nrow = n )
    if( length(sr) ) result[sr, ] <- gamma

    return( result )

  }

  ##  ##############################################################################

  ## begin of main body of function

  RFoptions(newAniso=FALSE)

  ## match arguments

  type <- match.arg( type )

  ## mandatory argument

  if( is.null(newdata) && missing(pred.coords) ) stop(
    "either 'newdata' or 'pred.coords' must be specified"
  )


#### -- coordinates and generalized covariance matrix of support data

  ## if object is missing

  if( is.null(object) ){

    ## compute required items because object is missing

    if( missing(support.coords) || missing(locations) || missing(variogram.object) ) stop(
      "some mandatory arguments are missing"
    )

    ## signal variance, xi, nugget

    tmp <- f.reparam.fwd( variogram.object )

    var.signal <- attr(tmp, "var.signal" )
    xi <- sapply( tmp, function(x) x[["param"]]["variance"] )
    nugget <- variogram.object[[1L]][["param"]]["nugget"]

    ## lag vectors

    isotropic <- all( sapply(variogram.object, function(x) x[["isotropic"]] ) )

    if( !isotropic ){
      i.pairs <- combn( NROW( support.coords ), 2L )
      lag.vectors <- support.coords[ i.pairs[2L,], ] - support.coords[ i.pairs[1L,], ]
    } else {
      lag.vectors <- as.vector( dist( support.coords ) )
    }
    attr( lag.vectors, "ndim.coords" ) <- NCOL(support.coords)

    ## compute generalized covariance matrix of support data

    Valpha <- f.aux.gcr(
      lag.vectors = lag.vectors,
      variogram.object = variogram.object,
      control.pcmp = control[["pcmp"]],
      verbose = verbose
    )

    if( any( sapply( Valpha, function(x) x[["error"]] ) ) ) stop(
      "an error occurred when computing semivariances between prediction points"
    )

    gcvmat <- list(
      diag = var.signal * rowSums(
        sapply(
          1L:length(Valpha),
          function( i, x, xi ){
            xi[i] * x[[i]][["Valpha"]][["diag"]]
          }, x = Valpha, xi = xi
        )
      ) + ( 1. - sum(xi) ) + nugget,
      tri = var.signal * rowSums(
        sapply(
          1L:length(Valpha),
          function( i, x, xi ){
            xi[i] * x[[i]][["Valpha"]][["tri"]]
          }, x = Valpha, xi = xi
        )
      )
    )
    attr( gcvmat, "struc" ) <- "sym"

    gcvmat <- expand( gcvmat )

    ## extract gcr.constant

    gcr.constant <- lapply(
      Valpha,
      function(x) x[["gcr.constant"]]
    )

  } else {

    ## get required items from object

    ## coordinates and formula for locations

    support.coords <- object[["locations.objects"]][["coordinates"]][!duplicated( object[["Tmat"]] ), , drop = FALSE]

    if( missing( locations ) ){
      locations <-  object[["locations.objects"]][["locations"]]
    }

    ## signal variance, xi, nugget and gcr.constant

    variogram.object  <- object[["variogram.object"]]

    tmp <- f.reparam.fwd( variogram.object )

    var.signal <- attr(tmp, "var.signal" )
    xi <- sapply( tmp, function(x) x[["param"]]["variance"] )

    nugget <- variogram.object[[1L]][["param"]]["nugget"]

    ## generalized covariance matrix of support data

    Valphaxi.objects  <- expand( object[["Valphaxi.objects"]] )

    gcvmat <- var.signal * Valphaxi.objects[["Valphaxi"]]
    diag( gcvmat ) <- diag( gcvmat ) + unname(nugget)
    attr( gcvmat, "struc" ) <- "sym"

    ## gcr.constant

    gcr.constant <- lapply(
      Valphaxi.objects[["Valpha"]],
      function(x) x[["gcr.constant"]]
    )

  }

#### -- inverse generalized covariance matrix of support data

  gcvmat.inverse <- try( chol2inv( chol( gcvmat ) ) )
  if( identical( class(gcvmat.inverse), "try-error" ) ) stop(
    "an error occurred when computing the inverse generalized covariance matrix of the data"
  )
  attr( gcvmat.inverse, "struc" ) <- "sym"


#### -- prepare items for block kriging

  if( !missing( newdata ) && class( newdata )[1] == "SpatialPolygonsDataFrame" ){

    ## check whether pwidth and pheight were provided

    if( is.null( control[["pwidth"]] ) || is.null( control[["pheight"]] ) ) stop(
      "'pwidth' and 'pheight' must be provided for block kriging"
    )

    ## map names of variogram models of RandomFields version 3 to version 2

    variogram.object <- lapply(
      variogram.object,
      function(x){

        variogram.model <- x[["variogram.model"]]
        isotropic <- x[["isotropic"]]
        param <- x[["param"]]

        if( variogram.model %in% control.georob()[["irf.models"]] ) stop(
          "block kriging not yet implemented for unbounded variogram models"
        )
        if( !isotropic ) stop(
          "block kriging not yet implemented for anisotropic variograms"
        )

        variogram.model.v2 <- gsub("^RM", "", variogram.model )

        variogram.model.v2 <- switch(
          variogram.model.v2[1],
          askey = stop(
            "variogram model 'RMaskey' not implemented in package constrainedKriging"
          ),
          dagum = stop(
            "variogram model 'RMdagum' not implemented in package constrainedKriging"
          ),
          dewijsian = stop(
            "variogram model 'RMdewijsian' not implemented in package constrainedKriging"
          ),
          fbm = stop(
            "variogram model 'RMfbm' not implemented in package constrainedKriging"
          ),
          genfbm = stop(
            "variogram model 'RMgenfbm' not implemented in package constrainedKriging"
          ),
          dampedcos = "dampedcosine",
          exp = "exponential",
          lgd = "lgd1",
          qexp = "qexponential",
          spheric = "spherical",
          variogram.model.v2
        )

        if( identical( variogram.model.v2, "gengneiting" ) ) param[6L] <- sum( param[5L:6L] ) + 0.5

        x[["variogram.model.v2"]] <- variogram.model.v2
        x[["param"]] <- param

        x

      }
    )

  }


#### -- compute generalized covariance matrix between prediction and support sites

  Terms.loc <- terms( locations )
  attr( Terms.loc, "intercept" ) <- 0

  ## get coordinates of prediction sites from newdata if missing

  if( missing(pred.coords) ){

    pred.coords <- switch(
      class( newdata )[1],
      "data.frame" = model.matrix(
        Terms.loc,
        model.frame(
          Terms.loc, newdata, na.action = na.pass
        )
      ),
      "SpatialPoints" = ,
      "SpatialPixels" = ,
      "SpatialGrid" = ,
      "SpatialPointsDataFrame" = ,
      "SpatialPixelsDataFrame" = ,
      "SpatialGridDataFrame" = coordinates( newdata ),
      "SpatialPolygonsDataFrame" = NULL
    )

  } else {

    if(!(
        NCOL( support.coords ) == NCOL( pred.coords ) &&
        all( colnames( pred.coords ) ==
          colnames( support.coords ) )
      )
    ) stop(
      "inconsistent number and/or names of coordinates in 'object' and in 'newdata' or 'pred.coords'"
    )

  }


  ## number of items to predict

  if( !is.null( pred.coords ) ){
    m <- NROW( pred.coords )
  } else {
    m <- NROW( newdata )
  }


  ## determine number of prediction parts

  n.part <- ceiling( m / control[["mmax"]] )
  rs <- ( 0L:(n.part-1L)) * control[["mmax"]] + 1L
  re <- ( 1L:(n.part  )) * control[["mmax"]]; re[n.part] <- m

  ncores <- min( n.part, control[["ncores"]] )

  ncores.available <- control[["pcmp"]][["max.ncores"]]
  if( sfIsRunning() ) sfStop()

  control.pcmp <- control[["pcmp"]]
  control.pcmp[["pmm.ncores"]] <- min(
    control.pcmp[["pmm.ncores"]],
    max( 1L, floor( (ncores.available - ncores) / ncores ) )
  )
  if( ncores > 1L && !control.pcmp[["allow.recursive"]] ) control.pcmp[["pmm.ncores"]] <- 1L

  if( control[["full.covmat"]] && n.part > 1L ) stop(
    "full covariance matrix of prediction errors cannot ",
    "be computed\n  if prediction problem is split into several parts\n",
    "-> increase 'mmax' to avoid splitting"
  )


  ## handle parallel processing

  ## auxiliary function to compute the predictions for one part

  f.aux <- function( i ){

    ## objects
    ##
    ##       rs, re, n.part,
    ##       type,
    ##       support.coords, pred.coords, newdata,
    ##       variogram.object, var.signal, xi, nugget, gcr.constant,
    ##       pwidth, pheight, napp,
    ##       control.pcmp, verbose
    ##
    ## are taken from parent enviroment

    if( verbose > 0. )
    cat( "  predicting part ", i, " of ", n.part, "\n" )

    ## select the data for the current part

    if( !is.null( pred.coords ) ) {
      pred.coords <- pred.coords[ rs[i]:re[i], , drop = FALSE]
    } else {
      newdata <- newdata[ rs[i]:re[i], ]
    }

    ## compute simple kriging weights for the current part

    result <- f.gamma(
      type = type,
      support.coords = support.coords, pred.coords = pred.coords, newdata = newdata,
      variogram.object = variogram.object, var.signal = var.signal,
      xi = xi, nugget = nugget, gcr.constant = gcr.constant,
      pwidth = pwidth, pheight = pheight, napp = napp,
      control.pcmp = control.pcmp,
      verbose = verbose
    )

    return( result )

  }

  ## prepare items to pass to function f.aux

  pwidth            <- control[["pwidth"]]
  pheight           <- control[["pheight"]]
  napp              <- control[["napp"]]

  ## set default value for control of forking if missing (required for backward compatibility)

  if( is.null( control[["pcmp"]][["fork"]] ) ){
    control[["pcmp"]][["fork"]] <- !identical( .Platform[["OS.type"]], "windows" )
  }
  ## compute the predictions for all the parts

  if( ncores > 1L && !control[["pcmp"]][["fork"]] ){

    ## create a SNOW cluster on windows OS

    clstr <- makeCluster( ncores, type = "SOCK")
    save( clstr, file = "SOCKcluster.RData" )
    options( error = f.stop.cluster )

    ## export required items to workers

    junk <- clusterEvalQ( clstr, require( georob, quietly = TRUE ) )

    junk <- clusterExport(
      clstr,
      c(
        "rs", "re", "n.part",
        "type",
        "support.coords",
        "pred.coords", "newdata",
        "variogram.object", "var.signal",
        "xi", "nugget", "gcr.constant",
        "pwidth", "pheight", "napp",
        "control.pcmp",  "verbose"
      ),
      envir = environment()
    )


    t.result <- try(
      parLapply(
        clstr,
        1L:n.part,
        f.aux
      )
    )

    f.stop.cluster( clstr )

  } else {

    ## fork child processes on non-windows OS

    t.result <- try(
      mclapply(
        1L:n.part,
        f.aux,
        mc.cores = ncores,
        mc.allow.recursive = control.pcmp[["allow.recursive"]]
      )
    )

  }

  has.error <- sapply(
    t.result, function( x ) identical( class(x), "try-error" )
  )

  if( any( has.error ) ){
    cat( "\nerror(s) occurred when computing kriging predictions in parallel:\n\n" )
    sapply( t.result[has.error], cat)
    cat( "\nuse 'ncores=1' and 'verbose = 1' to avoid parallel computations and to see where problem occurs\n\n" )
    stop()
  }

  ## collect results of the various parts into a matrix

  gcvmat.pred <- t.result[[1L]]
  if( length( t.result ) > 1L ){
    for( i in 2L:length( t.result ) ) {
      gcvmat.pred <- rbind( gcvmat.pred, t.result[[i]] )
    }
  }

#### -- prepare output

  result <- pmm( gcvmat.pred, gcvmat.inverse, control = control[["pcmp"]] )

  if( covariances ){
    result <- list(
      gcvmat = gcvmat,
      gcvmat.inverse = gcvmat.inverse,
      gcvmat.pred = gcvmat.pred,
      skw = result
    )
    result <- compress( result )
  }

  invisible( result )

}




# ##  ###########################################################################
#
# # functions to evaluate
# #
# #   1) cdf,
# #   2) pdf,
# #   3) quantiles,
# #   4) mean and variance,
# #   5) continuous ranked probability score of the predictive distribution
# #      of the response given the observations
# #
# # and
# #
# #   6) to simulate from these distributions.
# #
# # the predictive distribution is modelled by the convolution of the
# # Gaussian predictive distribution of the signal given the observations
# # (parametrized by m and s) and the empirical distribution of the
# # scaled residuals residuals/resscl of the fitted model object.  this results in a
# # mixture of equally weighted gaussian distributions with parameters
# # m + residuals/resscl and s
#
# # the functions use functions of the package nor1mix
#
# # common arguments
#
# # q      vector with quantiles for which pdf and cdf should be evaluated
# # p      vector with probabilities for which quantiles should be evaluated
# # y      vector with values of observations for which crps should be evaluated
# # m      vector with kriging predictions of signal
# # s      vector with kriging standard error of signal
# # r      vector with unscaled resiuals (estimated epsilons)
# # resscl scaling factor of residuals
#
# ##  ###########################################################################
#
# ppd.resp.rob <- function( q, m, s, r, resscl = 1., lower.tail = TRUE, log.p = FALSE ){
#
#   # cdf of robust predictive distribution of response
#
#   # 2015-06-24 A. Paprritz
#
#   param <- na.exclude( cbind( m, s ) )
#   r <- r[!is.na(r)] / resscl
#   q <- q[!is.na(q)]
#
#   # using function pnorMix{nor1mix}
#
#   result <- sapply(
#     1L:NROW(param),
#     function( i, m, s, r, q, lower.tail, log.p ){
#       pnorMix( q, norMix( m[i] + r, sigma = rep( s[i], length( r ) ) ),
#         lower.tail = lower.tail, log.p = log.p
#       )
#     },
#     m = param[, "m"], s = param[, "s"],
#     r = r, q = q, lower.tail = lower.tail, log.p = log.p
#   )
#
#   if( is.matrix( result ) ) result <- t( result )
#   napredict( attr( param, "na.action" ), result )
#
# }
#
#
# ##  ###########################################################################
#
# dpd.resp.rob <- function( q, m, s, r, resscl = 1., log = FALSE ){
#
#   # pdf of robust predictive distribution of response
#
#   # 2015-06-24 A. Paprritz
#
#   param <- na.exclude( cbind( m, s ) )
#   r <- r[!is.na(r)] / resscl
#   q <- q[!is.na(q)]
#
#   # using function pnorMix{nor1mix}
#
#   result <- sapply(
#     1L:NROW(param),
#     function( i, m, s, r, q, lower.tail, log.p ){
#       dnorMix( q, norMix( m[i] + r, sigma = rep( s[i], length( r ) ) ), log = log )
#     },
#     m = param[, "m"], s = param[, "s"],
#     r = r, q = q, lower.tail = lower.tail, log.p = log.p
#   )
#
#   if( is.matrix( result ) ) result <- t( result )
#   napredict( attr( param, "na.action" ), result )
#
# }
#
#
# ##  ###########################################################################
#
# qpd.resp.rob <- function( p, m, s, r, resscl = 1., lower.tail = TRUE ){
#
#   # quantiles of robust predictive distribution of response
#
#   # 2015-06-24 A. Paprritz
#
#   param <- na.exclude( cbind( m, s ) )
#   r <- r[!is.na(r)] / resscl
#   p <- p[!is.na(p)]
#
#   # using function qnorMix{nor1mix}
#
#   result <- sapply(
#     1L:NROW(param),
#     function( i, m, s, r, p, lower.tail ){
#       qnorMix( p, norMix( m[i] + r, sigma = rep( s[i], length( r ) ) ),
#         lower.tail = lower.tail
#       )
#     },
#     m = param[, "m"], s = param[, "s"],
#     r = r, p = p, lower.tail = lower.tail
#   )
#
#   if( is.matrix( result ) ) result <- t( result )
#   napredict( attr( param, "na.action" ), result )
#
# }


# ##  ###########################################################################
#
# mean.pd.resp.rob <- function( m, s, r, resscl = 1. ){
#
#   # means of robust predictive distribution of response
#
#   # 2015-06-24 A. Paprritz
#
#   param <- na.exclude( cbind( m, s ) )
#   r  <- r[!is.na(r)] / resscl
#
#   # using function mean.norMix{nor1mix}
#
#   result <- sapply(
#     1L:NROW(param),
#     function( i, m, s, r ){
#       mean( norMix( m[i] + r, sigma = rep( s[i], length( r ) ) ) )
#     },
#     m = param[, "m"], s = param[, "s"], r = r
#   )
#
#   napredict( attr( param, "na.action" ), result )
#
# }
#
#
# ##  ###########################################################################
#
# var.pd.resp.rob <- function( m, s, r, resscl = 1. ){
#
#   # variances of robust predictive distribution of response
#
#   # 2015-06-24 A. Paprritz
#
#   param <- na.exclude( cbind( m, s ) )
#   r  <- r[!is.na(r)] / resscl
#
#   # using function var.norMix{nor1mix}
#
#   result <- sapply(
#     1L:NROW(param),
#     function( i, m, s, r ){
#       var.norMix( norMix( m[i] + r, sigma = rep( s[i], length( r ) ) ) )
#     },
#     m = param[, "m"], s = param[, "s"], r = r
#   )
#
#   napredict( attr( param, "na.action" ), result )
#
# }


##  ###########################################################################

# crps

# cf.  equations 5 & 6 of
#
# @Article{Grimit-etal-2006,
#   Title                    = {The continuous ranked probability score for circular variables and its application to mesoscale forecast ensemble verification},
#   Author                   = {Grimit, E.P. and Gneiting, T. and Berrocal, V.J. and Johnson, N.A.},
#   Journal                  = {Quarterly Journal of the Royal Meteorological Society},
#   Year                     = {2006},
#   Number                   = {621 C},
#   Pages                    = {2925--2942},
#   Volume                   = {132},
#
#   Doi                      = {10.1256/qj.05.235},
#   File                     = {Grimit-etal-2006.pdf:PDFs/G/Grimit-etal-2006.pdf:PDF},
#   Separatanr               = {523},
#   Url                      = {http://www.scopus.com/inward/record.url?eid=2-s2.0-33947236600&partnerID=40&md5=24b826e4f3805ce68268993f3f0e208a}
# }


##  ###########################################################################

f.aux.crpsnorm <- function( m, s ){

  # auxiliary function for computing crps of a normal random variate with
  # mean m and standard deviation s (cf. equation 6, Grimit et al., 2006

  # 2015-06-24 A. Paprritz

  x <- m / s
  2. * s * dnorm(x) + m * ( 2.* pnorm(x) - 1. )

}

##  ###########################################################################

crpsnorm <- function( y, m, s ){

  # function to compute continuous ranked probability score for a normal
  # distribution with mean m and standard deviation s

  # 2015-06-24 A. Paprritz

  param <- na.exclude( cbind( m, s, y ) )

  result <- f.aux.crpsnorm( param[, "y"] - param[, "m"], param[, "s"] ) -
    0.5 * f.aux.crpsnorm( 0., sqrt(2.) * param[, "s"] )

  napredict( attr( param, "na.action" ), result )

}


# ##  ###########################################################################
#
# crpspd.resp.rob <- function( y, m, s, r, resscl = 1. ){
#
#   # computinng crps of robust predictive distribution of response
#
#   # 2015-06-24 A. Paprritz
#
#   # auxiliary function for computing crps of a normal random variate with
#   # mean m and standard deviation s (cf. equation 6, Grimit et al., 2006
#
#   f.aux.crpsnorm <- function( m, s ){
#     x <- m / s
#     2. * s * dnorm(x) + m * ( 2.* pnorm(x) - 1. )
#   }
#
#   param <- na.exclude( cbind( m, s, y ) )
#   r  <- r[!is.na(r)] / resscl
#
#   result <- sapply(
#     1L:NROW(param),
#     function( i, m, s, r, y ){
#
#       mm <- m[i] + r
#       ss <- rep( s[i], length( r ) )
#       yy <- rep( y[i], length( r ) )
#
#       mmm <- as.vector( outer( mm, mm, FUN = "-" ) )
#       sss <- rep( sqrt(2.) * s[i], length = length( mmm ) )
#
#       mean( f.aux.crpsnorm( yy - mm, ss ) ) - 0.5 * mean( f.aux.crpsnorm( mmm, sss ) )
#
#     },
#     m = param[, "m"], s = param[, "s"], r = r, y = param[, "y"]
#   )
#
#   napredict( attr( param, "na.action" ), result )
#
# }
#
#
# ##  ###########################################################################
#
# rpd.resp.rob <- function( n, m, s, r, resscl = 1. ){
#
#   # simulating random deviates from robust predictive distribution of response
#
#   # 2015-06-24 A. Paprritz
#
#   param <- na.exclude( cbind( m, s ) )
#   r  <- r[!is.na(r)] / resscl
#
#   # using function rnorMix{nor1mix}
#
#   result <- sapply(
#     1L:NROW(param),
#     function( i, m, s, r, n ){
#       rnorMix( n, norMix( m[i] + r, sigma = rep( s[i], length( r ) ) ) )
#     },
#     m = param[, "m"], s = param[, "s"], r = r, n = n
#   )
#
#   if( is.matrix( result ) ) result <- t( result )
#   napredict( attr( param, "na.action" ), result )
#
# }


##  ###########################################################################

check.newdata <- function(newdata){

  ## function checks whether newdata is a valid object to be passed to
  ## predict.georob or lgnpp

  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  stopifnot(
    all(
      class(newdata)[1] %in% c(
        "data.frame", "SpatialPointsDataFrame", "SpatialPixelsDataFrame",
        "SpatialGridDataFrame", "SpatialPolygonsDataFrame",
        "SpatialPoints", "SpatialPixels", "SpatialGrid"
      )
    )
  )

  return(NULL)

}
