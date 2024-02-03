##  ###########################################################################

### control.condsim -------------

control.condsim <- function(
  use.grid = FALSE,
  grid.refinement = 2.,
  condsim = TRUE,
  ce.method = c( "standard", "approximate" ),
  ce.grid.expansion = 1.,
#   ce.method = c( "standard", "approximate", "cutoff", "intrinsic" ),
  include.data.sites = FALSE,
  means = FALSE,
  trend.covariates = FALSE,
  covariances = FALSE,
  ncores = 1,
  mmax = 10000,
  pcmp = control.pcmp()
){

  ## auxiliary function to set meaningful default values for
  ## condsim

  ## 2018-01-12 A. Papritz
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()
  ## 2024-01-24 AP new arguments ce.method, cr.grid.expansion, mmax

  ## check arguments

  stopifnot(identical(length(use.grid), 1L)           && is.logical(use.grid))
  stopifnot(identical(length(condsim), 1L)            && is.logical(condsim))
  stopifnot(identical(length(include.data.sites), 1L) && is.logical(include.data.sites))
  stopifnot(identical(length(means), 1L)              && is.logical(means))
  stopifnot(identical(length(trend.covariates), 1L)   && is.logical(trend.covariates))
  stopifnot(identical(length(covariances), 1L)        && is.logical(covariances))

#   stopifnot(identical(length(grid.refinement), 1L) && is.numeric(grid.refinement) )
  stopifnot(
    identical(length(grid.refinement), 1L) &&
    is.numeric(grid.refinement) && grid.refinement > 1
  )
  stopifnot(
    identical(length(ce.grid.expansion), 1L) &&
    is.numeric(ce.grid.expansion) && ce.grid.expansion >= 1
  )
  stopifnot(
    identical(length(ncores), 1L) && is.numeric(ncores) && ncores >= 1
  )

  stopifnot(is.list(pcmp))

  if( grid.refinement < 1. ) warning( "grid refinement factor < 1")

  ce.method <- match.arg( ce.method )
  ## prepare list

  list(
    use.grid = use.grid,
    grid.refinement = grid.refinement,
    condsim = condsim,
    ce.method = ce.method,
    ce.grid.expansion = ce.grid.expansion,
    include.data.sites = include.data.sites,
    means = means,
    trend.covariates = trend.covariates,
    covariances = covariances,
    ncores = ncores, mmax = mmax,
    pcmp = pcmp
  )

}


##  ###########################################################################

### condsim -------------

condsim <- function(
  object, newdata, nsim, seed,
  type = c("response", "signal"),
  locations,
  trend.coef = NULL,
  variogram.model = NULL, param = NULL, aniso = NULL,
  variogram.object = NULL,
  control = control.condsim(),
  verbose = 0
){

  ## simulate (un)conditional realizations of Gaussian processes with the
  ## parameters estimated by georob for a spatial linear model

  ## 2018-01-27 A. Papritz
  ## 2018-11-25 A. Papritz small change in function f.aux.sim.1
  ## 2019-03-29 A. Papritz restore default for RFoption spConform
  ## 2019-03-29 A. Papritz conditioning data is passed as SpatialPointsDataFrame to
  ##                       RFsimulate (only applies for exact coordinates)
  ## 2019-05-24 A. Papritz correction of error that occured when preparing output consisting
  ##                       of a single realization
  ## 2019-12-17 A. Papritz correction of errors in output from RFsimulate
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()
  ## 2023-01-28 AP Forcing socket clusters for parallelized computations
  ## 2023-12-14 AP checking class by inherits()
  ## 2023-12-20 AP added on.exit(options(old.opt)), deleted options(error = NULL)
  ## 2023-12-21 AP replacement of identical(class(...), ...) by inherits(..., ...)
  ## 2024-01-21 AP new function sim.chol.decomp for simulation with exact coordinates
  ## 2024-01-21 AP new function sim.circulant.embedding for simulation on rectangular grid

  ## auxiliary function for aggregating a response variable x for the
  ## levels of the variable by.one

  f.aggregate.by.one <- function(x, by.one, FUN = mean, ...){
    by.one <- as.numeric( factor( by.one, levels = unique( by.one ) ) )
    aggregate( x, list( by = by.one ), FUN = FUN, ...  )[, -1]
  }

#### -- check arguments -------------

  ## match arguments

  type <- match.arg( type )

  ## mandatory argument

  if( missing( object ) || missing( newdata ) || missing( nsim ) || missing( seed ) ) stop(
    "some mandatory arguments are missing"
  )

  if(!missing(newdata)) check.newdata(newdata)

  stopifnot(identical(length(nsim), 1L)    && is.numeric(nsim)    && nsim >= 1)
  stopifnot(identical(length(seed), 1L)    && is.numeric(seed))
  stopifnot(identical(length(verbose), 1L) && is.numeric(verbose) && verbose >= 0)

  stopifnot(is.null(trend.coef) || is.numeric(trend.coef))
  stopifnot(is.null(param)      || is.numeric(param))
  stopifnot(is.null(aniso)      || is.numeric(aniso))

  stopifnot(is.null(variogram.model)  || (identical(length(variogram.model), 1L) && is.character(variogram.model)))

  stopifnot(is.list(control))
  stopifnot(is.null(variogram.object) || is.list(variogram.object))

  if( verbose > 0 ) cat( "  prepare data ...\n" )


  ## check consistency of provided arguments

  #   ## conditionally simulating signal
  #
  #   if( identical(type, "signal") && control[["condsim"]] && !control[["use.grid"]] ) stop(
  #     "conditional simulation of signal requires control argument 'use.grid = TRUE'"
  #   )

  ## class

  if( !inherits( object, "georob") ) stop(
    "object must be of class georob"
  )

  if( inherits( newdata, "SpatialPolygonsDataFrame") ) stop(
    "newdata is a SpatialPolygonsDataFrame; simulation for such targetst is not yet implemented"
  )

  ## coordinates only as covariates for trend model if newdata is a Spatialxx object

  formula.fixef <- as.formula(
    paste( as.character( formula( object ) )[-2L], collapse = " " )
  )
  if( class( newdata )[1] %in% c( "SpatialPoints", "SpatialPixels", "SpatialGrid" ) ){
    tmp <- try(
      get_all_vars( formula.fixef, as.data.frame( coordinates( newdata ) ) ),
      silent = TRUE
    )
    if( inherits( tmp, "try-error" ) ) stop(
      "'newdata' is a SpatialPoints, SpatialPixels or SpatialGrid object\n but drift covariates are not functions of coordinates"
    )
  }

  ## mean parameters

  if( !is.null( trend.coef) && !identical(length(trend.coef), length(coef(object)) ) ) stop(
    "'trend.coef' inconsistent with trend model of 'object'"
  )


#### -- re-fit model for new variogram -------------

  ## re-fit model if variogram parameters have been specified

  ## setup or check contents of variogram.object

  if( !all(
      is.null(variogram.model), is.null(param),  is.null(aniso),
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

      ## match names of param, aniso, fit.param, fit.aniso

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


#### -- support data sites: means and coordinates -------------

  ## extract required items about data sites (support data) from object

  ## coordinates

  if( missing( locations ) ){
    locations <- object[["locations.objects"]][["locations"]]
  }

  Terms.loc <- terms( locations )
  attr( Terms.loc, "intercept" ) <- 0

  coords.d <- as.data.frame(
    object[["locations.objects"]][["coordinates"]][!duplicated( object[["Tmat"]] ), , drop = FALSE]
  )

  ## terms for fixed effects

  tt <- terms( object )
  Terms <- delete.response( tt )

  ## model.frame

  mf.d <- model.frame( object )

  ## observations

  y.d <- f.aggregate.by.one(
    model.response( mf.d ), object[["Tmat"]], na.rm = TRUE
  )

  ## (un-)conditional mean

  if( !is.null( trend.coef ) ){

    ## specified fixed effects

    ## design matrix

    X.d <- model.matrix( object )[!duplicated( object[["Tmat"]] ), , drop = FALSE]

    ## deal with non-NULL offset

    if( !is.null( attr( tt, "offset" ) ) || !is.null( object[["call"]][["offset"]] ) ){
      offset <- model.offset( mf.d )[!duplicated( object[["Tmat"]] )]
    } else offset <- rep( 0., NROW(X.d) )

    ## unconditional mean

    mn.d <- drop( X.d %*% trend.coef ) + offset

    ## conditional mean

    if( control[["condsim"]] ){

      if( identical( type, "response" ) ){

        ## observations (kriging is an exact predictor for type ==
        ## "response"

        mn.d <- y.d

      } else {

        ## simple kriging weights for prediction of signal

        skw <- simple.kriging.weights(
          pred.coords = coords.d,
          object = object,
          type = type,
          covariances = FALSE,
          control = control.predict.georob(
            mmax = control[["mmax"]],
            ncores = control[["ncores"]]
          )
        )

        ## conditional mean

        mn.d <- mn.d + drop( skw %*% ( y.d - mn.d ) )

      }

    }


  } else {

    ## estimated fixed effects

    if( !control[["condsim"]] ){

      ## unconditional mean

      mn.d <- fitted( object )[!duplicated( object[["Tmat"]] )]

    } else {

      ## conditional mean

      if( identical( type, "response" ) ){

        ## observations (kriging is an exact predictor for type ==
        ## "response"

        mn.d <- y.d

      } else {

        ## kriging prediction of signal

        tmp <- predict(
          object,
          newdata = cbind(
            model.frame( Terms.loc, mf.d, na.action = na.pass ),
            get_all_vars( formula.fixef, data = eval( getCall(object)[["data"]] ) )
          ),
          locations = locations,
          type = type,
          control = control.predict.georob(
            mmax = control[["mmax"]],
            ncores = control[["ncores"]]
          )
        )

        ## conditional mean

        mn.d <- f.aggregate.by.one(
          tmp[, "pred"], object[["Tmat"]], na.rm = TRUE
        )

      }

    }

  }


#### -- simulation sites: means and coordinates -------------

  ## prepare required items for simulation sites

  if( verbose > 0 ) cat( "  compute (conditional) means ...\n" )


#### ---- using specified trend coefficients -------------

  if( !is.null( trend.coef ) ){

    ## convert SpatialGridDataFrame to SpatialPixelDataFrame

    if( inherits( newdata, "SpatialGridDataFrame" ) ){
      fullgrid.newdata <- TRUE
      fullgrid( newdata ) <- FALSE
    } else {
      fullgrid.newdata <- FALSE
    }

    ## get modelframe for newdata

    mf.s <- switch(
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
        class( newdata ), "'"
      )
    )

    ## check whether variables that will be used to compute the
    ## predictions agree with those in object

    if( !is.null( cl <- attr(Terms, "dataClasses" ) ) )
      .checkMFClasses( cl, mf.s )

    ## get fixed effects design matrix for newdata

    X.s <- model.matrix( Terms, mf.s, contrasts.arg = object[["contrasts"]] )

    ## deal with non-NULL offset

    offset <- rep( 0., NROW(X.s) )
    if( !is.null( off.num <- attr( tt, "offset" ) ) ){
      for( i in off.num ) {
        offset <- offset + eval( attr( tt, "variables" )[[i + 1L]], newdata )
      }
    }
    if( !is.null( object[["call"]][["offset"]] ) ){
      offset <- offset + eval( object[["call"]][["offset"]], newdata )
    }

    ## get matrix of coordinates of newdata for point kriging

    coords.s <- as.data.frame(switch(
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
    ))

    if( !is.null( coords.s ) &&
      !(
        NCOL( coords.d ) == NCOL( coords.s ) &&
        all( colnames( coords.s ) ==
          colnames(object[["locations.objects"]][["coordinates"]])
        )
      )
    ) stop(
      "inconsistent number and/or names of coordinates in 'object' and in 'newdata'"
    )

    ## compute unconditional mean

    mn.s <- drop( X.s %*% trend.coef ) + offset

    ## compute conditional mean

    if( control[["condsim"]] ){

      ## simple kriging weights

      skw <- simple.kriging.weights(
        pred.coords = coords.s,
        object = object,
        type = type,
        covariances = FALSE,
        control = control.predict.georob(
          mmax = control[["mmax"]],
          ncores = control[["ncores"]]
        )
      )

      ## conditional mean

      mn.s <- mn.s + drop( skw %*% ( y.d - mn.d ) )

    }

    ## re-convert newdata to SpatialGridDataFrame

    if( fullgrid.newdata ) fullgrid( newdata ) <- TRUE


  } else {


#### ---- using fitted trend coefficients -------------

    if( !control[["condsim"]] ){

      ## compute unconditional mean by predict.georob

      tmp <- predict(
        object, newdata = newdata, locations = locations,
        type = "trend"
      )

    } else {

      ## compute conditional mean by predict.georob

      tmp <- predict(
        object, newdata = newdata, locations = locations,
        type = type
      )

    }

    ## convert to dataframe if newdata is a Spatial... object

    if( !identical( class(tmp)[1], "data.frame") ) tmp <- as.data.frame( tmp )

    ## get coordinates, (un-)conditional mean

    coords.s <- tmp[, attr(terms(locations), "term.labels"), drop = FALSE]
    mn.s     <- tmp[, "pred"]

  }

  ## omit any duplicated simulation sites

  key.s <- apply( coords.s, 1, paste, collapse = " " )
  dup.s <- duplicated( key.s )

  if( any( dup.s ) ){

    warning( "duplicated simulation sites eliminated" )

    key.s <- key.s[!dup.s]
    coords.s <- coords.s[!dup.s, , drop = FALSE]
    mn.s <- mn.s[!dup.s]
    if( !is.null( covariates.s ) ) covariates.s <- covariates.s[!dup.s, , drop = FALSE]

  }


#### -- covariates for data and simulation sites -------------

  ## covariates for trend model if required

  covariates.d <- NULL
  covariates.s <- NULL

  if( control[["trend.covariates"]] ){

    covariates.d <- get_all_vars(
      formula.fixef,  data = eval( getCall(object)[["data"]] )
    )
    covariates.d <- covariates.d[!duplicated( object[["Tmat"]] ), , drop = FALSE]

    covariates.s <- get_all_vars( formula.fixef,  data = as.data.frame( newdata ) )

  }

  ## further preparations

  include.data.sites <- control[["include.data.sites"]]
  if( control[["condsim"]] ) include.data.sites <- TRUE


#### -- simulation by cholesky decomposition using exact coordinates of sites -------------

  if( !control[["use.grid"]] ){

    if( verbose > 0 ) cat(
      "  simulate by Cholesky decomposion of covariace matrix (exact coordinates) ...\n"
    )


#### ---- prepare indices of sites for which simulations are computed -------------

    ## indices of sites of support data for which simulations are generated

    if( include.data.sites ){
      i.d <- 1:NROW( coords.d )
    } else {
      i.d <- integer(0L)
    }

   ## indices of simulation sites for which simulation results are reported

    i.s <- 1:NROW( coords.s )
    i.s.d <- integer(0L)  ## indices of eliminated simulation sites in coords.s
    i.d.s <- integer(0L)  ## indices of eliminated simulation sites in coords.d

    if( include.data.sites ){

      ## indices of simulation sites if some simulation sites coincide
      ## with sites of support data

      key.d <- apply( coords.d, 1, paste, collapse = " " )
      i.s   <- which( !key.s %in% key.d )  ## not coinciding, will be kept
      i.s.d <- which(  key.s %in% key.d )  ## coinciding, will be omitted

      ## indices of sites of support data that coincide with simulation
      ## sites

      i.d.s <- which( key.d %in% key.s )

      if( length( i.s.d ) && control[["include.data.sites"]] ) warning(
        "some simulations sites coincide with sites of support data",
        " and are therefore omitted"
      )



    }


#### ---- simulate unconditional zero mean realizations

    sim.values <- sim.chol.decomp(
      coords = rbind(
        coords.d[i.d, , drop = FALSE],
        coords.s[i.s, , drop = FALSE]
      ),
      type = type, variogram.object = object[["variogram.object"]],
      nsim = nsim, seed = seed,
      control.pcmp = control[["pcmp"]], verbose = verbose
    )

    colnames( sim.values ) <- paste0( "sim.", seq_len( nsim) )


#### ---- add mean or condition to support data

    if( control[["include.data.sites"]] ){
      ii.d <- i.d
    } else {
      ii.d <- i.d.s
    }

    ii.s <- length( i.d ) + ( 1:length( i.s ) )

    if( !control[["condsim"]] ){

      ## unconditional simulation: add unconditional mean

      sim.values <- sim.values[c(ii.d, ii.s), , drop = FALSE] +
        c( mn.d[ii.d], mn.s[i.s] )

    } else {

      ## conditional simulation: compute conditional error and add
      ## conditional mean

      ## simple kriging weights

      skw <- simple.kriging.weights(
        pred.coords = rbind(
          coords.d[ii.d, , drop = FALSE],
          coords.s[i.s, , drop = FALSE]
        ),
        object = object,
        type = type,
        covariances = control[["covariances"]],
        control = control.predict.georob(
          mmax = control[["mmax"]],
          ncores = control[["ncores"]]
        )
      )

      ## store covariance matrices

      if( control[["covariances"]] ){
        t.gcvmat <- skw[["gcvmat"]]
        t.gcvmat.pred <- skw[["gcvmat.pred"]]
        skw <- skw[["skw"]]
      }

      ## conditional errors

      tmp <- sim.values[c(ii.d, ii.s), , drop = FALSE] -
        skw %*% sim.values[i.d, , drop = FALSE]

      ## conditional realizations

      sim.values <- c( mn.d[ii.d], mn.s[i.s] ) + tmp

    }


#### ---- prepend coordinates and optionally covariates and (un-)conditional mean

    tmp <- rbind(
      coords.d[ii.d, , drop = FALSE],
      coords.s[i.s, , drop = FALSE]
    )
    if( control[["means"]] ) tmp <- cbind(
      tmp, expct = c( mn.d[ii.d], mn.s[i.s] )
    )
    if( control[["trend.covariates"]] ) tmp <- cbind(
      tmp, rbind(
        covariates.d[ii.d, , drop = FALSE],
        covariates.s[i.s, , drop = FALSE]
      )
    )
    result <- cbind( tmp, sim.values )

    ## rearrange rows if simulation sites coincided with support data sites

    if(
      control[["condsim"]] && !control[["include.data.sites"]] &&
      length( i.s.d )
    ){
      result <- result[order( c(i.s.d, i.s) ), , drop = FALSE]
      n.d <- 0L
      n.s <- length( c(i.s.d, ii.s) )
    } else {
      n.d <- length( ii.d )
      n.s <- length( i.s )
    }

    ## ---------- end !control[["use.grid"]]

  } else {


#### -- simulation by circular embedding using simulation grid -------------

    if( verbose > 0 ) cat( "  simulate with coordinates assigned to grid  ...\n" )

    if( verbose > 1 ) cat( "   assign points to grid nodes ...\n" )


#### ---- define grid for simulating (conditional) realizations

    ## (refined) grid increment

    incr.grid <- sapply(
      1:NCOL( coords.s ),
      function( i ){
        x.s <- sort( unique( coords.s[, i] ) )
        if( length( x.s ) > 2L ){
          incr <- min( diff( x.s ))
        } else {
          x.d <- sort( unique( coords.d[, i] ) )
          if( length( x.d ) > 2L ){
            incr <- min( diff( x.d ))
          } else incr <- NA_real_
        }
        incr / control[["grid.refinement"]]
      }
    )

    if( any( is.na( incr.grid ) ) ) stop(
      "error when generating simulation grid: some coordinates are constant ",
       "for both simulation and support data sites;\n  omit these coordinates in the data objects"
    )

    ## extremal coordinates of support data and simulation sites

    range.coords.d <- sapply( coords.d, range, na.rm = TRUE )
    range.coords.s <- sapply( coords.s, range, na.rm = TRUE )

    ## nodes of simulation grid covering support data and simulation sites

    grid.nodes <- lapply(
      1:NCOL( coords.d ),
      function(i, incr, rc.d, rc.s){
        less <- min(0, floor( ( rc.d[1, i] - rc.s[1, i] ) / incr[i] ) )
        more <- max(0, ceiling( ( rc.d[2, i] - rc.s[2, i] ) / incr[i] ) )
        rc.s[, i] <- rc.s[, i] + c(less, more) * incr[i]
        seq(rc.s[1, i], rc.s[2, i], by = incr[i] )
      }, incr = incr.grid, rc.d = range.coords.d, rc.s = range.coords.s
    )
    names( grid.nodes ) <- colnames( coords.d )

    ## generate all nodes simulation grid

    sim.grid <- expand.grid( grid.nodes )


#### ---- assign sites to grid nodes

    ## convert to SpatialPixelsDataFrame and SpatialPoints for overlay

    sg <- sim.grid
    coordinates( sg )  <- locations
    gridded( sg )      <- TRUE

    c.d <- coords.d
    c.s <- coords.s
    coordinates( c.d ) <- locations
    coordinates( c.s ) <- locations

    ## determine indices of grid nodes closest to coords.d and coords.s

    i.d <- over( c.d, sg )
    i.s <- over( c.s, sg )

    ## check for data sites that could not be assigned to grid nodes

    if( sum( sel <- is.na( i.d ) ) ){
      stop(
        "some support data sites (rows ",
        paste( which( sel ), collapse = ", "),
        ") could not be assigned grid nodes; refit spatial model without these sites"
      )
    }

    ## omit simulation sites that could not be assigned to grid nodes

    if( sum( sel <- is.na( i.s) ) ){
      warning(
        "some simulation sites could not be assigned grid nodes and will be omitted"
      )
      mn.s      <- mn.s[!sel]
      coords.s  <- coords.s[!sel, ]
      i.s <- i.s[!sel]
      if( !is.null( covariates.s ) ) covariates.s <- covariates.s[!sel, , drop = FALSE]
    }

    ## redefine i.d if simulations are not required for support data sites

    if( !include.data.sites ) i.d <- integer(0L)

    ## grid indices of simulation sites for which simulation results are reported

    ii.s  <- i.s
    i.s.d <- integer(0L)  ## grid indices of eliminated simulation sites in coords.s
    i.d.s <- integer(0L)  ## grid indices of eliminated simulation sites in coords.d

    if( include.data.sites ){

      ## handle grid indices of sites of support data that coincide with
      ## simulation sites

      i.d.s <- i.d[i.d %in% i.s]

      ## handle grid indices of simulation sites that coincide with sites
      ## of support data

      i.s.d <- i.s[ i.s %in% i.d] ## coinciding, will be omitted
      ii.s  <- i.s[!i.s %in% i.d] ## not coinciding, will be kept

      if( length( i.s.d ) ){

        if( control[["include.data.sites"]] ) warning(
          "some simulations and support data sites are assigned to the same grid nodes;",
          " these simulation sites will be omitted"
        ) else if( control[["condsim"]] && identical( type, "response" ) ) warning(
          "some simulations and support data sites are assigned to the same grid nodes;",
          " conditionally simulated values at these sites are equal to observations"
        )
      }
    }


#### ---- simulate unconditional zero mean realizations

    sim.values <- sim.circulant.embedding(
      grid.nodes = grid.nodes,
      type = type, variogram.object = object[["variogram.object"]],
      nsim = nsim, seed = seed, ce.method = control[["ce.method"]],
      ce.grid.expansion = control[["ce.grid.expansion"]],
      ncores = control[["ncores"]], control.pcmp = control[["pcmp"]],
      verbose = verbose
    )

    colnames( sim.values ) <- paste0( "sim.", seq_len( nsim) )


#### ---- add mean or condition to support data

    ## grid indices of included sites of support data

    if( control[["include.data.sites"]] ){
      ii.d <- i.d
    } else {
      ii.d <- i.d.s
    }

    ## indices of included elements in coords, mn and covariates for
    ## support data and simulation sites

    jj.d <- which( i.d %in% ii.d )
    jj.s <- which( i.s %in% ii.s )

    if( !control[["condsim"]] ){

      ## unconditional simulation: add unconditional mean

      sim.values <- sim.values[c(ii.d, ii.s), , drop = FALSE] +
      c( mn.d[jj.d], mn.s[jj.s] )

    } else {

      ## conditional simulation: compute conditional error and add
      ## conditional mean

      ## simple kriging weights

      skw <- simple.kriging.weights(
        pred.coords = rbind(
          coords.d[jj.d, , drop = FALSE],
          coords.s[jj.s, , drop = FALSE]
        ),
        object = object,
        type = type,
        covariances = control[["covariances"]],
        control = control.predict.georob(
          mmax = control[["mmax"]],
          ncores = control[["ncores"]]
        )
      )

      ## store covariance matrices

      if( control[["covariances"]] ){
        t.gcvmat <- skw[["gcvmat"]]
        t.gcvmat.pred <- skw[["gcvmat.pred"]]
        skw <- skw[["skw"]]
      }

      ## conditional errors

      tmp <- sim.values[c(ii.d, ii.s), , drop = FALSE] -
      skw %*% sim.values[i.d, , drop = FALSE]

      ## conditional realizations

      sim.values <- c( mn.d[jj.d], mn.s[jj.s] ) + tmp

    }


#### ---- prepend coordinates and optionally covariates and (un-)conditional mean

    tmp <- sim.grid[c(ii.d, ii.s), , drop = FALSE]
    if( control[["means"]] ) tmp <- cbind(
      tmp, expct = c( mn.d[jj.d], mn.s[jj.s] )
    )
    if( control[["trend.covariates"]] ) tmp <- cbind(
      tmp, rbind(
        covariates.d[jj.d, , drop = FALSE],
        covariates.s[jj.s, , drop = FALSE]
      )
    )
    result <- cbind( tmp, sim.values )

    ## rearrange rows if simulation sites coincided with support data sites

    if(
      control[["condsim"]] && !control[["include.data.sites"]] &&
      length( i.s.d )
    ){
      k <- match( i.s, c(i.s.d, ii.s) )
      result <- result[k, , drop = FALSE]
      n.d <- 0L
      n.s <- length( i.s )
    } else {
      n.d <- length( jj.d )
      n.s <- length( jj.s )
    }

  }     ## ---------- end control[["use.grid"]]


#### -- finalizing output -------------

  ## convert result to same class as newdata

  result <- switch(
    class( newdata )[1],
    "data.frame" = result,
    "SpatialPoints" = ,
    "SpatialPointsDataFrame" = {
      coordinates( result ) <- locations
      result
    },
    "SpatialPixels" = ,
    "SpatialPixelsDataFrame" = {
      coordinates( result ) <- locations
      if( !control[["include.data.sites"]] ){
        gridded( result ) <- TRUE
      }
      result
    },
    "SpatialGrid" = ,
    "SpatialGridDataFrame" = {
      coordinates( result ) <- locations
      if( !control[["include.data.sites"]] ){
        gridded( result ) <- TRUE
        fullgrid( result ) <- TRUE
      }
      result
    }
  )

  attr( result, "n" ) <- c( n.d = n.d, n.s = n.s )

  if( control[["covariances"]] & control[["condsim"]] ){

    attr( result, "gcvmat.d.d" ) <- t.gcvmat
    attr( result, "gcvmat.s.d" ) <- t.gcvmat.pred

  }

  result

}


##  ###########################################################################

### sim.chol.decomp -------------

sim.chol.decomp <- function(
  coords,
  type, variogram.object, nsim, seed,
  control.pcmp, verbose
){

  ## function computes unconditional realizations of Gaussian random
  ## fields by the Cholesky matrix decomposition method, cf.  Davies (1987)

  ## 2024-01-06 A. Papritz


#### -- prepare required objects

  n.coords <- NROW( coords )


#### -- generalized covariance matrix

  ## lag vectors for all distinct pairs of points in coords

  if( !all(
      sapply( variogram.object, function(x) x[["isotropic"]] )
    )
  ){
    lag.vectors <- sapply(
      coords,
      function( x ){
        nx <- length( x )
        tmp <- matrix( rep( x, nx ), ncol = nx)
        sel <- lower.tri( tmp )
        tmp[sel] - (t( tmp ))[sel]
      }
    )
  } else {
    lag.vectors <- as.vector( dist( coords ) )
  }
  attr( lag.vectors, "ndim.coords" ) <- NCOL( coords )

  ## list with generalized correlation matrices of signal

  Valpha <- f.aux.gcr(
    lag.vectors = lag.vectors,
    variogram.object = variogram.object,
    gcr.constant = NULL,
    symmetric = TRUE,
    control.pcmp = control.pcmp,
    verbose = verbose
  )

  if( any( sapply( Valpha, function(x) x[["error"]] ) ) ) stop(
    "error in computing generalized correlation matrix"
  )

  ## initialize generalized covariance matrix

  Sigma <- list(
    diag = rep( 0., length( Valpha[[1]][["Valpha"]][["diag"]] ) ),
    tri =  rep( 0., length( Valpha[[1]][["Valpha"]][["tri"]] ) )
  )
  attr( Sigma, "struc" ) <- "sym"

  ## loop over all list components of Valpha

  for( i in 1:length( variogram.object ) ){

    ## check whether variogram is stationary or has a locally stationary
    ## representation

    if(
      variogram.object[[i]][["variogram.model"]] %in%
      control.georob()[["irf.models"]][!control.georob()[["irf.models"]] %in% "RMfbm"] &&
      verbose > 0
    ){
      warning(
        "simulation of an intrinsic random field without a known ",
        "locally stationary representation"
      )
    }

    ## sill parameters

    param <- variogram.object[[i]][["param"]]

    psill <- tsill <- param[["variance"]]
    if( "snugget" %in% names( param ) ){
      tsill <- tsill + param["snugget"]
    }

    ## convert correlations of signal to covariances

    Valpha[[i]][["Valpha"]][["diag"]] <- tsill *
      Valpha[[i]][["Valpha"]][["diag"]]
    Valpha[[i]][["Valpha"]][["tri"]]  <- psill *
      Valpha[[i]][["Valpha"]][["tri"]]

    ## handle nugget (measurement error variance)

    if( "nugget" %in% names( param ) && identical( type, "response" ) ){

      ## simulation of response: add nugget to diagonal elements
      ## corresponding to simulation sites

      Valpha[[i]][["Valpha"]][["diag"]] <-
        Valpha[[i]][["Valpha"]][["diag"]] + param["nugget"]

    }

    ## add up generalized covariances

    Sigma[["diag"]] <- Sigma[["diag"]] + Valpha[[i]][["Valpha"]][["diag"]]
    Sigma[["tri"]]  <- Sigma[["tri"]]  + Valpha[[i]][["Valpha"]][["tri"]]

  }

  ## expand Sigma to matrix

  Sigma <- expand( Sigma )


#### -- simulate realizations

  set.seed( seed )

  ## Cholesky decomposition of Sigma

  U <- try( chol( Sigma ), silent = TRUE )

  if( inherits( U, "try-error" ) ) stop(
    "generalized covariance matrix not positive definite"
  )

  ## simulation of vector with iid N(0, 1) random values for simulation sites

  w <- matrix( rnorm( nsim * n.coords, mean = 0., sd = 1.), ncol = nsim )

  ## compute L_22 %*% w, cf. Davis, 1987, p. 96

  result <- t( U ) %*% w

  ## return simulated values

  return( result )

}

##  ###########################################################################

### sim.circulant.embedding -------------

sim.circulant.embedding <- function(
  grid.nodes,
  type, variogram.object, nsim, seed,
  ce.method, ce.grid.expansion,
  ncores, control.pcmp, verbose
){

  ## function computes unconditional realizations of Gaussian random fields
  ## by the circular embedding method, cf.  Davies and Bryant (2013), Gneiting
  ## et al.  (2006a), Stein (2002), Wood and Chan, (1994)


  ## 2024-01-21 A. Papritz
  ## 2024-01-29 AP minor changes in handling negative eigenvalues of base matrix
  ## 2024-02-01 AP setting random seeds for parallel processing
  ## 2024-02-01 AP saving SOCKcluster.RData to tempdir()

#### -- auxiliary function to simulate realizations

  f.aux.sim <- function(
    i
  ){

    ## objects taken from parent environment
    ## rs, re,
    ## ce.n.nodes, ce.n.nodes.total, snt,
    ## n.nodes, slambda


    ## number of relizations to simulate

    nsim <- re[i] - rs[i] + 1L

    ## dataframe with realizations of standard normal random numbers

    stn <- as.data.frame(
      matrix(
        rnorm( nsim * ce.n.nodes.total ),
        nrow = ce.n.nodes.total
      )
    )

    ## simulate realizations

    sapply(
      stn,
      function( x ){

        ## convert x to array if necessary
        if( length( ce.n.nodes ) > 1L ) dim( x ) <- ce.n.nodes

        ## simulate realization
        tmp <- slambda * ( fft( x ) / snt )
        tmp <- Re( fft( tmp, inverse = TRUE ) / snt )

        ## select realizations for simulation sites on target grid and
        ## convert to vector

        as.vector(
          switch(
            length( n.nodes ),
            tmp[1:n.nodes[1]],
            tmp[1:n.nodes[1], 1:n.nodes[2]],
            tmp[1:n.nodes[1], 1:n.nodes[2], 1:n.nodes[3]]
          )
        )

      }

    )

  }

#### -- extend simulation grid such that number of nodes are highly composite numbers

  ## highly composite numbers, cf.
  ## https://en.wikipedia.org/wiki/Highly_composite_number,
  ## https://gist.github.com/dario2994/fb4713f252ca86c1254d

  hcn <- c(
    4L, 6L, 12L, 24L, 36L, 48L, 60L, 120L, 180L, 240L, 360L, 720L, 840L,
    1260L, 1680L, 2520L, 5040L, 7560L, 10080L, 15120L, 20160L, 25200L, 27720L, 45360L,
    50400L, 55440L, 83160L, 110880L, 166320L, 221760L, 277200L, 332640L, 498960L,
    554400L, 665280L, 720720L, 1081080L, 1441440L, 2162160L, 2882880L, 3603600L,
    4324320L, 6486480L, 7207200L, 8648640L, 10810800L, 14414400L, 17297280L,
    21621600L, 32432400L, 36756720L, 43243200L, 61261200L, 73513440L
  )

  ## grid increments

  grid.incr <- sapply(
    grid.nodes,
    function( x ){
      dx <- unique( diff( x ) )
      stopifnot( identical( length(dx), 1L ) )
      dx
    }
  )

  ## highly composite number of grid nodes

  n.nodes <- sapply( grid.nodes, length ) * ce.grid.expansion

  ex.n.nodes <- sapply(
    n.nodes,
    function( x, hcn ){
      stopifnot( x < max( hcn ) )
      hcn[ which( hcn / x >= 1. )[1] ]
    }, hcn = hcn
  )

  ## nodes of extended grid

  ex.grid.nodes <- lapply(
    1:length( grid.nodes ),
    function( i, grid.nodes, grid.incr, n.nodes ){
      seq(
        from = grid.nodes[[i]][1], by = grid.incr[i],
        length.out = n.nodes[i]
      )
    },
    grid.nodes = grid.nodes, grid.incr = grid.incr,
    n.nodes = ex.n.nodes
  )
  names( ex.grid.nodes ) <- names( grid.nodes )


#### -- construct circular embedding grid

  ## number of grid nodes

  ce.n.nodes <- ceiling( ex.n.nodes * 2L )
  ce.n.nodes.total <- prod( ce.n.nodes )

  ## coordinates of nodes of circular embedding grid

  ce.grid.nodes <- lapply(
    1:length( grid.nodes ),
    function( i, grid.nodes, grid.incr, n.nodes ){
      seq(
        from = grid.nodes[[i]][1], by = grid.incr[i],
        length.out = n.nodes[i]
      )
    },
    grid.nodes = grid.nodes, grid.incr = grid.incr,
    n.nodes = ce.n.nodes
  )
  names( ce.grid.nodes ) <- names( grid.nodes )

  ## extent of circular embedding grid

  ce.grid.extent <- sapply(
    ce.grid.nodes,
    function( x ) diff( range( x ) )
  ) + grid.incr


  ## dataframe with coordinates of all circular embedding grid nodes

  ce.grid <- as.matrix( expand.grid( ce.grid.nodes ) )
  colnames( ce.grid ) <- names( ce.grid.nodes )


#### -- lag distances
  ## (first row of distance matrix on circular embedding grid )

  ## differences of coordinates between lower left node of ce.grid with
  ## all nodes of ce.grid, cf. Davies and Bryant (2013) pp. 7

  tmp1 <- t( ce.grid ) - ce.grid[1, ]

  ## differences of coordinates between upper right node of ce.grid with
  ## all nodes plus grid increment

  tmp2 <- ce.grid.extent - tmp1
  #     bla <- ce.grid[172800, ] - t( ce.grid ) + grid.incr
  #     range(bla - tmp2)

  ## select shorter of the 2 differences of coordinates

  sel2 <- tmp2 < tmp1
  #   bla <- colSums(sel2)

  lag.vectors <- sapply(
    1:NROW( tmp1 ),
    function( i, dist1, dist2, sel2 ){
      dist1[i, sel2[i, ]] <- - dist2[i, sel2[i, ]]
      dist1[i, ]
    },
    dist1 = tmp1, dist2 = tmp2, sel2 = sel2
  )

  #     bla <- sqrt( rowSums( tmp^2 ) )
  #
  #     mygrid <- mygrid.prep(ex.grid.nodes, 180, 240)
  #     Rx <- mygrid$M.ext * mygrid$cell.width
  #     Ry <- mygrid$N.ext * mygrid$cell.height
  #     m.abs.diff.row1 <- abs(mygrid$mcens.ext[1] - mygrid$mcens.ext)
  #     m.diff.row1 <- pmin(m.abs.diff.row1, Rx - m.abs.diff.row1)
  #     n.abs.diff.row1 <- abs(mygrid$ncens.ext[1] - mygrid$ncens.ext)
  #     n.diff.row1 <- pmin(n.abs.diff.row1, Ry - n.abs.diff.row1)
  #     cent.ext.row1 <- expand.grid(m.diff.row1, n.diff.row1)
  #     D.ext.row1 <- sqrt(cent.ext.row1[, 1]^2 + cent.ext.row1[, 2]^2)
  #     range(bla - D.ext.row1)

  if( all(
      sapply( variogram.object, function(x) x[["isotropic"]] )
    )
  ){
    ## isotropic autocorelations
    lag.vectors <- sqrt( rowSums( lag.vectors^2 ) )
  }
  attr( lag.vectors, "ndim.coords" ) <- length( ce.n.nodes )


#### -- generalized covariance matrix (standard embedding)
  ## (base matrix, first row of generalized covariance matrix for circulant
  ## embedding grid)

  ## list with generalized correlation matrices of signal

  Valpha <- f.aux.gcr(
    lag.vectors = lag.vectors,
    variogram.object = variogram.object,
    gcr.constant = NULL,
    symmetric = FALSE,
    control.pcmp = control.pcmp,
    verbose = verbose
  )

  if( any( sapply( Valpha, function(x) x[["error"]] ) ) ) stop(
    "error in computing generalized correlation matrix"
  )

  ## initialize generalized covariance matrix

  Sigma <- rep( 0., ce.n.nodes.total )

  ## loop over all list components of Valpha

  for( i in 1:length( variogram.object ) ){

    ## check whether variogram is stationary or has a locally stationary
    ## representation

    if(
      variogram.object[[i]][["variogram.model"]] %in%
      control.georob()[["irf.models"]][!control.georob()[["irf.models"]] %in% "RMfbm"] &&
      verbose > 0
    ){
      warning(
        "simulation of an intrinsic random field without a known ",
        "locally stationary representation"
      )
    }

    ## sill parameters

    param <- variogram.object[[i]][["param"]]

    psill <- tsill <- param[["variance"]]
    if( "snugget" %in% names( param ) ){
      tsill <- tsill + param["snugget"]
    }

    ## convert correlations of signal to covariances

    Valpha[[i]][["Valpha"]][ 1] <- tsill * Valpha[[i]][["Valpha"]][ 1]
    Valpha[[i]][["Valpha"]][-1] <- psill * Valpha[[i]][["Valpha"]][-1]

    ## handle nugget (measurement error variance)

    if( "nugget" %in% names( param ) && identical( type, "response" ) ){

      ## simulation of response: add nugget to diagonal elements
      ## corresponding to simulation sites

      Valpha[[i]][["Valpha"]][1] <-
        Valpha[[i]][["Valpha"]][1] + param["nugget"]

    }

    ## add up generalized covariances

    Sigma <- Sigma + Valpha[[i]][["Valpha"]]

  }

  ## convert to matrix or array

  if( length( ce.n.nodes ) > 1L ) dim( Sigma ) <- ce.n.nodes


#### -- eigenvalues of base matrix

  Lambda <- as.vector( Re( fft( Sigma, inverse = TRUE ) ) )

  ## handle negative eigenvalues

  sel.neg   <- Lambda < 0.
  sel.small <- abs( Lambda ) < sqrt( .Machine[["double.eps"]] )

  ## set small negative eigenvalues equal to zero

  Lambda[sel.neg & sel.small] <- 0.

  ## handle large negative eigenvalues

  if( any( sel <- ( sel.neg & !sel.small ) ) ){

    ## some eigenvalues are negative

    if( identical( ce.method, "approximate" ) ){

      ## approximate circulant embedding: set negative eigenvalues to zero
      ## and make variance adjustment, cf.  Chan & Wood, 1994, eq.  4.2

      tmp <- Lambda
      tmp[sel] <- 0.

      Lambda <- tmp * sum( Lambda ) / sum( tmp )

      if( length( ce.n.nodes ) > 1L ) dim( Lambda ) <- ce.n.nodes

    } else {

      if(verbose > 0 ){
        cat( "\n  summary of eigenvalues of base matrix:\n" )
        print( summary( c( Lambda ) ) )
        cat( "\n" )
      }

      cat(
        "  some eigenvalues of base matrix are negative,",
        "standard circulant embedding is not possible;\n",
        " enlarge simulation domain by increasing the parameter 'ce.grid.expansion' or\n",
        " use another circulant embedding method ('ce.method' %in% c(approximate))\n"
#         " use another circulant embedding method ('ce.method' %in% c(approximate, cutoff, intrinsic))\n"
      )

      stop()

    }
  }

#### -- simulate realizations

  ## initialize random number generation

  old.kind <- RNGkind( "L'Ecuyer-CMRG" )
  on.exit( RNGkind( old.kind[1] ) )

  ## sqrt of total number of nodes of circulant embedding grid

  snt <- sqrt( ce.n.nodes.total )

  ## sqrt of eigenvalues of base matrix

  slambda <- sqrt( Lambda )

  ## prepare parallel execution

  ncores <- min( nsim, ncores )
  n.part <- ceiling( nsim / ncores )
  rs <- ( 0L:(ncores - 1L) ) * n.part + 1L
  re <- ( 1L:ncores ) * n.part; re[ncores] <- nsim

  ## compute the realizations for all the parts

  if( ncores > 1L && !control.pcmp[["fork"]] ){

    ## create a PSOCK cluster on windows OS

    fname <- file.path( tempdir(), "SOCKcluster.RData" )

    clstr <- makePSOCKcluster( ncores )
    save( clstr, file = fname )
    old.opt <- options( error = f.stop.cluster )
    on.exit( options( old.opt ) )

    ## export required items to workers

    junk <- clusterExport(
      clstr,
      c(
        "rs", "re",
        "ce.n.nodes", "ce.n.nodes.total", "snt",
        "n.nodes", "slambda"
      ),
      envir = environment()
    )

    ## set random seeds for workers

    junk <- clusterSetRNGStream(cl = clstr, seed)

    ## compute simulations on workers

    result <- try(
      parLapply(
        clstr,
        1L:ncores,
        f.aux.sim
      )
    )

    f.stop.cluster( clstr, fname )

  } else {

    ## fork child processes on non-windows OS

    ## set random seed

    set.seed( seed )

    ## compute simulations on children

    result <- try(
      mclapply(
        1L:ncores,
        f.aux.sim,
        mc.cores = ncores, mc.set.seed = TRUE,
        mc.allow.recursive = control.pcmp[["allow.recursive"]]
      )
    )

  }

  has.error <- sapply(
    result, function( x ) inherits( x, "try-error" )
  )

  if( any( has.error ) ){

    cat( "\nerror(s) occurred when computing realizations in parallel:\n\n" )
    sapply( result[has.error], cat)
    cat( "\nuse 'ncores=1' to avoid parallel computations and to see where problem occurs\n\n" )
    stop()

  } else {

    ## convert list to matrix and return result

    matrix( unlist( result ), ncol = nsim )

  }

}
