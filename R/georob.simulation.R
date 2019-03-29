##  ###########################################################################

control.condsim <- function(
  use.grid = FALSE,
  grid.refinement = 2.,
  condsim = TRUE,
  include.data.sites = FALSE,
  means = FALSE,
  trend.covariates = FALSE,
  covariances = FALSE,
  ncores = detectCores(),
  pcmp = control.pcmp()
){
  
  ## auxiliary function to set meaningful default values for 
  ## condsim
  
  ## 2018-01-12 A. Papritz
  
  if( grid.refinement < 1. ) warning( "grid refinement factor < 1")
  
  list(
    use.grid = use.grid,
    grid.refinement = grid.refinement,
    condsim = condsim,
    include.data.sites = include.data.sites,
    means = means,
    trend.covariates = trend.covariates,
    covariances = covariances,
    ncores = ncores, pcmp = pcmp
  )
  
}


##  ###########################################################################

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
  ## 2018-11-25	A. Papritz small change in function f.aux.sim.1
  ## 2019-03-29 A. Papritz restore default for RFoption spConform
  ## 2019-03-29 A. Papritz conditioning data is passed as SpatialPointsDataFrame to
  ##                       RFsimulate (only applies for exact coordinates)
  
  
  ## auxiliary function for aggregating a response variable x for the
  ## levels of the variable by.one
  
  f.aggregate.by.one <- function(x, by.one, FUN = mean, ...){ 
    by.one <- as.numeric( factor( by.one, levels = unique( by.one ) ) ) 
    aggregate( x, list( by = by.one ), FUN = FUN, ...  )[, -1]
  }
  
  
  ## match arguments
  
  type <- match.arg( type )
  
  if( verbose > 0 ) cat( "  prepare data ...\n" )
  
  
  ## check consistency of provided arguments
  
  ## conditionally simulating signal
  
  if( identical(type, "signal") && control[["condsim"]] && !control[["use.grid"]] ) stop(
    "conditional simulation of signal requires control argument 'use.grid = TRUE'"
  )
  
  ## mandatory argument
  
  if( missing( object ) || missing( newdata ) || missing( nsim ) || missing( seed ) ) stop(  
    "some mandatory arguments are missing" 
  )
  
  ## class
  
  if( !identical( class( object ) , "georob") ) stop(
    "object must be of class georob"
  )
  
  if( identical( class( newdata ) , "SpatialPolygonsDataFrame") ) stop(
    "newdata is a SpatialPolygonsDataFrame; simulation for such targetst is not yet implemented"
  )
  
  ## coordinates only as covariates for trend model if newdata is a Spatialxx object
  
  if( class( newdata ) %in% c( "SpatialPoints", "SpatialPixels", "SpatialGrid" ) ){
    t.formula <- as.formula( paste( as.character( formula( object ) )[-2L], collapse = "" ) )
    tmp <- try(
      get_all_vars( t.formula, as.data.frame( coordinates( newdata ) ) ),
      silent = TRUE
    )
    if( identical( class( tmp ), "try-error" ) ) stop(
      "'newdata' is a SpatialPoints, SpatialPixels or SpatialGrid object\n but drift covariates are not functions of coordinates"
    )
  }
  
  ## mean parameters
  
  if( !is.null( trend.coef) && !identical(length(trend.coef), length(coef(object)) ) ) stop(
    "'trend.coef' inconsistent with trend model of 'object'"
  )
  
    
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
  

  ## extract required items about data sites (support data) from object
  
  ## terms
  
  tt <- terms( object )
  Terms <- delete.response( tt )
  
  ## model.frame
  
  mf.d <- model.frame( object )
  
  ## design matrix
  
  X.d <- model.matrix( object )[!duplicated( object[["Tmat"]] ), , drop = FALSE]
  
  ## unconditional mean (fitted values)
  
  if( !is.null( trend.coef ) ){
    
    ## deal with non-NULL offset    
    
    if( !is.null( attr( tt, "offset" ) ) || !is.null( object[["call"]][["offset"]] ) ){
      offset <- model.offset( mf.d )[!duplicated( object[["Tmat"]] )]
    } else offset <- rep( 0., NROW(X.d) )
    
    mean.d <- drop( X.d %*% trend.coef ) + offset
    
  } else {
    
    mean.d <- fitted( object )[!duplicated( object[["Tmat"]] )]
    
  }
  
  ## response (compute arith. average for replicated locations) or trend
  
  if( control[["condsim"]] ){
    
    y.d <- f.aggregate.by.one( model.response( mf.d ), object[["Tmat"]], na.rm = TRUE )
    
  } else {
    
    y.d <- mean.d
    
  }  
  
  ## coordinates
  
  if( missing( locations ) ){
    locations <- object[["locations.objects"]][["locations"]]
  }
  
  Terms.loc <- terms( locations )
  attr( Terms.loc, "intercept" ) <- 0
  
  coords.d <- as.data.frame(
    object[["locations.objects"]][["coordinates"]][!duplicated( object[["Tmat"]] ), , drop = FALSE]
  )
  
  
  ## extract signal variance, xi, nugget and gcr.constant
  
  tmp <- f.reparam.fwd( object[["variogram.object"]] )
  
  var.signal <- attr(tmp, "var.signal" )
  xi <- sapply( tmp, function(x) x[["param"]]["variance"] )
  nugget <- object[["variogram.object"]][[1L]][["param"]]["nugget"]
  
  isotropic <- all( sapply( object[["variogram.object"]], function(x) x[["isotropic"]] ) )
  
  gcr.constant <- lapply( 
    object[["Valphaxi.objects"]][["Valpha"]],
    function(x) x[["gcr.constant"]]
  )
  
  
  ## create variogram object for simulation by RFsimulate
  
  variogram <- lapply(
    object[["variogram.object"]],
    function( x, type ){
      
      param           <- x[["param"]]
      aniso           <- x[c("aniso", "sclmat", "rotmat")]
      
      ## continous part of signal
      
      A <- aniso[["sclmat"]] * aniso[["rotmat"]] / param["scale"]
      
      signal <- list( x[["variogram.model"]] )
      signal <- c( 
        signal, 
        as.list( param[!names(param) %in% c("variance", "snugget", "nugget", "scale")] ) 
      )
      signal <- list( "$", var = param["variance"], A = A, signal )
      
      ## non-continuous part of signal and independent error
      
      nugget <- 0.
      if( "snugget" %in% names(param) ) nugget <- nugget + param["snugget"]
      if( identical( type, "response" ) && "nugget" %in% names(param) ){
        nugget <- nugget + param["nugget"]
      }
      
      list( signal = signal, nugget = nugget )
      
    }, type = type
  )
  
  
  ## covariates for trend model if required
  
  covariates.d <- NULL
  covariates.s <- NULL
  
  if( control[["trend.covariates"]] ){
    
    t.formula <- as.formula( paste( as.character( formula(object))[-2L], collapse = " " ) )
    
    covariates.d <- get_all_vars( t.formula,  data = eval( getCall(object)[["data"]] ) )
    covariates.d <- covariates.d[!duplicated( object[["Tmat"]] ), , drop = FALSE]
    
    covariates.s <- get_all_vars( t.formula,  data = as.data.frame( newdata ) )
    
  }
  
  
  ## prepare required items for simulation sites
  
  if( verbose > 0 ) cat( "  compute (conditional) mean ...\n" )
  
  ## use provided trend function parameters
  
  if( !is.null( trend.coef ) ){
    
    ## get modelframe for newdata
    
    mf.s <- switch(
      class( newdata ),
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
      class( newdata ),
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
          colnames(object[["locations.objects"]][["coordinates"]]) )
      )
    ) stop(
      "inconsistent number and/or names of coordinates in 'object' and in 'newdata'"
    )
    
    ## compute unconditional mean
    
    y.s <- drop( X.s %*% trend.coef ) + offset
    
    ## compute conditional mean
    
    if( control[["condsim"]] ){
      
      ## simple kriging weights
      
      skw <- simple.kriging.weights( 
        pred.coords = coords.s,
        object = object, 
        type = type,
        covariances = FALSE,
        control = control.predict.georob(ncores = control[["ncores"]])
      )
            
      ## conditional mean
      
      y.s <- y.s + drop( skw %*% ( y.d - mean.d ) )
      
    }
  
  } else {
    
    ## compute (conditional) mean by predict.georob
    
    tmp <- predict( 
      object, newdata = newdata, locations = locations,
      type = if( control[["condsim"]] ) "response" else "trend"
    )
    
    ## convert to  dataframe if newdata is a Spatial... object
    
    if( !identical( class(tmp), "data.frame") )  tmp <- as.data.frame(tmp)
    
    ## get coordinates and (conditional) mean
    
    coords.s <- tmp[, attr(terms(locations), "term.labels"), drop = FALSE]
    y.s <- tmp[, "pred"]
    
    rm( tmp ); gc()    
  
  }
  
  ## compute (conditional) realizations ...
  
  if( !control[["use.grid"]] ){
    
    ## ... using exact coordinates of sites
    
    if( verbose > 0 ) cat( "  simulate with exact coordinates ...\n" )
    
    if( control[["include.data.sites"]] && !control[["condsim"]] ){
      
      ## include also data sites
      
      ## omit simulation sites that coincide with data sites
      
      key <-  apply( coords.d, 1, paste, collapse = " " )
      sel <- !apply( coords.s, 1, paste, collapse = " " ) %in% key
      
      coords.s <- coords.s[sel, , drop = FALSE]
      y.s <- y.s[sel]
      if( !is.null( covariates.s ) ) covariates.s <- covariates.s[sel, , drop = FALSE]
      
      ## prepare coordinates and data
      
      n.d <- length( y.d )
      n.s <- length( y.s )
      
      coords <- rbind( coords.d, coords.s )
      y <- c( y.d, y.s )
      if( !is.null( covariates.s ) ) covariates <- rbind( covariates.d, covariates.s )      
      
    } else {
      
      # only simulation sites
      
      n.d <- 0L
      n.s <- length( y.s )
      
      coords <- coords.s
      y <- y.s
      if( !is.null( covariates.s ) ) covariates <- covariates.s    
      
    }
    
    ## prepare conditioning data
    
    if( control[["condsim"]] ){
      data <- cbind( coords.d, values = rep( 0., NROW(coords.d) ) )
      coordinates(data) <- colnames(coords.d)
    } else {
      data <- NULL
    }
    
    
    ## auxiliary function for simulating (conditional) realizations in parallel
    
    f.aux.sim.1 <- function( i ){
      
      ## s, e, sim.seed, variogram, coords, data are used from .GlobalEnv
      
      RFoptions( spConform = FALSE )
      
      set.seed( sim.seed[i] )
      
      nsim <- e[i] - s[i] + 1L
      
      tmp <- lapply(
        variogram,
        
        function( x, nsim ){
          
          signal <- RFsimulate(
            model = x[["signal"]], n = nsim,
            x = coords[, 1], 
            y = if( NCOL(coords) > 1 ) coords[, 2] else NULL,
            z = if( NCOL(coords) > 2 ) coords[, 3] else NULL,
            data = data
          )
          
          if( x[["nugget"]] > 0. ){
            nugget <- rnorm( prod( dim( signal ) ), mean = 0., sd = sqrt( x[["nugget"]] ) )
            dim( nugget ) <- dim( signal )
            nugget + signal
          } else {
            signal
          }
                    
        }, nsim = nsim
      )
      
      res <- tmp[[1]]
      if( length(tmp) - 1L ){
        for( i in 2L:length(tmp) ){
          res <- res + tmp[[i]]
        }
      }
      
      RFoptions( spConform = TRUE )      

      res
      
    }
    
    ## determine number of cores
    
    ncores <- control[["ncores"]]
    
    ## definition of junks to be evaluated in parallel
    
    k <- ncores
    n <- nsim
    dn <- floor( n / k )
    s <- ( (0L:(k-1L)) * dn ) + 1L
    e <- (1L:k) * dn
    e[k] <- n
    
    ## random number generator seed used in junks
    
    set.seed( seed )
    sim.seed <- sample.int( 100000, k)
        
    ## simulate (unconditional) realizations in parallel
    
    if( ncores > 1L && !control[["pcmp"]][["fork"]] ){
      
      if( !sfIsRunning() ){
        options( error = f.stop.cluster )
        junk <- sfInit( parallel = TRUE, cpus = ncores )
      }
      
      #       junk <- sfLibrary( RandomFields, verbose = FALSE )
      junk <- sfExport( "s", "e", "sim.seed", "variogram", "coords", "data" ) 

      res <- sfLapply( 1L:k, f.aux.sim.1)
      
      if( control[["pcmp"]][["sfstop"]] ){
        junk <- sfStop()
        options( error = NULL )
      }
      
    } else {

      res <- mclapply( 1L:k, f.aux.sim.1, mc.cores = ncores )

    }
    
    ## store realizations in a matrix 
    
    sim.values <- res[[1]]
    if( length(res) - 1L ){
      for( i in 2L:length(res) ){
        sim.values <- cbind( sim.values, res[[i]] )
      }
    }
    colnames( sim.values ) <- paste0( "sim.", seq_len( nsim) )
    
    #     sv.exact <<- sim.values
    
    junk <- gc()
    
    tmp <- coords
    if( control[["means"]] ) tmp <- cbind( tmp, expct = y )
    if( control[["trend.covariates"]] ) tmp <- cbind( tmp, covariates )
    result <- cbind( tmp, sim.values + y )
    
    ## add conditioning data
    
    if( control[["include.data.sites"]] & control[["condsim"]] ){
      
      n.d <- length( y.d )
      
      tmp <- coords.d
      if( control[["means"]] ) tmp <- cbind( tmp, expct = y.d )
      if( control[["trend.covariates"]] ) tmp <- cbind( tmp, covariates.d )
      tmp <- cbind( tmp, matrix( rep( y.d, nsim ), ncol = nsim ) )
      colnames( tmp ) <- colnames( result )
      
      result <- rbind( tmp, result )
      
    }
    
    ## ---------- end !control[["use.grid"]]
    
  } else {
    
    ## ... using simulation grid    

    if( verbose > 0 ) cat( "  simulate with coordinates assigned to grid  ...\n" )
    
    if( verbose > 1 ) cat( "   assign grid nodes to points  ...\n" )
    
    ## ... assigning sites to grid
    
    ## define (refined) grid for simulating conditional realizations
    
    ## (refined) grid increment
    
    incr.grid <- sapply(
      coords.s,
      function( x ) min( diff( sort( unique( x ) ) ) )
    ) / control[["grid.refinement"]]
    
    ## extremal coordinates of data and simulation sites
    
    range.coords.d <- sapply( coords.d, range, na.rm = TRUE )
    range.coords.s <- sapply( coords.s, range, na.rm = TRUE )
    
    ## nodes of simulation grid covering data and simulation sites
    
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
    
    ## simulation grid
    
    sim.grid <- expand.grid( grid.nodes )
    
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
    
    ## omit data sites that could not be assigned to grid nodes
    
    if( sum( sel <- is.na( i.d ) ) ){
      warning( "some data sites could not be assigned to grid" )
      y.d 			 <- y.d[!sel]
      coords.d  <- coords.d[!sel, ]
      i.d <- i.d[!sel]
      if( !is.null( covariates.d ) ) covariates.d <- covariates.d[!sel, , drop = FALSE]
    }
    
    ## omit simulation sites that could not be assigned to grid nodes
    
    if( sum( sel <- is.na( i.s) ) ){
      warning( "some simulation sites could not be assigned to grid" )
      y.s 			 <- y.s[!sel]
      coords.s  <- coords.s[!sel, ]
      i.s <- i.s[!sel]
      if( !is.null( covariates.s ) ) covariates.s <- covariates.s[!sel, , drop = FALSE]
    }
    
    ## auxiliary function for simulating unconditional realizations in parallel
    
    f.aux.sim.2 <- function( i ){
      
      ## s, e, sim.seed, variogram, grid.nodes are used from .GlobalEnv
      
      RFoptions( spConform = FALSE )
      
      set.seed( sim.seed[i] )
      
      nsim <- e[i] - s[i] + 1L
      
      tmp <- lapply(
        variogram,
        function( x, grid.nodes, nsim ){
          
          signal <- RFsimulate(
            model = x[["signal"]], grid = TRUE, n = nsim,
            x = grid.nodes[[1]], 
            y = if( length(grid.nodes) > 1 ) grid.nodes[[2]] else NULL,
            z = if( length(grid.nodes) > 2 ) grid.nodes[[3]] else NULL
          )
          
          if( x[["nugget"]] > 0. ){
            nugget <- rnorm( prod( dim( signal ) ), mean = 0., sd = sqrt( x[["nugget"]] ) )
            dim( nugget ) <- dim( signal )
            nugget + signal
          } else {
            signal
          }
          
        }, grid.nodes = grid.nodes, nsim = nsim
      )
      
      res <- tmp[[1]]
      if( length(tmp) - 1L ){
        for( i in 2L:length(tmp) ){
          res <- res + tmp[[i]]
        }
      }
      
      RFoptions( spConform = TRUE )

      res
      
    }
    
    if( verbose > 1 ) cat( "    simulate unconditional realizations  ...\n" )
    
    ## determine number of cores
    
    ncores <- control[["ncores"]]
    
    ## definition of junks to be evaluated in parallel
    
    k <- ncores
    n <- nsim
    dn <- floor( n / k )
    s <- ( (0L:(k-1L)) * dn ) + 1L
    e <- (1L:k) * dn
    e[k] <- n
    
    ## random number generator seed used in junks
    
    set.seed( seed )
    sim.seed <- sample.int( 100000, k)
        
    ## simulate (unconditional) realizations in parallel
    
    if( ncores > 1L && !control[["pcmp"]][["fork"]] ){
            
      if( !sfIsRunning() ){
        options( error = f.stop.cluster )
        junk <- sfInit( parallel = TRUE, cpus = ncores )
      }
      
      #       junk <- sfLibrary( RandomFields, verbose = FALSE )
      junk <- sfExport( "s", "e", "sim.seed", "variogram", "grid.nodes" ) 
      
      res <- sfLapply( 1L:k, f.aux.sim.2 )
      
      if( control[["pcmp"]][["sfstop"]] ){
        junk <- sfStop()
        options( error = NULL )
      }
      
    } else {
      
      res <- mclapply( 1L:k, f.aux.sim.2, mc.cores = ncores )
      
    }
    
    ## store unconditional realizations in a matrix consistent with
    ## dimensions of simulation grid
    
    sim.values <- res[[1]]
    if( length(res) - 1L ){
      for( i in 2L:length(res) ){
        sim.values <- abind( sim.values, res[[i]] )
      }
    }
        
    dim( sim.values ) <- c( NROW( sim.grid ), nsim )
    colnames( sim.values ) <- paste0( "sim.", seq_len( nsim) )
    
    #     sv.grid <<- sim.values[i.s, , drop = FALSE]
    
    junk <- gc()
    
    if( !control[["condsim"]] ){
      
      ## extract unconditionally simulated values for simulation sites
      
      n.s <- length(i.s)
      n.d <- 0L
      
      tmp <- sim.grid[i.s, , drop = FALSE]
      if( control[["means"]] ) tmp <- cbind( tmp, expct = y.s )
      if( control[["trend.covariates"]] ) tmp <- cbind( tmp, covariates.s )
      
      result <- cbind( tmp, sim.values[i.s, , drop = FALSE] + y.s )
      
      if( control[["include.data.sites"]] ){
        
        ## include also data sites
        
        ## omit simulation sites that coincide with data sites
        
        i.s <- i.s[!i.s %in% i.d]
        n.s <- length(i.s)
        n.d <- length(i.d)
        
        tmp <- sim.grid[i.d, , drop = FALSE]
        if( control[["means"]] ) tmp <- cbind(tmp, expct = y.d )
        if( control[["trend.covariates"]] ) tmp <- cbind(tmp, covariates.d )
        
        result <- rbind(
          cbind( tmp, sim.values[ i.d, , drop = FALSE] + y.d ),
          result[i.s, ]
        )
        
      }
      
      #     library(lattice)
      #     levelplot(sim.1 ~ x + y, result, aspect  ="iso")
      
    } else {
      
       ## condition to data by kriging method
      
      if( verbose > 1 ) cat( "    condition to data by kriging method ...\n" )
      
      ## assign grid nodes that were assigned both to simulation and data
      ## locations only to the latter
      
      sel <- !i.s %in% i.d
      
      y.s 			 <- y.s[sel]
      coords.s  <- coords.s[sel, ]
      i.s <- i.s[sel]
      if( !is.null( covariates.s ) ) covariates.s <- covariates.s[sel, , drop = FALSE]
      
      ## average data for data sites that were assigned to same grid node
      
      if( sum( duplicated( i.d ) ) ){
        y.d 				<- f.aggregate.by.one( y.d, i.d, na.rm = TRUE)
        #       X.d 				<- f.aggregate.by.one( X.d, i.d, na.rm = TRUE)
        coords.d 	<- f.aggregate.by.one( coords.d, i.d, na.rm = TRUE)
        i.d	<- i.d[!duplicated(i.d)]
        if( !is.null( covariates.d ) ) covariates.d <- covariates.s[!duplicated(i.d), , drop = FALSE]
      }
      
      ## average data for simulation sites that were assigned to same grid node
      
      if( sum( duplicated( i.s ) ) ){
        y.s 				<- f.aggregate.by.one( y.s, i.s, na.rm = TRUE)
        #       X.s 				<- f.aggregate.by.one( X.s, i.s, na.rm = TRUE)
        coords.s 	<- f.aggregate.by.one( coords.s, i.s, na.rm = TRUE)
        i.s	<- i.s[!duplicated(i.s)]
        if( !is.null( covariates.s ) ) covariates.s <- covariates.s[!duplicated(i.s), , drop = FALSE]

      }
      
      n.d <- length( i.d )
      n.s <- length( i.s )
      
      c.d <- sim.grid[i.d, , drop = FALSE]
      c.s <- sim.grid[i.s, , drop = FALSE]
      
       
      ## compute simple kriging weights and covariances
      
      if( verbose > 2 ) cat( "      compute simple kriging weights ...\n" )
      
      skw <- simple.kriging.weights( 
        pred.coords = c.s,
        support.coords = c.d,
        locations = locations,
        variogram.object = object[["variogram.object"]],
        type = type,
        covariances = control[["covariances"]],
        control = control.predict.georob(ncores = control[["ncores"]])
      )
      
      if( control[["covariances"]] ){
        lambdaT <- skw[["skw"]]      
      } else {
        lambdaT <- skw
      }
            
      ## compute conditional realizations at data and simulation sites
      
      result <- y.s + sim.values[i.s, , drop = FALSE] - pmm(
        lambdaT,
        sim.values[i.d, , drop = FALSE],
        control = control[["pcmp"]]
      )
      
      tmp <- c.s
      if( control[["means"]] ) tmp <- cbind( tmp, expct = y.s )
      if( control[["trend.covariates"]] ) tmp <- cbind( tmp, covariates.s )
      
      result <- cbind( tmp, result )
      
      if( control[["include.data.sites"]] ){
        tmp <- matrix( y.d, nrow = n.d, ncol = nsim )
        colnames( tmp ) <- paste0( "sim.", seq_len( nsim) )
        
        tmp1 <- c.d
        if( control[["means"]] ) tmp1 <- cbind( tmp1, expct = y.d )
        if( control[["trend.covariates"]] ) tmp1 <- cbind( tmp1, covariates.d )
        tmp <- cbind( tmp1, tmp)
        colnames( tmp ) <- colnames( result )
        
        result <- rbind( tmp, result )
      }
      
    }
    
    
  }
  
  ## convert result to same class as newdata
  
  result <- switch(
    class( newdata ),
    "data.frame" = result,
    "SpatialPoints" = ,
    "SpatialPointsDataFrame" = {
      coordinates( result ) <- locations
      result    
    },
    "SpatialPixels" = ,
    "SpatialPixelsDataFrame" = {
      coordinates( result ) <- locations
      gridded( result ) <- TRUE
      result    
    },
    "SpatialGrid" = ,
    "SpatialGridDataFrame" = {
      coordinates( result ) <- locations
      gridded( result ) <- TRUE
      fullgrid( result ) <- TRUE
      result    
    } 
  )
  
  attr( result, "n" ) <- c(n.d = n.d, n.s = n.s )
  
  if( control[["use.grid"]] & control[["covariances"]] & control[["condsim"]] ){
    
    attr( result, "gcvmat.d.d" ) <- skw[["gcvmat"]]
    attr( result, "gcvmat.s.d" ) <- skw[["gcvmat.pred"]]
    
  }
  
  result
  
}
