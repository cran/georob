##  ############################################################################

sample.variogram <- function( object, ... ) UseMethod( "sample.variogram" )


## ##############################################################################

sample.variogram.georob <- function(
  object,
  lag.dist.def,
  xy.angle.def = c( 0., 180. ),
  xz.angle.def = c( 0., 180. ),
  max.lag = Inf,
  estimator = c( "qn", "mad", "matheron", "ch" ),
  mean.angle = TRUE, ...
)
{
  
  ## purpose:          georob method for generic sample.variogram
  ##                   
  ## author:           A. Papritz
  ## date:             2015-11-27
  
  ## checking mandatory arguments
  
  if( missing( lag.dist.def ) ) stop(
    "some mandatory arguments are missing"
  )
  
  ## preparing response vector and matrix of coordinates
  
  x <- residuals( object, level = 0L )
  locations <- object[["locations.objects"]][["coordinates"]]
  
  ##computing sample.variogram
  
  sample.variogram.default(
    x, locations, 
    lag.dist.def = lag.dist.def,
    xy.angle.def = xy.angle.def, xz.angle.def = xz.angle.def,
    max.lag = max.lag, estimator = estimator, mean.angle = mean.angle
  )
  
}


## ##############################################################################

sample.variogram.formula <- function(
  object,
  data, subset, na.action, 
  locations,
  lag.dist.def,
  xy.angle.def = c( 0., 180. ),
  xz.angle.def = c( 0., 180. ),
  max.lag = Inf,
  estimator = c( "qn", "mad", "matheron", "ch" ),
  mean.angle = TRUE, ...
)
{
  
  ## purpose:          formula method for generic sample.variogram
  ##                   
  ## author:           A. Papritz
  ## date:             2015-11-27
  
  ## checking mandatory arguments
  
  if( missing( locations ) || missing( lag.dist.def ) ) stop(
    "some mandatory arguments are missing"
  )
  
  # get model frame, response vector, matrix of coordinates
  
  ### build combined formula for response and locations
  
  extended.formula <- update( 
    object,
    paste( as.character( object )[2L], as.character( locations )[2L], sep = " ~ " )
  )
  
  ## setting-up model frame
  
  cl <- match.call()
  mf <- match.call( expand.dots = FALSE )
  m <- match( 
    c( "data", "subset", "na.action" ),
    names(mf), 0L 
  )
  mf <- mf[c(1L, m)]
  mf[["formula"]] <- extended.formula
  mf[["drop.unused.levels"]] <- TRUE
  mf[[1L]] <- as.name( "model.frame" )
  
  mf <- eval( mf, parent.frame() )
  
  attr( attr( mf, "terms" ), "intercept" ) <- 0
  
  ## preparing response vector and matrix of coordinates
  
  x <- model.response( mf )
  locations <- model.matrix( terms( mf ), mf )
  
  ## computing sample.variogram
  
  sample.variogram.default(
    x, locations, 
    lag.dist.def = lag.dist.def,
    xy.angle.def = xy.angle.def, xz.angle.def = xz.angle.def,
    max.lag = max.lag, estimator = estimator, mean.angle = mean.angle
  )
  
}


## ##############################################################################

sample.variogram.default <- 
  function(
    object,
    locations,
    lag.dist.def,
    xy.angle.def = c( 0., 180. ),
    xz.angle.def = c( 0., 180. ),
    max.lag = Inf,
    estimator = c( "qn", "mad", "matheron", "ch" ),
    mean.angle = TRUE, ...
  )
{
  
  # purpose:          function computes the sample variogram of response
  #                   by various (non-)robust estimators
  #                   
  # author:           A. Papritz
  # date:             2012-04-13
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-04-07 AP correcting error when computing directional variograms
  ## 2015-11-27 default method for generic sample.variogram, checking mandatory arguments
  
  # checking mandatory arguments
  
  if( missing( object ) || missing( locations ) || missing( lag.dist.def ) ) stop(
	"some mandatory arguments are missing"
  )
    
  d2r <- pi / 180.

  estimator <- match.arg( estimator )
  
  # pad missing coordinates
  
  if( ( ndim <- NCOL( locations ) ) < 3L ){
    locations <- cbind( 
      locations, 
      matrix( 0., nrow = NROW( locations ), ncol = 3L - NCOL( locations ) )
    )
  }
  colnames( locations ) <- c( "x", "y", "z" )
  
  # compute lag vectors for all pairs of coordinates
  
  indices.pairs <- combn( NROW( locations ), 2L )
  lag.vectors <- locations[ indices.pairs[2L,], ] - locations[ indices.pairs[1L,], ]
  
  # reflect lag vectors onto half circle
  
  neg.x <- lag.vectors[, 1L] < 0.
  
  lag.vectors[neg.x, 1L] <- -lag.vectors[neg.x, 1L]
  lag.vectors[neg.x, 2L] <- -lag.vectors[neg.x, 2L]
  lag.vectors[neg.x, 3L] <- -lag.vectors[neg.x, 3L]
  
  # compute pairwise differences of responses
  
  t.diff <- object[indices.pairs[2L,]] - object[indices.pairs[1L,]]
  
  # compute Euclidean distances
  
  t.dist <- sqrt( rowSums( lag.vectors^2 ) )
  
  # compute angles in xy- and xz-planes
  
  xy.angle <- -( atan2( lag.vectors[,2L], lag.vectors[,1L] ) - pi/2. ) / d2r
  xz.angle <- -( atan2( lag.vectors[,3L], lag.vectors[,1L] ) - pi/2. ) / d2r

  # define lag class upper limits
  
  if( length( lag.dist.def ) == 1L ) {
    t.lag.limits <- seq( 0., max( c( t.dist ) ) + lag.dist.def, by = lag.dist.def )
  } else {
    t.lag.limits <- lag.dist.def
  }
  
  # group the lag vectors into classes
  
  t.lag.class <- cut( t.dist, breaks = t.lag.limits, include.lowest = TRUE )
  xy.angle.class <- cut( xy.angle, breaks = xy.angle.def, include.lowest = TRUE )
  xz.angle.class <- cut( xz.angle, breaks = xz.angle.def, include.lowest = TRUE )
  
  # join first and last angle classes if end points match and there is
  # more than one angle
  
  # xy-plane
  
  n <- length( xy.angle.def )
  d <- diff( xy.angle.def )
  xy.angle.mid.class <- 0.5 * ( xy.angle.def[-1L] + xy.angle.def[-n] )
  if( 
    n > 2L &&
    identical( xy.angle.def[1L], 0. ) && 
    identical( xy.angle.def[n], 180. ) &&
    !all( d[1L] == d[-1L] )
  ){
    
    nb <- 2L
    right <- TRUE
    include.lowest <- TRUE
    dig.lab <- 3L
    breaks <- c( -diff( xy.angle.def[(n-1L):n] ), xy.angle.def[2L] )
    
    for (dig in dig.lab:max(12L, dig.lab)) {
      ch.br <- formatC(breaks, digits = dig, width = 1L)
      if (ok <- all(ch.br[-1L] != ch.br[-nb])) 
      break
    }
    labels <- if(ok){
      paste( 
        if(right){
          "("  
        } else {
          "["
        }, ch.br[-nb], ",", ch.br[-1L], 
        if (right){ 
          "]"
        } else {
          ")"
        }, sep = ""
      )
    } else {
      paste( "Range", seq_len(nb - 1L), sep = "_")
      if( ok && include.lowest ){
        if (right)  substr(labels[1L], 1L, 1L) <- "["
        else substring( labels[nb - 1L], nchar(labels[nb - 1L], "c")) <- "]"
      }
      
    }
    sel <- as.integer( xy.angle.class ) == nlevels( xy.angle.class )
    xy.angle[sel] <- xy.angle[sel] - 180.
    levels( xy.angle.class )[c( 1L, nlevels( xy.angle.class )) ] <- labels
    xy.angle.mid.class[1L] <- xy.angle.mid.class[1L] - (180. - xy.angle.mid.class[n-1L])
    xy.angle.mid.class <- xy.angle.mid.class[-(n-1L)]
  }
  
  # xz-plane
  
  n <- length(xz.angle.def)
  d <- diff( xz.angle.def )
  xz.angle.mid.class <- 0.5 * ( xz.angle.def[-1L] + xz.angle.def[-n] )
  if( 
    n > 2L &&
    identical( xz.angle.def[1L], 0. ) && 
    identical( xz.angle.def[n], 180. ) &&
    !all( d[1L] == d[-1L] )
  ){
    
    nb <- 2L
    right <- TRUE
    include.lowest <- TRUE
    dig.lab <- 3L
    breaks <- c( -diff( xz.angle.def[(n-1L):n] ), xz.angle.def[2L] )
    for (dig in dig.lab:max(12L, dig.lab)) {
      ch.br <- formatC(breaks, digits = dig, width = 1)
      if (ok <- all(ch.br[-1L] != ch.br[-nb])) 
      break
    }
    labels <- if(ok){
      paste( 
        if(right){
          "("  
        } else {
          "["
        }, ch.br[-nb], ",", ch.br[-1L], 
        if (right){ 
          "]"
        } else {
          ")"
        }, sep = ""
      )
    } else {
      paste( "Range", seq_len(nb - 1L), sep = "_")
      if( ok && include.lowest ){
        if (right)  substr(labels[1L], 1L, 1L) <- "["
        else substring( labels[nb - 1L], nchar(labels[nb - 1L], "c")) <- "]"
      }
      
    }
    sel <- as.integer( xz.angle.class ) == nlevels( xz.angle.class )
    xz.angle[sel] <- xz.angle[sel] - 180.
    levels( xz.angle.class )[c( 1L, nlevels( xz.angle.class )) ] <- labels
    xz.angle.mid.class[1L] <- xz.angle.mid.class[1L] - (180. - xz.angle.mid.class[n-1L])
    xz.angle.mid.class <- xz.angle.mid.class[-(n-1L)]
  }
  
  t.classes <- list( t.lag.class, xy.angle.class, xz.angle.class )
  
  # compute mean distance, mean lag vector, angle classes and and number of
  # pairs per class
  
  lag.mean <- as.vector( tapply(
      t.dist,
      t.classes,
      mean
    ))
  xy.angle.mean <- as.vector( tapply(
      xy.angle,
      t.classes,
      mean
    ))
  xz.angle.mean <- as.vector( tapply(
      xz.angle,
      t.classes,
      mean
    ))
  xy.angle.centre <- as.vector( tapply(
      xy.angle.class,
      t.classes,
      function( x ) unique( x )
    ))
  xz.angle.centre <- as.vector( tapply(
      xz.angle.class,
      t.classes,
      function( x ) unique( x )
    ))
  
  xy.angle.centre[!is.na(xy.angle.centre)] <- 
    levels( xy.angle.class )[xy.angle.centre[!is.na(xy.angle.centre)]]
  xz.angle.centre[!is.na(xz.angle.centre)] <- 
    levels( xz.angle.class )[xz.angle.centre[!is.na(xz.angle.centre)]]
  
  t.lag.npairs <- as.vector( table( t.classes ) )
  t.lag.select <- lag.mean <= max.lag
  
  # compute semivariance per class
  
  if( estimator %in% c( "mad", "qn" ) ) {
    
    if( estimator == "mad" ) { 
      
      t.gamma <- as.vector( tapply( 
          t.diff, 
          t.classes, 
          mad, center = 0. 
        ))
      
    } else {
      
      t.gamma <- as.vector( tapply( 
          t.diff, 
          t.classes, 
          Qn, finite.corr = TRUE
        ))
    } 
    
    t.gamma <- 0.5 * c( t.gamma )^2
    
  } else if ( estimator == "ch" ) {
    
    t.gamma <- as.vector( tapply( 
        t.diff, 
        t.classes, 
        function( x ) 0.5 * mean( sqrt(abs(x)) )^4 / ( 0.457+0.494/length(x) )
      ))
    
  } else if( estimator == "matheron" ) {
    
    t.gamma <- as.vector( tapply( 
        t.diff, 
        t.classes, 
        function( x ) 0.5 * mean( x^2 )
      ))
    
  }
  
  # collect results
  
  t.sel <- t.lag.npairs > 0L & lag.mean <= max.lag
  r.result <- data.frame(
    lag.dist = lag.mean[ t.sel],
    xy.angle = factor( xy.angle.centre[ t.sel ], levels = levels( xy.angle.class ) ),
    xz.angle = factor( xz.angle.centre[ t.sel ], levels = levels( xz.angle.class ) ),
    gamma = t.gamma[ t.sel ],
    npairs = t.lag.npairs[ t.sel ]
  )
    
  # compute lag vectors
  
  # center of circle sectors
  
  if( mean.angle ){
    
    # mean angle of lag pairs
    
    t.aux <- r.result[["lag.dist"]] * sin( xz.angle.mean[t.sel] * d2r )
    r.result[["lag.x"]] <- t.aux * sin( xy.angle.mean[t.sel] * d2r )
    r.result[["lag.y"]] <- t.aux * cos( xy.angle.mean[t.sel] * d2r )
    r.result[["lag.z"]] <- r.result[["lag.dist"]] * cos( xz.angle.mean[t.sel] * d2r )
    
  } else {
    
    t.aux <- r.result[["lag.dist"]] * sin( xz.angle.mid.class[ as.integer( r.result[["xz.angle"]] ) ] * d2r )
    r.result[["lag.x"]] <- t.aux * sin( xy.angle.mid.class[ as.integer( r.result[["xy.angle"]] ) ] * d2r )
    r.result[["lag.y"]] <- t.aux * cos( xy.angle.mid.class[ as.integer( r.result[["xy.angle"]] ) ] * d2r )
    r.result[["lag.z"]] <- r.result[["lag.dist"]] * 
      cos( xz.angle.mid.class[ as.integer( r.result[["xz.angle"]] ) ] * d2r )
    
  }
  
  class( r.result) <-  c( "sample.variogram", "data.frame" )
  
  attr( r.result, "ndim")                <- ndim
  attr( r.result, "lag.dist.def")       <- lag.dist.def
  attr( r.result, "xy.angle.mid.class")  <- xy.angle.mid.class
  attr( r.result, "xz.angle.mid.class")  <- xz.angle.mid.class
  attr( r.result, "estimator" )          <- estimator
    
  return( r.result )  
  
}

## ##############################################################################

summary.sample.variogram <- 
  function( object, ... )
{
  
  ## Summary method for class summary.sample.variogram
  
  ## 2012-11-12 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  ans <- object[c( "lag.dist", "npairs", "xy.angle", "xz.angle" ) ]
  attr( ans, "estimator" ) <- attr( object, "estimator" )
  class( ans ) <- "summary.sample.variogram"
  return( ans )
  
}


## ##############################################################################

print.summary.sample.variogram <- 
  function( 
    x, digits = max(3, getOption("digits") - 3), ...
  )
{
  
  ## Print method for class summary.sample.variogram
  
  ## 2012-11-12 A. Papritz
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  cat( "\nSample variogram estimator:", attr( x, "estimator" ), "\n" )
  
  cat( "\nSummary of lag distances\n" )
  print( summary( x[["lag.dist"]] ) )
  
  cat( "\nSummary of number of pairs per lag and distance classes\n" )
  print( summary( x[["npairs"]] ) )
  
  cat( "\nAngle classes in xy-plane:", levels( x[["xy.angle"]] ), "\n" )
  cat( "Angle classes in xz-plane:", levels( x[["xz.angle"]] ), "\n" )
  
  invisible( x )
  
}

## ##############################################################################

plot.sample.variogram <- 
  function(
    x,
    type = "p", add = FALSE, 
    xlim = c( 0., max( x[["lag.dist"]] ) ),
    ylim = c( 0., 1.1 * max( x[["gamma"]] ) ),
    col,
    pch,
    lty,
    cex = 0.8,
    xlab = "lag distance", ylab = "semivariance",
    annotate.npairs = FALSE,
    npairs.pos = 3, npairs.cex = 0.7,
    legend = nlevels( x[["xy.angle"]] ) > 1 || nlevels( x[["xz.angle"]] ) > 1,
    legend.pos = "topleft",
    ...
  )
{
  
  ## Plot method for class sample.variogram
  
  ## 2012-12-12 A. Papritz
  ## 2012-12-21 AP correction for using col and pch
  ## 2013-05-12 AP correction for using ...
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2016-08-24 AP new argument lty
  
  if( !add ) plot( 
    gamma ~ lag.dist, x, type = "n",
    xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...
  )
  
  if( missing( col ) ){
    col <- 1L:nlevels( x[["xy.angle"]] )
  } else if( length( col ) < nlevels( x[["xy.angle"]] ) ) stop(
    "number of colors less than number of directions in x-y-plane for which semivariances were computed"
  )
  if( missing( pch ) ){
    pch <- 1L:nlevels( x[["xz.angle"]] )
  } else if( length( pch ) < nlevels( x[["xz.angle"]] ) ) stop(
    "number of colors less than number of directions in x-z-plane for which semivariances were computed"
  )
  if( missing( lty ) ) lty <- "solid"
  
  tapply( 
    1L:NROW( x ),
    list( x[["xy.angle"]], x[["xz.angle"]] ),
    function( i, x, type, col, pch, cex ){
      points( 
        gamma ~ lag.dist, x, subset = i, 
        type = type,
        col = col[as.numeric(x[i, "xy.angle"])],
        pch = pch[as.numeric(x[i, "xz.angle"])], 
        lty = lty, cex = cex
      )
      
    },
    x = x, type = type, col = col, pch = pch, cex = cex
  )
  
  if( annotate.npairs ){
    with(
      x, 
      text( lag.dist, gamma, npairs, pos = npairs.pos, cex = npairs.cex, col = col ) 
    )
  }
  
  if( legend ){
    legend(
      x = legend.pos, bty = "n", 
      col = c( 
        1L:nlevels( x[["xy.angle"]] ), 
        if( nlevels( x[["xz.angle"]] ) > 1L ) rep( 1L, nlevels( x[["xz.angle"]] ) ) 
      ),
      pch = c( 
        rep( 1L, nlevels( x[["xy.angle"]] ) ), 
        if( nlevels( x[["xz.angle"]] ) > 1L ) 1L:nlevels( x[["xz.angle"]] ) 
      ),
      legend = c( 
        paste( "xy.angle:", levels( x[["xy.angle"]] ) ), 
        if( nlevels( x[["xz.angle"]] ) > 1L ) paste( "xz.angle:", levels( x[["xz.angle"]] ) )
      ),
      pt.cex = cex
    )
    
  }
  
  invisible( x )
  
}

## ##############################################################################

control.fit.variogram.model <- function(
  maximizer = c( "nlminb", "optim" ),
  param.tf = param.transf(),
  fwd.tf = fwd.transf(), 
  bwd.tf = bwd.transf(),
  hessian = TRUE,
  optim = control.optim(),
  nlminb = control.nlminb()
){

  ## auxiliary function to set meaningful default values for georob

  ## 2019-04-07 A. Papritz

  list(
    maximizer = match.arg( maximizer ),
    param.tf = param.tf, fwd.tf = fwd.tf, bwd.tf = bwd.tf,
    hessian = hessian,
    optim = optim, 
    nlminb = nlminb
  )

}


## ##############################################################################

fit.variogram.model <-
  function(
    sv,
    variogram.model = c( "RMexp", "RMaskey", "RMbessel", "RMcauchy", 
      "RMcircular", "RMcubic", "RMdagum", "RMdampedcos", "RMdewijsian", "RMfbm",
      "RMgauss", "RMgencauchy", "RMgenfbm", "RMgengneiting", "RMgneiting", "RMlgd",
      "RMmatern", "RMpenta", "RMqexp", "RMspheric", "RMstable",
      "RMwave", "RMwhittle"
    ), 
    param, fit.param = default.fit.param()[names(param)],
	  aniso = default.aniso(), fit.aniso = default.fit.aniso(),
    variogram.object = NULL,
    max.lag = max( sv[["lag.dist"]] ),
    min.npairs = 30,
    weighting.method = c( "cressie", "equal", "npairs" ),
    control = control.fit.variogram.model(),
    verbose = 0
  )
{
  
  ## Function to fit a variogram model to a sample variogram generated by
  ## sample.variogram
  
  ## 2012-12-10 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  ## 2015-11-27 AP checking mandatory arguments, issuing warnings
  ## 2016-02-08 AP correcting error in setting default values for fit.param
  ## 2016-08-18 AP changes for nested variogram models
  ## 2016-11-14 AP correcting error in 3d rotation matrix for geometrically anisotropic variograms
  ## 2019-04-07 AP new control argument and possibility to choose between 
  ##               optim and nlminb for non-linear least squares 
  
  ## auxiliary function called by optim to compute objective function
  
  f.aux <- function( 
      adjustable.param.aniso, 
      envir, 
      fixed.param.aniso, name.param.aniso, tf.param.aniso,
      bwd.tf, lag.vectors, gamma, npairs, 
      weighting.method, d2r, verbose
    )
  {
    
    ## load variogram.item object
    
    variogram.item <- get( "variogram.item", pos = as.environment( envir ) )
    
    #   print(str(variogram.item[c("variogram.object", "eta", "xi")]))
    
    ##  transform variogram parameters back to original scale
    
    param.aniso <- c( adjustable.param.aniso, fixed.param.aniso )[name.param.aniso]
    
    param.aniso <- sapply(
      name.param.aniso,
      function( x, bwd.tf, param.tf, param ) bwd.tf[[param.tf[x]]]( param[x] ),
      bwd.tf = bwd.tf,
      param.tf = tf.param.aniso,
      param = param.aniso
    )
    names( param.aniso ) <- name.param.aniso
    
    ## correct parameters names and determine model component
    
    tmp <- strsplit( names(param.aniso), control.georob()[["sepstr"]], fixed = TRUE )
    
    name.param.aniso <- names( param.aniso ) <- names( tf.param.aniso ) <- sapply(
      tmp, function(x) x[1L]
    )
    
    cmp <- sapply( tmp, function(x) as.integer(x[2L]) )
    
    ## build temporary variogram.object with new parameter values
    
    variogram.object <- tapply(
      param.aniso,
      factor( cmp ),
      function( x ) list(
        param = x[-((length(x)-4L):length(x))],
        aniso = x[((length(x)-4L):length(x))]
      ), simplify = FALSE
    )
    
    variogram.object <- lapply(
      1L:length(variogram.object),
      function( i, x, variogram.item ){
        c(
          variogram.item[["variogram.object"]][[i]][c("variogram.model", "isotropic")],
          x[[i]]
        )
      }, x = variogram.object, variogram.item = variogram.item
    )
    
    ## check whether variogram parameters are within reasonable bounds and
    ## return an error otherwise
    
    if( length( param.aniso ) && any( param.aniso > control.georob()[["safe.param"]] ) ){
      
      lapply(
        1L:length(variogram.object),
        function( i, x, d2r ){
          
          x <- x[[i]]
          
          t.param <- x[["param"]]
          
          if( !x[["isotropic"]] ) t.param <- c(
            t.param, x[["aniso"]] / c( rep( 1., 2L), rep( d2r, 3L ) )
          )
          cat( "\n\n                      ",
            format( names( t.param ), width = 14L, justify = "right" ),
            "\n", sep = ""
          )
          cat( "  Variogram parameters",
            format(
              signif( t.param, digits = 7L ),
              scientific = TRUE, width = 14L
            ), "\n" , sep = ""
          )
          
        }, x = variogram.object, d2r = d2r
      )
      
      return( NA_real_ )
      
    }
    
    ## check whether extra variogram parameters are within allowed bounds and
    ## return an error otherwise
    
    lapply(
      1L:length(variogram.object),
      function( i, x, lag.vectors ){
        x <- x[[i]]
        ep <- param.names( model = x[["variogram.model"]] )
        param.bounds <- param.bounds( x[["variogram.model"]], attr( lag.vectors, "ndim.coords" ) )
        ep.param <- x[["param"]][ep]
        if( !is.null( param.bounds ) ) t.bla <- sapply(
          1L:length( ep.param ),
          function( i, param, bounds ){
            if( param[i] < bounds[[i]][1L] || param[i] > bounds[[i]][2L] ) cat(
              "value of parameter '", names( param[i] ), "' outside of allowed range\n\n", sep = ""
            )
            return( NA_real_ )
          },
          param = ep.param,
          bounds = param.bounds
        )
      }, x = variogram.object, lag.vectors = lag.vectors
    )
    
    ##  update variogram and parameters
    
    variogram.item[["variogram.object"]] <- lapply(
      1L:length(variogram.object),
      function( i, x, vo, n ){
        vo <- vo[[i]]
        vo[["param"]] <- x[[i]][["param"]]
        vo[["aniso"]] <- aniso <- x[[i]][["aniso"]]
        vo[["sincos"]] <- list(
          co = unname( cos( aniso["omega"] ) ),
          so = unname( sin( aniso["omega"] ) ),
          cp = unname( cos( aniso["phi"] ) ),
          sp = unname( sin( aniso["phi"] ) ),
          cz = unname( cos( aniso["zeta"] ) ),
          sz = unname( sin( aniso["zeta"] ) )
        )
        
        if( n <= 3L ){
          
          vo[["rotmat"]] <- with(
            vo[["sincos"]],
            rbind(
              c(             so*sp,             co*sp,       cp ),
              c( -co*cz - so*cp*sz,  so*cz - co*cp*sz,    sp*sz ),
              c(  co*sz - so*cp*cz, -so*sz - co*cp*cz,    sp*cz )
            )[ 1L:n, 1L:n, drop = FALSE ]
          )
          
          vo[["sclmat"]] <- 1. / c( 1., aniso[ c("f1", "f2") ] )[ 1L:n ]
          
        } else {  # only isotropic case for n > 3
          
          vo[["rotmat"]] <- diag( n )
          vo[["sclmat"]] <- rep( 1., n )
          
        }
        
        vo
      }, 
      x = variogram.object, vo = variogram.item[["variogram.object"]], 
      n = attr( lag.vectors, "ndim.coords" )
    )
    
    ##  print updated variogram parameters
    
    if( verbose > 1. ) {
      
      lapply(
        1L:length(variogram.object),
        function( i, x, d2r, reparam ){
          
          x <- x[[i]]
          
          t.param <- x[["param"]]
                    
          if( !x[["isotropic"]] ) t.param <- c(
            t.param, x[["aniso"]] / c( rep( 1., 2L), rep( d2r, 3L ) )
          )
          tmp <- names( t.param )
          if( identical( i, length(variogram.object) ) ) tmp <- c( tmp, "SSE" )
          cat( "\n\n                      ",
            format( 
              tmp, 
              width = 14L, justify = "right" ),
            "\n", sep = ""
          )
          cat( "  Variogram parameters",
            format(
              signif( t.param, digits = 7L ),
              scientific = TRUE, width = 14L
            ), 
            if( identical( i, length(variogram.object) ) ) "" else "\n" , sep = ""
          )
          
        }, x = variogram.object, d2r = d2r, reparam = FALSE
      )
      
    }
    
    ## compute the generalized covariance of signal
    
    Valpha <- f.aux.gcr( 
      lag.vectors = lag.vectors,
      variogram.object = variogram.item[["variogram.object"]], symmetric = FALSE,
      control.pcmp = control.pcmp(max.ncores=1L), verbose = verbose
    )
    
    if( any( sapply( Valpha, function(x) x[["error"]] ) ) ){
      warning( "there were errors: call function with argument 'verbose' > 1" )
      if( verbose > 1. ) cat( "\nan error occurred when computing semivariances\n" )
      return( NA_real_ )
    }
    
    ## compute semivariance of signal and sum up
    
    t.model <- rowSums( sapply(
        1L:length(Valpha), 
        function( i, x, variogram.object ){
          x <- Valpha[[i]]
          -( x[["Valpha"]] - x[["gcr.constant"]] ) * variogram.object[[i]][["param"]]["variance"]
        }, x = Valpha, variogram.object = variogram.item[["variogram.object"]]
      ))
    
    ## add nugget
    
    sel <- sqrt( rowSums( lag.vectors^2 ) ) > 0.
    t.model[sel] <- t.model[sel] + sum( 
      variogram.item[["variogram.object"]][[1L]][["param"]][c("nugget", "snugget")] 
    )
    
        
    ## compute weights
    
    t.weights <- rep( 1., NROW( lag.vectors ) )
    if( weighting.method != "equal" ) t.weights <- npairs
    if( weighting.method == "cressie" ) t.weights <- t.weights / t.model^2
    
    ##  ... and compute the objective function that is minimized
    
    t.res <- gamma - t.model
    sse <- sum( t.weights * t.res^2 )
    
    if( verbose > 0. ) cat(
      format(
        signif( sse, digits = 7L ),
        scientific = TRUE, width = 14L
      ), "\n", sep = ""
    )
    
    ## store copies of model semivariance, residuals and weights
    
    variogram.item[["fitted"]] <- t.model
    variogram.item[["residuals"]] <- t.res
    variogram.item[["weights"]] <- t.weights
    
    ## store variogram.item
    
    assign( "variogram.item", variogram.item, pos = as.environment( envir ) )
    
    return( sse )
    
  }
  
  ## begin of main body of function
  
  ## check whether all mandatory arguments have been provided
  
  if( missing( sv ) ) stop( 
    "some mandatory arguments are missing" 
  )
  
  d2r <- pi / 180.
      
  ## match arguments
  
  cl <- match.call()

  weighting.method = match.arg( weighting.method )
  
  ## setup or check contents of variogram.object
  
  if( is.null( variogram.object ) ){
    
    ## match variogram model
    
    variogram.model <- match.arg( variogram.model )  
        
    ## match names of param, aniso, fit.param, fit.aniso
    ## !! CARE: code works only if fit.param has not been used before (lazy evaluation)
  
    if( !missing( param ) ){
      tmp <- names( param )
      tmp <- sapply(tmp, function(x, choices){
          match.arg(x, choices)
        },
        choices = names( default.fit.param() )
      )
      names( param ) <- tmp
    }
    
    if( !missing( fit.param ) ){
      tmp <- names( fit.param )
      tmp <- sapply(tmp, function(x, choices){
          match.arg(x, choices)
        },
        choices = names( default.fit.param() )
      )
      names( fit.param ) <- tmp
      fit.param <- fit.param[names( fit.param ) %in% names( param )]
    }
    
    if( !missing( aniso ) ){
      tmp <- names( aniso )
      tmp <- sapply(tmp, function(x, choices){
          match.arg(x, choices)
        },
        choices = names( default.aniso() )
      )
      names( aniso ) <- tmp
    }
    
    if( !missing( fit.aniso ) ){
      tmp <- names( fit.aniso )
      tmp <- sapply(tmp, function(x, choices){
          match.arg(x, choices)
        },
        choices = names( default.aniso() )
      )
      names( fit.aniso ) <- tmp
    }
    
    ## create variogram.object
    
    variogram.object <- list(
      list( 
        variogram.model = variogram.model,
        param = param, fit.param = fit.param,
        aniso = aniso, fit.aniso = fit.aniso
      )
    )
  
  } else {
    
    if( !missing( param ) ) cat( 
      "\n information on initial parameters in 'param' ignored because 'variogram.object' specified\n\n"
    )
    if( !missing( fit.param ) ) cat( 
      "\n information on initial parameters in 'fit.param' ignored because 'variogram.object' specified\n\n"
    )
    if( !missing( aniso ) ) cat( 
      "\n information on initial parameters in 'aniso' ignored because 'variogram.object' specified\n\n"
    )
    if( !missing( fit.aniso ) ) cat( 
      "\n information on initial parameters in 'fit.aniso' ignored because 'variogram.object' specified\n\n"
    )
    
    variogram.object <- lapply(
      variogram.object,
      function( y, vm ){
            
        ## match names of components
        
        tmp <- names( y )
        tmp <- sapply(tmp, function(x, choices){
            match.arg(x, choices)
          },
          choices = c( "variogram.model", "param", "fit.param", "aniso", "fit.aniso" )
        )
        names( y ) <- tmp
        
        ## check whether component variogram.model and param are present
        
        if( is.null( y[["variogram.model"]] ) ) stop(
          "'variogram.model' missing in some component of variogram.object"       
        )
        if( is.null( y[["param"]] ) ) stop(
          "'param' missing in some component of variogram.object"       
        )
        
        ## match variogram.model
        
        y[["variogram.model"]] <- match.arg( y[["variogram.model"]], vm )
        
                
        ## match names of param, aniso, fit.param, fit.aniso and set
        ## default values if missing
        
        tmp <- names( y[["param"]] )
        tmp <- sapply(tmp, function(x, choices){
            match.arg(x, choices)
          },
          choices = names( default.fit.param() )
        )
        names( y[["param"]] ) <- tmp
        
        if( is.null( y[["fit.param"]] ) ){
          y[["fit.param"]] <- default.fit.param()
        } else {
          tmp <- names( y[["fit.param"]] )
          tmp <- sapply(tmp, function(x, choices){
              match.arg(x, choices)
            },
            choices = names( default.fit.param() )
          )
          names( y[["fit.param"]] ) <- tmp
        }
        y[["fit.param"]] <- y[["fit.param"]][names(y[["fit.param"]] ) %in% names( y[["param"]] )]
        
        if( is.null( y[["aniso"]] ) ){
          y[["aniso"]] <- default.aniso()
        } else {
          tmp <- names( y[["aniso"]] )
          tmp <- sapply(tmp, function(x, choices){
              match.arg(x, choices)
            },
            choices = names( default.aniso() )
          )
          names( y[["aniso"]] ) <- tmp
        }
        
        if( is.null( y[["fit.aniso"]] ) ){
          y[["fit.aniso"]] <- default.fit.aniso()
        } else {
          tmp <- names( y[["fit.aniso"]] )
          tmp <- sapply(tmp, function(x, choices){
              match.arg(x, choices)
            },
            choices = names( default.aniso() )
          )
          names( y[["fit.aniso"]] ) <- tmp
        }
        y
      }, vm = variogram.model
    )
  }
  
  ## process contents of variogram.object

   variogram.object <- lapply(
    variogram.object,
    function( x, TT, d2r, n ){

      ## create local copies of objects

      variogram.model <- x[["variogram.model"]]
      param <- x[["param"]]
      fit.param <- x[["fit.param"]]
      aniso <- x[["aniso"]]
      fit.aniso <- x[["fit.aniso"]]

      ##  check whether fitting of chosen variogram model is implemented and
      ##  return names of extra parameters (if any)

      ep <- param.names( model = variogram.model )

      ## check names of initial variogram parameters and flags for fitting

      param.name <- c( "variance", "snugget", "nugget", "scale", ep )

      if( !all( param.name[-(2L:3L)]  %in% names( param ) ) ) stop(
        "no initial values provided for parameter(s) '",
        paste( (param.name[-(2L:3L)])[ !(param.name[-(2L:3L)]) %in% names( param ) ], collapse= ", "), "'"
      )

      if( !all( param.name[-(2L:3L)]  %in% names( fit.param ) ) ) stop(
        "no fitting flagss provided for parameter(s) '",
        paste( (param.name[-(2L:3L)])[ !(param.name[-(2L:3L)]) %in% names( fit.param ) ], collapse= ", "), "'"
      )

      if( length( param ) != length( fit.param ) ||
        !all( names( fit.param ) %in% names( param ) )
      ) stop(
        "names of variogram parameters and control flags for fitting do not match"
      )

      if( !all( is.numeric( param ) ) ) stop(
        "initial values of variogram parameters must be of mode 'numeric'"
      )
      if( !all( is.logical( fit.param ) ) ) stop(
        "fitting control flags of variogram parameters must be of mode 'logical'"
      )

      ##  rearrange initial variogram parameters

      param <- param[param.name[param.name %in% names(param)]]
      fit.param <- fit.param[param.name[param.name %in% names(param)]]

      ## check whether intitial values of variogram parameters are valid

      if( param["variance"] < 0. ) stop("initial value of 'variance' must be >= 0" )
      if( !is.na(param["snugget"]) && param["snugget"] < 0. )  stop("initial value of 'snugget' must be >= 0" )
      if( !is.na(param["nugget"]) && param["nugget"] <= 0. ) stop("initial value of 'nugget' must be > 0" )
      if( param["scale"] < 0. ) stop("initial value of 'scale' must be >= 0" )

      param.bounds <- param.bounds( variogram.model, n )
      ep.param <- param[ep]

      if( !is.null( param.bounds ) ) t.bla <- sapply(
        1L:length( ep.param ),
        function( i, param, bounds ){
          if( param[i] < bounds[[i]][1L] || param[i] > bounds[[i]][2L] ) stop(
            "initial value of parameter '", names( param[i] ), "' outside of allowed range"
          )
        },
        param = ep.param,
        bounds = param.bounds
      )

      ##  rearrange and check flags controlling variogram parameter fitting

      if(
        variogram.model %in% (t.models <- c( "RMfbm" ) ) && all(
          fit.param[c( "variance", "scale" ) ]
        )

      ) stop(
        "'variance', 'scale' cannot be fitted simultaneously for variograms ",
        paste( t.models, collapse = " or "), "; \n  'scale' parameter must be fixed"
      )

      ## check names of initial anisotropy parameters and flags for fitting

      aniso.name <- c( "f1", "f2", "omega", "phi", "zeta" )

      #       if( !all( names( aniso ) %in% aniso.name ) ) stop(
      #         "error in names of initial values of anisotropy parameters"
      #       )

      if( !all( aniso.name  %in% names( aniso ) ) ) stop(
        "no initial values provided for parameter(s) '",
        aniso.name[ !aniso.name %in% names( aniso ) ], "'"
      )

      if( length( aniso ) != length( fit.aniso ) ||
        !all( names( fit.aniso ) %in% names( aniso ) )
      ) stop(
        "names of anisotropy parameters and control flags for fitting do not match"
      )

      if( !all( is.numeric( aniso ) ) ) stop(
        "initial values of anisotropy parameters must be of mode 'numeric'"
      )
      if( !all( is.logical( fit.aniso ) ) ) stop(
        "fitting control flags of anisotropy parameters must be of mode 'logical'"
      )

      ##  rearrange initial anisotropy parameters

      aniso <- aniso[aniso.name]
      fit.aniso <- fit.aniso[aniso.name]

      ## check whether intitial values of anisotropy parameters are valid

      if( aniso["f1"] < 0. ||  aniso["f1"] > 1. ) stop(
        "initial value of parameter 'f1' must be in [0, 1]"
      )
      if( aniso["f2"] < 0. ||  aniso["f1"] > 1. ) stop(
        "initial value of parameter 'f2' must be in [0, 1]"
      )
      if( aniso["omega"] < 0. ||  aniso["omega"] > 180. ) stop(
        "initial value of parameter 'omega' must be in [0, 180]"
      )
      if( aniso["phi"] < 0. ||  aniso["phi"] > 180. ) stop(
        "initial value of parameter 'phi' must be in [0, 180]"
      )
      if( aniso["zeta"] < -90. ||  aniso["zeta"] > 90. ) stop(
        "initial value of parameter 'zeta' must be in [-90, 90]"
      )


      ## check whether variogram is isotropic

      if( identical( aniso, default.aniso() ) ){
        isotropic <- TRUE
      } else {
        isotropic <- FALSE
      }      

      ## adjust default initial values of anisotropy parameters if these are
      ## fitted

      if( fit.aniso["omega"] && identical( aniso["f1"], 1. ) ){
        aniso["f1"] <- aniso["f1"] - sqrt( .Machine$double.eps )
      }

      if( fit.aniso["phi"] ){
        if( identical( aniso["f1"], 1. ) ) aniso["f1"] <- aniso["f1"] - sqrt( .Machine$double.eps )
        if( identical( aniso["f2"], 1. ) ) aniso["f2"] <- aniso["f2"] - sqrt( .Machine$double.eps )
      }
      if( fit.aniso["zeta"] && identical( aniso["f2"], 1. ) ){
        aniso["f2"] <- aniso["f2"] - sqrt( .Machine$double.eps )
      }


      ##  convert angles to radian

      aniso[c("omega", "phi", "zeta" )] <- aniso[c("omega", "phi", "zeta" )] * d2r

      ## complement aniso components with sin/cos terms, rotation and scaling matrices

      sincos <- list(
        co = unname( cos( aniso["omega"] ) ),
        so = unname( sin( aniso["omega"] ) ),
        cp = unname( cos( aniso["phi"] ) ),
        sp = unname( sin( aniso["phi"] ) ),
        cz = unname( cos( aniso["zeta"] ) ),
        sz = unname( sin( aniso["zeta"] ) )
      )

      if( n <= 3L ){

        rotmat <- with(
          sincos,
          rbind(
            c(             so*sp,             co*sp,       cp ),
            c( -co*cz - so*cp*sz,  so*cz - co*cp*sz,    sp*sz ),
            c(  co*sz - so*cp*cz, -so*sz - co*cp*cz,    sp*cz )
          )[ 1L:n, 1L:n, drop = FALSE ]
        )

        sclmat <- 1. / c( 1., aniso[ c("f1", "f2") ] )[ 1L:n ]

      } else {  # only isotropic case for n > 3

        rotmat <- diag( n )
        sclmat <- rep( 1., n )

      }

      ## return all items

      list(
        variogram.model = variogram.model,
        param = param, fit.param = fit.param,
        isotropic = isotropic,
        aniso = aniso, fit.aniso = fit.aniso,
        sincos = sincos, rotmat = rotmat, sclmat = sclmat
      )

    }, d2r = d2r, n = attr( sv, "ndim" )
  )

  #   print(str(variogram.object))

  ## set consistent values for nugget and snugget of nested variogram models

  ## no nugget and snugget in any model component

 is.na.nugget <- sapply( variogram.object, function( x ) is.na( x[["param"]]["nugget"] ) )
  if( all( is.na.nugget ) ) stop(
    "one of the variogram components must contain a 'nugget' effect"
  )

  is.na.snugget <- sapply( variogram.object, function( x ) is.na( x[["param"]]["snugget"] ) )
  if( all( is.na.snugget ) ){
    param <- variogram.object[[1L]][["param"]]
    fit.param <- variogram.object[[1L]][["fit.param"]]
    param <- c( param[1L], snugget = 0., param[-1L] )
    fit.param <- c( fit.param[1L], snugget = FALSE, fit.param[-1L] )
    variogram.object[[1L]][["param"]] <- param
    variogram.object[[1L]][["fit.param"]] <- fit.param
  }

  ## nugget and snugget are combined and shifted to first model component

  tmp <- names(variogram.object[[1L]][["param"]])
  tmp <- tmp[!tmp %in% c("variance", "snugget", "nugget")]

  variogram.object[[1L]][["param"]]["nugget"] <- sum(
    sapply( variogram.object, function(x) x[["param"]]["nugget"] ),
    na.rm = TRUE
  )
  variogram.object[[1L]][["fit.param"]]["nugget"] <- any(
    sapply( variogram.object, function(x) x[["fit.param"]]["nugget"] ),
    na.rm = TRUE
  )
  variogram.object[[1L]][["param"]]["snugget"] <- sum(
    sapply( variogram.object, function(x) x[["param"]]["snugget"] ),
    na.rm = TRUE
  )
  variogram.object[[1L]][["fit.param"]]["snugget"] <- any(
    sapply( variogram.object, function(x) x[["fit.param"]]["snugget"] ),
    na.rm = TRUE
  )

  ## rearrage order of parameters in first component

  variogram.object[[1L]][["param"]] <-  variogram.object[[1L]][["param"]][c("variance", "snugget", "nugget", tmp)]
  variogram.object[[1L]][["fit.param"]] <-  variogram.object[[1L]][["fit.param"]][c("variance", "snugget", "nugget", tmp)]

  ## set snugget to zero if snugget has not been specified or if there are
  ## no replicated observations

  variogram.object[[1L]][["param"]]["nugget"] <- sum(  variogram.object[[1L]][["param"]][c("nugget", "snugget") ])
  variogram.object[[1L]][["fit.param"]]["nugget"] <- any(  variogram.object[[1L]][["fit.param"]][c("nugget", "snugget") ])
  variogram.object[[1L]][["param"]]["snugget"] <- 0.
  variogram.object[[1L]][["fit.param"]]["snugget"] <- FALSE
  
  ## eliminate nugget and snugget from second and following components

  if( length( variogram.object ) > 1L ){

    variogram.object <- c(
      variogram.object[1L],
      lapply( variogram.object[-1L], function(x){
          sel <- names( x[["param"]] )[!names(x[["param"]]) %in% c( "snugget", "nugget" )]
          x[["param"]] <- x[["param"]][sel]
          x[["fit.param"]] <- x[["fit.param"]][sel]
          x
        }
      )
    )

  }
  
  ## adjust zero initial values of variance and scale parameters if they
  ## are fitted

  variogram.object <- lapply(
    variogram.object,
    function(x){
      sel <- names(x[["param"]]) %in% c("variance", "snugget", "nugget", "scale") &
        x[["fit.param"]] & x[["param"]] <= 0.
      x[["param"]][sel] <- sqrt( .Machine$double.eps )
      x
    }
  )
  
  ## transform variogram parameters

  param.tf <-  control[["param.tf"]]
  fwd.tf <- control[["fwd.tf"]]
  bwd.tf <- control[["bwd.tf"]]
  
  tmp <- f.aux.tf.param.fwd( variogram.object, param.tf, fwd.tf )

  transformed.param.aniso <- tmp[["transformed.param.aniso"]]
  tf.param.aniso <- tmp[["tf.param.aniso"]]
  fit.param.aniso <- tmp[["fit.param.aniso"]]

  ## select lag distances that are used for fitting
  
  t.lag.select <- sv[["lag.dist"]] <= max.lag & sv[["npairs"]] >= min.npairs
  
  ##  create environment to store items required to compute likelihood and
  ##  estimating equations that are provided by
  ##  likelihood.calculations

  envir <- new.env()
  variogram.item <- list()

  ##  initialize values of variogram object item stored in the environment

  variogram.item[["variogram.object"]] <- lapply(
    variogram.object,
    function(x) x[c("variogram.model", "param", "isotropic",
      "aniso", "sincos", "sclmat", "rotmat")]
  )

  assign( "variogram.item", variogram.item, pos = as.environment( envir ) )
  
  ## prepare lag vectors
  
  lag.vectors <- as.matrix( 
    sv[t.lag.select, c( "lag.x", "lag.y", "lag.z" )[1L:attr( sv, "ndim")]] 
  )
  attr( lag.vectors, "ndim.coords" ) <- attr( sv, "ndim")
  
  ## fit the model
    
  if(identical( control[["maximizer"]], "optim" )){
  
    r.fit <- optim(
      par = transformed.param.aniso[fit.param.aniso],
      fn = f.aux,
      method = control[["optim"]][["method"]],
      lower = control[["optim"]][["lower"]],
      upper = control[["optim"]][["upper"]],
      control = control[["optim"]][["control"]],
      hessian = control[["hessian"]],
      envir = envir,
      fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
      name.param.aniso = names(transformed.param.aniso),
      tf.param.aniso = tf.param.aniso,
      bwd.tf = bwd.tf,
      lag.vectors = lag.vectors,
      gamma = sv[["gamma"]][t.lag.select],
      npairs = sv[["npairs"]][t.lag.select],
      weighting.method = weighting.method,
      d2r = d2r,
      verbose = verbose
    )
    
  } else {
    
    r.fit <- nlminb(
      start = transformed.param.aniso[fit.param.aniso],
      objective = f.aux,
      lower = control[["nlminb"]][["lower"]],
      upper = control[["nlminb"]][["upper"]],
      control = control[["nlminb"]][["control"]],
      envir = envir,
      fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
      name.param.aniso = names(transformed.param.aniso),
      tf.param.aniso = tf.param.aniso,
      bwd.tf = bwd.tf,
      lag.vectors = lag.vectors,
      gamma = sv[["gamma"]][t.lag.select],
      npairs = sv[["npairs"]][t.lag.select],
      weighting.method = weighting.method,
      d2r = d2r,
      verbose = verbose
    )
    
    ## compute hessian
    
    if(control[["hessian"]]){
      r.fit[["hessian"]] <- optimHess(
        par = r.fit[["par"]],
        fn = f.aux,
        control = control[["optim"]][["control"]],
        envir = envir,
        fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
        name.param.aniso = names(transformed.param.aniso),
        tf.param.aniso = tf.param.aniso,
        bwd.tf = bwd.tf,
        lag.vectors = lag.vectors,
        gamma = sv[["gamma"]][t.lag.select],
        npairs = sv[["npairs"]][t.lag.select],
        weighting.method = weighting.method,
        d2r = d2r,
        verbose = verbose
      )
    }   
    
    ## rename components 
    
    tmp <- names( r.fit )
    tmp[match(c("objective", "evaluations"), tmp)] <- c("value", "counts")
    names( r.fit ) <- tmp
    
  }
  
  ## get variogram.item
  
  variogram.item <- get( "variogram.item", pos = as.environment( envir ) )
  
  ## convert angles to degree and add fit.param and fit.aniso
  
  variogram.item[["variogram.object"]] <- lapply(
    1L:length(variogram.item[["variogram.object"]]),
    function( i, x, vo, d2r ){
      x <- x[[i]]
      x[["aniso"]] <- x[["aniso"]]  / c( rep( 1., 2L ), rep( d2r, 3L ) )
      x <- c( x, 
        fit.param = list( vo[[i]][["fit.param"]] ), 
        fit.aniso = list( vo[[i]][["fit.aniso"]] )
      )
      x[c(
        "variogram.model", "param", "fit.param", 
        "isotropic", "aniso", "fit.aniso", "sincos", "rotmat", "sclmat"
      )]
    }, x = variogram.item[["variogram.object"]], vo = variogram.object, d2r = d2r
  )
    
  ## collect results
  
  r.result <- list(
    sse = r.fit[["value"]],
    variogram.object = variogram.item[["variogram.object"]],
    param.tf = param.tf,
    fwd.tf = fwd.tf,
    bwd.tf = bwd.tf,
    converged = if( sum( fit.param.aniso ) == 0L ){ 
      NA
    } else {
      r.fit[["convergence"]] == 0L
    },
    convergence.code = r.fit[["convergence"]],      
    iter = r.fit[["counts"]],
    call = cl,
    residuals = variogram.item[["residuals"]],
    fitted = variogram.item[["fitted"]],
    weights = variogram.item[["weights"]]
  )
  if( control[["hessian"]] ) r.result[["hessian"]] <- r.fit[["hessian"]]
  
  class( r.result ) <- "fitted.variogram"
  
  attr( r.result, "ndim" )               <- attr( sv, "ndim" )
  attr( r.result, "lag.dist.def" )       <- attr( sv, "lag.dist.def" )
  attr( r.result, "xy.angle.mid.class" ) <- attr( sv, "xy.angle.mid.class" )
  attr( r.result, "xz.angle.mid.class" ) <- attr( sv, "xz.angle.mid.class" )
  
  return( r.result )
}

## ##############################################################################

print.fitted.variogram <- 
  function( 
    x, digits = max(3, getOption("digits") - 3), ...
  )
{
  
  ## print method for fitted.variogram
  
  ## 2012-04-13 A. Papritz
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2016-08-18 AP changes for nested variogram models
  
  ## print variogram parameters
  
  cat("\n")
  
  lapply( 
    x[["variogram.object"]], 
    function(x){
      
      cat( "Variogram: ", x[["variogram.model"]], "\n" )
      
      param <- x[["param"]]
      names( param ) <- ifelse(
        x[["fit.param"]],
        names( param ),
        paste( names( param ), "(fixed)", sep = "" )
      )
      print( 
        format( param, digits = digits, 
          width = 16L, justify = "right" ), print.gap = 2L, 
        quote = FALSE
      )
      
      cat("\n")
      
      if( !x[["isotropic"]] ){
        
        aniso <- x[["aniso"]]
        names( aniso ) <- ifelse(
          x[["fit.aniso"]],
          names( aniso ),
          paste( names( aniso ), "(fixed)", sep = "" )
        )
        
        print( 
          format( aniso, digits = digits, 
            width = 16L, justify = "right" ), print.gap = 2L, 
          quote = FALSE
        )
        cat("\n")
      }
    }
  )
  
  invisible( x )
  
}

## ##############################################################################

summary.fitted.variogram <- 
  function (
    object, correlation = FALSE,
    signif = 0.95,
    ...
  )
{
    
  ## summary method for fitted.variogram
  
  ## 2012-12-10 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  ## 2016-08-18 AP changes for nested variogram models

  ans <- object[c(
    "call", "residuals", "weights", "converged", "convergence.code", 
    "iter", "sse", "variogram.object"
  )]
  
  ans[["param.aniso"]] <- lapply(
    ans[["variogram.object"]],
    function(object){
      res <- as.matrix( object[["param"]], ncol = 1L )
      if( !object[["isotropic"]] ) res <- rbind( 
        res, as.matrix( object[["aniso"]], ncol = 1L )
      )
      colnames( res ) <- "Estimate"
      res
    }
  )
  
  ## compute confidence intervals of variogram parameters from observed
  ## Fisher information matrix (Gaussian REML only)
  
  if( !is.null( object[["hessian"]] ) ){
    
    ## initialization
    
    cor.tf.param <- cov.tf.param <- matrix( 
      NA, nrow = NROW( object[["hessian"]] ), ncol = NROW( object[["hessian"]] ),
      dimnames = dimnames( object[["hessian"]] )
    )
    
#     se <- rep( NA_real_, NROW( object[["hessian"]] ) )
#     names( se ) <- rownames( object[["hessian"]])
#     
#     ci <- matrix( NA_real_, nrow = NROW( object[["hessian"]] ), ncol = 2L )
#     colnames( ci ) <- c( "Lower", "Upper" )
#     rownames( ci ) <- rownames( object[["hessian"]] )
    
    ## select parameters that are not on boundary of parameter space
    
    sr  <- !apply( object[["hessian"]], 1L, function( x ) all( is.na( x ) ) )
    
    if( sum( sr ) > 0L ){
      
      t.chol <- try( chol( object[["hessian"]][sr, sr, drop = FALSE] ), silent = TRUE )
      
      if( !identical( class( t.chol ), "try-error" ) ){
        
        ## compute covariance matrix of fitted transformed parameters
        
        cov.tf.param <- chol2inv( t.chol )
        dimnames( cov.tf.param ) <- dimnames( t.chol )
        
        ## correlation matrix and standard errors of fitted transformed
        ## parameters
        
        if( correlation ){
          ans[["cor.tf.param"]] <- cov2cor( cov.tf.param )
          colnames( ans[["cor.tf.param"]] ) <- rownames( ans[["cor.tf.param"]] ) <-
            gsub( ".__...__.", ".", colnames( ans[["cor.tf.param"]] ), fixed = TRUE )
        }
          
        se <- sqrt( diag( cov.tf.param ) )
                
        tmp <- f.aux.tf.param.fwd( 
          ans[["variogram.object"]], object[["param.tf"]],
          object[["fwd.tf"]]
        )
        
        param <- tmp[["transformed.param.aniso"]][tmp[["fit.param.aniso"]]]
        param <- param[names(param) %in% names(se)]
        
        param.tf <- tmp[["tf.param.aniso"]][names(param)]

        se <- se[match( names(se), names(param))]
        
        ## confidence intervals
        
        ci <- t( sapply(
            1L:length(se),
            function( i, m, se, signif ){
              m[i] + c(-1., 1.) * se[i] * qnorm( (1.-signif)/2., lower.tail = FALSE ) 
            }, m = param, se = se, signif = signif
          ))
        colnames( ci ) <- c("Lower", "Upper")
        rownames( ci ) <- names( se )

        ci <- apply(
          ci, 2L,
          function( x, nme, bwd.tf, param.tf ){
            sapply( nme,
              function( x, bwd.tf, param.tf, param ) bwd.tf[[param.tf[x]]]( param[x] ),
              bwd.tf = bwd.tf, param.tf = param.tf, param = x
            )
          }, nme = names(se), bwd.tf = object[["bwd.tf"]],
          param.tf = param.tf
        )
        
        colnames( ci ) <- c("Lower", "Upper")
        rownames( ci ) <- names( se )
        
        ## convert to list
        
        tmp <- strsplit( rownames(ci), ".__...__.", fixed = TRUE )
        name.tmp <- rownames(ci) <- sapply( tmp, function(x) x[1L] )
        cmp <- sapply( tmp, function(x) x[2L] )
        
        ci <- tapply( 
          1L:NROW(ci), factor( cmp ), 
          function( i, ci ){
            ci[i, , drop = FALSE]
          }, ci = ci, simplify = FALSE
        ) 
        
        ## merge into ans[["param.aniso"]]
        
        ans[["param.aniso"]] <- lapply(
          1L:length(ans[["param.aniso"]]),
          function( i, pa, ci ){
            pa <- pa[[i]]
            if( i <= length(ci) ){
              ci <- ci[[i]]
            } else {
              ci <- matrix( rep( NA_real_, 2. * NROW( pa ) ), ncol = 2L )
              rownames( ci ) <- rownames(pa)
            }
            
            pa <- cbind( pa, Lower = rep( NA_real_, NROW(pa) ),
              Upper = rep( NA_real_, NROW(pa) ) )
            i <- match( rownames( ci ), rownames( pa ) )
            pa[i, 2L:3L] <- ci
            tmp <- rownames( pa )
            tmp[is.na(pa[, 2L])] <- paste( tmp[is.na(pa[, 2L])], "(fixed)", sep="" )
            rownames(pa) <- tmp
            pa
          }, pa = ans[["param.aniso"]], ci = ci
        )
        
      } else {
        warning(
          "Hessian not positive definite:",
          "\nconfidence intervals of variogram parameters cannot be computed"        
        )
      }
    } 
  }
    
  class( ans ) <- c( "summary.fitted.variogram" )
  
  ans
}

## ##############################################################################

print.summary.fitted.variogram <- 
  function (
    x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"),
    ...
  ) 
{
    
  ## print.summary method for fitted.variogram
  
  ## 2012-04-13 A. Papritz
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  ## 2016-08-18 AP changes for nested variogram models
  
  cat("\nCall:")
  cat( paste( deparse(x[["call"]]), sep = "\n", collapse = "\n"),  "\n", sep = "" )
  
  if( is.na( x[["converged"]] ) ){
    cat( "\nEstimation with fixed variogram parameters\n" )
    
  } else {
    
    if(!(x[["converged"]])) {
      cat( 
        "\nAlgorithm did not converge, diagnostic code: ", 
        x[["convergence.code"]], "\n"
      )
    } else {
      cat(
        "\nConvergence in", x[["iter"]][1L], "function and", 
        x[["iter"]][2L], "Jacobian/gradient evaluations\n"
      )
    }
    
    cat(
      "\nResidual Sum of Squares:", 
      x[["sse"]], "\n"
    )
    
  }
    
  resid <- x[["residuals"]]
  cat( "\nResiduals (epsilon):\n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure( quantile(resid), names = nam )
  print( rq, digits = digits, ...)
  
  lapply(
    1L:length(x[["param.aniso"]]),
    function(i, x, vo){
      x <- x[[i]]
      tmp <- x[ ,1L]
      tmp <- tmp[!is.na(tmp) & tmp > 0. ]
      t.digits <- -floor( log10( min( tmp ) ) )
      cat( "\nVariogram: ", vo[[i]][["variogram.model"]], "\n" )
      printCoefmat(
        x, digits = max( digits, t.digits) , signif.stars = FALSE, ...
      )
      #       print(format(
#         signif( x, digits = 7L ), width = 16L, scientific = TRUE, justify = "left"), quote = FALSE, ...
#       )
    }, x = x[["param.aniso"]], vo = x[["variogram.object"]]
  )  
  
  if( !is.null( x[["cor.tf.param"]] ) ){
    
    correl <- x[["cor.tf.param"]]
    p <- NCOL(correl)
    if( p > 1L ){
      cat("\nCorrelation of (transformed) variogram parameters:\n")
      correl <- format(round(correl, 2L), nsmall = 2L, 
        digits = digits)
      correl[!lower.tri(correl)] <- ""
      print(correl[-1L, -p, drop = FALSE], quote = FALSE)
    }
    
  }
  
    
  invisible( x )
}

## ##############################################################################
plot.georob <- 
  function(
    x, what = c( "variogram", "covariance", "correlation", 
      "ta", "sl", "qq.res", "qq.ranef" ),
    add = FALSE,
    lag.dist.def, 
    xy.angle.def = c( 0., 180. ),
    xz.angle.def = c( 0., 180. ),
    max.lag = Inf,
    estimator = c( "mad", "qn", "ch", "matheron" ),
    mean.angle = TRUE,
    level = what != "ta", 
    smooth = what == "ta" || what == "sl",
    id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
    label.pos = c(4,2),    
    col, pch, xlab, ylab, main, lty = "solid", 
    ...    
  )
{
  
  ## Function plots the graph of a variogram fitted by f.georob
  
  ## 2012-12-11 A. Papritz
  ## 2012-12-21 AP correction for using col and pch 
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-05-08 AP changes for plotting covariances and correlations
  ## 2015-11-27 AP changes for Tukey-Anscombe, QQnorm plots and for plotting variograms
  ## 2016-08-19 AP changes for nested variogram models
  
  x[["na.action"]] <- NULL
  
  estimator <- match.arg( estimator )
  what <- match.arg( what )
  
  ## labelling of points in residual diagnostic plots (taken from plot.lmrob())
  
  n <- length(x[["fitted.values"]])
  if (is.null(id.n)){
    id.n <- 0L
  } else {
    id.n <- as.integer(id.n)
    if(id.n < 0L || id.n > n) stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
  }
  if(id.n > 0L){ ## label the largest residuals
    if(is.null(labels.id))
    labels.id <- paste(1L:n)
    iid <- 1L:id.n
    ##    show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
    ## if(any(show[2L:3L]))
    ##     show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
    text.id <- function(x, y, ind, adj.x = TRUE) {
      labpos <-
      if(adj.x) label.pos[1L+as.numeric(x > mean(range(x)))] else 3L
      text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
        pos = labpos, offset = 0.25)
    }
  }
  yh <- x[["fitted.values"]]
  
  switch(
    what,
    ta = {
      
      ## Tukey-Anscombe plot
      
      if( missing( col ) ) col <- 1L; if( missing( pch ) ) pch <- 1L
      if( missing( xlab ) ) xlab <- "Fitted values"
      if( missing( ylab ) ) ylab <- "Residuals"
      if( missing( main ) ) main <- "Residuals vs. Fitted"
      r <- residuals( x, level = level )
      if( !add ){
        plot( yh, r, col = col, pch = pch,
          xlab = xlab, ylab = ylab, main = main, ...
        )
      } else {
        points( yh, r, col = col, pch = pch, ... )
      }
      if( smooth ){
        lines( loess.smooth( yh, r, ... ), col = "red", lty = 1L, ... )
      }
      if(id.n > 0L){   ## adapted from plot.lmrob()
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        y.id <- r[show.r]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(yh[show.r], y.id, show.r)
      }
      
    },
    sl = {
      
      ## scale.location plot
      
      if( missing( col ) ) col <- 1L; if( missing( pch ) ) pch <- 1L
      if( missing( xlab ) ) xlab <- "Fitted values"
      if( missing( ylab ) ) ylab <- "Sqrt of abs(Residuals)"
      if( missing( main ) ) main <- "Scale-Location"
      r <- sqrt( abs( residuals( x, level = level ) ) )
      if( !add ){
        plot( yh, r, col = col, pch = pch,
          xlab = xlab, ylab = ylab, main = main, ...
        )
      } else {
        points( yh, r, col = col, pch = pch, ... )
      }
      if( smooth ){
        lines( loess.smooth( yh, r, ... ), col = "red", lty = 1L )
      }
      if(id.n > 0L) {   ## adapted from plot.lmrob()
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        y.id <- r[show.r]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(yh[show.r], y.id, show.r)
      }
      
    },
    qq.res = {
      
      ## qqnorm of standardized residuals (eps or eps + b, depending of
      ## value of level passed in ...)
      
      if( missing( col ) ) col <- 1L; if( missing( pch ) ) pch <- 1L
      if( missing( xlab ) ) xlab <- "Theoretial quantiles"
      if( missing( ylab ) ) ylab <- "Standardized residuals"
      if( missing( main ) ) main <- "Normal Q-Q residuals"
      r <- rstandard( x, level = level )
      tmp <- qqnorm( r, col = col, pch = pch,
        xlab = xlab, ylab = ylab, main = main, 
        plot.it = !add, ...
      )
      if( add ) points( y~x, tmp, col = col, pch = pch, ... )
      if(id.n > 0L){   ## adapted from plot.lmrob()
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        text.id(tmp$x[show.r], tmp$y[show.r], show.r)
      }
    },
    qq.ranef = {
      
      ## qqnorm of standardized random effects
      
      if( missing( col ) ) col <- 1L; if( missing( pch ) ) pch <- 1L
      if( missing( xlab ) ) xlab <- "Theoretial quantiles"
      if( missing( ylab ) ) ylab <- "Standardized random effects"
      if( missing( main ) ) main <- "Normal Q-Q random effects"
      tmp <- qqnorm( ranef( x, standard = TRUE ), 
        col = col, pch = pch,
        xlab = xlab, ylab = ylab, main = main, 
        plot.it = !add, ...
      )
      if( add ) points( y~x, tmp, col = col, pch = pch, ... )
      
    },
    {
      
      ## plotting variogram, covariance or correlation
      
      ## compute and plot sample variogram
      
      if( !missing( lag.dist.def ) ){
        
        ## compute and plot sample variogram of regression residuals
        
        r.sv <- sample.variogram(
          x, 
          lag.dist.def = lag.dist.def,
          xy.angle.def = xy.angle.def, 
          xz.angle.def = xz.angle.def,
          max.lag = max.lag,
          estimator = estimator,
          mean.angle = mean.angle
        )
        
        if( missing( col ) ){
          col <- 1L:nlevels( r.sv[["xy.angle"]] )
        } else if( length( col ) < nlevels( r.sv[["xy.angle"]] ) ) stop(
          "number of colors less than number of directions in x-y-plane for which semivariances are computed"
        )
        if( missing( pch ) ){
          pch <- 1L:nlevels( r.sv[["xz.angle"]] )
        } else if( length( pch ) < nlevels( r.sv[["xz.angle"]] ) ) stop(
          "number of colors less than number of directions in x-z-plane for which semivariances are computed"
        )
        if( missing( xlab ) ) xlab <- "lag distance"
        if( missing( ylab ) ) ylab <- "semivariance"
        
        plot( 
          r.sv, add = add, col = col, pch = pch, xlab = xlab, ylab = ylab, ... 
        )
        
        xmax <- max( r.sv[["lag.dist"]] )
        xy.angle.mid.class <- attr( r.sv, "xy.angle.mid.class" )
        xz.angle.mid.class <- attr( r.sv, "xz.angle.mid.class" )
        
      } else {
        
        ## setup window for plotting variogram, covariance, correlation model
        
        if( is.finite( max.lag ) ){
          xmax <- max.lag
        } else {
          xmax <- sqrt( max( rowSums( x[["locations.objects"]][["lag.vectors"]]^2 ) ) )
        }
        
        sill <- sum( unlist( lapply( 
            x[["variogram.object"]], 
            function(x) x[["param"]][names(x[["param"]]) %in% c("variance", "snugget", "nugget")]
            )))
        
        ymax <- 1.1 * switch(
          what,
          correlation = 1.,
          sill
        )
        
        ## see sample.variogram.default
        
        # join first and last angle classes if end points match and there is
        # more than one angle
        
        # xy-plane
        
        n <- length( xy.angle.def )
        d <- diff( xy.angle.def )
        xy.angle.mid.class <- 0.5 * ( xy.angle.def[-1L] + xy.angle.def[-n] )
        if( 
          n > 2L &&
          identical( xy.angle.def[1L], 0. ) && 
          identical( xy.angle.def[n], 180. ) &&
          !all( d[1L] == d[-1L] )
        ){
          
          xy.angle.mid.class[1L] <- xy.angle.mid.class[1L] - (180. - xy.angle.mid.class[n-1L])
          xy.angle.mid.class <- xy.angle.mid.class[-(n-1L)]
        }
        
        # xz-plane
        
        n <- length(xz.angle.def)
        d <- diff( xz.angle.def )
        xz.angle.mid.class <- 0.5 * ( xz.angle.def[-1L] + xz.angle.def[-n] )
        if( 
          n > 2L &&
          identical( xz.angle.def[1L], 0. ) && 
          identical( xz.angle.def[n], 180. ) &&
          !all( d[1L] == d[-1L] )
        ){
          
          xz.angle.mid.class[1L] <- xz.angle.mid.class[1L] - (180. - xz.angle.mid.class[n-1L])
          xz.angle.mid.class <- xz.angle.mid.class[-(n-1L)]
        }
                
        if( missing( col ) ) col <- 1L:length( xy.angle.mid.class )
        if( missing( pch ) ) pch <- 1L:length( xz.angle.mid.class )
        if( missing( xlab ) ) xlab <- "lag distance"
        if( missing( ylab ) ) ylab <- what
        
        if( !add ){
          plot( c(0., xmax), c(0., ymax ), xlab = xlab, 
            ylab = ylab, type = "n", ... )
        }
        
      }
      
      ## add graph of fitted variogram/covariance/correlation model
      
      lines( 
        x, 
        what,
        to = xmax,
        xy.angle = xy.angle.mid.class,
        xz.angle = xz.angle.mid.class,
        col = col, pch = pch, lty = lty, ...
      )
      
    }
    
  )
  
  invisible()
  
}

 ##############################################################################

lines.georob <- 
  function( 
    x, 
    what = c("variogram", "covariance", "correlation"),
    from = 1.e-6, to, n = 501, 
    xy.angle = 90.,
    xz.angle = 90.,
    col = 1:length( xy.angle ), pch = 1:length( xz.angle ), lty = "solid", ...
  )
{
  
  ## Function plots the graph of a variogram fitted by f.georob
  
  ## 2012-12-12 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-05-08 AP changes for plotting covariances and correlations
  ## 2016-08-19 AP changes for nested variogram models
    
  d2r <- pi / 180.
  
  what <- match.arg( what )
  
  tmp <- sapply( x[["variogram.object"]], function(x) x[["variogram.model"]] )
  if( 
    any( tmp  %in% control.georob()[["irf.models"]] ) && 
    what != "variogram" 
  ) stop(
    "stationary covariance and correlation does not exist for intrinsic variogram models"  
  )
  
  ## generate grid of angle classes
  
  angle <- expand.grid(
    xy = xy.angle * d2r,
    xz = xz.angle * d2r
  )
  
  ## set up lag distances
  
  if( missing( to ) ){
    u <- par( "usr" )[1L:2L]
    to <- (u[2L] + 0.04/1.04*u[1L]) / (1.04 - 0.04^2/1.04)
  }
  
  lag.class <- seq( from, to, length = n )
  
  ## determine number of dimensions
  
  if( identical( class( x ), "fitted.variogram" ) ){
    ndim <- attr( x, "ndim" )
  } else {
    ndim <- NCOL( x[["locations.objects"]][["coordinates"]] )
  }
  
  ## loop over all angles
  
  t.bla <- lapply(
    1L:NROW( angle ),
    function( 
      i, what, angle, lag.class, 
      variogram.object, 
      nxy, nxz, ndim, col, pch, lty, ...
    ){
      
      ## generate lag vectors
      
      t.aux <- lag.class * sin( angle[i, "xz"] )
      lag.vectors <- cbind(
        t.aux * sin( angle[i, "xy"] ),
        t.aux * cos( angle[i, "xy"] ),
        lag.class * cos( angle[i, "xz"] )
      )
      
      ## drop unneeded components
      
      lag.vectors <- lag.vectors[, 1L:ndim, drop = FALSE ]
      attr( lag.vectors, "ndim.coords" ) <- NCOL( lag.vectors )
      
      ## compute the generalized covariance
      
      Valpha <- f.aux.gcr( 
        lag.vectors = lag.vectors,
        variogram.object = x[["variogram.object"]], symmetric = FALSE,
        control.pcmp = control.pcmp(max.ncores=1L), ...
      )
      
      if( any( sapply( Valpha, function(x) x[["error"]] ) ) ) stop( 
        "\nan error occurred when computing semivariances\n" 
      )
      
      ## compute semivariance of signal and sum up
      
      r.gamma <- rowSums( sapply(
          1L:length(Valpha), 
          function( i, x, variogram.object ){
            x <- Valpha[[i]]
            -( x[["Valpha"]] - x[["gcr.constant"]] ) * variogram.object[[i]][["param"]]["variance"]
          }, x = Valpha, variogram.object = x[["variogram.object"]]
        ))
      
      ## add nugget
      
      sel <- sqrt( rowSums( lag.vectors^2 ) ) > 0.
      r.gamma[sel] <- r.gamma[sel] + sum( 
        x[["variogram.object"]][[1L]][["param"]][c("nugget", "snugget")] 
      )
            
      ## plot semivariance
      
      sel.pch <- ((i-1L) %/% nxy) + 1L
      sel.col <- i - (sel.pch-1L) * nxy
      type <- if( nxz > 1L ) "o" else "l"
      
      sill <- sum( unlist( lapply( 
            x[["variogram.object"]], 
            function(x) x[["param"]][names(x[["param"]]) %in% c("variance", "snugget", "nugget")]
          )))
      
      r.gamma <- switch(
        what,
        variogram = r.gamma,
        covariance = sill - r.gamma,
        correlation = (sill - r.gamma ) / sill
      )

      lines( 
        lag.class, r.gamma, 
        col = col[sel.col], pch = pch[sel.pch], 
        lty = lty, type = type, ... )
    },
    what = what,
    angle = angle,
    lag.class = lag.class,
    variogram.object = x[["variogram.object"]],
    nxy = length( xy.angle ),
    nxz = length( xz.angle ),
    ndim = ndim,
    col = col, pch = pch, lty = lty
  )
  
  invisible( NULL )
  
}

lines.fitted.variogram <- lines.georob
