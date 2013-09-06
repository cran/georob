##  ############################################################################

cv <- function( object, ... ) UseMethod( "cv" )

##  ############################################################################

cv.georob <-
  function(
    object,
    formula = NULL, subset = NULL,
    nset = 10, seed = NULL, sets = NULL,
    duplicates.in.same.set = TRUE,
    re.estimate = TRUE, param = object[["param"]], 
    fit.param = object[["initial.objects"]][["fit.param"]],
    aniso = object[["aniso"]][["aniso"]], 
    fit.aniso = object[["initial.objects"]][["fit.aniso"]],
    return.fit = FALSE, reduced.output = TRUE,
    lgn = FALSE,
    ncores = min( nset, detectCores() ),
    verbose = 0,
    ...
  )
{
  
#  \$([[:alnum:]\.]+)([\^\r,$\[\] \(\)])          [["\1"]]\2
  ## Function computes nset-fold cross-validation predictions from a
  ## fitted georob object
  
  ## Arguments:
  
  ## object        fitted georob object
  ## formula       a formula passed by update to georob
  ## nset          integer scalar for the number of cross-validation subsets
  ## seed integer scalar passed to set.seed before selecting the 
  ##               cross-valdation subsets by a call to runif()
  ## sets          an integer vector with length nrow(data) defining the 
  ##               cross-validation sets and over-riding the values provided 
  ##               for nset and seed
  ## duplicates.in.same.set   logical flag controlling whether replicated observations
  ##               at a given location are assigned to the same cross-validation set
  ## re.estimate   logical flag controlling whether the variogram parameters should 
  ##               be re-estimated for each cross-validation subset
  ## param         initial values of variogram parameters when the variogram is 
  ##               re-estimated for each cross-validation subset
  ## return.fit    logical flag to control whether the information about the fit are
  ##               should be returned for each cross-valdiation subset when re-estimating the
  ##               model
  ## reduced.output    logical flag controlling whether for each cross-valdiation subset the
  ##               the full fitted object or just a selection (information about convergence,
  ##               variogram and fixed-effects parameter estimates) should be returned when 
  ##               re-estimating the model
  ## lgn   logical flag controlling whether lognormal kriging predictions should be computed 
  ## ncores        integer scalar with the number of cores to used in parallel processing
  ## verbose       integer scalar, controlling verbosity of the information sent to standard output
  ## ...           further arguments passed by update to georob or to 
  ##               mclapply on non-windows platforms
  
  ## ToDos:
  
  ## - Klasse und Methoden definieren fuer cv (kompatibel mit geoR)
  
  ## History:
  
  ## 2011-10-24 Korrektur Ausschluss von nichtbenoetigten Variablen fuer lognormal kriging
  ## 2011-12-23 AP modified for replicated observations and for parallel computing
  ## 2012-03-02 AP eliminated possibility for logging to file in parallel processing
  ## 2012-03-19 AP correction of error in parallel processing on Windows
  ## 2012-05-01 AP correct handling of NAs
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2012-05-09 AP correction of error if a new formula is passed via update to georob
  ## 2012-05-22 AP correction of error in passing param and fit.param to georob
  ## 2012-06-05 AP correction of error in handling optional subset argument
  ## 2012-11-04 AP handling compressed cov.betahat
  ## 2012-12-04 AP modifiction for changes in predict.georob
  ## 2013-04-24 AP changes for parallelization on windows os
  ## 2013-05-23 AP correct handling of missing observations
  ## 2013-05-24 AP separate initial variogram parameters for each cross-validation set
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-02 AP passing initial values of aniso and fit.aniso to georob via update
  ## 2013-07-05 AP return "variogram.model" as part of fit componnent
    
  ## auxiliary function that fits the model and computes the predictions of
  ## a cross-validation set
  
  f.aux <- function( 
    ..i.., object, formula, data, sets, re.estimate, 
    param, fit.param, aniso, fit.aniso, lgn, verbose, ...
  ){  ## cv function
    
    if (verbose) cat( "\n\n  processing cross-validation set", ..i.., "\n" ) 
    
    ## fit model to complement of current set
    
    if( !re.estimate ){
      fit.param <- c( 
        variance = FALSE, snugget = FALSE, nugget = FALSE, scale = FALSE, 
        a = FALSE, alpha = FALSE, beta = FALSE, delta = FALSE, 
        gamma = FALSE, lambda = FALSE, n = FALSE, nu = FALSE,
        f1 = FALSE, f2  =FALSE, omega = FALSE, phi = FALSE, zeta = FALSE      
      )[names( param )]
      fit.aniso <- c( f1 = FALSE, f2 = FALSE, omega = FALSE, phi = FALSE, zeta = FALSE )
    }
    
    ## change environment of terms and formula so that subset selection works for update            
    
    environment( formula ) <- environment()
    environment( object[["terms"]] ) <- environment()
    
    ## read-off initial values of variogram parameters
    
    if( ( is.matrix( param ) || is.data.frame( param ) ) )  param <- param[..i..,]
    if( ( is.matrix( fit.param ) || is.data.frame( fit.param ) ) )  fit.param <- fit.param[..i..,]
    if( ( is.matrix( aniso ) || is.data.frame( param ) ) )  aniso <- aniso[..i..,]
    if( ( is.matrix( fit.aniso ) || is.data.frame( fit.aniso ) ) )  fit.aniso <- fit.aniso[..i..,]

    t.georob <- update( 
      object, 
      formula = formula,
      data = data,
      subset = -sets[[..i..]] ,
      param = param, fit.param = fit.param,
      aniso = aniso, fit.aniso = fit.aniso,
      verbose = verbose,
      ...
    )
    
    if( verbose > 0 ){
      cat( "\n\n" )
      print( summary( t.georob ) )
    }
    
    ## compute predictions for current set
    
    t.predict <- predict( 
      t.georob, newdata = data[sets[[..i..]], ], type = "response",
      mmax = length( sets[[..i..]] ),
      extended.output = lgn,
      ncores = 1
    )
    
    ## backtransformation for log-normal kriging
    
    if( lgn ){
      t.predict <- lgnpp( t.predict )
      t.predict <- t.predict[, -match( 
        c( "trend", "var.pred", "cov.pred.target", "var.target" ), names( t.predict ) 
      )]
    }
    
    t.predict <- data.frame( i = sets[[..i..]], t.predict )
    
    t.ex <- c( 
      grep( "lower", colnames( t.predict ), fixed = TRUE ),
      grep( "upper", colnames( t.predict ), fixed = TRUE )
    )
    
    t.predict <- t.predict[, -t.ex]
    
    if( reduced.output ){
      
      if( !is.null( t.georob[["cov"]][["cov.betahat"]] ) ){
        t.se.coef <- sqrt( diag( expand( t.georob[["cov"]][["cov.betahat"]] ) ) )
      } else {
        t.se.coef <- NULL
      }
      
      t.georob <- t.georob[c(  
        "tuning.psi", "converged", "convergence.code",
        "gradient", "variogram.model", "param", "aniso",
        "coefficients"
      )]
      
      t.georob[["aniso"]] <- t.georob[["aniso"]][["aniso"]]
      
      if( !is.null( t.se.coef ) ) t.georob[["se.coefficients"]] <- t.se.coef
      
    }
    
    return( list( pred = t.predict, fit = t.georob ) )
    ## end cv function
  }
  
  ## redefine na.action component of object
  
  if( identical( class( object[["na.action"]] ), "exclude" ) ){
    class( object[["na.action"]] ) <- "omit"
  }
  
  ## update terms of object is formula is provided
  
  if( !is.null( formula ) ){
    formula <- update( formula( object ), formula )
    object[["terms"]] <- terms( formula )
  } else {
    formula <- formula( object )
  }
    
  ## get data.frame with required variables (note that the data.frame passed
  ## as data argument to georob must exist in GlobalEnv)
  
  data <- cbind(
    get_all_vars( formula( object ), eval( getCall(object)[["data"]] ) ),
    get_all_vars( object[["locations.objects"]][["locations"]], eval( getCall(object)[["data"]] ) )
  )
  
  ## select subset if appropriate
  
  if( !is.null( subset ) ){
    data <- data[subset, ]
    object[["Tmat"]] <- object[["Tmat"]][subset]
  } else if( !is.null( getCall(object)[["subset"]] ) ){
   data <- data[eval( getCall(object)[["subset"]] ), ]
  }
  
#   if( !is.null( getCall(object)[["subset"]] ) )
   
  ## define cross-validation sets
  
  if( is.null( sets ) ){
    
    if( !is.null( seed ) ) set.seed( seed )
    sets <- runif( NROW( data ) )
    sets <- cut( 
      sets, 
      breaks = c( -0.1, quantile( sets, probs = ( 1:(nset-1)/nset ) ), 1.1 )
    )
    sets <- factor( as.integer( sets ) )

  } else {
    
    nset <- length( unique( sets ) )
    if( length( sets ) != NROW( data ) ) stop(
      "sets must be an integer vector with length equal to the number of observations"      
    )
    
  }
  
  if( duplicates.in.same.set ){
    dups <- duplicated( object[["Tmat"]] )
    idups <- match( object[["Tmat"]][dups], object[["Tmat"]][!dups] )
    sets[dups] <- (sets[!dups])[idups]
  }
  
  sets <- tapply(
    1:NROW( data ),
    sets,
    function( x ) x
  )
  
  ## check dimension of param and fit.param
  
  if( ( is.matrix( param ) || is.data.frame( param ) ) && nrow( param )!= nset ) stop(
    "'param' must have 'nset' rows if it is a matrix or data frame"  
  )
    
  if( ( is.matrix( fit.param ) || is.data.frame( fit.param ) ) && nrow( param )!= nset ) stop(
    "'fit.param' must have 'nset' rows if it is a matrix or data frame"  
  )
    
  if( ( is.matrix( aniso ) || is.data.frame( aniso ) ) && nrow( aniso )!= nset ) stop(
    "'aniso' must have 'nset' rows if it is a matrix or data frame"  
  )
    
  if( ( is.matrix( fit.aniso ) || is.data.frame( fit.aniso ) ) && nrow( aniso )!= nset ) stop(
    "'fit.aniso' must have 'nset' rows if it is a matrix or data frame"  
  )
    
  ## loop over all cross-validation sets
  
  if( .Platform[["OS.type"]] == "windows" ){
    
    ## create a SNOW cluster on windows OS
    
    cl <- makePSOCKcluster( ncores, outfile = "")
    
    ## export required items to workers
    
    junk <- clusterEvalQ( cl, require( georob, quietly = TRUE ) )
    
    t.result <- parLapply(
      cl, 
      1:length( sets ),
      f.aux, 
      object = object,
      formula = formula,
      data = data,
      sets = sets,
      re.estimate = re.estimate,
      param = param, fit.param = fit.param,
      aniso = aniso, fit.aniso = fit.aniso,
      lgn = lgn,
      verbose = verbose,
      ...
    )
    
    stopCluster(cl)
    
  } else {
        
    ## fork child processes on non-windows OS
    
    t.result <- mclapply(
      1:length( sets ),
      f.aux, 
      object = object,
      formula = formula,
      data = data,
      sets = sets,
      re.estimate = re.estimate,
      param = param, fit.param = fit.param,
      aniso = aniso, fit.aniso = fit.aniso,
      lgn = lgn,
      verbose = verbose,
      mc.cores = ncores,
      mc.allow.recursive = FALSE,
      ...
    )
    
  }
    
  ## create single data frame with cross-validation results 
  
  result <- t.result[[1]][["pred"]]
  result[["subset"]] <- rep( 1, nrow( t.result[[1]][["pred"]] ) )
  
  for( t.i in 2:length( t.result ) ) {
    result <- rbind( 
      result, 
      data.frame( 
        t.result[[t.i]][["pred"]],
        subset = rep( t.i, nrow( t.result[[t.i]][["pred"]] ) )
      )
    )
  }
  t.ix <- sort( result[["i"]], index.return = T )[["ix"]]
  result <- result[t.ix, ]
  result[["data"]] <- model.response( 
    model.frame( formula( object), data, na.action = na.pass ) 
  )
  
  if( lgn ) result[["lgn.data"]] <- exp( result[["data"]] )
  
  result <- result[, -match("i", colnames( result) )]
  
  isubset <- match( "subset", colnames( result ) )
  idata <- grep( "data", colnames( result ), fixed = TRUE )
  ipred<-  grep( "pred", colnames( result ), fixed = TRUE )
  ise <-   grep( "se", colnames( result ), fixed = TRUE )
  ise <- ise[ise != isubset]
  
  result <- cbind(
    result[, -c(isubset, idata, ipred, ise)],
    result[,  c(isubset, idata, ipred, ise)]
  )
  
  t.fit <- lapply( t.result, function( x ) return( x[["fit"]] ) )
  
  if( re.estimate && !all( sapply( t.fit, function(x) x[["converged"]] ) ) ) warning(
    "lack of covergence for  ", 
    sum( !sapply( t.fit, function(x) x[["converged"]] ) ), " cross-validation sets"
  )
  
  result <- list( 
    pred = result, 
    fit = if( return.fit ) t.fit else NULL
  )
  
  class( result ) <- "cv.georob"
  
  invisible( result )
  
}

##  ###########################################################################

plot.cv.georob <-
  function( 
    x, type = c( "sc", "lgn.sc", "ta", "qq", "pit", "mc", "bs" ), 
    ncutoff = NULL, 
    add = FALSE, 
    col, pch, lty,
    main, xlab, ylab, 
    ... 
  )
{
  
  ## plot method for class "cv.georob"  
  
  ## 2011-12-21 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  x <- x[["pred"]]
  
  type = match.arg( type )
  
  if( type == "sc.lgn" && !"lgn.pred" %in% names( x ) ) stop(
    "lognormal kriging results missing, use 'lgn = TRUE' for cross-validation"
  )
  
  if( type %in% c( "pit", "mc", "bs" ) ){
    
    result <- validate.predictions( 
      data = x[["data"]],
      pred = x[["pred"]],
      se.pred = x[["se"]],
      statistic = type, ncutoff = ncutoff
    )
    
  }
  
  if( missing( col ) ) col <- 1
  if( missing( pch ) ) pch <- 1
  if( missing( lty ) ) lty <- 1

  
    
  switch(
    type,
    sc = {
      
      ##  scatterplot of (transformed) measurements vs. predictions
      
      if( missing( main ) ) main <- "data vs. predictions"
      if( missing( xlab ) ) xlab <- "predictions"
      if( missing( ylab ) ) ylab <- "data"
      
      if( add ){
        points( data ~ pred, x, col = col, pch = pch, ... )
      } else {        
        plot( 
          data ~ pred, x, col = col, pch = pch, 
          main = main, xlab = xlab, ylab = ylab, ... 
        )
      }
      
      
    },        
    lgn.sc = {
      
      ##  scatterplot of original measurements vs.  back-transformded
      ##  lognormal predictions
      
      if( missing( main ) ) main <- "data vs. back-transformed predictions"
      if( missing( xlab ) ) xlab <- "back-transformed predictions"
      if( missing( ylab ) ) ylab <- "data"
      
      if( add ){
        points( lgn.data ~ lgn.pred, x, col = col, pch = pch, ... )
      } else {        
        plot( 
          lgn.data ~ lgn.pred, x, col = col, pch = pch, 
          main = main, xlab = xlab, ylab = ylab, ... 
        )
      }
      
      
      
    }, 
    ta = {
      
      ##  Tukey-Anscombe plot
      
      if( missing( main ) ) main <- "Tukey-Anscombe plot"
      if( missing( xlab ) ) xlab <- "predictions"
      if( missing( ylab ) ) ylab <- "standardized prediction errors"
      
      if( add ){
        points( I((data-pred)/se) ~ pred, x, col = col, pch = pch, ... )
      } else {        
        plot( 
          I((data-pred)/se) ~ pred, x, col = col, pch = pch, 
          main = main, xlab = xlab, ylab = ylab, ... 
        )
      }
      
    },
    qq = {
      
      ##  normal QQ-Plot of standardized prediction errors
      
      if( missing( main ) ) main <- "normal-QQ-plot of standardized prediction errors"
      if( missing( xlab ) ) xlab <- "quantile N(0,1)"
      if( missing( ylab ) ) ylab <- "quantiles of standardized prediction errors"
      
      r.qq <- with( x, qqnorm( ( data - pred ) / se, plot.it = FALSE ) )
      
      if( add ){
        points( r.qq, col = col, pch = pch, ... )
      } else {
        plot( r.qq, col = col, pch = pch, main = main, xlab = xlab, ylab = ylab, ... )
      }
    },
    pit = {
      
      ##  histogramm of probability-integral-transformation
      
      if( missing( main ) ) main <- "histogramm PIT-values"
      if( missing( xlab ) ) xlab <- "PIT"
      if( missing( ylab ) ) ylab <- "density"
      
      r.hist <- hist( 
        result, 
        col = col, lty = lty, 
        main = main, xlab = xlab, ylab = ylab, freq = FALSE, ... )
    },
    mc = {
      
      ##  narginal calibration plots: ecdf of measurements and mean
      ##  predictive cdf
      
      if( missing( main ) ) main <- "empirical cdf of data and mean predictive cdfs"
      if( missing( xlab ) ) xlab <- "data or predicitons"
      if( missing( ylab ) ) ylab <- "probability"
      
      matplot( 
        result[["y"]], 
        result[, c( "ghat", "fbar" )], type = "l",
        col = c( "black", "red" ),
        lty = c( "solid", "dashed" ),
        main = main, xlab = xlab, ylab = ylab,
        ...
      )
      
      t.usr <- par( "usr" )
      t.usr[3:4] <- with( result, range( fbar - ghat ) ) *1.04
      par( usr = t.usr )
      with( result, lines( y, fbar-ghat, col= "blue", lty = "dotted" ) )
      axis(2, pos = t.usr[2], col.axis = "blue", col.ticks = "blue" )
      legend( 
        "topleft",
        lty = c("solid", "dashed", "dotted" ), 
        col = c( "black", "red", "blue" ), bty = "n", cex = 1,
        legend = c(
          expression( paste( "empirical cdf ", hat(G) ) ),
          expression( paste( "mean predictive cdf ", bar(F) ) ),
          expression( bar(F)-hat(G) )
        )
      )
    },
    bs  ={
      
      # plot of brier score vs. cutoff
      
      if( missing( main ) ) main <- "Brier score vs. cutoff"
      if( missing( xlab ) ) xlab <- "cutoff"
      if( missing( ylab ) ) ylab <- "Brier score"
      
      if( add ){
        lines( result[["y"]], result[["bs"]], col = col, lty = lty, ... )
      } else {
        plot( result[["y"]], result[["bs"]], type = "l", col = col, lty = lty, 
          main = main, xlab = xlab, ylab = ylab, ...
        )
      }
    }
  )
  invisible( NULL )
}

  

##  ###########################################################################

print.cv.georob <-
  function( 
    x, digits = max(3, getOption("digits") - 3), ...
  )
{   ## print method for class "cv.georob"  
  
  ## 2011-10-13 A. Papritz
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  x <- x[["pred"]]
  
  st <- validate.predictions( 
    data = x[["data"]],
    pred = x[["pred"]],
    se.pred = x[["se"]],
    statistic = "st",
    ...
  )
  
  print(
    format( st, digits = digits ), print.gap = 2, 
    quote = FALSE
  )
  
  invisible( x )
}


##  ###########################################################################

summary.cv.georob <-
  function( object, se = FALSE, ... )
{
  
  ## summary method for class "cv.georob" 
  
  ## function computes statistics of the cross-validation errors
  
  ## 2011-10-13 A. Papritz
  ## 2012-05-21 ap
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  
  object <- object[["pred"]]
  
  bs <- validate.predictions( 
    data = object[["data"]],
    pred = object[["pred"]],
    se.pred = object[["se"]],
    ncutoff = length( object[["data"]] ),
    statistic = "bs"
  )
  
  t.d <- diff( bs[["y"]] )
  crps <- sum( bs[["bs"]] * 0.5 * ( c( 0., t.d ) + c( t.d, 0. ) ) )
  
  st <- validate.predictions( 
    data = object[["data"]],
    pred = object[["pred"]],
    se.pred = object[["se"]],
    statistic = "st"
  )
  
  if( !is.null( object[["lgn.pred"]] ) ){
    st.lgn <- validate.predictions( 
      data = object[["lgn.data"]],
      pred = object[["lgn.pred"]],
      se.pred = object[["lgn.se"]],
      statistic = "st"
    )
  } else {
    st.lgn <- NULL
  }
  
  ## collect results
  
  result <- list( st = st, crps = crps, st.lgn = st.lgn )

  ## compute standard errors of criteria across cross-validation sets
  
  if( se && !is.null( object[["subset"]] ) ){
    
    criteria <- t( sapply(
        tapply(
          1:nrow( object ),
          factor( object[["subset"]] ),
          function( i, data, pred, se.pred, lgn.data, lgn.pred, lgn.se.pred ){
            
            bs <- validate.predictions( 
              data = data[i],
              pred = pred[i],
              se.pred = se.pred[i],
              ncutoff = length( data[i] ),
              statistic = "bs"
            )
            
            t.d <- diff( bs[["y"]] )
            crps <- c( crps = sum( bs[["bs"]] * 0.5 * ( c( 0., t.d ) + c( t.d, 0. ) ) ) )
            
            st <- validate.predictions( 
              data = data[i],
              pred = pred[i],
              se.pred = se.pred[i],
              statistic = "st"
            )
            
            if( !is.null( lgn.pred ) ){
              st.lgn <- validate.predictions( 
                data = lgn.data[i],
                pred = lgn.pred[i],
                se.pred = lgn.se.pred[i],
                statistic = "st"
              )
              names( st.lgn ) <- paste( names( st.lgn ), "lgn", sep = "." )
            } else {
              st.lgn <- NULL
            }
            
            
            return( c( st, st.lgn, crps ) )
            
          },
          data = object[["data"]],
          pred = object[["pred"]],
          se.pred = object[["se"]],
          lgn.data = object[["lgn.data"]],
          lgn.pred = object[["lgn.pred"]],
          lgn.se.pred = object[["lgn.se"]]
        ),
        function( x ) x
      ))
    
    se.criteria <- apply(
      criteria, 2,
      function( x ) sd( x ) / sqrt( length( x ) )
    )
    
    result[["se.st"]] <- se.criteria[c( "me", "mede", "rmse", "made", "qne", "msse", "medsse")]
    result[["se.crps"]] <- se.criteria["crps"]
    if( !is.null( st.lgn ) ){
      result[["se.st.lgn"]] <- se.criteria[
        c( "me.lgn", "mede.lgn", "rmse.lgn", "made.lgn", "qne.lgn", "msse.lgn", "medsse.lgn")
      ]
      names( result[["se.st.lgn"]] ) <- gsub( ".lgn", "", names( result[["se.st.lgn"]] ) )
    }
    
  }
  
  class( result ) <- "summary.cv.georob"
 
  
  return( result )
  
}

##  ###########################################################################

print.summary.cv.georob <-
  function( 
    x, digits = max(3, getOption("digits") - 3), ...
  )
{
  
  ## print method for class "summary.cv.georob"  
  
  ## 2011-12-20 A. Papritz
  ## 2012-05-21 ap
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  

  result <- c( x[["st"]], crps = x[["crps"]] )
  if( !is.null( x[["se.st"]] ) ){
    result <- rbind( result, c( x[["se.st"]], crps = x[["se.crps"]] ) )
    rownames( result ) <- c( "", "se" )
  }
  
  cat( "\nStatistics of cross-validation prediction errors\n" )
  print(
    format( result, digits = digits ), print.gap = 2, 
    quote = FALSE
  )
  
  if( !is.null( x[["st.lgn"]] ) ){
    result <- x[["st.lgn"]]
    if( !is.null( x[["se.st.lgn"]] ) ){
      result <- rbind( x[["st.lgn"]], x[["se.st.lgn"]] )
      rownames( result ) <- c( "", "se" )
    }
    
    cat( "\nStatistics of back-transformed cross-validation prediction errors\n" )
    print(
      format( result, digits = digits ), print.gap = 2, 
      quote = FALSE
    )
  }
  
  invisible( x )
  
}

##  ###########################################################################

rstudent.cv.georob <-
  function( model, ... )
{
  
  ## Function extracts studentized residuals from cv.georob object
  
  ## Arguments:
  
  ## model     cv.georob object
  ## ...       further arguments (currently not used)
  
  ## 2011-10-13 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  if( !identical( class( model )[1], "cv.georob" ) ) stop(
    "model is not of class 'cv.georob'" 
  )
  
  model <- model[["pred"]]
  
  ( model[["data"]] - model[["pred"]] ) / model[["se"]]
  
}

# ##  ###########################################################################
# 
# cv.variomodel <-
#   function( object, geodata, ... )
# {
#   
#   ## Wrapper function for cross-validation of object of class variomodel{geoR}
#   ## by function xvalid{geoR}
#   
#   ## Arguments:
#   
#   ## model     an object of class "variomodel{geoR}
#   ## ...       further arguements passed to xvalid{geoR), cf. respective help page
#   
#   ## 2012-11-22 A. Papritz
#   
#   call.fc <- match.call()
# 
#   res <- geoR::xvalid( model = object, ... )
#   
#   if( !is.null( attr( res, "geodata.xvalid" ) ) ){
#     attr( res, "geodata.xvalid" ) <- call.fc[["geodata"]]
#   }
#   if( !is.null( attr( res, "locations.xvalid" ) ) ){
#     attr( res, "locations.xvalid" ) <- call.fc[["locations.xvalid"]]
#   }
#   
#   return(res)
#   
# }
# 
# cv.likGRF <-
#   function( object, geodata, ... )
# {
#   
#   ## Wrapper function for cross-validation of object of class variomodel{geoR}
#   ## by function xvalid{geoR}
#   
#   ## Arguments:
#   
#   ## model     an object of class "likGRF{geoR}
#   ## ...       further arguements passed to xvalid{geoR), cf. respective help page
#   
#   ## 2012-11-22 A. Papritz
#   
#   call.fc <- match.call()
# 
#   res <- geoR::xvalid( model = object, geodata = geodata, ... )
#   
#   if( !is.null( attr( res, "geodata.xvalid" ) ) ){
#     attr( res, "geodata.xvalid" ) <- call.fc[["geodata"]]
#   }
#   if( !is.null( attr( res, "locations.xvalid" ) ) ){
#     attr( res, "locations.xvalid" ) <- call.fc[["locations.xvalid"]]
#   }
#   
#   return(res)
#   
# }

## ======================================================================
validate.predictions <- 
  function( 
    data,
    pred,
    se.pred,
    statistic = c( "pit", "mc", "bs", "st" ),	
    ncutoff = NULL
  )
{
  
  ## function computes several statistics to validate probabilistic
  ## predictions, cf.  Gneiting et al., 2007, JRSSB
  
  # 2011-20-21 A. Papritz
  # 2012-05-04 AP coping with NAs
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  statistic = match.arg( statistic )
  
  ## exclude item with NAs
  
  t.sel <- complete.cases( data, pred, se.pred )
  
  if( sum( t.sel ) < length( t.sel ) ) warnings(
    "missing values encountered when validating predictions"
  )
  
  data    <- data[t.sel]
  pred    <- pred[t.sel]
  se.pred <- se.pred[t.sel]
  
  if( missing( ncutoff ) || is.null( ncutoff ) ) ncutoff <- min( 500, length( data ) )
  
  result <- switch(
    statistic,
    pit = {
      
      ## probability integral transformation
      
      pnorm( data, mean = pred, sd = se.pred )
    },
    mc = ,
    bs = {
      
      ## marginal calibration and brier score
      
      margin.calib <- data.frame(
        y = t.x <- unique( t.y <- sort( c( data ) ) ),
        ghat = cumsum( tabulate( match( t.y, t.x ) ) ) / length(t.y)
      )
      t.sel <- trunc( 
        seq( 
          from = as.integer(1), 
          to = nrow( margin.calib ), 
          length.out = min( nrow( margin.calib ), ncutoff ) 
        ) 
      )
      margin.calib <- margin.calib[t.sel,]
      
      t.bla <- t(
        sapply(
          margin.calib[["y"]],
          function( q, m, s, y ){
            t.p <- pnorm( q, mean = m, sd = s )
            c( 
              fbar = mean( t.p ), 
              bs = mean( ( t.p - as.numeric( y <= q ) )^2 )
            )
          },
          m = pred,
          s = se.pred,
          y = data
        )
      )
      cbind(
        margin.calib, as.data.frame( t.bla ) 
      )
    },
    st = {
      
      ## statistics of (standardized) prediction errors
      
      error <- data - pred
      std.error <- error / se.pred
      
      statistics <- c( 
        me = mean( error ),
        mede = median( error ),
        rmse = sqrt( mean( error^2 ) ),
        made = mad( error, center = 0 ),
        qne = Qn( error, finite.corr = TRUE ),
        msse = mean( std.error^2 ),
        medsse = median( std.error^2 )
      )
      
    }
  )
  
  return( result )   
  
}

