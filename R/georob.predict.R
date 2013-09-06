##  ###########################################################################

predict.georob <-
function( 
  object, newdata, 
  type = c( "signal", "response", "trend", "terms" ),
  terms = NULL, se.fit = TRUE,
  signif = 0.95,
  mmax = 10000,
  locations,
  full.covmat = FALSE,
  pwidth = NULL, pheight = NULL, napp = 1,
  extended.output = FALSE,
  ncores = detectCores(),
  verbose = 0, 
  ...
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
  
  ## Arguments:
  
  ## object            object of class inheriting from "georob"
  ## newdata           an optional data frame in which to look for variables 
  ##                   for point kriging
  ##                   or optionally an object of class SpatialPointsDataFrame{sp} 
  ##                   for block kriging
  ##                   with which to predict. If omitted, the fitted values are used
  ## type              character string indicating what quantity should be predicted,
  ##                   possible values are:
  ##                   "trend"         trend, i.e., X %*% beta
  ##                   "terms"         terms of trend, i.e. X[,terms] %*% beta[terms]
  ##                   "signal"        trend plus autocorrelated component, 
  ##                                   i.e., X %*% beta + Z
  ##                   "response"   signal component plus independent error, 
  ##                                   i.e. X %*% beta + Z + epsilon
  ## terms             optional character vector, only used if type == "terms",
  ##                   see predict.lm{stats}
  ## se.fit            logical, only used if type == "terms", see predict.lm{stats}
  ## signif             tolerance/confidence signif
  ## mmax              maximum number of prediction items, computed in a single 
  ##                   computation step, ceiling( m / mmax ) such step are executed if
  ##                   m items are predicted
  ## locations         a one-sided formula specifying what variables of newdata are the
  ##                   coordinates of the prediction points (applies only if newdata is 
  ##                   data frame)
  ## full.covmat       logical, indicating whether the full covariance matrix of the
  ##                   prediction errors should be computed
  ## pwidth, pheight,  see preCKrige{constrainedKriging}
  ## napp              see preCKrige{constrainedKriging}
  ## extended.output   logical, flag controlling whether the covariance matrix of the
  ##                   kriging predictions and the covariance matrix of the kriging predictions 
  ##                   and the data should be computed
  ## ncores            integer scalar with the number of cores to used in parallel processing
   ##           
  
  
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

  
  ##  ##############################################################################
  
  ## auxiliary function for computinge robust kriging predictions for a part
  ## of prediction targets
  
  f.robust.uk <-
    function(
      type, terms,
      locations.coords, betahat, bhat,
      pred.X, pred.coords, offset, newdata,
      variogram.model, param, aniso,
      cov.delta.bhat.betahat.l, cov.betahat.l, cov.bhat.betahat, cov.p.t, Valpha.objects,
      pwidth, pheight, napp,
      signif,
      extended.output, full.covmat
    )
  { ## f.robust.uk
    
    ## function computes robust universal point or block kriging
    ## predictions 
    
    ## 2011-07-29 A. Papritz
    ## 2012-05-04 AP modifications for lognormal block kriging
    
    n <- length( bhat )
    
    ## exclude prediction items with missing information
    
    ex <- attr( na.omit( pred.X ), "na.action" )
    
    if( !is.null( pred.coords ) ){
      ex <- unique( c( ex, attr( na.omit( pred.coords ), "na.action" ) ) )
    }
    
    if( !is.null( ex ) ) {
      ex <- ( 1:NROW(pred.X) ) %in% sort( ex )
    } else {
      ex <- rep( FALSE, NROW(pred.X) )
    }
    
    ## compute trend surface prediction
    
    if( any( !ex ) ){
      
      t.pred <- t.trend <- drop( pred.X[!ex, , drop = FALSE ] %*% betahat ) + offset[!ex]
      
      if( !identical( type, "trend" ) ){
        
        ## compute point or block kriging predictions
        
        ## get covariance matrix (cov.target) of z at predictons locations and
        ## covariance matrix (gamma) between z at prediction and support
        ## locations
        
        
        sill <- Valpha.objects[["gcr.constant"]] * sum( param[c("variance", "snugget")] )
        
        if( !is.null( pred.coords ) ){ 
          
          ## point kriging
          
          ## matrix for coordinate transformation
          
          A <- aniso[["sclmat"]] * aniso[["rotmat"]] / param["scale"]
          
          ## variogram model
          
          model.list <- list( variogram.model )
          model.list <- c( model.list, as.list( param[-(1:4)] ) )
          
          ## generalized (co-)variance (matrix) of prediction
          ## points
          
          if( full.covmat && NROW( pred.coords[!ex, , drop = FALSE ] ) > 1 ){
            
            ## lag vectors for all distinct pairs
            
            indices.pairs <- combn( NROW( pred.coords[!ex, , drop = FALSE ] ), 2 )
            lag.vectors <- (pred.coords[!ex, , drop = FALSE ])[ indices.pairs[2,], ] - 
            (pred.coords[!ex, , drop = FALSE ])[ indices.pairs[1,], ]
            
            ##  negative semivariance matrix
            
            ## functions of version 3 of RandomFields
  
            RFoptions(newAniso=FALSE)
            
            t.var.target <- try(
              -RFvariogram(
                lag.vectors,
                model = list( "+",
                  list( "$", var = param["variance"], A = A, model.list ),
                  list( "$", var = param["snugget"], list( "nugget" ) )
                ),
                dim = NCOL( lag.vectors ), grid = FALSE
              ),
              silent = TRUE
            )

            ## functions of version 2 of RandomFields
  
            ##             RFoldstyle()
            ##             RFoptions(newAniso=FALSE)
            ##             t.var.target <- try(
            ##               -Variogram(
            ##                 lag.vectors,
            ##                 model = list( "+",
            ##                   list( "$", var = param["variance"], A = A, model.list ),
            ##                   list( "$", var = param["snugget"], list( "nugget" ) )
            ##                 )
            ##               ),
            ##               silent = TRUE
            ##             )
            
            if( 
              identical( class( t.var.target ), "try-error" ) || 
              any( is.na( t.var.target ) ) 
            ) stop(
              "an error occurred when computing semivariances between ",
              "prediction locations"
            )
            
            t.var.target <- sill + t.var.target
            
            ## convert to symmetric matrix
            
            t.var.target <- list(
              diag = rep( sill, 0.5 * ( 1 + sqrt( 1 + 8 * length( t.var.target ) ) ) ),
              tri = t.var.target
            )
            attr( t.var.target, "struc" ) <- "sym"
            
            t.var.target <- expand( t.var.target )
            
            #           pred.dist <- as.matrix( dist( pred.coords[!ex, ] ) )
            #           
            #           alt <- try( 
            #             Variogram(
            #               c( pred.dist ), model = variogram.model, 
            #               param = c( NA, param["variance"], param["snugget"], param[-(1:3)] ), 
            #               dim = NCOL( pred.coords )
            #             )
            #           )
            #           if( 
            #             identical( class( alt ), "try-error" ) || 
            #             any( is.na( alt ) ) 
            #           ) stop(
            #             "error when computing semivariances between ",
            #             "prediction locations"
            #           )
            #           alt <- sill - alt
            #           dim( alt ) <- dim( pred.dist )
            
          } else {
            
            t.var.target <- sill
            
          }
          
          ## generalized covariance matrix between prediction and
          ## support points
          
          
          indices.pairs <- expand.grid( 
            1:NROW( pred.coords[!ex, , drop = FALSE ] ),
            1:NROW( locations.coords )
          )
          lag.vectors <- locations.coords[ indices.pairs[, 2], ] - 
          (pred.coords[!ex, , drop = FALSE ])[ indices.pairs[, 1], ]
          
          ## functions of version 3 of RandomFields
          
          RFoptions(newAniso=FALSE)
          
          gamma <- try(
            -RFvariogram(
              lag.vectors,
              model = list( "+",
                list( "$", var = param["variance"], A = A, model.list ),
                list( "$", var = param["snugget"], list( "nugget" ) )
              ), 
              dim = NCOL( lag.vectors ), grid = FALSE
            ),
            silent = TRUE
          )
          ##           RFoldstyle()
          ##           RFoptions(newAniso=FALSE)
          ##           gamma <- try(
          ##             -Variogram(
          ##               lag.vectors,
          ##               model = list( "+",
          ##                 list( "$", var = param["variance"], A = A, model.list ),
          ##                 list( "$", var = param["snugget"], list( "nugget" ) )
          ##               )
          ##             ),
          ##             silent = TRUE
          ##           )
          
          if( 
            identical( class( gamma ), "try-error" ) || 
            any( is.na( gamma ) ) 
          ) stop(
            "an error occurred when computing semivariances between support ",
            "and prediction points"
          )
          
          gamma <- sill + gamma
          
          gamma <- matrix( gamma, nrow = NROW( pred.coords[!ex, , drop = FALSE ] ) )
          
          #         xdist <- fields::rdist( pred.coords[!ex, , drop = FALSE ], locations.coords )
          #         
          #         alt <- try( 
          #           Variogram(
          #             c( xdist ), model = variogram.model, 
          #             param = c( NA, param["variance"], param["snugget"], param[-(1:3)] ), 
          #             dim = NCOL( pred.coords )
          #           )
          #         )
          #         
          #         if( 
          #           identical( class( alt ), "try-error" ) || 
          #           any( is.na( alt ) ) 
          #         ) stop(
          #           "error when computing semivariances between support ",
          #           "and prediction points"
          #         )
          #         
          #         alt <- sill - alt
          #         dim( alt ) <- dim( xdist ) ## end of point kriging
          
        } else {
          
          ## block kriging
          
          ## setup covariance model list
          
          t.covmodel <- covmodel(
            modelname = variogram.model,
            mev = switch(
              type,
              response = 0,
              signal = param["nugget"]
            ),
            nugget = switch(
              type,
              response = sum( param[c("snugget", "nugget")] ),
              signal = param["snugget"]
            ),
            variance = param["variance"],
            scale = param["scale"]
          )
          if( length( param ) > 4 ) t.covmodel[[3]][["parameter"]] <- param[-(1:4)]
          
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
            
            t.neighbours <- list()
            for( i in 1:length( newdata[!ex, , drop = FALSE ] ) ){
              t.neighbours[[i]] <- integer(0)
            }
            t.neighbours[[1]] <- 2:length( newdata[!ex, , drop = FALSE ] )
            
            t.preCKrige.aux <- preCKrige(
              newdata = newdata[!ex, , drop = FALSE ], 
              neighbours = t.neighbours,
              model = t.covmodel,
              pwidth = pwidth, pheight = pheight, napp = napp
            )
            
            #           cat( "\n\n!!!!!!BLOCK-BLOCK KOVARIANZMATRIX WIRD EINGELESEN!!!!!!\n\n")
            #           
            #           load( "r.preCK_3" )
            #           t.preCKrige.aux <- r.preCK_3
            
            t.se <- sqrt( t.var.target )
            t.var.target <- t.se * t( t.se *
              cov2cor( t.preCKrige.aux@covmat[[1]] ) )
            
          } 
          
          ## get rid of mev component in covariance model list
          
          t.covmodel <- t.preCKrige@model[
          unlist( 
            lapply( 
              1:length(t.preCKrige@model), 
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
        
        ## compute uk predictions
        
        gammaValphai <- gamma %*% Valpha.objects[["Valpha.inverse"]] / sum( param[c("variance", "snugget")] )
        t.pred <- t.pred + drop( gammaValphai %*% bhat )
        
        ## compute uk variance (= (co-)variance of prediction errors)
        
        aux <- cbind(
          gammaValphai %*% cov.delta.bhat.betahat.l[1:n, 1:n] - pred.X[!ex, , drop = FALSE ] %*% cov.delta.bhat.betahat.l[-(1:n), 1:n],
          - pred.X[!ex, , drop = FALSE ] %*% cov.delta.bhat.betahat.l[-(1:n), -(1:n)]
        )
        
        if( full.covmat ){
          t.mse.pred <- tcrossprod( aux ) + t.var.target - gammaValphai %*% t( gamma )
        } else {
          t.mse.pred <- rowSums( aux^2 ) + t.var.target - rowSums( gammaValphai * gamma )
        }
        
        if( identical( type, "response" ) && !is.null( pred.coords ) ){
          
          if( full.covmat ){
            diag( t.mse.pred ) <- diag( t.mse.pred ) + unname( param["nugget"] )
            diag( t.var.target ) <- diag( t.var.target ) + unname( param["nugget"] )
          } else {  
            t.mse.pred <- t.mse.pred + unname( param["nugget"] )
            t.var.target <- t.var.target + unname( param["nugget"] )
          }
          
        }
        
        if( extended.output ){
          
          ## compute covariance matrix of uk predictions and
          ## covariance matrix of uk predictions and prediction targets
          ## (needed for lognormal kriging)
          
          aux <- cbind( gammaValphai, pred.X[!ex, , drop = FALSE ] )
          t.var.pred <- aux %*% cov.bhat.betahat %*% t( aux )
          t.cov.pred.target <- aux %*% cov.p.t %*% t( gamma )
          if( !full.covmat ){
            t.var.pred <- diag( t.var.pred )
            t.cov.pred.target <- diag( t.cov.pred.target )
          }
          
        }
        
        ## end compute kriging predictions
        
      } else {
        
        ## compute variance of trend surface prediction
        
        if( full.covmat ){
          t.var.pred <- tcrossprod( pred.X[!ex, , drop = FALSE ] %*% cov.betahat.l )
          t.mse.pred <- matrix( NA_real_, nrow = nrow( t.var.pred ), ncol = ncol( t.var.pred ) )
          t.var.target <- matrix( 0., nrow = nrow( t.var.pred ), ncol = ncol( t.var.pred ) )
          t.cov.pred.target <- matrix( 0., nrow = nrow( t.var.pred ), ncol = ncol( t.var.pred ) )
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
    
    ## add items with missing information back
    
    sr <- (1:NROW(pred.X))[!ex]
    
    pred <- rep( NA_real_, NROW(pred.X) )
    pred[!ex] <- t.pred
    if( extended.output ){
      trend <- rep( NA_real_, NROW(pred.X) )
      if( any( !ex ) ) trend[!ex] <- t.trend            
    }
    
    if( full.covmat ){
      
      mse.pred <- matrix( rep( NA_real_, NROW(pred.X)^2 ), nrow = NROW(pred.X) )
      mse.pred[sr, sr ] <- t.mse.pred
      
      if( identical( type, "trend" ) || extended.output ){
        var.pred <- matrix( rep( NA_real_, NROW(pred.X)^2 ), nrow = NROW(pred.X) )
        var.pred[sr, sr ] <- t.var.pred
      }
      if( extended.output ){
        cov.pred.target <- matrix( rep( NA_real_, NROW(pred.X)^2 ), nrow = NROW(pred.X) )
        cov.pred.target[sr, sr ] <- t.cov.pred.target
        var.target <- matrix( rep( NA_real_, NROW(pred.X)^2 ), nrow = NROW(pred.X) )
        var.target[sr, sr ] <- t.var.target
      }
      
    } else {
      
      mse.pred <- rep( NA_real_, NROW(pred.X) )
      mse.pred[sr] <- t.mse.pred
      if( identical( type, "trend" ) || extended.output ){
        var.pred <- rep( NA_real_, NROW(pred.X) )
        var.pred[sr] <- t.var.pred
      }
      if( extended.output ){
        cov.pred.target <- rep( NA_real_, NROW(pred.X) )
        cov.pred.target[sr] <- t.cov.pred.target
        var.target <- rep( NA_real_, NROW(pred.X) )
        var.target[sr] <- t.var.target
      }
      
    }
    
    ## collect results
    
    pred.se <-  sqrt( 
      if( full.covmat ){ 
        diag( mse.pred )
      } else {
        mse.pred
      }
    )
    
    result <- data.frame( 
      pred = pred,
      se = pred.se,
      lower = pred + qnorm( 0.5 * ( 1-signif[1] ) ) * pred.se,
      upper = pred + qnorm( 0.5 * ( 1+signif[1] ) ) * pred.se
    )
    
    if( extended.output ) result[["trend"]] <- trend
    
    if( identical( type, "trend" ) || extended.output ){
      if( full.covmat ){
        result[["var.pred"]] <- diag( var.pred )
        if( extended.output ){
          result[["cov.pred.target"]] <- diag( cov.pred.target )
          result[["var.target"]]      <- diag( var.target )
        }  
      } else {
        result[["var.pred"]] <- var.pred
        if( extended.output ){
          result[["cov.pred.target"]] <- cov.pred.target
          result[["var.target"]]      <- var.target
        }
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
  
  ## begin of body of main function
  
  ## expand matrices
  
  object[["Valpha.objects"]] <- expand( object[["Valpha.objects"]] )
  object[["cov"]]            <- expand( object[["cov"]]  )
  
  if( missing( locations ) ){
    locations <- object[["locations.objects"]][["locations"]]
  } else {
    locations <- as.formula( 
      paste( deparse( locations ), "-1" ), env = parent.frame() 
    )  
  }
  
  ## check the consistency of the provided arguments
  
  if( !missing( newdata ) && class( newdata ) == "SpatialPolygonsDataFrame" ){
    
    if( is.null( pwidth ) || is.null( pheight ) )
    stop( 
      "'pwidth' and 'pheight' must be provided for block kriging"
    )
    
    if( object[["variogram.model"]] %in% georob.control()[["irf.models"]] )
    stop(
      "block kriging not yet implemented for unbounded variogram models"
    )
    
    if( !object[["aniso"]][["isotropic"]]) stop(
      "block kriging not yet implemented for anisotropic variograms"          
    )
    
    
  }
  
  if( full.covmat ){
    if( verbose > 0 ){
      cat(
        "\ncomputing full covariance matrix of prediction errors\n" 
      )
      if( 
        !missing( newdata ) && class( newdata ) == "SpatialPolygonsDataFrame" 
      ) cat( 
        "requires some computing time for block kriging, be patient ...\n"
      )
    }
  }
  
  if( identical( type, "terms" ) && 
    !( missing( newdata ) || is.null( newdata )) ) stop(
    "predicting terms for newdata not yet implemented"        
  )
  
  type <- match.arg( type )
  
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
  
  ## if needed compute missing covariance matrices
  
  cov.betahat    <- is.null( object[["cov"]][["cov.betahat"]] )
  cov.delta.bhat   <- is.null( object[["cov"]][["cov.delta.bhat"]] ) ||
    !is.matrix( object[["cov"]][["cov.delta.bhat"]] )
  cov.delta.bhat.betahat <- is.null( object[["cov"]][["cov.delta.bhat.betahat"]] )
  cov.bhat    <- extended.output & (
    is.null( object[["cov"]][["cov.bhat"]] ) || !is.matrix( object[["cov"]][["cov.bhat"]] )
  )
  cov.bhat.betahat  <-  extended.output & is.null( object[["cov"]][["cov.bhat.betahat"]] )
  cov.p.t  <-  extended.output & is.null( object[["cov"]][["cov.pred.target"]] )
  
  if( any( c( cov.betahat, cov.delta.bhat, cov.delta.bhat.betahat, 
        extended.output & ( cov.bhat || cov.bhat.betahat || cov.p.t )
      )
    )
  ){ ## cov
    
    r.cov <- compute.covariances(
      Valpha.objects = object[["Valpha.objects"]],
      Aalpha = object[["Aalpha"]], Palpha = object[["Palpha"]],
      rweights = object[["rweights"]],
      XX = X, TT = object[["Tmat"]], names.yy = rownames( model.frame( object ) ),
      nugget = object[["param"]]["nugget"],
      eta = sum( object[["param"]][c( "variance", "snugget")] ) / object[["param"]]["nugget"],
      expectations = object[["expectations"]],
      cov.bhat = cov.bhat, full.cov.bhat = cov.bhat,
      cov.betahat = cov.betahat, 
      cov.bhat.betahat = cov.bhat.betahat,
      cov.delta.bhat = cov.delta.bhat, full.cov.delta.bhat = cov.delta.bhat,
      cov.delta.bhat.betahat = cov.delta.bhat.betahat,
      cov.ehat = FALSE, full.cov.ehat = FALSE,
      cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
      aux.cov.pred.target = cov.p.t,
      verbose = verbose
    )
    
    if( r.cov[["error"]] ) stop( 
      "an error occurred when computing the covariances of fixed and random effects",
    )
    
    if( is.null( object[["cov"]] ) ) object[["cov"]] <- list()
    
    if( cov.betahat )    object[["cov"]][["cov.betahat"]]    <- r.cov[["cov.betahat"]]
    if( cov.delta.bhat )   object[["cov"]][["cov.delta.bhat"]] <- r.cov[["cov.delta.bhat"]]
    if( cov.delta.bhat.betahat ) object[["cov"]][["cov.delta.bhat.betahat"]] <- 
      r.cov[["cov.delta.bhat.betahat"]]
    if( extended.output && cov.bhat )   object[["cov"]][["cov.bhat"]] <- r.cov[["cov.bhat"]]
    if( extended.output && cov.bhat.betahat ) object[["cov"]][["cov.bhat.betahat"]] <- 
      r.cov[["cov.bhat.betahat"]]
    if( extended.output && cov.p.t ) object[["cov"]][["cov.pred.target"]] <- 
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
  
  if( extended.output ){
    
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
  
  ## compute predictions
  
  if( missing( newdata ) || is.null( newdata ) ){ 
    
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
      
      if( nterms > 0 ){
        
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
        predictor <- ip <- matrix(0, length( object[["Tmat"]] ), 0L)
      }
      
      attr( predictor, "constant" ) <- 
      if( hasintercept )  termsconst  else  0
      
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
          residual.scale = unname( sqrt( object[["param"]]["nugget"] ) )
        )
      } else {
        predictor
      }
      ## end "terms"
      
    } else {  
      
      ## no newdata object: compute predictions for support locations
            
      pred <- switch(
        type,
        "response" = model.response( model.frame( object ) ),
        "signal" = object[["fitted.values"]] + object[["bhat"]][object[["Tmat"]]],
        "trend" = object[["fitted.values"]] 
      )
      
      var.pred <- NULL
      var.target <- NULL
      cov.pred.target <- NULL
      
      if( extended.output ){
        V <- sum( object[["param"]][c("variance", "snugget")] ) * object[["Valpha.objects"]][["Valpha"]]
      }
      
      t.result <- switch(
        type,
        response = {  ## response
          mse.pred <- rep( 0., length( object[["Tmat"]] ) )
          if( extended.output ) var.pred <- var.target <- cov.pred.target <- rep( 
            sum( object[["param"]][c("variance", "nugget", "snugget")] ), 
            length( object[["Tmat"]] ) 
          )
          if( full.covmat ){
            mse.pred <- diag( mse.pred )
            if( extended.output ){
              var.pred <- V
              diag( var.pred ) <- diag( var.pred ) + object[["param"]][c("nugget")]
              var.pred <- var.target <- cov.pred.target <- var.pred[object[["Tmat"]], object[["Tmat"]]]
            }
          }
          c( 
            list( mse.pred = mse.pred ),
            if( extended.output ) list( 
              var.pred = var.pred, var.target = var.target, cov.pred.target = cov.pred.target
            )
          )
        },
        signal = {    ## signal
          aux <- cbind( 
            cov.delta.bhat.betahat.l[1:n,1:n] - X %*% cov.delta.bhat.betahat.l[-(1:n),1:n],
            - X  %*% cov.delta.bhat.betahat.l[-(1:n),-(1:n)]
          )
          aux <- aux[object[["Tmat"]], , drop = FALSE]
          if( full.covmat ){
            mse.pred <- tcrossprod( aux )
          } else {
            mse.pred <- rowSums( aux^2 )
          }
          if( extended.output ){
            aux <- cov.bhat.betahat[1:n, -(1:n), drop = FALSE] %*% t(X)
            var.pred <- cov.bhat.betahat[1:n, 1:n, drop = FALSE] + aux + t(aux) + 
              X %*% cov.bhat.betahat[-(1:n), -(1:n), drop = FALSE] %*% t(X)
            var.target <- V[object[["Tmat"]], object[["Tmat"]]]
            cov.pred.target <- (cov.p.t[1:n,] + X %*% cov.p.t[-(1:n),]) %*% V
            cov.pred.target <- cov.pred.target[object[["Tmat"]], object[["Tmat"]]]
            if( !full.covmat ){
              var.pred <- diag( var.pred )
              var.target <- diag( var.target )
              cov.pred.target <- diag( cov.pred.target )
            }
          }
          c( 
            list( mse.pred = mse.pred ),
            if( extended.output ) list( 
              var.pred = var.pred, var.target = var.target, cov.pred.target = cov.pred.target
            )
          )
        },
        trend = {     ## trend
          aux <- X %*% cov.betahat.l
          aux <- aux[object[["Tmat"]], , drop = FALSE]
          if( full.covmat ){
            mse.pred <- matrix( NA_real_, length( object[["Tmat"]] ), length( object[["Tmat"]] ) )
            var.pred <- tcrossprod( aux )
            if( extended.output ){
              var.target <- matrix( 0., length( object[["Tmat"]] ), length( object[["Tmat"]] ) )
              cov.pred.target <- matrix( 0., length( object[["Tmat"]] ), length( object[["Tmat"]] ) )
            }
          } else {
            mse.pred <- rep( NA_real_, length( object[["Tmat"]] ) )
            var.pred <- rowSums( aux^2 )
            if( extended.output ){
              var.target <- rep( 0., length( object[["Tmat"]] ) )
              cov.pred.target <- rep( 0., length( object[["Tmat"]] ) )
            }
          }
          list( mse.pred = mse.pred, var.pred = var.pred )
        }
      )
      
      ## collect results
      
      pred.se <-  sqrt( 
        if( full.covmat ){ 
          diag( t.result[["mse.pred"]] )
        } else {
          t.result[["mse.pred"]]
        }
      )
      
      result <- data.frame( 
        pred = pred,
        se = pred.se,
        lower = pred + qnorm( 0.5 * ( 1-signif[1] ) ) * pred.se,
        upper = pred + qnorm( 0.5 * ( 1+signif[1] ) ) * pred.se
      )
      
      if( !is.null(t.result[["var.pred"]]) ){
        result[["var.pred"]] <- if( full.covmat ){
          diag( t.result[["var.pred"]] ) 
        } else {
          t.result[["var.pred"]] 
        }
      }
      if( !is.null(t.result[["var.target"]]) ){
        result[["var.target"]]<- if( full.covmat ){
          diag( t.result[["var.target"]] ) 
        } else {
          t.result[["var.target"]] 
        }
      }
      if( !is.null(t.result[["cov.pred.target"]]) ){
        result[["cov.pred.target"]]<- if( full.covmat ){
          diag( t.result[["cov.pred.target"]] ) 
        } else {
          t.result[["cov.pred.target"]] 
        }
      }
      
      result <- as.data.frame( napredict( object[["na.action"]], as.matrix( result ) ) )
            
      if( full.covmat ){
        
        result <- c(
          list( pred = result, 
            mse.pred = napredict( object[["na.action"]], t( napredict( object[["na.action"]], mse.pred ) ) ) 
          ), 
          if( !is.null(var.pred) ) list( 
            var.pred = napredict( object[["na.action"]], t( napredict( object[["na.action"]], var.pred ) ) ) 
          ),
          if( !is.null(var.target) ) list( 
            var.target = napredict( object[["na.action"]], t( napredict( object[["na.action"]], var.target ) ) )
          ),
          if( !is.null(cov.pred.target) ) list( 
            cov.pred.target = napredict( 
              object[["na.action"]], t( napredict( object[["na.action"]], cov.pred.target ) ) 
            )
          )
        )
      }
      
        
    }
    ## end no newdata object
    
  } else { 
    
    ## compute predictions for newdata
    
    Terms <- delete.response( tt )
    Terms.loc <- locations
    
    ## get the model frame for newdata
    
    mf.newdata <- switch( 
      class( newdata ),
      "data.frame" = model.frame( 
        Terms, newdata, na.action = na.pass, xlev = object[["xlevels"]] 
      ),
      "SpatialPointsDataFrame" = model.frame(
        Terms, slot( newdata, "data" ), na.action = na.pass, 
        xlev = object[["xlevels"]] 
      ),
      "SpatialPixelsDataFrame" = model.frame(
        Terms, slot( newdata, "data" ), na.action = na.pass, 
        xlev = object[["xlevels"]] 
      ),
      "SpatialGridDataFrame" = model.frame(
        Terms, slot( newdata, "data" ), na.action = na.pass, 
        xlev = object[["xlevels"]] 
      ),
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
    
    offset <- rep( 0, NROW(pred.X) )
    if( !is.null( off.num <- attr( tt, "offset" ) ) ){
      for( i in off.num ) {
        offset <- offset + eval( attr( tt, "variables" )[[i + 1]], newdata )
      }
    }
    if( !is.null( object[["call"]][["offset"]] ) ){
      offset <- offset + eval( object[["call"]][["offset"]], newdata )
    }
    
    ## get matrix of coordinates of newdata for point kriging
    
    pred.coords <- switch( 
      class( newdata ),
      "data.frame" = model.matrix(
        Terms.loc,
        model.frame( 
          Terms.loc, newdata, na.action = na.pass 
        )
      ),
      "SpatialPointsDataFrame" = model.matrix(
        Terms.loc,
        model.frame(  
          Terms.loc, as.data.frame( coordinates( newdata ) ),
          na.action = na.pass
        )
      ),
      "SpatialPixelsDataFrame" = model.matrix(
        Terms.loc,
        model.frame(  
          Terms.loc, as.data.frame( coordinates( newdata ) ),
          na.action = na.pass
        )
      ),
      "SpatialGridDataFrame" = model.matrix(
        Terms.loc,
        model.frame(  
          Terms.loc, as.data.frame( coordinates( newdata ) ),
          na.action = na.pass
        )
      ),
      "SpatialPolygonsDataFrame" = NULL
    )
    
    if( !is.null( pred.coords ) &&
      NCOL( locations.coords ) != NCOL( pred.coords ) 
    ) stop(
      "inconsistent number of coordinates in object and in newdata"
    )   
    
    ## number of items to predict
    
    m <- NROW( newdata )
    
    ## determine number of prediction parts
    
    n.part <- ceiling( m / mmax )
    rs <- ( 0:(n.part-1)) * mmax + 1
    re <- ( 1:(n.part  )) * mmax; re[n.part] <- m 
    
    ncores <- min( n.part, ncores )
    
    parallel <- ncores > 1
    
    if( full.covmat && n.part > 1 ) stop(
      "full covariance matrix of prediction errors cannot ",
      "be computed\n  if prediction problem is split into several parts\n",
      "-> increase 'mmax' to avoid splitting"
    )
    
    ## handle parallel processing
    
    ## auxiliary function to compute the predictions for one part
    
    f.aux <- function(
      i, 
      rs, re, n.part,
      type,
      locations.coords, betahat, bhat,
      pred.X, pred.coords, offset, newdata, 
      variogram.model, param, aniso,
      cov.delta.bhat.betahat.l, cov.betahat.l, cov.bhat.betahat, cov.p.t, Valpha.objects,
      pwidth, pheight, napp,
      signif,
      extended.output, full.covmat,
      verbose
    ){
      
      if( verbose > 0 )
      cat( "  predicting part ", i, " of ", n.part, "\n" )
      
      ## select the data for the current part
      
      pred.X <- pred.X[ rs[i]:re[i], , drop = FALSE]
      offset <- offset[ rs[i]:re[i] ]
      newdata <- newdata[ rs[i]:re[i], ]
      if( !is.null( pred.coords ) ) {
        pred.coords <- pred.coords[ rs[i]:re[i], , drop = FALSE]
      }
      
      ## compute the predictions for the current part
      
      result <- f.robust.uk(
        type = type, terms = terms,
        locations.coords = locations.coords, 
        betahat = betahat,
        bhat = bhat,
        pred.X = pred.X, pred.coords = pred.coords, offset = offset, newdata = newdata,
        variogram.model = variogram.model, param = param, aniso = aniso,
        cov.delta.bhat.betahat.l = cov.delta.bhat.betahat.l,
        cov.betahat.l = cov.betahat.l, 
        cov.bhat.betahat = cov.bhat.betahat,
        cov.p.t = cov.p.t,
        Valpha.objects = Valpha.objects,
        pwidth = pwidth, pheight = pheight, napp = napp,
        signif = signif,
        extended.output = extended.output,
        full.covmat = full.covmat
      )
      
      return( result )              
    }
    
    ## compute the predictions for all the parts 
    
    if( parallel ){
      
      if( .Platform[["OS.type"]] == "windows" ){
        
        ## create a SNOW cluster on windows OS
        
        cl <- makePSOCKcluster( ncores, outfile =  "" )
        
        ## export required items to workers
        
        junk <- clusterEvalQ( cl, require( georob, quietly = TRUE ) )
        
        t.result <- parLapply(
          cl,
          1:n.part,
          f.aux,
          rs = rs, re = re, n.part = n.part,
          type = type,
          locations.coords = locations.coords,
          betahat = object[["coefficients"]],
          bhat = object[["bhat"]],
          pred.X = pred.X, pred.coords = pred.coords, offset = offset, newdata = newdata,
          variogram.model = object[["variogram.model"]],
          param = object[["param"]],
          aniso = object[["aniso"]],
          cov.delta.bhat.betahat.l = cov.delta.bhat.betahat.l,
          cov.betahat.l = cov.betahat.l,
          cov.bhat.betahat = cov.bhat.betahat,
          cov.p.t = cov.p.t,
          Valpha.objects = object[["Valpha.objects"]],
          pwidth = pwidth, pheight = pheight, napp = napp,
          signif = signif,
          extended.output = extended.output, 
          full.covmat = full.covmat,
          verbose = verbose
        )
        
        stopCluster(cl)
        
      } else {
        
        ## fork child processes on non-windows OS
        
        t.result <- mclapply(
          1:n.part,
          f.aux,
          rs = rs, re = re, n.part = n.part,
          type = type,
          locations.coords = locations.coords,
          betahat = object[["coefficients"]],
          bhat = object[["bhat"]],
          pred.X = pred.X, pred.coords = pred.coords, offset = offset, newdata = newdata,
          variogram.model = object[["variogram.model"]],
          param = object[["param"]],
          aniso = object[["aniso"]],
          cov.delta.bhat.betahat.l = cov.delta.bhat.betahat.l,
          cov.betahat.l = cov.betahat.l,
          cov.bhat.betahat = cov.bhat.betahat,
          cov.p.t = cov.p.t,
          Valpha.objects = object[["Valpha.objects"]],
          pwidth = pwidth, pheight = pheight, napp = napp,
          signif = signif,
          extended.output = extended.output, 
          full.covmat = full.covmat,
          verbose = verbose,
          mc.cores = ncores
        )
        
      }
      
    } else {
      
      t.result <- lapply(
        1:n.part,
        f.aux,
        rs = rs, re = re, n.part = n.part,
        type = type,
        locations.coords = locations.coords,
        betahat = object[["coefficients"]],
        bhat = object[["bhat"]],
        pred.X = pred.X, pred.coords = pred.coords, offset = offset, newdata = newdata,
        variogram.model = object[["variogram.model"]],
        param = object[["param"]],
        aniso = object[["aniso"]],
        cov.delta.bhat.betahat.l = cov.delta.bhat.betahat.l,
        cov.betahat.l = cov.betahat.l,
        cov.bhat.betahat = cov.bhat.betahat,
        cov.p.t = cov.p.t,
        Valpha.objects = object[["Valpha.objects"]],
        pwidth = pwidth, pheight = pheight, napp = napp,
        signif = signif,
        extended.output = extended.output, 
        full.covmat = full.covmat,
        verbose = verbose
      )
      
    }
    
    ## collect results of the various parts into a single list
    
    result <- t.result[[1]]
    if( length( t.result ) > 1 ){
      for( i in 2:length( t.result ) ) {
        result <- rbind( result, t.result[[i]] )                
      } 
    }
    
    ## end compute predictions for newdata
    
  }
  
  ## complement kriging result with coordinate information
  
  if( missing( newdata ) || is.null( newdata ) ){
    
    coords <- napredict( object[["na.action"]], object[["locations.objects"]][["coordinates"]] )
    
    if( !identical( type, "terms" ) ){
      if( full.covmat ){
        result[["pred"]] <- data.frame( coords, result[["pred"]] )
      } else {
        result <- data.frame( coords, result )
      }
    }
    
  } else {
    
    result <- switch(
      class( newdata ),
      data.frame = {
        if( full.covmat ){
          result[["pred"]] <- data.frame( pred.coords, result[["pred"]] )
        } else {
          result <- data.frame( pred.coords, result )
        }
        result
      },
      SpatialPointsDataFrame = {
        if( full.covmat ){
          result[["pred"]] <- SpatialPointsDataFrame( 
            coords = coordinates( newdata ), 
            data = result[["pred"]] 
          )        
        } else {
          result <- SpatialPointsDataFrame( 
            coords = 
            coordinates( newdata ), 
            data = result 
          )        
        }
        result
      },
      SpatialPixelsDataFrame = {
        if( full.covmat ){
          result[["pred"]] <- SpatialPixelsDataFrame( 
            points = coordinates( newdata ), 
            data = result[["pred"]]
          )        
        } else {
          result <- SpatialPixelsDataFrame( 
            points = coordinates( newdata ), 
            data = result 
          )        
        }
        result
      },
      SpatialGridDataFrame = {
        aux <- newdata
        if( full.covmat ){
          aux@data <- result[["pred"]]
          result[["pred"]] <- aux
        } else {
          aux@data <- result
          result <- aux
        }
        result
      },
      SpatialPolygonsDataFrame = {
        if( full.covmat ){
          result[["pred"]] <- SpatialPolygonsDataFrame( 
            Sr = SpatialPolygons( newdata@polygons ), 
            data = result[["pred"]]
          )        
        } else {
          result <- SpatialPolygonsDataFrame( 
            Sr = SpatialPolygons( newdata@polygons ), 
            data = result
          )        
        }
        result
      }
    )
  }
  
  ## set attributes required for back-transformation by lgnpp
  
  if( !identical( type, "terms" ) ){
    if( full.covmat ){
      if( is.data.frame( result[["pred"]] ) ){
        attr( result[["pred"]], "variogram.model" )      <- object[["variogram.model"]]
        attr( result[["pred"]], "param" )                <- object[["param"]]
        attr( result[["pred"]], "type" )                 <- type
      } else {
        attr( result[["pred"]]@data, "variogram.model" ) <- object[["variogram.model"]]
        attr( result[["pred"]]@data, "param" )           <- object[["param"]]
        attr( result[["pred"]]@data, "type" )            <- type
        if( class( result[["pred"]] ) == "SpatialPolygonsDataFrame" ){
          attr( result[["pred"]]@data, "coefficients" )  <- object[["coefficients"]]
          attr( result[["pred"]]@data, "terms" )         <- object[["terms"]]
          attr( result[["pred"]]@data, "locations" )     <- object[["locations.objects"]][["locations"]]
        }
      }
    } else {
      if( is.data.frame( result ) ){
        attr( result, "variogram.model" )      <- object[["variogram.model"]]
        attr( result, "param" )                <- object[["param"]]
        attr( result, "type" )                 <- type
      } else {
        attr( result@data, "variogram.model" ) <- object[["variogram.model"]]
        attr( result@data, "param" )           <- object[["param"]]
        attr( result@data, "type" )            <- type
        if( class( result ) == "SpatialPolygonsDataFrame" ){
          attr( result@data, "coefficients" )  <- object[["coefficients"]]
          attr( result@data, "terms" )         <- object[["terms"]]
          attr( result@data, "locations" )     <- object[["locations.objects"]][["locations"]]
        }
      }
    }
  }

  invisible( result )
  
}


