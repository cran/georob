####################################
#                                  #
#   Hilfsfunktionen fuer georob    #
#                                  #
####################################

##  ##############################################################################

compute.covariances <- 
  function(
    Valpha.objects, 
    Aalpha, Palpha,
    rweights, XX, TT, names.yy,
    nugget, eta, 
    expectations,
    cov.bhat, full.cov.bhat,
    cov.betahat, 
    cov.bhat.betahat,
    cov.delta.bhat, full.cov.delta.bhat,
    cov.delta.bhat.betahat,
    cov.ehat, full.cov.ehat,
    cov.ehat.p.bhat, full.cov.ehat.p.bhat,
    aux.cov.pred.target,
    control.parallel,
    verbose
  )
{
  
  ##  ToDos:
  
  ##  function computes the covariance matrices of 
  ##  - bhat
  ##  - betahat
  ##  - bhat and betahat
  ##  - delta.b = b - bhat
  ##  - delta.b and betahat
  ##  - residuals ehat = y - X betahat - bhat
  ##  - residuals ehat.p.bhat = y - X betahat = ehat + bhat
  ##  - auxiliary matrix to compute covariance between kriging predictions of
  ##    y and y
  
  ## 2011-10-13 A. Papritz
  ## 2011-12-14 AP modified for replicated observations
  ## 2012-02-23 AP checking new variant to compute covariances of betahat and bhat
  ## 2012-04-27 AP scaled psi-function
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2012-11-04 AP unscaled psi-function
  ## 2013-02-05 AP covariance matrix of xihat
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-05-06 AP changes for solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-19 AP correcting error when computing covariance of regression residuals
  
  ## adjust flags for computing covariance matrices
  
  cov.bhat.b       <- FALSE
  cov.bhat.e       <- FALSE
  cov.betahat.b    <- FALSE
  cov.betahat.e    <- FALSE
  
  if( any( c( cov.delta.bhat, aux.cov.pred.target )))                          cov.bhat.b <- TRUE
  if( any( c( cov.ehat, aux.cov.pred.target )))                                cov.bhat.e <- TRUE
  if( any( c( cov.delta.bhat.betahat, cov.ehat.p.bhat, aux.cov.pred.target ))) cov.betahat.b <- TRUE
  if( any( c( cov.ehat, cov.ehat.p.bhat, aux.cov.pred.target )))               cov.betahat.e <- TRUE
  if( any( c( cov.ehat, cov.ehat.p.bhat ) ) )                                  cov.betahat <- TRUE
  if( any( c( cov.delta.bhat, cov.delta.bhat.betahat ))){
    cov.bhat <- TRUE
    if( full.cov.delta.bhat )                                                  full.cov.bhat <- TRUE
  }
  if( cov.delta.bhat.betahat )                                                 cov.bhat.betahat <- TRUE
  if( cov.ehat ){
    cov.delta.bhat.betahat <- TRUE
    cov.delta.bhat <- TRUE
    if( full.cov.ehat )                                                        full.cov.delta.bhat <- TRUE
  }
  
  ## compute required auxiliary items 
  
  a <- expectations["psi2"] 
  b <- expectations["dpsi"]
  
  TtT <- as.vector( table( TT ) )

  V <- eta * nugget * Valpha.objects[["Valpha"]] 
  VTtT <- t( TtT * V )
  
  result.new <- list( error = FALSE )
  
  ## inverse of G_theta (note that Palpha is not symmetric!!)
  
  Gi <- try( solve( t(TtT * Valpha.objects[["Valpha"]]) + Palpha / ( b * eta ) ), silent = TRUE )
  if( identical( class( Gi ), "try-error" ) ){
    result.new[["error"]] <- TRUE
    return( result.new )            
  }
  
  ## factors to compute bhat and betahat from xihat
  
  if( any( c( 
        cov.bhat.b, cov.bhat.e,
        cov.bhat, cov.bhat.betahat, cov.betahat ) ) 
  ){
    aux <- pmm( Gi, Valpha.objects[["Valpha"]], control.parallel )
    PaGtiVa <- pmm( Palpha, aux, control.parallel )
  }
  
  if( any( c( 
        cov.betahat.b, cov.betahat.e,
        cov.betahat, cov.bhat.betahat
      ) ) 
  ){
    aux <- pmm( Gi, Valpha.objects[["Valpha"]], control.parallel )
    AaGtiVa <- pmm( Aalpha, aux, control.parallel )
  }
  
  ## covariance of huberized observations
  
  if( any( c( cov.bhat, cov.betahat, cov.bhat.betahat ) )
  ){
    cov.b.psi <- TtT * VTtT
    diag( cov.b.psi ) <- diag( cov.b.psi ) + (a * nugget / b^2) * TtT
    PaGtiVa.cov.b.psi <- pmm( PaGtiVa, cov.b.psi, control.parallel )
  }
  
  ## covariance of bhat and betahat with B and epsilon
  
  if( cov.bhat.b )    cov.bhat.b      <- pmm( PaGtiVa, t( VTtT ), control.parallel )
  if( cov.bhat.e )    cov.bhat.e      <- (nugget * PaGtiVa)[, TT]
  if( cov.betahat.b ){
    cov.betahat.b <- AaGtiVa %*% t( VTtT )
    TX.cov.betahat.bT <- (XX %*% cov.betahat.b)[TT,TT]
  }
  if( cov.betahat.e ){
    cov.betahat.e <- (nugget * AaGtiVa)[, TT]
    TX.cov.betahat.e <- (XX %*% cov.betahat.e)[TT,]
  }
  
  ## compute now the requested covariances ...
  
  ## ... of bhat (debugging status ok)
  
  if( cov.bhat ){
    result.new[["cov.bhat"]] <- if( full.cov.bhat )
    {
      #       aux <- tcrossprod( aux, PaGtiVa )
      aux <- pmm( PaGtiVa.cov.b.psi, t(PaGtiVa ), control.parallel )
      attr( aux, "struc" ) <- "sym"
      aux
    } else {
      aux <- rowSums( aux * PaGtiVa )
      names( aux ) <- rownames( XX )
      aux
    }
  }
  
  ## ... of betahat (debugging status ok)
  
  if( cov.betahat ){
    result.new[["cov.betahat"]] <- tcrossprod( tcrossprod( AaGtiVa, cov.b.psi ), AaGtiVa )
    attr( result.new[["cov.betahat"]], "struc" ) <- "sym"
  }
  
  ##  ... of bhat and betahat (debugging status ok)
    
  if( cov.bhat.betahat ){
    result.new[["cov.bhat.betahat"]] <- tcrossprod( PaGtiVa.cov.b.psi, AaGtiVa )
  }
  
  ## ... of (b - bhat) (debugging status ok)
  
  if( cov.delta.bhat ){
    result.new[["cov.delta.bhat"]] <- if( full.cov.delta.bhat )
    {
      aux <- V + result.new[["cov.bhat"]] - cov.bhat.b - t( cov.bhat.b )
      attr( aux, "struc" ) <- "sym"
      dimnames( aux ) <- list( rownames( XX ), rownames( XX ) )
      aux
    } else {
      aux <- diag( V ) - 2 * diag( cov.bhat.b ) + (if( full.cov.bhat ){
        diag( result.new[["cov.bhat"]] )
      } else {
        result.new[["cov.bhat"]]
      })
      names( aux ) <- rownames( XX )
      aux
    }
  }
  
  ## ... of (b - bhat) and betahat (debugging status ok)
  
  if( cov.delta.bhat.betahat ){
    result.new[["cov.delta.bhat.betahat"]] <- t( cov.betahat.b ) - result.new[["cov.bhat.betahat"]]
    dimnames( result.new[["cov.delta.bhat.betahat"]] ) <- dimnames( XX )
  }
  
  ## ... of ehat (debugging status ok)
 
  if( cov.ehat ){
    aux1 <- tcrossprod( result.new[["cov.delta.bhat.betahat"]], XX )[TT,TT]
    result.new[["cov.ehat"]] <- if( full.cov.ehat )
    {
      aux <- bla <- result.new[["cov.delta.bhat"]][TT,TT] + 
        tcrossprod( tcrossprod( XX, result.new[["cov.betahat"]] ), XX )[TT,TT] -
        aux1 - t(aux1) - cov.bhat.e[TT,] - t(cov.bhat.e)[,TT] - 
        TX.cov.betahat.e - t(TX.cov.betahat.e)
      diag( aux ) <- diag( aux ) + nugget
      attr( aux, "struc" ) <- "sym"
      dimnames( aux ) <- list( names.yy, names.yy )
      aux   
    } else {
      aux <- (if( full.cov.delta.bhat ){
        diag( result.new[["cov.delta.bhat"]] )[TT] 
      } else {
        result.new[["cov.delta.bhat"]][TT]
      }) + rowSums( XX * (XX %*% result.new[["cov.betahat"]]) )[TT] -
        2 * diag( aux1 ) - 2 * diag( cov.bhat.e[TT,] ) - 2 * diag( TX.cov.betahat.e ) + 
        nugget
      names( aux ) <- names.yy
      aux
    }
  }
  
  ## ... of ehat + bhat (debugging status ok)
  
  if( cov.ehat.p.bhat ){
    result.new[["cov.ehat.p.bhat"]] <- if( full.cov.ehat.p.bhat )
    {
      aux <- tcrossprod( tcrossprod( XX, result.new[["cov.betahat"]] ), XX )[TT,TT] - 
        TX.cov.betahat.bT - t(TX.cov.betahat.bT) -
        TX.cov.betahat.e - t(TX.cov.betahat.e) + V[TT,TT]
      diag( aux ) <- diag( aux ) + nugget
      attr( aux, "struc" ) <- "sym"
      dimnames( aux ) <- list( names.yy, names.yy )
      aux   
    } else {
      aux <- rowSums( XX * (XX %*% result.new[["cov.betahat"]]) )[TT] - 
        2 * diag( TX.cov.betahat.bT ) - 
        2 * diag( TX.cov.betahat.e ) + diag( V )[TT] + nugget
      names( aux ) <- names.yy
      aux
    }
  }
  
  ## ...  auxiliary item to compute covariance of kriging predictions
  ## and observations
  
  if( aux.cov.pred.target ){
    result.new[["cov.pred.target"]] <- pmm(
      rbind( cov.bhat.b, cov.betahat.b ),
      Valpha.objects[["Valpha.inverse"]] / eta / nugget,
      control.parallel
    )
    
  }
    
  return( result.new )
  
}

##   ##############################################################################

update.xihat <- 
  function( 
    XX, yy, res, TT, 
    nugget, eta, 
    Valpha.inverse.Palpha,
    psi.function, tuning.psi, 
    verbose
  )
{
  
  ## 2013-02-04 AP solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  ## function computes (1) updated IRWLS estimates xihat of linearized
  ## normal equations, (2) the associated rweights,
  ## (3) the unstandardized residuals (= estimated epsilons); the results
  ## are returned as a list
  
  ## compute rweights (cf. p. 7, proposal HRK of 26 July 2010)

  ## 2013-04-23 AP new names for robustness weights

  std.res <- res / sqrt( nugget )
  
  ##  construct left-hand side matrix M and right-hand side vector of
  ##  linear equation system
  
  Wi <- ifelse( 
    abs( std.res ) < sqrt( .Machine[["double.eps"]] ),
    1.,
    psi.function( std.res, tuning.psi ) / std.res
  )
  
  ##  aggregate rweights for replicated observations
  
  if( sum( duplicated( TT ) ) > 0 ){
    
    TtWiT  <- as.vector( tapply( Wi, factor( TT ), sum ) )
    TtWiyy <- as.vector( tapply( Wi * yy, factor( TT ), sum ) )
    
  } else {
    
    TtWiT <- Wi
    TtWiyy <- Wi * yy
    
  }
  
  ##  construct left-hand side matrix M and right-hand side vector b of
  ##  linearized system of equations
    
  M <- Valpha.inverse.Palpha / eta
  diag( M ) <- diag( M ) + TtWiT
  
  b <- TtWiyy
  
  ##  solve linear system
  
  result <- list( error = TRUE )
  
  r.solve <- try( solve( M, b ), silent = TRUE ) 
  
  if( !identical( class( r.solve ), "try-error" ) ) {
    
    ##  collect output
    
    result[["error"]]      <- FALSE
    result[["xihat"]]      <- r.solve
    result[["residuals"]]  <- yy - result[["xihat"]][TT]
    result[["rweights"]]   <- Wi
    
  }
  
  return( result )
  
}


##    ##############################################################################

estimating.eqations.xihat <- function( 
  res, TT, xihat, nugget, eta, Valpha.inverse.Palpha, 
  psi.function, tuning.psi
){
  
  ## auxiliary function to compute estimating equations for xihat

  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)

  Ttpsi <- psi.function( res / sqrt( nugget ), tuning.psi )
  
  if( sum( duplicated( TT ) > 0 ) ){
    Ttpsi <- as.vector( tapply( Ttpsi, factor( TT ), sum ) )
  }
  
  Ttpsi - drop( Valpha.inverse.Palpha %*% xihat ) / sqrt( nugget ) / eta
}

##    ##############################################################################

estimate.xihat <- 
  function(
    compute.xihat,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, tuning.psi, tuning.psi.nr, 
    maxit, reltol,
    nugget, eta, Valpha.inverse,
    control.parallel,
    verbose
  )
{
  
  ## 2013-02-04 AP solving estimating equations for xi
  ## 2013-06-03 AP handling design matrices with rank < ncol(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
   
  ## function computes (1) estimates xihat, bhat, betahat by
  ## solving robustified estimating equations by IRWLS,
  ## (2) the weights of the IRWLS, (3) the unstandardized residuals
  ## (= estimated epsilons); the results are returned as a list
  
  ##  compute projection matrix Palpha and related items
  
  result <- list( error = FALSE )
  
  aux <- t( XX ) %*% Valpha.inverse
  
  if( rankdef.x ){
    s <- svd( aux %*% XX )
    s[["d"]] <- ifelse( s[["d"]] / max( s[["d"]] ) <= min.condnum, 0., 1. / s[["d"]] )
    Palpha <- s[["v"]] %*% ( s[["d"]] * t( s[["u"]] ) )                # Moore-Penrose inverse
  } else {
    t.chol <- try( chol( aux %*% XX ), silent = TRUE )
    if( !identical( class( t.chol ), "try-error" ) ){
      Palpha <- chol2inv( t.chol )    
    } else {
      result[["error"]] <- TRUE
      return( result )    
    }
  }
  
  result[["Aalpha"]]             <- Palpha %*% aux
  dimnames( result[["Aalpha"]] ) <- dimnames( t(XX) )
  
  result[["Palpha"]]             <- -XX %*% result[["Aalpha"]]
  diag( result[["Palpha"]] )     <- diag( result[["Palpha"]] ) + 1.
  rownames( result[["Palpha"]] ) <- rownames( XX )
  colnames( result[["Palpha"]] ) <- rownames( XX )
  
  result[["Valpha.inverse.Palpha"]] <- pmm( 
    Valpha.inverse, result[["Palpha"]], control.parallel 
  )
  rownames( result[["Valpha.inverse.Palpha"]] ) <- rownames( XX )
  colnames( result[["Valpha.inverse.Palpha"]] ) <- rownames( XX )
  
  if( compute.xihat ){
  
    ##  initialization
    
    res <- yy - xihat[TT]
    
    eeq.old <- estimating.eqations.xihat(     
      res, TT, xihat, nugget, eta, result[["Valpha.inverse.Palpha"]], 
      psi.function, tuning.psi
    )
    eeq.old.l2 <- sum( eeq.old^2 )
    
    if( !is.finite( eeq.old.l2 ) ) {
      result[["error"]] <- TRUE
      return( result )
    }
    
    converged <- FALSE
    
    if( verbose > 2 ) cat(
      "\n  IRWLS\n",
      "      it        L2.old        L2.new      delta.L2\n", sep = ""
    )
    
    ##  IRWLS
    
    for( i in 1:maxit ){
      
      ##  compute new estimates 
      
      new <- update.xihat(
        XX, yy, res, TT, 
        nugget, eta, 
        result[["Valpha.inverse.Palpha"]],
        psi.function, tuning.psi, 
        verbose
      )
      
      if( new[["error"]] ) {
        result[["error"]] <- TRUE
        return( result )
      }
      
      
      ##  evaluate estimating equations for xi and compute its l2 norm
      
      eeq.new <- estimating.eqations.xihat(       
        new[["residuals"]], TT, new[["xihat"]], nugget, eta, result[["Valpha.inverse.Palpha"]], 
        psi.function, tuning.psi
      )
      eeq.new.l2 <- sum( eeq.new^2 )
      
      if( !is.finite( eeq.new.l2 ) ) {
        result[["error"]] <- TRUE
        return( result )
      }
      
      if( verbose > 2 ) cat( 
        format( i, width = 8 ),
        format( 
          signif( 
            c( eeq.old.l2, eeq.new.l2, eeq.old.l2 - eeq.new.l2 ), digits = 7 
          ), scientific = TRUE, width = 14 
        ), "\n", sep = ""
      )
      
      ##  check for convergence (cf. help( optim ) )
      
      if( max( abs( res - new[["residuals"]] ) ) < sqrt(  reltol ) * sqrt( nugget ) ) {
        converged <- TRUE
        break
      }
      
      ##  update xihat, residuals and eeq.old.l2
      
      eeq.old.l2 <- eeq.new.l2
      xihat      <- new[["xihat"]]
      res        <- new[["residuals"]]
      
    }
    
    ##  collect output
    
    result[["xihat"]]            <- new[["xihat"]]
    names( result[["xihat"]] )   <- rownames( XX )
    
    result[["residuals"]]        <- new[["residuals"]]
    result[["rweights"]]         <- new[["rweights"]]
    result[["converged"]]        <- converged
    result[["nit"]]              <- i
    
  } else {
    
    result[["xihat"]]            <- xihat
    names( result[["xihat"]] )   <- rownames( XX )
    
    result[["residuals"]]        <- yy - xihat[TT]
    
    result[["rweights"]]         <- ifelse( 
      abs( std.res <- result[["residuals"]] / sqrt( nugget ) ) < sqrt( .Machine[["double.eps"]] ),
      1.,
      psi.function( std.res, tuning.psi ) / std.res
    )
    result[["converged"]]        <- NA
    result[["nit"]]              <- NA_integer_
    
  }
  
  result[["bhat"]]             <- drop( result[["Palpha"]] %*% result[["xihat"]] )
  names( result[["bhat"]] )    <- rownames( XX )
  
  result[["betahat"]]          <- drop( result[["Aalpha"]] %*% result[["xihat"]] )
  names( result[["betahat"]] ) <- colnames( XX )
  
  result[["z.star"]]           <- drop( Valpha.inverse %*% result[["bhat"]] )

  return( result )
  
}


##    ##############################################################################

gcr <- 
  function( 
    lag.vectors, variogram.model, param, aniso, 
    irf.models = georob.control()[["irf.models"]],
    control.parallel, verbose
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
  
  result <- list( error = TRUE )
  
  ## matrix for coordinate transformation
  
  A <- aniso[["sclmat"]] * aniso[["rotmat"]] / param["scale"]
  
  ## prepare model
  
  model.list <- list( variogram.model )
  model.list <- c( model.list, as.list( param[-(1:4)] ) )
  model.list <- list( "$", var = 1., A = A, model.list )

  ##  negative semivariance matrix
  
  ## functions of version 3 of RandomFields
  
  RFoptions(newAniso=FALSE)
  
  ## auxiliary function to compute generalized correlations in parallel
  
  f.aux <- function(i, s, e, lag.vectors, model.list ){
    result <- try(
      -RFvariogram(
        x = lag.vectors[s[i]:e[i], ], model = model.list, dim = NCOL( lag.vectors ), grid = FALSE
      ),
      silent = TRUE
    )
    if( !(identical( class( result ), "try-error" ) || any( is.na( result ) )) ){
      result
    } else {
      "RFvariogram.error"
    }
  }
  
  ## definition of junks to be evaluated in parallel
  
  k <- control.parallel[["f"]] * control.parallel[["pmm.ncores"]]
  n <- NROW(lag.vectors)
  dn <- floor( n / k )
  s <- ( (0:(k-1)) * dn ) + 1
  e <- (1:k) * dn
  e[k] <- n
  
  ## compute generalized correlations in parallel
  
  if( control.parallel[["pmm.ncores"]] > 1L ){
    
    if( identical( .Platform[["OS.type"]], "windows") ){
      
      if( !sfIsRunning() ){
        options( error = f.stop.cluster )
        junk <- sfInit( parallel = TRUE, cpus = control.parallel[["pmm.ncores"]] )
        #         junk <- sfLibrary( RandomFields, verbose = FALSE )
        junk <- sfLibrary( georob, verbose = FALSE )
      }
      
        
      Valpha0 <- sfLapply( 
        1:k, f.aux, s = s, e = e, lag.vectors = lag.vectors, model.list = model.list 
      )
      
      if( control.parallel[["sfstop"]] ){
        junk <- sfStop()
        options( error = NULL )
      }
      
    } else {
      
      Valpha0 <- mclapply( 
        1:k, f.aux, s = s, e = e, lag.vectors = lag.vectors, model.list = model.list,  
        mc.cores = control.parallel[["pmm.ncores"]] 
      )
      
    }
    
    not.ok <- any( sapply( Valpha0, function( x ) identical( x, "RFvariogram.error" ) ) )
    
  } else {
    
    Valpha0 <- try(
      -RFvariogram(
        x = lag.vectors, model = model.list, dim = NCOL( lag.vectors ), grid = FALSE
      ),
      silent = TRUE
    )
    
    not.ok <- identical( class( Valpha0 ), "try-error" ) || any( is.na( Valpha0 ) )
    
  }
  
  if( !not.ok ){
    
    Valpha0 <- unlist( Valpha0 )
    
    ## partial sill to total variance of z process
    
    ps <- unname( param["variance"] / sum( param[c( "variance", "snugget" )] ) )
    
    ## convert semivariance vectors to symmetric matrices
    
    Valpha0 <- list( 
      diag = rep( 0., 0.5 * ( 1 + sqrt( 1 + 8 * length( Valpha0 ) ) ) ),
      tri = Valpha0
    )
    attr( Valpha0, "struc" ) <- "sym"
    
    Valpha <- Valpha0
    Valpha[["tri"]] <- ps * ( Valpha[["tri"]] + 1. ) - 1.
    
    Valpha0 <- expand( Valpha0 )
    Valpha  <- expand( Valpha )
    
    ##  compute additive constant for positive definiteness and
    
    if( variogram.model %in% irf.models ){
      gcr.constant.Valpha0 <- max( -Valpha0 ) * 2.                    
      gcr.constant.Valpha  <- max( -Valpha )  * 2.                    
    } else {
      gcr.constant.Valpha0 <- gcr.constant.Valpha <- 1.
    }
    
    ##  collect results
    
    result[["error"]]        <- FALSE
    result[["gcr.constant"]] <- gcr.constant.Valpha
    result[["Valpha0"]]      <- Valpha0 + gcr.constant.Valpha0  # correlation matrix that does not contain spatial nugget
    result[["Valpha"]]       <- Valpha  + gcr.constant.Valpha   # correlation matrix that includes spatial nugget
    
  } else {
    
    if( verbose > 3 ) cat(
      "\n an error occurred when computing the negative semivariance matrix\n"
    )
    
  }
  
  return( result )
    
}




##    ##############################################################################

prepare.likelihood.calculations <- 
  function(
    envir,
    adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, dpsi.function, tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.reltol,
    compute.xihat = TRUE,
    compute.Q,
    control.parallel, 
    verbose
  )
{
  
  ## 2011-12-10 AP modified for replicated observations
  ## 2012-02-03 AP modified for geometrically anisotropic variograms
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-02 AP modification computing ilcf
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2012-11-27 AP changes in check allowed parameter range
  ## 2013-02-04 AP solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-02 AP new transformation of rotation angles
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  
  ##  function transforms (1) the variogram parameters back to their
  ##  original scale; computes (2) the correlation matrix, its inverse
  ##  and its inverse lower cholesky factor; (3) computes betahat,
  ##  bhat and further associates items; and (4) computes the
  ##  matrices A and the cholesky factor of the matrix Q
  
  ##  transform variogram and anisotropy parameters back to original scale
  
  param <- c( adjustable.param, fixed.param )[param.name]
  
  param <- sapply(
    param.name,
    function( x, param.tf, param ) bwd.tf[[param.tf[x]]]( param[x] ),
    param.tf = param.tf,
    param = param
  )
  names( param ) <- param.name
  
  fwd.tf.aniso <- aniso<- c( adjustable.param, fixed.param )[aniso.name]
  
  aniso <- sapply(
    aniso.name,
    function( x, param.tf, param ) bwd.tf[[param.tf[x]]]( param[x] ),
    param.tf = param.tf,
    param = aniso
  )
  names( aniso ) <- aniso.name
  
  ##  check whether the current variogram parameters and the variogram
  ##  parameters that were used in the previous call to
  ##  prepare.likelihood.calculations are the same
  
  lik.item <- get( "lik.item", pos = as.environment( envir ) )
  
  same.param <- isTRUE( all.equal( 
      c( param, aniso ), 
      c( lik.item[["param"]], lik.item[["aniso"]][["aniso"]] )
    ))
  
  if( same.param && !is.null( lik.item[["effects"]] ) ){
    
    if( verbose > 4 ) cat(
      "\n     prepare.likelihood.calculations: exit without computing any objects\n"    
    )
    return( lik.item )
  }
  
    
  ##  compute updates of required likelihood items if the
  ##  variogram parameters differ
  
  if( !( same.param && length( lik.item[["Valpha"]] ) && !length( adjustable.param ) ) ){
    
    if( verbose > 4 ) cat(
      "\n     prepare.likelihood.calculations: computing 'Valpha' object\n"    
    )
    
    lik.item[["Valpha"]][["error"]] <- TRUE

    ## check whether variogram parameters are within reasonalble bounds and
    ## return an error otherwise
    
    if( length( c( param, aniso ) ) && any( c( param, aniso ) > safe.param ) ){
      if( verbose > 1 ){
        if( !lik.item[["aniso"]][["isotropic"]] ) param <- c( param, aniso )
        cat( "\n" )
        print( signif( param ) )
      }
      return( lik.item )  
    }
    
    ## check whether extra variogram parameters are within allowed bounds and
    ## return an error otherwise
    
    ep <- param.names( model = variogram.model )
    param.bounds <- param.bounds( variogram.model, NCOL( lag.vectors ) )
    ep.param <- param[ep]
    
    if( !is.null( param.bounds ) ) t.bla <- sapply(
      1:length( ep.param ),
      function( i, param, bounds ){
        if( param[i] < bounds[[i]][1] || param[i] > bounds[[i]][2] ) cat(
          "value of parameter '", names( param[i] ), "' outside of allowed range", sep = "" 
        )
        return( lik.item )
      }, 
      param = ep.param,
      bounds = param.bounds
    )
    
    ##  update variogram and parameters and compute eta
    
    lik.item[["param"]] <- param
    lik.item[["eta"]] <- sum( param[c( "variance", "snugget" )] ) / param["nugget"] 
    
    ##  update anisotropy parameters and compute the coordinate
    ##  transformation matrices
    
    lik.item[["aniso"]][["aniso"]] <- aniso
    lik.item[["aniso"]][["sincos"]] <- list(
      co = unname( cos( fwd.tf.aniso["omega"] ) ),
      so = unname( sin( fwd.tf.aniso["omega"] ) ),
      cp = unname( cos( fwd.tf.aniso["phi"] ) ),
      sp = unname( sin( fwd.tf.aniso["phi"] ) ),
      cz = unname( cos( fwd.tf.aniso["zeta"] ) ),
      sz = unname( sin( fwd.tf.aniso["zeta"] ) )
    )
    
    n <- NCOL( lag.vectors)
    
    if( n <= 3 ){
      
      lik.item[["aniso"]][["rotmat"]] <- with( 
        lik.item[["aniso"]][["sincos"]],
        rbind(
          c(             sp*so,             sp*co,       cp ),
          c( -cz*co + sz*cp*so,  co*sz*cp + cz*so,   -sp*sz ),
          c( -co*sz - cz*cp*so, -cz*co*cp + sz*so,    cz*sp )
        )[ 1:n, 1:n, drop = FALSE ]
      )
      
      
      lik.item[["aniso"]][["sclmat"]] <- 1. / c( 1., aniso[ c("f1", "f2") ] )[ 1:n ]
      
    } else {  # only isotropic case for n > 3
      
      lik.item[["aniso"]][["rotmat"]] <- diag( n )
      lik.item[["aniso"]][["sclmat"]] <- rep( 1., n )
      
    }
    
    ##  calculate generalized correlation matrix, its inverse and its
    ##  inverse cholesky factor
    
    t.Valpha <- f.aux.Valpha(
      lag.vectors = lag.vectors, variogram.model = variogram.model,
      param = param, aniso = lik.item[["aniso"]], control.parallel = control.parallel,
      verbose = verbose
    )
    
    if( !t.Valpha[["error"]] ){
      lik.item[["Valpha"]] <- t.Valpha
    } else {
      return( lik.item )
    }
    
  } else if( is.null( lik.item[["aniso"]][["sincos"]] ) ){
    
    ## add missing anisotropy items
    
    lik.item[["aniso"]][["sincos"]] <- list(
      co = unname( cos( fwd.tf.aniso["omega"] ) ),
      so = unname( sin( fwd.tf.aniso["omega"] ) ),
      cp = unname( cos( fwd.tf.aniso["phi"] ) ),
      sp = unname( sin( fwd.tf.aniso["phi"] ) ),
      cz = unname( cos( fwd.tf.aniso["zeta"] ) ),
      sz = unname( sin( fwd.tf.aniso["zeta"] ) )
    )
    
    n <- NCOL( lag.vectors)
    
    if( n <= 3 ){
      
      lik.item[["aniso"]][["rotmat"]] <- with( 
        lik.item[["aniso"]][["sincos"]],
        rbind(
          c(             sp*so,             sp*co,       cp ),
          c( -cz*co + sz*cp*so,  co*sz*cp + cz*so,   -sp*sz ),
          c( -co*sz - cz*cp*so, -cz*co*cp + sz*so,    cz*sp )
        )[ 1:n, 1:n, drop = FALSE ]
      )
      
      
      lik.item[["aniso"]][["sclmat"]] <- 1. / c( 1., aniso[ c("f1", "f2") ] )[ 1:n ]
      
    } else {  # only isotropic case for n > 3
      
      lik.item[["aniso"]][["rotmat"]] <- diag( n )
      lik.item[["aniso"]][["sclmat"]] <- rep( 1., n )
      
    }
    
  }
  
  t.param <- lik.item[["param"]]
  if( !lik.item[["aniso"]][["isotropic"]] ) t.param <- c( t.param, lik.item[["aniso"]][["aniso"]] )
  if( verbose > 1 ) {
    cat( "\n" )
    print( signif( t.param ) )
  }
  
  ##  estimate fixed and random effects (xihat, betahat, bhat,
  ##  residuals )
  
  ##  either take initial guess of betahat and bhat for the current
  ##  irwls iteration from initial.object or from previous iteration
  
  if( is.null( lik.item[["effects"]] ) || !same.param ){
    
    if( verbose > 4 ) cat(
      "\n     prepare.likelihood.calculations: computing 'effects' object\n"    
    )
    
    if( 
      !irwls.initial && !is.null( lik.item[["effects"]][["xihat"]] ) 
    ){
      xihat <- lik.item[["effects"]][["xihat"]]
    }
    
    lik.item[["effects"]] <- estimate.xihat( 
      compute.xihat,
      XX, min.condnum, rankdef.x, yy, TT, xihat, 
      psi.function, tuning.psi, tuning.psi.nr, 
      irwls.maxiter, irwls.reltol,
      lik.item[["param"]]["nugget"], lik.item[["eta"]], 
      lik.item[["Valpha"]][["Valpha.inverse"]],
      control.parallel,
      verbose
    )      
    
    if( lik.item[["effects"]][["error"]] ) return( lik.item )     ##  an error occurred
  }
  
  
  ##  compute Q matrix and its Cholesky factor (required for
  ##  non-robust (RE)ML estimation)
  
  if( compute.Q && (is.null( lik.item[["Q"]] ) || !same.param ) ) {
    
    if( verbose > 4 ) cat(
      "\n     prepare.likelihood.calculations: computing 'Q' object\n"    
    )
    
    t.Q <- f.aux.Q( 
      TT = TT, XX = XX,  rankdef.x =  rankdef.x, 
      lik.item = lik.item, min.condnum = min.condnum,
      ml.method = ml.method, control.pmm = control.parallel
    )
    
    if( !t.Q[["error"]] ){
      lik.item[["Q"]] <- t.Q
    } else {
      return( lik.item )
    }
    
  }
  
  ##  store updated lik.item object
  
  assign( "lik.item", lik.item, pos = as.environment( envir ) )
  
#   print( str( lik.item ) ); stop()
  
  return( lik.item )
    
}


##   ##############################################################################

dcorr.dparam <- 
  function(
    x, variogram.model, param, d.param, aniso, verbose 
  )
{
  
  ##  Function to compute partial derivatives of generalized
  ##  correlation matrix with respect to scale and extra parameters
  
  ##  Arguments:
  ##  x             lag vectors for all pairs of distinct locations
  ##  variogram.model         Covariance Model as in Variogram{RandomFields}
  ##  param         Vector with variogram parameters
  ##  d.param       String, Parameter for which to determine the derivative
  
  ##  Value:
  ##  Vector or Matrix with partial derivative of Valpha for scale and extra parameters
  ##                named a, b, c, ... as in Variogram{RandomFields}
  
  ##  References:
  ##  help(Variogram)
  ##  Chiles and Delfiner, Section 2.5
  
  ##  06 Apr 2011  C.Schwierz
  ##  2011-07-17 ap
  ##  2012-01-24 ap RMcauchytbm and RMlgd models added
  ##  2012-01-25 ap extra model parameter with same names as in Variogram{RandomFields}
  ##  2012-02-07 AP modified for geometrically anisotropic variograms
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-05-15 AP changes for version 3 of RandomFields
  
  aniso.name <- names( aniso[["aniso"]] )
  alpha <- unname( param["scale"] )
  n = NCOL( x )
  aux <- aniso[["rotmat"]] %*% t(x)
  
  ## scaled lag distance
  
  hs <- sqrt( colSums( ( aniso[["sclmat"]] * aux )^2 ) ) / alpha
  
  ## partial derivatives of scaled lag distance with respect to
  ## anisotropy parameters
  
  dhs.daniso <- switch(
    
    d.param,
    
    f1 = {
      colSums(
        ( c( 0., -1. / aniso[["aniso"]]["f1"]^2, 0. )[1:n] * aniso[["sclmat"]] ) * aux^2 
      )
    },
    
    f2 = { 
      colSums(
        ( c( 0., 0., -1. / aniso[["aniso"]]["f2"]^2 )[1:n] * aniso[["sclmat"]] ) * aux^2 
      )
    },
    omega = {
      drotmat <- with(
        aniso[["sincos"]],
        rbind(
          c(             sp*co,            -sp*so, 0. ),
          c(  co*sz*cp + cz*so,  cz*co - sz*cp*so, 0. ),
          c( -cz*co*cp + sz*so,  co*sz + cz*cp*so, 0. )
        )[ 1:n, 1:n, drop = FALSE ]
      )
      colSums( 
        ( aniso[["sclmat"]] * drotmat %*% t(x) ) * ( aniso[["sclmat"]] * aux ) 
      )
    },
    
    phi = {
      drotmat <- with(
        aniso[["sincos"]],
        rbind(
          c(     cp*so,     cp*co,    -sp ),
          c( -sz*sp*so, -co*sz*sp, -cp*sz ),
          c(  cz*sp*so,  cz*co*sp,  cz*cp )
        )[ 1:n, 1:n, drop = FALSE ]
      )
      colSums( 
        ( aniso[["sclmat"]] * drotmat %*% t(x) ) * ( aniso[["sclmat"]] * aux ) 
      )
    },
    
    zeta = {
      drotmat <- with(
        aniso[["sincos"]],
        rbind(
          c(                0.,               0.,     0. ),
          c(  co*sz + cz*cp*so, cz*co*cp - sz*so, -cz*sp ),
          c( -cz*co + sz*cp*so, co*sz*cp + cz*so, -sp*sz )
        )[ 1:n, 1:n, drop = FALSE ]
      )
      colSums( 
        ( aniso[["sclmat"]] * drotmat %*% t(x) ) * ( aniso[["sclmat"]] * aux ) 
      )
    },
    
    NA
  ) / ( hs * alpha^2 )
  
  ##  partial derivative of scaled lag distance with respect to scale
  ##  parameter
  
  dhs.dscale <- -hs / alpha
  
  ##  compute derivative of generalized correlation matrix with
  ##  respect to scale and extra parameters
  
  result <- switch(
    variogram.model,
    
    RMbessel = {
      
      A <- unname( param["nu"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( 2^A * besselJ( hs, 1+A ) * gamma( 1+A ) ) / hs^A
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( 2^A * besselJ( hs, 1+A ) * gamma(1 + A) ) / hs^A,
      #   0.
      # )
      
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        nu = {
          myenv <- new.env()
          assign( "hs", hs, envir = myenv )
          assign( "nu", param["nu"], envir = myenv )
          as.vector( 
            attr( 
              numericDeriv( 
                expr = quote( 
                  2^nu * gamma( nu+1 ) * besselJ( hs, nu ) / hs^nu 
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
    
    RMcauchy = {
      
      A <- unname( param["gamma"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -2 * A * hs * ( 1+hs^2 )^(-1-A)        
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        gamma = {
          -( 1 + hs^2 )^(-A) * log( 1 + hs^2 )
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
#         B * hs^(-1+A) * (1+hs^A)^(-2-B/A) * ( A + C - (-B + C) * hs^A ) 
#       ) / C
#       # dgc.dhs <- ifelse(
#       #   hs > 0.,
#       #   -( 
#       #     B * hs^(-1+A) * (1+hs^A)^(-2-B/A) * ( A + C - B * hs^A + C * hs^A) 
#       #   ) / C,
#       #   if( A > 1. ){
#       #     0.
#       #   } else if( identical( A, 1. ) ){
#       #     -B * (1+C) / C
#       #   } else {
#       #     -Inf
#       #   }
#       # )
#       
#       
#       switch(
#         d.param,
#         scale = dgc.dhs * dhs.dscale,
#         # scale = {
#         #   ( B * hs^A * (1+hs^A)^(-2-B/A) * (A + C + (-B+C) * hs^A ) ) / ( C * scale )
#         # },
#         alpha = {
#           ( B * (1+hs^A)^(-2 - B/A) * (
#               -( A * hs^A * ( A + C + (-B+C) * hs^A ) * log(hs) ) + 
#               ( 1 + hs^A) * (C + (-B+C) * hs^A ) * log( 1+hs^A ) 
#             ) 
#           ) / (A^2 * C )
#         },
#         # alpha = {
#         #   ifelse(
#         #     hs > 0.,
#         #     ( B * (1+hs^A)^(-2 - B/A) * (
#         #         -( A * hs^A * ( A + C + (-B+C) * hs^A ) * log(hs) ) + 
#         #         ( 1 + hs^A) * (C + (-B+C) * hs^A ) * log( 1+hs^A ) 
#         #       ) 
#         #     ) / (A^2 * C ),
#         #     0.
#         #   )
#         # },
#         beta = {
#           ( -( A * hs^A) - (C + (-B+C) * hs^A ) * log( 1+hs^A ) ) / 
#           ( A*C * (1+hs^A)^( (A+B)/A ) )
#         },
#         gamma = {
#           ( B * hs^A ) / ( C^2 * (1+hs^A)^( (A+B)/A) )
#         },
#         dgc.dhs * dhs.daniso
#       )
#     }, ##  end case RMcauchytbm
    
    RMcircular = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- ( -4 * sqrt( 1-hs[sel]^2 ) ) / pi
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMcircular
    
    RMcubic = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- hs[sel] * ( -14. + 26.25*hs[sel] - 17.5*hs[sel]^3 + 5.25*hs[sel]^5 )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMcubic
    
    RMdagum = {
      
      A <- unname( param["beta"] )
      B <- unname( param["gamma"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -B / ( hs * ( 1+hs^(-A) )^(B/A) * ( 1+hs^A ) )
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( B / ( hs * ( 1+ hs^(-A) )^(B/A) * (1 + hs^A ) ) ),
      #   -Inf
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse(
        #     hs > 0.,
        #     B / ( ( 1 + hs^(-A) )^(B/A) * (1 + hs^A) * scale ),
        #     0.
        #   )
        # },
        beta = {
          -( B * ( A * log(hs) + (1+hs^A) * log( 1+hs^(-A) ) ) ) /
          ( A^2 * ( 1+hs^(-A) )^(B/A) * ( 1+hs^A ) )
        },
        # beta = {
        #   ifelse(
        #     hs > 0.,
        #     -( B * ( A * log(hs) + (1+hs^A) * log( 1+hs^(-A) ) ) ) /
        #     ( A^2 * ( 1+hs^(-A) )^(B/A) * ( 1+hs^A ) ),
        #     0.
        #   )
        # },
        gamma = {
          log( 1 + hs^(-A) ) / ( A * (1 + hs^(-A) )^(B/A) )
        },
        # gamma = {
        #   ifelse(
        #     hs > 0.,
        #     log( 1 + hs^(-A) ) / ( A * (1 + hs^(-A) )^(B/A) ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMdagum
    
    RMdampedcos = {
      
      A <- unname( param["lambda"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( ( A * cos(hs) + sin(hs) ) / exp( A*hs ) )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        lambda = {
          -exp( -A * hs ) * hs * cos( hs )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMdampedcos
    
    RMdewijsian = {
      
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -A / ( hs + hs^(1-A) )
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( A / ( hs + hs^(1-A) ) ),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( A * hs^A )/( scale + hs^A * scale )
        # },
        alpha = {
          -( ( hs^A * log( hs ) ) / ( 1 + hs^A ) )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -( ( hs^A * log( hs ) ) / ( 1 + hs^A ) ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMdewijsian
    
    
    RMexp = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -exp( -hs )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case exponential
    
    RMfbm = {
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -A * hs^(-1+A)
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( A * hs^(-1+A) ),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
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
    
    RMgauss = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -2 * hs / exp( hs^2 )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgauss
    
    RMgenfbm = {
      
      A <- unname( param["alpha"] )
      B <- unname( param["delta"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( A * B * hs^(-1+A) * (1+hs^A)^(-1+B))
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( A * B * hs^(-1+A) * (1+hs^A)^(-1+B)),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -B.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( A * B * hs^A * (1+hs^A)^(-1+B) ) / scale
        # },
        alpha = {
          -( B * hs^A * (1+hs^A)^(-1+B) * log(hs) )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     -( B * hs^A * (1+hs^A)^(-1+B) * log(hs) ),
        #     0.
        #   )
        # },
        delta = {
          -( (1 + hs^A )^B * log( 1 + hs^A ) )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgenfbm
    
    RMgencauchy = {
      
      A <- unname( param["alpha"] )
      B <- unname( param["beta"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( ( B * hs^(-1+A)) / (1+hs^A)^((A+B)/A))
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( ( B * hs^(-1+A)) / (1+hs^A)^((A+B)/A)),
      #   if( A < 1. ){
      #     -Inf
      #   } else if( identical( A, 1 ) ){
      #     -B.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( B * hs^A ) / ( (  1+hs^A)^((A+B)/A) * scale )
        # },
        alpha = {
          B * ( 1 + hs^A )^(-(A+B)/A) * (
            -A * hs^A * log( hs ) +
            ( 1 + hs^A ) * log( 1 + hs^A )
          ) / A^2
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     B * ( 1 + hs^A )^(-(A+B)/A) * (
        #       -A * hs^A * log( hs ) +
        #       ( 1 + hs^A ) * log( 1 + hs^A )
        #     ) / A^2,
        #     0.
        #   )
        # },
        beta = {
          -( log( 1+hs^A ) / ( A * (1+hs^A)^(B/A) ) )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgencauchy
    
    #     gengneiting = { Version 2 of package RandomFields
    #       
    #       
    #       A <- unname( param["n"] )
    #       B <- unname( param["alpha"] )
    #       
    #       ## derivative with of generalized covariance with respect to
    #       ## scaled lag distance
    #       
    #       dgc.dhs <- rep( 0., length( hs ) )
    #       sel <- hs < 1.
    #       dgc.dhs[sel] <- if( identical( A, 1 ) ){
    #         -( (1+B) * (2+B) * (1-hs[sel])^B * hs[sel] )
    #       } else if( identical( A, 2 ) ){
    #         -( (3+B) * (4+B) * (1-hs[sel])^(1+B) * hs[sel] * ( 1 + hs[sel] + B*hs[sel]) ) / 3.
    #       } else if( identical( A, 3 ) ){
    #         -( 
    #           (5+B) * (6+B) * (1-hs[sel])^(2+B) * hs[sel] * ( 3 + 3 * (2+B) * hs[sel] + (1+B) * (3+B) * hs[sel]^2 ) 
    #         ) / 15.
    #       } else {
    #         stop( "gengneiting model undefined for 'n' != 1:3" )
    #       }
    #       
    #       result <- rep( 0., length( hs ) )
    #       
    #       switch(
    #         d.param,
    #         scale = dgc.dhs * dhs.dscale,
    #         alpha = {
    #           result[sel] <- if( identical( A, 1 ) ){
    #             (1-hs[sel])^(1+B) * ( hs[sel] + (1 + hs[sel] + B*hs[sel]) * log( 1-hs[sel]) )
    #             
    #           } else if( identical( A, 2 ) ){
    #             (
    #               (1-hs[sel])^(2+B) * ( 
    #                 hs[sel] * ( 3 + 2 * (2+B) *hs[sel] ) + 
    #                 ( 3 + 3 * ( 2+B) * hs[sel] + ( 1+B) * (3+B) * hs[sel]^2 ) * log( 1-hs[sel] )
    #               )
    #             ) / 3.
    #           } else if( identical( A, 3 ) ){
    #             ( 
    #               (1-hs[sel])^(3+B) * ( 
    #                 hs[sel] * ( 15 + hs[sel] * ( 36 + 23*hs[sel] + 3 * B * ( 4 + (6+B)*hs[sel] ) ) ) + 
    #                 ( 15 + 15 * (3+B) * hs[sel] + ( 45 + 6 * B * (6+B) ) * hs[sel]^2 + (1+B) * (3+B) * (5+B) * hs[sel]^3 ) * 
    #                 log( 1-hs[sel]) 
    #               ) 
    #             ) / 15.
    #           }
    #           result
    #         },
    #         dgc.dhs * dhs.daniso
    #       )
    #       
    #       
    #     }, ##  end case Gengneiting
    
    RMgengneiting = {
      
      
      A <- unname( param["kappa"] )
      B <- unname( param["mu"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- if( identical( A, 1 ) ){
        (2.5+B) * (1-hs[sel])^(2.5+B) - (2.5+B) * (1-hs[sel])^(1.5+B) * (1+(2.5+B) * hs[sel])
      } else if( identical( A, 2 ) ){
        (1 - hs[sel])^(4.5+B) * (4.5+B + 2/3 * (3.5+B) * (5.5+B) * hs[sel] ) - 
        (4.5+B) * (1 - hs[sel])^(3.5+B) * (
          1 + hs[sel]*(4.5 + B + 6.416666666666666*hs[sel] + B/3. * (9.+B) * hs[sel] ) 
        )
      } else if( identical( A, 3 ) ){
        (1 - hs[sel])^(6.5+B) * (6.5 + B + 0.8 * (5.275255128608411+B) * (7.724744871391589+B) * hs[sel] + 
          0.2 * (4.5+B) * (6.5+B) * (8.5+B) * hs[sel]^2) - 
        (6.5+B) * (1 - hs[sel])^(5.5+B) * (1 + (6.5+B) * hs[sel] + 0.4 * (5.275255128608411+B) * 
          (7.724744871391589+B) * hs[sel]^2 + 0.2/3 * (4.5+B) * (6.5+B) * (8.5+B) * hs[sel]^3
        )
      } else {
        stop( "RMgengneiting model undefined for 'n' != 1:3" )
      }
      
      result <- rep( 0., length( hs ) )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        mu = {
          result[sel] <- if( identical( A, 1 ) ){
            (1 - hs[sel])^(2.5+B) * (hs[sel] + (1 + (2.5+B) * hs[sel]) * log(1 - hs[sel]))
          } else if( identical( A, 2 ) ){
            (1 - hs[sel])^(4.5+B) * (hs[sel] + 2/3 * (4.5+B) * hs[sel]^2 + (1 + hs[sel] * (
                  4.5 + B + 6.416666666666666*hs[sel] +  B/3. * (9.+B) * hs[sel]) ) * log(1 - hs[sel])
            )
          } else if( identical( A, 3 ) ){
            (1 - hs[sel])^(6.5 + B)*
            (hs[sel] + (5.2 + 0.8*B)*hs[sel]^2 + 
              0.2*(5.345299461620754 + B)*
              (7.654700538379246 + B)*hs[sel]^3 + 
              (1 + hs[sel]*(6.5 + 1.*B + 
                  0.4*(5.275255128608411 + B)*
                  (7.724744871391589 + B)*hs[sel] + 
                  0.06666666666666667*(4.5 + B)*

                  (6.5 + B)*(8.5 + B)*hs[sel]^2))*
              log(1 - hs[sel]))          
          }
          result
        },
        dgc.dhs * dhs.daniso
      )
      
      
    }, ##  end case Gengneiting
    

    RMgneiting = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < (1. / 0.301187465825)
      dgc.dhs[sel] <- (1. - 0.301187465825*hs[sel])^7 * (
        -1.9957055705418814*hs[sel] -  4.207570523270417*hs[sel]^2 - 2.896611435848653*hs[sel]^3
      )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMgneiting
    
    RMlgd = {
      
      A <- unname( param["alpha"] )
      B <- unname( param["beta"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( A * B * hs^(-1-B) ) / (A+B)
      sel <- hs <= 1.
      dgc.dhs[sel] <- -( A * B * hs[sel]^(-1+A) ) / (A+B)
      
      # dgc.dhs <- ifelse(
      #   hs > 0.
      #   ifelse(
      #     hs <= 1.,
      #     -( A * B * hs^(-1+A) ) / (A+B),
      #     -( A * B * hs^(-1-B) ) / (A+B)
      #   ),
      #   if( identical( A, 1. ) ){
      #     -B / ( B + 1 )
      #   } else if( A < 1. ){
      #     -Inf
      #   } 
      # )
      
      switch(
        d.param,
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
          result[sel] <- -( B * hs[sel]^A * ( -1 + (A+B ) * log( hs[sel] ) ) ) / (A+B)^2
          result
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     ifelse(
        #       hs <= 1.,
        #       -( B * hs^A * ( -1 + (A+B ) * log( hs ) ) ) / (A+B)^2,
        #       B / ( (A+B)^2 * hs^B )              
        #     ),
        #     0.
        #   )
        # },
        beta = {
          result <- -A * ( 1 + (A+B) * log( hs ) ) / ( (A+B)^2 * hs^B )
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
        #       -A * ( 1 + (A+B) * log( hs ) ) / ( (A+B)^2 * hs^B )
        #     ),
        #     0.
        #   )
        # },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMlgd
    
    RMmatern = {
      
      A <- unname( param["nu"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( 
        2^(1.5 - A/2.) * sqrt(A) * ( sqrt(A) * hs )^A * besselK( sqrt(2*A)*hs, -1+A )
      ) / gamma(A)
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -(
      #     ( 2^(1.5 - A/2.) * sqrt(A) * ( sqrt(A) * hs )^A * besselK( sqrt(2) * sqrt(A) * hs , -1+A )
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
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse( 
        #     hs > 0.,
        #     ( 2^(1.5 - A/2.) * ( sqrt(A) * hs )^(1+A) * 
        #       besselK( sqrt(2*A) * hs, A-1) 
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
                  2^(1.-nu) / gamma(nu) * 
                  ( sqrt( 2*nu ) * hs )^nu * besselK( sqrt( 2*nu ) * hs, nu )
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
    
    RMpenta = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- ( 11 * (-1+hs[sel])^5 * hs[sel] * (2+hs[sel]) * ( 4 + hs[sel] * ( 18 + 5 * hs[sel] * (3+hs[sel]) ) ) ) / 6.
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMpenta
    
    RMaskey = {
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- -(A * (1-hs[sel])^(-1+A))         
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        alpha = {
          result <- rep( 0., length( hs ) )
          sel <- hs < 1.
          result[sel] <- ( 1 - hs[sel] )^A * log( 1 - hs[sel] )
          result
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMaskey
    
    RMqexp = {
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- 2 * (-A + exp(hs) ) / ( (-2+A ) * exp(2*hs) )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        alpha = {
          ( 2 * exp( -2*hs ) * ( -1 + exp( hs ) ) ) / (-2+A)^2
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMqexp
    
    
    RMspheric = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- -1.5 + 1.5 * hs[sel]^2          
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case RMspheric
    
    RMstable = {
      
      A <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( ( A * hs^(-1+A) ) / exp(hs^A) )
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( ( A * hs^(-1+A) ) / exp(hs^A) ),
      #   if( A > 1. ){
      #     0.
      #   } else if( identical( A, 1. ) ){
      #     -1.
      #   } else {
      #     -Inf            
      #   }
      # )
      
      switch(
        d.param,
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
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case wave
    
    RMwhittle = {
      
      A <- unname( param["nu"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( 2^(1-A) * hs^A * besselK( hs, -1+A ) ) / gamma(A)
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( 2^(1-A) * hs^A * besselK( hs, -1+A ) ) / gamma(A),
      #   if( A < 0.5 ){
      #     -Inf
      #   } else if( identical( A, 0.5 ) ){
      #     -1.
      #   } else {
      #     0.
      #   }
      # )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ifelse(
        #     hs > 0.,
        #     ( 
        #       2^(1-A) * h * hs^A * besselK( hs, -1+A ) 
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
                  2^(1.-nu) / gamma(nu) * hs^nu * besselK( hs, nu )
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
  
  ## account for spatial nugget
  
  result <- result * unname( param["variance"] / sum( param[c( "variance", "snugget" )] ) )
  
  ##  convert to matrix
  
  result <- list(
    diag = rep( 0., 0.5 * ( 1 + sqrt( 1 + 8 * length( result ) ) ) ),
    tri = result
  )
  attr( result, "struc" ) <- "sym"
  result <- expand( result )
  
  return( result )
  
}

##   ##############################################################################

compute.estimating.equations <- 
  function(
    adjustable.param,
    #     slv,
    envir,
    variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, dpsi.function, 
    tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.reltol,
    force.gradient,
    expectations,
    control.parallel, 
    verbose
  )
{
  
  ## function evaluates the robustified estimating equations of
  ## variogram parameters derived from the Gaussian log-likelihood
  
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-05-06 AP changes for solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  
  ##  get lik.item
  
  lik.item <- prepare.likelihood.calculations(
    envir,
    adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, dpsi.function, tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.reltol,
    compute.xihat = TRUE, compute.Q = FALSE,
    control.parallel = control.parallel,
    verbose = verbose
  )
  
  ##  check whether generalized covariance matrix is positive definite
  
  if( lik.item[["Valpha"]][["error"]] ) {
    if( verbose > 0 ) cat(
      "\n(generalized) correlation matrix Valpha is not positive definite\n"
    )
    t.result <- rep( Inf, length( adjustable.param ) )
    names( t.result ) <- names( adjustable.param )
    return( t.result )
  }
  
  ##  check whether computation of betahat and bhat failed
  
  if( lik.item[["effects"]][["error"]] ) {
    if( verbose > 0 ) cat(
      "\nan error occurred when estimating the fixed and random effects\n"
    )
    t.result <- rep( Inf, length( adjustable.param ) )
    names( t.result ) <- names( adjustable.param )
    return( t.result )
  }
  
  ##  check whether estimating equations should be computed for fixed parameters
  
  if( length( adjustable.param ) == 0 && force.gradient ){
    adjustable.param <- fixed.param
  }
  
  ##  evaluate estimating equations
  
  if( length( adjustable.param ) > 0 ){
    
    ##  compute auxiliary items
    
    TtT <- as.vector( table( TT ) )
    
    ##  compute Cov[bhat]
    
    r.cov <- compute.covariances(
      Valpha.objects = lik.item[["Valpha"]],
      Aalpha = lik.item[["effects"]][["Aalpha"]],
      Palpha = lik.item[["effects"]][["Palpha"]],
      rweights = lik.item[["effects"]][["rweights"]],
      XX = XX, TT = TT, names.yy = names( yy ),
      nugget = lik.item[["param"]]["nugget"],
      eta = lik.item[["eta"]],
      expectations = expectations,
      cov.bhat = TRUE, full.cov.bhat = TRUE,
      cov.betahat = FALSE,
      cov.bhat.betahat = FALSE,
      cov.delta.bhat = FALSE, full.cov.delta.bhat = FALSE,
      cov.delta.bhat.betahat = FALSE,
      cov.ehat = FALSE, full.cov.ehat = FALSE,
      cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
      aux.cov.pred.target = FALSE,
      control.parallel = control.parallel,
      verbose = verbose
    )
    
    if( r.cov[["error"]] ) {
      if( verbose > 0 ) cat(
        "\nan error occurred when computing the covariances of fixed and random effects\n"
      )
      t.result <- rep( Inf, length( adjustable.param ) )
      names( t.result ) <- names( adjustable.param )
      return( t.result )
    }
    
    ##  parallelized computation of estimating equations for all elements
    ##  of adjustable.param
    
    ncores <- min( length( adjustable.param ), control.parallel[["gradient.ncores"]] )
    
    ncores.available <- control.parallel[["max.ncores"]]
    if( sfIsRunning() ) ncores.available <- ncores.available - sfCpus()
    
    control.pmm <- control.parallel
    control.pmm[["pmm.ncores"]] <- min(
      control.pmm[["pmm.ncores"]],
      max( 1L, floor( (ncores.available - ncores) / length( adjustable.param ) ) )
    )
    
    if( ncores > 1L ){
      
      if( !control.parallel[["allow.recursive"]] ) control.pmm[["pmm.ncores"]] <- 1L

      if( identical( .Platform[["OS.type"]], "windows") ){
        
        ## start SNOW cluster if needed
        
        if( length( lik.item[["defaultCluster"]] ) == 0 ){
        
          ## create a SNOW cluster on windows OS
          
          cl <- makePSOCKcluster( ncores, outfile =  "" )
          
          ## export required items to workers
          
#           junk <- clusterExport( cl, c( "dcorr.dparam", "pmm", "expand" ) )
#           junk <- clusterEvalQ( cl, require( snowfall, quietly = TRUE ) )
          junk <- clusterEvalQ( cl, require( georob, quietly = TRUE ) )
          
          lik.item[["defaultCluster"]] <- cl
          
          assign( "lik.item", lik.item, pos = as.environment( envir ) )
          
          save( cl, file = "SOCKcluster.RData" )
          
          options( error = f.stop.cluster )

        } else {
        
          cl <- lik.item[["defaultCluster"]]
          
        }

        t.eeq <- parLapply( 
          cl,
          names( adjustable.param ),
          f.aux.eeq,
          lik.item = lik.item, TtT = TtT, r.cov = r.cov,
          lag.vectors = lag.vectors, variogram.model = variogram.model, 
          control.pmm = control.pmm, verbose = verbose
        )
                
      } else {
        
        t.eeq <- mclapply(
          names( adjustable.param ),
          f.aux.eeq,
          lik.item = lik.item, TtT = TtT, r.cov = r.cov,
          lag.vectors = lag.vectors, variogram.model = variogram.model, 
          control.pmm = control.pmm, verbose = verbose,
          mc.cores = ncores
        )
        
      }
      
    } else {
      
      t.eeq <- lapply(
        names( adjustable.param ),
        f.aux.eeq,
        lik.item = lik.item, TtT = TtT, r.cov = r.cov,
        lag.vectors = lag.vectors, variogram.model = variogram.model, 
        control.pmm = control.pmm, verbose = verbose
      )
      
    }
    
    t.eeq <- simplify2array( t.eeq )
     
    colnames( t.eeq ) <- names( adjustable.param )
    
    eeq.exp <- t.eeq["eeq.exp", ]
    eeq.emp <- t.eeq["eeq.emp", ]
    
      
    
    if( verbose > 1 ) {
      cat( "\n                      ",
        format( names( eeq.emp), width = 14, justify = "right" ), 
        "\n", sep =""
      )
      cat( "  EEQ                :", 
        format( 
          signif( eeq.emp / eeq.exp - 1, digits = 7 ), 
          scientific = TRUE, width = 14
        ), "\n", sep = "" 
      )
      if( verbose > 2 ){
        cat( "      empirical terms:", 
          format( 
            signif( eeq.emp, digits = 7 ), 
            scientific = TRUE, width = 14
          ), "\n", sep = "" 
        )
        cat( "      expected  terms:", 
          format( 
            signif( eeq.exp, digits = 7 ), 
            scientific = TRUE, width = 14
          ), "\n", sep = ""
        )
      }
      cat("\n")
    }
    
    ##  store terms in lik.item object
    
    lik.item[["eeq"]] <- list(
      eeq.emp = eeq.emp,
      eeq.exp = eeq.exp
    )
    
    assign( "lik.item", lik.item, pos = as.environment( envir ) )
    
    #     if( slv ){
    
    return( eeq.emp / eeq.exp - 1. )
    
    #     } else {
    #       res <- sum( (eeq.emp / eeq.exp - 1.)^2 )
    #       if( verbose > 1 ) cat( 
    #         "  sum(EEQ^2)         :",
    #         format( 
    #           signif( res, digits = 7 ), 
    #           scientific = TRUE, width = 14
    #         ), "\n", sep = "" 
    #       )
    #       
    #       return( res )
    #     }
    
  } else {
    
    ##  all parameters are fixed
    
    return( NA_real_ )
    
  }
  
}


##   ##############################################################################

## compute.expanded.estimating.equations <- 
##   function(
##     allpar,
##     slv,
##     envir,
##     variogram.model, fixed.param, param.name, aniso.name,
##     param.tf, bwd.tf, safe.param,
##     lag.vectors,
##     XX, min.condnum, rankdef.x, yy, TT, 
##     psi.function, dpsi.function, 
##     tuning.psi, tuning.psi.nr, 
##     irwls.initial, irwls.maxiter, irwls.reltol,
##     force.gradient,
##     expectations,
##     verbose
##   )
## {
##   
##   ## function evaluates the robustified estimating equations of
##   ## variogram parameters derived from the Gaussian log-likelihood
##   
##   ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
##   
##   ## select xihat and variogram parameters
##   
##   xihat <- allpar[ 1:NROW(XX) ]
##   adjustable.param <- allpar[ -(1:NROW(XX)) ]
## 
##   ##  get lik.item
##   
##   lik.item <- prepare.likelihood.calculations(
##     envir,
##     adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
##     param.tf, bwd.tf, safe.param,
##     lag.vectors,
##     XX, min.condnum, rankdef.x, yy, TT, xihat, 
##     psi.function, dpsi.function, tuning.psi, tuning.psi.nr, ml.method, 
##     irwls.initial, irwls.maxiter, irwls.reltol,
##     compute.xihat = FALSE, compute.Q = FALSE,
##     control.parallel = control.parallel,
##     verbose
##   )
##   
##   ##  check whether generalized covariance matrix is positive definite
##   
##   if( lik.item[["Valpha"]][["error"]] ) {
##     if( verbose > 0 ) cat(
##       "\n(generalized) correlation matrix Valpha is not positive definite\n"
##     )
##     t.result <- rep( Inf, length( adjustable.param ) )
##     names( t.result ) <- names( adjustable.param )
##     return( t.result )
##   }
##     
##   ##  check whether estimating equations should be computed for fixed parameters
##   
##   if( length( adjustable.param ) == 0 && force.gradient ){
##     adjustable.param <- fixed.param
##   }
##   
##   ##  evaluate estimating equations
##   
##   ##  compute auxiliary items
##   
##   TtT <- as.vector( table( TT ) )
##   
##   ##  compute Cov[bhat]
##   
##   r.cov <- compute.covariances(
##     Valpha.objects = lik.item[["Valpha"]],
##     Aalpha = lik.item[["effects"]][["Aalpha"]],
##     Palpha = lik.item[["effects"]][["Palpha"]],
##     rweights = lik.item[["effects"]][["rweights"]],
##     XX = XX, TT = TT, names.yy = names( yy ),
##     nugget = lik.item[["param"]]["nugget"],
##     eta = lik.item[["eta"]],
##     expectations = expectations,
##     cov.bhat = TRUE, full.cov.bhat = TRUE,
##     cov.betahat = FALSE,
##     cov.bhat.betahat = FALSE,
##     cov.delta.bhat = FALSE, full.cov.delta.bhat = FALSE,
##     cov.delta.bhat.betahat = FALSE,
##     cov.ehat = FALSE, full.cov.ehat = FALSE,
##     cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
##     aux.cov.pred.target = FALSE,
##     control.parallel = control.parallel,
##     verbose = verbose
##   )
##   
##   if( r.cov[["error"]] ) {
##     if( verbose > 0 ) cat(
##       "\nan error occurred when computing the covariances of fixed and random effects\n"
##     )
##     t.result <- rep( Inf, length( adjustable.param ) )
##     names( t.result ) <- names( adjustable.param )
##     return( t.result )
##   }
##   
##   ## estimating equations for xihat
##   
##   eeq.xihat <- estimating.eqations.xihat(
##     res = lik.item[["effects"]][["residuals"]],
##     TT = TT, xihat = xihat, 
##     nugget = lik.item[["param"]]["nugget"],
##     eta = lik.item[["eta"]],
##     Valpha.inverse.Palpha = lik.item[["effects"]][["Valpha.inverse.Palpha"]],
##     psi.function = psi.function, 
##     tuning.psi = tuning.psi  
##   )
##   
##   ##  initialize estimating equations for variogram parameters
##   
##   eeq.emp <- rep( NA, length( adjustable.param ) )
##   names( eeq.emp ) <- names( adjustable.param )
##   
##   eeq.exp <- rep( NA, length( adjustable.param ) )
##   names( eeq.exp ) <- names( adjustable.param )
##   
##   ##  estimation equation for nugget
##   
##   if( "nugget" %in% names( adjustable.param ) ) {
##     
##     ##  compute trace of Cov[ psi( residuals/sqrt(nugget) ) ]
##     
##     eeq.exp["nugget"] <- sum( 
##       diag( 
##         lik.item[["Valpha"]][["Valpha.inverse"]] %*%             
##         ( 1/TtT * lik.item[["Valpha"]][["Valpha.inverse"]] ) %*% 
##         r.cov[["cov.bhat"]] 
##       ) 
##     )
##     eeq.emp["nugget"] <- sum( 
##       ( lik.item[["effects"]][["z.star"]] )^2 / TtT
##     )
##     
##   }
##   
##   ##  estimation equation for spatial nugget
##   
##   if( "snugget" %in% names( adjustable.param ) ) {
##     
##     ##  compute trace( Valpha^-1 Cov[bhat] )
##     
##     eeq.exp["snugget"] <- sum(
##       rowSums( 
##         (lik.item[["Valpha"]][["Valpha.inverse"]] %*% lik.item[["Valpha"]][["Valpha.inverse"]] ) * 
##         r.cov[["cov.bhat"]]
##       )
##     )
##     eeq.emp["snugget"] <- sum( lik.item[["effects"]][["z.star"]]^2 )
##     
##   }
##   
##   ##  estimation equation for variance
##   
##   if( "variance" %in% names( adjustable.param ) ) {
##     
##     ##  compute trace( Valpha^-1 Cov[bhat] )
##     
##     eeq.exp["variance"] <- sum(
##       rowSums( 
##         ( lik.item[["Valpha"]][["Valpha.inverse"]] %*% lik.item[["Valpha"]][["Valpha0"]] %*% lik.item[["Valpha"]][["Valpha.inverse"]] ) * 
##         r.cov[["cov.bhat"]]
##       )
##     )
##     eeq.emp["variance"] <- sum( 
##       lik.item[["effects"]][["z.star"]] * drop( lik.item[["Valpha"]][["Valpha0"]] %*% lik.item[["effects"]][["z.star"]] )
##     )
##     
##   }
##   
##   ##  estimation equations for scale, extra variogram and anisotropy
##   ##  parameters
##   
##   extra.par <- names( adjustable.param )[ !( 
##     names( adjustable.param ) %in% c( "variance", "snugget", "nugget" )
##   )]
##   
##   for( t.i in extra.par ){
##     
##     ##  compute trace( Valpha^-1 * dValpha/dalpha * Valpha^-1 * Cov[bhat] )
##     
##     dValpha <- dcorr.dparam(
##       x = lag.vectors, variogram.model = variogram.model, param = lik.item[["param"]], 
##       d.param = t.i,
##       aniso = lik.item[["aniso"]],
##       verbose = verbose
##     )
##     ##       if( identical( class( dValpha ), "try-error" ) ){
##     ##         if( verbose > 0 ) cat( "error in dcorr.dparam\n\n" )
##     ##         t.result <- rep( Inf, length( adjustable.param ) )
##     ##         names( t.result ) <- names( adjustable.param )
##     ##         return( t.result )
##     ##       }
##     
##     eeq.exp[t.i] <- sum(
##       rowSums( 
##         (lik.item[["Valpha"]][["Valpha.inverse"]] %*% dValpha %*% lik.item[["Valpha"]][["Valpha.inverse"]]) * 
##         r.cov[["cov.bhat"]]
##       )
##     )
##     eeq.emp[t.i] <- sum( 
##       lik.item[["effects"]][["z.star"]] * drop( dValpha %*% lik.item[["effects"]][["z.star"]] )
##     )
##     
##   }
##   
##   if( verbose > 1 ) {
##     cat( "\n                      ",
##       format( c( "min(xihat)", "max(xihat)" ), width = 14, justify = "right" ), 
##       "\n", sep =""
##     )
##     cat( "  EEQ                :", 
##       format( 
##         signif( range(eeq.xihat), digits = 7 ), 
##         scientific = TRUE, width = 14
##       ), "\n", sep = "" 
##     )
##     cat( "\n                      ",
##       format( names( eeq.emp), width = 14, justify = "right" ), 
##       "\n", sep =""
##     )
##     cat( "  EEQ                :", 
##       format( 
##         signif( eeq.emp / eeq.exp - 1, digits = 7 ), 
##         scientific = TRUE, width = 14
##       ), "\n", sep = "" 
##     )
##     if( verbose > 2 ){
##       cat( "      empirical terms:", 
##         format( 
##           signif( eeq.emp, digits = 7 ), 
##           scientific = TRUE, width = 14
##         ), "\n", sep = "" 
##       )
##       cat( "      expected  terms:", 
##         format( 
##           signif( eeq.exp, digits = 7 ), 
##           scientific = TRUE, width = 14
##         ), "\n", sep = ""
##       )
##     }
##     cat("\n")
##   }
##   
##   ##  store terms in lik.item object
##   
##   lik.item[["eeq"]] <- list(
##     eeq.xihat = eeq.xihat,
##     eeq.emp = eeq.emp,
##     eeq.exp = eeq.exp
##   )
##   
##   assign( "lik.item", lik.item, pos = as.environment( envir ) )
##   
##   return( c( eeq.xihat, eeq.emp / eeq.exp - 1. ) )
##   
## }


##   ##############################################################################

negative.restr.loglikelihood <- 
  function(
    adjustable.param,
    envir,
    variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, dpsi.function, 
    tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.reltol,
    control.parallel, 
    verbose,
    ...
  )
{
  
  ## function computes laplace approximation to negative (un)restricted
  ## loglikelihood 
  
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2013-06-03 AP changes for estimating xihat
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  
  #     sel <- !c( param.name, aniso.name ) %in% names( fixed.param )
  #     names( adjustable.param ) <- c( param.name, aniso.name )[sel]
  
  ##  compute required items (param, eta, Valpha.inverse, Valpha.ilcf, 
  ##  betahat, bhat, residuals, etc.)
  
  lik.item <- prepare.likelihood.calculations(
    envir,
    adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, dpsi.function, tuning.psi, tuning.psi.nr, ml.method, 
    irwls.initial, irwls.maxiter, irwls.reltol,
    compute.xihat = TRUE, compute.Q = TRUE,
    control.parallel = control.parallel,
    verbose
  )
  
  ##  check whether generalized covariance matrix is positive definite
  
  if( lik.item[["Valpha"]][["error"]] ) {
    if( verbose > 0 ) cat(
      "\n(generalized) correlation matrix Valpha is not positive definite\n"
    )
    return( NA )
  }
  
  ##  check whether computation of betahat and bhat failed
  
  if( lik.item[["effects"]][["error"]] ) {
    if( verbose > 0 ) cat(
      "\nan error occurred when estimating the fixed and random effects\n"
    )
    return( NA )
  }
  
  ##  check whether Q matrix not positive definite
  
  if( lik.item[["Q"]][["error"]] ) {
    if( verbose > 0 ) cat(
      "\nan error occurred when determinants required for",
      "Gaussian log-likelihood were computed\n"
    )
    return( NA )
  }
  
  ##  compute laplace approximation of negative (un)-restricted
  ##  (profile-)loglikelihood
  
  t.dim <- dim( XX )
  t.c <- if( identical( ml.method, "REML" ) ) diff(t.dim) else -t.dim[1]
  
  term1 <- -0.5 * (
    t.c * log( 2 * pi ) + t.dim[1] * log( 1./lik.item[["eta"]] ) 
  ) + t.dim[1] * log( lik.item[["param"]]["nugget"] ) +
  sum( log( diag( lik.item[["Valpha"]][["Valpha.ucf"]] ) ) ) 
  
  Ttpsi <- lik.item[["effects"]][["residuals"]] / sqrt( lik.item[["param"]]["nugget"] ) *
    lik.item[["effects"]][["rweights"]]
  TtT   <- rep( 1., length( Ttpsi ) )
  if( sum( duplicated( TT ) > 0 ) ){
    Ttpsi <- as.vector( tapply( Ttpsi, factor( TT ), sum ) )
    TtT   <- as.vector( table( TT ) )
  }
  
  term2 <- 0.5 * ( sum( Ttpsi^2 / TtT ) + sum( 
      lik.item[["effects"]][["z.star"]] * lik.item[["effects"]][["bhat"]] 
    ) / lik.item[["param"]]["nugget"] / lik.item[["eta"]]
  )
  attributes( term2 ) <- NULL
  
  term3 <- 0.5 * lik.item[["Q"]][["log.det.Q"]]
  
  r.neg.restricted.loglik <- term1 + term2 + term3
  
  attributes( r.neg.restricted.loglik ) <- NULL
  
  if( verbose > 1 ) cat(
    "\n  Negative. restrict. loglikelihood:", 
    format( 
      signif( r.neg.restricted.loglik, digits = 7 ), 
      scientific = TRUE, width = 14
    ), "\n", sep = ""
  )
  
  return( r.neg.restricted.loglik )
  
}


##   ##############################################################################

gradient.negative.restricted.loglikelihood <- 
  function(
    adjustable.param,
    envir,
    variogram.model, fixed.param, param.name, aniso.name,
    param.tf, deriv.fwd.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, dpsi.function, d2psi.function, 
    tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.reltol,
    force.gradient,
    control.parallel,
    verbose
  )
{
  
  ##  function computes gradient of Laplace approximation of negative 
  ##  restricted log-likelihood with respect to covariance parameters
  
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-04 AP correction of values returned on error
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation

  ##   dtrafo.fct <- list(
  ##     log      = function( x ) 1/x,
  ##     identity = function( x ) rep( 1, length( x ) )
  ##   )
  
  ##  get lik.item
  
  lik.item <- prepare.likelihood.calculations(
    envir,
    adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, dpsi.function, tuning.psi, tuning.psi.nr, ml.method, 
    irwls.initial, irwls.maxiter, irwls.reltol,
    compute.xihat = TRUE, compute.Q = TRUE,
    control.parallel = control.parallel,
    verbose
  )
  
  ##  check whether generalized covariance matrix is positive definite
  
  if( lik.item[["Valpha"]][["error"]] ) {
    if( verbose > 0 ) cat(
      "\n(generalized) correlation matrix Valpha is not positive definite\n"
    )
    return( rep( NA, length( adjustable.param ) ) )
  }
  
  ##  check whether computation of betahat and bhat failed
  
  if( lik.item[["effects"]][["error"]] ) {
    if( verbose > 0 ) cat(
      "\nan error occurred when estimating the fixed and random effects\n"
    )
    return( rep( NA, length( adjustable.param ) ) )
  }
  
  ##  check whether Q matrix not positive definite
  
  if( lik.item[["Q"]][["error"]] ) {
    if( verbose > 0 ) cat(
      "\nan error occurred when determinants required for ",
      "Gaussian log-likelihood were computed\n"
    )
    return( rep( NA, length( adjustable.param ) ) )
  }
  
  ##  check whether gradient should be computed for fixed parameters
  
  if( length( adjustable.param ) == 0 && force.gradient ){
    adjustable.param <- fixed.param
  }
  
  ##  evaluate gradient
  
  if( length( adjustable.param ) > 0 ){
    
    ##  compute auxiliary items
    
    n <- nrow( XX )
    
    std.res <- lik.item[["effects"]][["residuals"]] / sqrt( lik.item[["param"]]["nugget"] )
    
    Qi <- lik.item[["Q"]][["Q.inverse"]]
    Vi <- lik.item[["Valpha"]][["Valpha.inverse"]] / 
      ( lik.item[["param"]]["snugget"] + lik.item[["param"]]["variance"] )
    TtT <- as.vector( table( TT ) )
    
    ## parallelized computation of gradient
    
    nme.adj.param <- names( adjustable.param )
    nme.var.snug <- nme.adj.param[nme.adj.param %in% c( "variance", "snugget" )]
    nme.adj.param[nme.adj.param %in% c( "snugget", "variance" )] <- "Vparam"
    nme.adj.param <- unique( nme.adj.param )
    
    ncores <- min( length( nme.adj.param ), control.parallel[["gradient.ncores"]] )

    ncores.available <- control.parallel[["max.ncores"]]
    if( sfIsRunning() ) ncores.available <- ncores.available - sfCpus()

    control.pmm <- control.parallel
    control.pmm[["pmm.ncores"]] <- min(
      control.pmm[["pmm.ncores"]],
      max( 1L, floor( (ncores.available - ncores) / length( nme.adj.param ) ) )
    )
    
    if( ncores > 1L && !control.parallel[["allow.recursive"]] ) control.pmm[["pmm.ncores"]] <- 1L
    
    if( ncores > 1L ){
      
      if( !control.parallel[["allow.recursive"]] ) control.pmm[["pmm.ncores"]] <- 1L
      
      if( identical( .Platform[["OS.type"]], "windows") ){
        
        ## start SNOW cluster if needed
        
        if( length( lik.item[["defaultCluster"]] ) == 0 ){
        
          ## create a SNOW cluster on windows OS
          
          cl <- makePSOCKcluster( ncores, outfile =  "" )
          
          ## export required items to workers
          
#           junk <- clusterExport( 
#             cl, 
#             c( "dcorr.dparam", "pmm", "expand", "f.stop.cluster" ) 
#           )
#           junk <- clusterEvalQ( cl, require( snowfall, quietly = TRUE ) )
          junk <- clusterEvalQ( cl, require( georob, quietly = TRUE ) )
          
          lik.item[["defaultCluster"]] <- cl
          
          assign( "lik.item", lik.item, pos = as.environment( envir ) )
          
          save( cl, file = "SOCKcluster.RData" )
          
          options( error = f.stop.cluster )
          
        } else {
        
          cl <- lik.item[["defaultCluster"]]
          
        }
        
        r.gradient <- parLapply( 
          cl,
          nme.adj.param,
          f.aux.gradient.nll,
          nme.var.snug = nme.var.snug,
          TT = TT, TtT = TtT, std.res = std.res, XX = XX, 
          Qi = Qi, Vi = Vi,
          n = n, 
          lag.vectors = lag.vectors, variogram.model = variogram.model,
          lik.item = lik.item,
          param.tf = param.tf, deriv.fwd.tf = deriv.fwd.tf, 
          ml.method = ml.method,
          control.pmm = control.pmm, verbose = verbose
        )
        
      } else {
        
        r.gradient <- mclapply(
          nme.adj.param,
          f.aux.gradient.nll,
          nme.var.snug = nme.var.snug,
          TT = TT, TtT = TtT, std.res = std.res, XX = XX, 
          Qi = Qi, Vi = Vi,
          n = n, 
          lag.vectors = lag.vectors, variogram.model = variogram.model,
          lik.item = lik.item,
          param.tf = param.tf, deriv.fwd.tf = deriv.fwd.tf, 
          ml.method = ml.method,
          control.pmm = control.pmm, verbose = verbose,
          mc.cores = ncores
        )
        
      }
      
    } else {
      
      r.gradient <- lapply(
        nme.adj.param,
        f.aux.gradient.nll,
        nme.var.snug = nme.var.snug,
        TT = TT, TtT = TtT, std.res = std.res, XX = XX, 
        Qi = Qi, Vi = Vi,
        n = n, 
        lag.vectors = lag.vectors, variogram.model = variogram.model,
        lik.item = lik.item,
        param.tf = param.tf, deriv.fwd.tf = deriv.fwd.tf, 
        ml.method = ml.method,
        control.pmm = control.pmm, verbose = verbose
      )
      
    }
    
    r.gradient <- unlist( r.gradient )
    names( r.gradient ) <- names( adjustable.param )
        
    ##  rearrange elements of gradient and change sign (for negative
    ##  log-likelihood)
    
    r.gradient <- -r.gradient[names( adjustable.param )]
    
    if( verbose > 1 ){
      cat( "\n                      ",
        format( names( r.gradient ), width = 14, justify = "right" ), 
        "\n", sep = ""
      )
      cat( "  Gradient           :", 
        format( 
          signif( r.gradient, digits = 7 ), 
          scientific = TRUE, width = 14
        ), "\n" , sep = ""
      )
    }
    
    return( r.gradient )
    
  } else {
    
    ##  all parameters are fixed
    
    return( NA_real_ )
    
  }
}


##  ##   ##############################################################################
##      
##      f.compute.df <- function( Valpha, XX, param ){
##          
##          ##  function computes three estimates of the degrees of freedom of
##          ##  the smoothing universal kriging predictor, cf.  Hastie &
##          ##  Tibshirani, 1990, Generalized additive models, pp.52
##          
##          ##  2011-07-05
##          ##  Andreas Papritz
##          
##          sigma <- param["variance"] * Valpha
##          diag( sigma ) <- diag( sigma ) + param["nugget"]
##          
##          ##  compute inverse lower cholesky factor of covariance matrix of
##          ##  data
##          
##          ilcf <- t( backsolve( chol( sigma ), diag( nrow( Valpha ) ), k = nrow( Valpha ) ) )
##          
##          ##  compute hat matrix
##          
##          q <- qr.Q( qr( xtilde <- ilcf %*% XX ) )
##          s <- -tcrossprod( q )
##          
##          diag( s ) <- diag( s ) + 1
##          s <- -param["nugget"] * t( ilcf ) %*% s %*% ilcf
##          diag( s ) <- diag( s ) + 1
##          
##          ##  compute degrees of freedom
##          
##          df.1 <- sum( diag( s ) )
##          df.3 <- sum( s^2 )
##          df.2 <- 2 * df.1 - df.3
##          
##          return( 
##              c( 
##                  df.SSt    = t.df.2 <- sum( s^2 ), 
##                  df.S      = t.df.1 <- sum( diag( s ) ), 
##                  df.2SmSSt = 2 * t.df.1 - t.df.2
##              ) 
##          )
##              
##      }

##   ##############################################################################

georob.fit <- 
  function(
    ## root.finding,
    #     slv,
    #     envir,
    initial.objects,
    variogram.model, param, fit.param,
    aniso, fit.aniso,
    param.tf, 
    fwd.tf, 
    deriv.fwd.tf, 
    bwd.tf,
    georob.object,
    safe.param,
    tuning.psi, 
    cov.bhat, full.cov.bhat,
    cov.betahat, 
    cov.bhat.betahat,
    cov.delta.bhat, full.cov.delta.bhat,
    cov.delta.bhat.betahat,
    cov.ehat, full.cov.ehat,
    cov.ehat.p.bhat, full.cov.ehat.p.bhat,
    aux.cov.pred.target,
    min.condnum, rankdef.x,
    psi.func,
    tuning.psi.nr,
    ml.method,
    irwls.initial,
    irwls.maxiter, 
    irwls.reltol, 
    force.gradient,
    zero.dist,
    control.nleqslv,
    ## bbsolve.method, bbsolve.control,
    control.optim, hessian,
    control.parallel,
    full.output,
    verbose
  )
{
  
  ## 2011-06-24 ap
  ## 2011-06-24 cs
  ## 2011-06-29 ap, cs
  ## 2011-07-22 ap
  ## 2011-07-28 ap
  ## 2011-08-12 ap
  ## 2011-10-14 ap
  ## 2011-12-19 ap
  ## 2011-12-22 ap
  ## 2011-12-23 AP modified for estimating variogram model with spatial
  ##               nugget (micro-scale variation)
  ## 2012-02-07 AP modified for geometrically anisotropic variograms
  ## 2012-02-20 AP replacement of ifelse
  ## 2012-02-27 AP rescaled rho-, psi-function etc.
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-21 AP arguments lower, upper passed to optim
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2012-11-27 AP changes in check allowed parameter range
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-06-03 AP handling design matrices with rank < ncol(x)
  ## 2013-05-06 AP changes for solving estimating equations for xi
  ## 2013-06-12 AP changes in stored items of Valpha object
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-03 AP new transformation of rotation angles
  ## 2013-07-09 AP catching errors occuring when fitting anisotropic
  ##               variograms with default anisotropy parameters
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-02-18 AP correcting error when fitting models with offset
  ## 2014-05-28 AP change in check for initial variogram parameter values
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  
  ##  ToDos:
  
  ##  main body of georob.fit
  
  ##  define rho-function and derivatives
  
  rho.psi.etc <- switch(
    psi.func,
    t.dist = list(
      rho.function = function( x, tuning.psi ){
        return( tuning.psi / 2 * log( ( 1 + ( x^2 / tuning.psi ) ) ) )
      },            
      psi.function = function( x, tuning.psi ){
        return( tuning.psi * x / ( tuning.psi + x^2 ) )
      },
      dpsi.function = function( x, tuning.psi ) {
        return( tuning.psi * ( tuning.psi - x^2 ) / ( tuning.psi + x^2 )^2 )
      },
      d2psi.function = function( x, tuning.psi ) {
        return( 
          2 * tuning.psi * x * ( x^2 - 3 * tuning.psi ) / 
          ( tuning.psi + x^2 )^3
        )
      }
    ),
    logistic = list(
      rho.function = function( x, tuning.psi ) {
        return( 
          tuning.psi * (-x + tuning.psi * 
            ( -log(2) + log( 1 + exp(( 2 * x ) / tuning.psi ) ) )
          )
        )
      }, 
      psi.function = function( x, tuning.psi ) {
        t.x <- exp(-(2*x)/tuning.psi)
        return( (2*tuning.psi / (1 + t.x) - tuning.psi) )
      }, 
      dpsi.function = function( x, tuning.psi ) {
        t.x <- exp(-(2*x)/tuning.psi)
        t.result <- ( 4 * t.x ) / ( 1 + t.x )^2
        t.result[is.nan(t.result)] <- 0.
        return( t.result )
      }, 
      d2psi.function = function( x, tuning.psi ) {
        t.x <- exp(-(2*x)/tuning.psi)
        t.result <- ( ( 16*t.x^2 / (1+t.x)^3 ) - ( 8*t.x / (1+t.x)^2 ) ) / tuning.psi
        t.result[is.nan(t.result)] <- 0.
        return( t.result )
      }
    ),
    huber = list(
      rho.function <- function( x, tuning.psi ) {
        ifelse( 
          abs( x ) <= tuning.psi, 
          0.5 * x^2, 
          tuning.psi * abs( x ) - 0.5 * tuning.psi^2 
        )
      },
      psi.function <- function( x, tuning.psi ) {
        ifelse( abs( x ) <= tuning.psi, x, sign(x) * tuning.psi )
      },
      dpsi.function <- function( x, tuning.psi ) {
        ifelse( abs( x ) <= tuning.psi, 1, 0 )
      },
      d2psi.function = function( x, tuning.psi ) {
        rep( 0, length( x ) )
      }
    )
  )
  
  ##  set number of IRWLS iterations for estimating bhat and betahat to
  ##  1 for non-robust REML case
  
  if( psi.func %in% c( "logistic", "huber" ) & tuning.psi >= tuning.psi.nr ){
    irwls.maxiter <- 1        
  }
  
  ##  copy items of initial.objects to local environment
  
  XX          <- initial.objects[["x"]]
  yy          <- initial.objects[["y"]]
  betahat     <- coefficients( initial.objects[["initial.fit"]] )
  bhat        <- initial.objects[["bhat"]]
  coordinates <- initial.objects[["locations.objects"]][["coordinates"]]
  
  ##  check for multiple observations at same location and generate
  ##  designmatrix of replicated observations
  
  dist0 <- as.matrix( dist( coordinates ) ) <= zero.dist
  first.dist0 <- unname( apply( dist0, 1, function( x ) ( (1:length(x))[x])[1] ) )
  
  
  TT <- matrix( 0, nrow = length( yy ), ncol = length( yy ) )
  TT[ cbind( 1:nrow(TT), first.dist0 ) ] <- 1
  rep.obs <- (1:ncol(TT))[ apply( TT, 2, function( x ) all( x == 0 ) ) ]
  if( length( rep.obs ) > 0 )  TT <- TT[, -rep.obs]
  
  ## check whether explanatory variables are the identical for the replicated
  ## observations and issue an error if not
  
  apply( 
    TT, 
    2, 
    function( i, XX ){
      XX <- XX[as.logical(i), , drop = FALSE]
      apply( 
        XX, 
        2, 
        function( x ){
          if( length(x) > 1 && any( x[-1] != x[1] ) ) warning(
            "explanatory variables differ for some replicated observations" 
          )
        }
      )    
    },
    XX = XX
  )
  
  ## store row indices of replicated observations only
  
  TT <- drop( TT %*% 1:ncol( TT ) )
  
  ##  omit elements corresponding to replicated observations in XX, bhat
  ##  and coordinates
  
  if( length( rep.obs ) > 0 ) {
    XX          <- XX[ -rep.obs, , drop = FALSE]
    bhat      <- bhat[ -rep.obs ]
    coordinates <- coordinates[ -rep.obs, , drop = FALSE]
    if( verbose > 0 ) cat( "\n", length(rep.obs), "replicated observations at", 
      length( unique( TT[rep.obs] ) ), "sampling locations\n" 
    )
  }
  
  ## compute lag vectors for all pairs of coordinates
  
  if( 
    !is.null( georob.object ) && 
    isTRUE( all.equal( 
        georob.object[["locations.objects"]][["coordinates"]],
        coordinates
      )
    )
  ){
    lag.vectors <- georob.object[["locations.objects"]][["lag.vectors"]]
  } else {
    indices.pairs <- combn( NROW( coordinates ), 2 )
    lag.vectors <- coordinates[ indices.pairs[2,], ] - coordinates[ indices.pairs[1,], ]
  }
  
  ## set snugget to zero if snugget has not been specified or if there are
  ## no replicated observations
  
  if( !"snugget" %in% names( param ) | sum( duplicated( TT ) ) == 0 ){
    param["snugget"] <- 0.
    fit.param["snugget"] <- FALSE
  }
  
  ##  check whether fitting of chosen variogram model is implemented and
  ##  return names of extra parameters (if any)
  
  ep <- param.names( model = variogram.model )
  
  ## check names of initial variogram parameters and flags for fitting
  
  param.name <- c( "variance", "snugget", "nugget", "scale", ep )
  
  if( !all( names( param ) %in% param.name ) ) stop( 
    "error in names of initial values of variogram parameters" 
  )
  
  if( !all( param.name  %in% names( param ) ) ) stop( 
    "no initial values provided for parameter(s) '", 
    paste( param.name[ !param.name %in% names( param ) ], collapse= ", "), "'"
  )
  
  if( !all( names( fit.param ) %in% param.name ) ) stop( 
    "error in names of control flags for fitting variogram parameters" 
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
  
  param <- param[param.name]
  
  ## check whether intitial values of variogram parameters are valid
  
  if( param["variance"] < 0. ) stop("initial value of 'variance' must be positive" )
  if( param["snugget"] < 0. )  stop("initial value of 'snugget' must be positive" )
  if( param["nugget"] < 0. ) stop("initial value of 'nugget' must be positive" )
  if( param["scale"] <= 0. ) stop("initial value of 'scale' must be positive" )
  
  param.bounds <- param.bounds( variogram.model, NCOL( coordinates ) )
  ep.param <- param[ep]
  
  if( !is.null( param.bounds ) ) t.bla <- sapply(
    1:length( ep.param ),
    function( i, param, bounds ){
      if( param[i] < bounds[[i]][1] || param[i] > bounds[[i]][2] ) stop(
        "initial value of parameter '", names( param[i] ), "' outside of allowed range" 
      )
    }, 
    param = ep.param,
    bounds = param.bounds
  )
  
  
  ##  rearrange and check flags controlling variogram parameter fitting 
  
  fit.param <- fit.param[param.name]
  
  if( 
    variogram.model %in% (t.models <- c( "RMfbm" ) ) && 
    ( 
      sum( duplicated( TT ) > 0 ) && all( 
        fit.param[c( "variance", "snugget", "scale" ) ] 
      ) ||
      sum( duplicated( TT ) == 0 ) && all( 
        fit.param[c( "variance", "scale" ) ] 
      ) 
    )
  ) stop( 
    "'variance', 'scale' (and 'snugget') cannot be fitted simultaneously for variograms ",
    paste( t.models, collapse = " or "), "; \n  'scale' parameter must be fixed"
  )
  
  ##  preparation for variogram parameter transformations
  
  all.param.tf <- param.tf
  
  t.sel <- match( param.name, names( all.param.tf ) )
  
  if( any( is.na( t.sel ) ) ){
    stop( "transformation undefined for some variogram parameters" )
  } else {
    param.tf <- all.param.tf[t.sel]
  }
  
  ##  transform initial variogram parameters
  
  transformed.param <- sapply(
    param.name,
    function( x, param.tf, param ) fwd.tf[[param.tf[x]]]( param[x] ),
    param.tf = param.tf,
    param = param
  )
  
  names( transformed.param ) <- param.name 
  
  ## check names of initial anisotropy parameters and flags for fitting
  
  aniso.name <- c( "f1", "f2", "omega", "phi", "zeta" )
  
  if( !all( names( aniso ) %in% aniso.name ) ) stop( 
    "error in names of initial values of anisotropy parameters" 
  )
  
  if( !all( aniso.name  %in% names( aniso ) ) ) stop( 
    "no initial values provided for parameter(s) '", 
    aniso.name[ !aniso.name %in% names( aniso ) ], "'"
  )
  
  if( !all( names( fit.aniso ) %in% aniso.name ) ) stop( 
    "error in names of control flags for fitting  anisotropy parameters"
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
  
  ## adjust default initial values of anisotropy parameters if these are
  ## fitted
  
  if( fit.aniso["omega"] && aniso["f1"] == 1. ) aniso["f1"] <- aniso["f1"] - sqrt( .Machine$double.eps )
  if( fit.aniso["phi"] ){
    if( aniso["f1"] == 1. ) aniso["f1"] <- aniso["f1"] - 0.0001
    if( aniso["f2"] == 1. ) aniso["f2"] <- aniso["f2"] - 0.0001
  }
  if( fit.aniso["zeta"] && aniso["f2"] == 1. ) aniso["f2"] <- aniso["f2"] - 0.0001
  
  ##  rearrange and check flags controlling anisotropy parameter fitting 
  
  fit.aniso <- fit.aniso[aniso.name]
  
  ##  preparation for anisotropy parameter transformations
  
  t.sel <- match( aniso.name, names( all.param.tf ) )
  
  if( any( is.na( t.sel ) ) ){
    stop( "transformation undefined for some anisotropy parameters" )
  } else {
    aniso.tf <- all.param.tf[t.sel]
  }
  
  #   if( !all( aniso.tf %in% c( "log", "identity" ) ) ) stop(
  #     "undefined transformation of anisotropy parameter"
  #   )
  
  ##  transform initial anisotropy parameters
  
  transformed.aniso <- sapply(
    aniso.name,
    function( x, param.tf, param ){
      fwd.tf[[param.tf[x]]]( param[x] )
    },
    param.tf = aniso.tf,
    param = aniso
  )
  names( transformed.aniso ) <- aniso.name 
  
  param.tf <- c( param.tf, aniso.tf )
  
  ##  create environment to store items required to compute likelihood and
  ##  estimating equations that are provided by
  ##  prepare.likelihood.calculations
  
  envir <- new.env()
  lik.item <- list()
  
  ##  initialize values of variogram parameters stored in the environment
  
  #   lik.item <- get( "lik.item", pos = as.environment( envir ) )
  
  #   lik.item[["param"]]   <-  rep( -1., length( param.name ) )
  #   names( lik.item[["param"]] ) <- param.name
  #   lik.item[["eta"]]     <- NA
  #   lik.item[["aniso"]]   <- list( 
  #     isotropic = initial.objects[["isotropic"]], 
  #     aniso = rep( -1., length( aniso.name ) )
  #   )
  #   names( lik.item[["aniso"]][["aniso"]] ) <- aniso.name
  lik.item[["param"]]   <- param
  lik.item[["eta"]]     <- sum( param[c( "variance", "snugget" )] ) / param["nugget"]
  lik.item[["aniso"]]   <- list( 
    isotropic = initial.objects[["isotropic"]], 
    aniso = aniso
  )
  names( lik.item[["aniso"]][["aniso"]] ) <- aniso.name
  if( is.null( georob.object ) ){
    lik.item[["Valpha"]]  <- list()
  } else {
    lik.item[["Valpha"]]  <- expand( georob.object[["Valpha.objects"]] )
  }
  #   lik.item[["effects"]] <- list()
  #   lik.item[["eeq"]]     <- list()
  #   lik.item[["defaultCluster"]] <- list()
  
  assign( "lik.item", lik.item, pos = as.environment( envir ) )
  
  ##  compute various expectations of psi, chi, etc.
  
  expectations <- numeric()
  
  ##  ... E[ Chi(x) ] (= E[ psi(x) * x ])
  
  t.exp <- integrate( 
    function( x, dpsi.function, tuning.psi ) {
      dnorm( x ) * dpsi.function( x, tuning.psi = tuning.psi )
    }, 
    lower = -Inf, upper = Inf, 
    dpsi.function = rho.psi.etc[["dpsi.function"]], 
    tuning.psi = tuning.psi
  )
  if( !identical( t.exp[["message"]], "OK" ) ) stop( t.exp[["message"]] )
  expectations["dpsi"] <- t.exp[["value"]]
  if( verbose > 1 ) cat( 
    "\nexpectation of psi'(epsilon/sigma)                    :", 
    signif( expectations["dpsi"] ), "\n" 
  )
  
  ##  ... E[ psi(x)^2 ] 
  
  t.exp <- integrate( 
    function( x, psi.function, tuning.psi ) {
      dnorm( x ) * ( psi.function( x, tuning.psi = tuning.psi ) )^2
    }, 
    lower = -Inf, upper = Inf, 
    psi.function = rho.psi.etc[["psi.function"]],
    tuning.psi = tuning.psi
  )
  if( !identical( t.exp[["message"]], "OK" ) ) stop( t.exp[["message"]] )
  expectations["psi2"] <- t.exp[["value"]]
  if( verbose > 1 ) cat( 
    "expectation of (psi(epsilon/sigma))^2                 :", 
    signif( t.exp[["value"]] ), "\n" 
  )
  
  ## xihat
  
  sel <- !is.na (betahat )
  xihat <- drop( XX[, sel, drop=FALSE] %*% betahat[sel] + bhat )
  names( xihat ) <- rownames( XX )
  
  r.hessian <- NULL
  
  if( tuning.psi < tuning.psi.nr ) {
    
    ## robust REML estimation
    
    if( any( c( fit.param, fit.aniso ) ) ){
      
      ##  find roots of estimating equations
      
      #       if( slv ){
      #         
      #         ##         if( identical( root.finding, "nleqslv" ) ){
      
      r.root <- nleqslv(
        x = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        fn = compute.estimating.equations,
        method = control.nleqslv[["method"]],
        global = control.nleqslv[["global"]],
        xscalm = control.nleqslv[["xscalm"]],
        control = control.nleqslv[["control"]],
        #         slv = slv,
        envir = envir,        
        variogram.model = variogram.model,
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = param.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
        yy = yy, TT = TT, xihat = xihat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        dpsi.function = rho.psi.etc[["dpsi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        force.gradient = force.gradient,
        expectations = expectations,
        control.parallel = control.parallel,
        verbose = verbose
      ) 
      
      #       r.param <- r.root[["x"]] names( r.param ) <- names(
      #       transformed.param[ fit.param ] )
      
      r.gradient <- r.root[["fvec"]]
      names( r.gradient ) <- c(
        names( transformed.param[ fit.param ] ),
        names( transformed.aniso[ fit.aniso ] )
      )
      
      r.converged <- r.root[["termcd"]] == 1
      r.convergence.code <- r.root[["termcd"]] 
      
      r.counts <- c( nfcnt = r.root[["nfcnt"]], njcnt = r.root[["njcnt"]] )
      
      ##         } 
      ##        
      ##         else if( identical( root.finding, "bbsolve" ) ) {
      ##           
      ##           r.root <- BBsolve(
      ##             par = c( 
      ##               xihat,
      ##               transformed.param[ fit.param ], 
      ##               transformed.aniso[ fit.aniso ] 
      ##             ),
      ##             fn = compute.expanded.estimating.equations,
      ##             method = bbsolve.method,
      ##             control = bbsolve.control,
      ##             quiet = verbose == 0,
      ##             slv = slv,
      ##             envir = envir,        
      ##             variogram.model = variogram.model,
      ##             fixed.param = c( 
      ##               transformed.param[ !fit.param ], 
      ##               transformed.aniso[ !fit.aniso ]
      ##             ),
      ##             param.name = param.name, 
      ##             aniso.name = aniso.name,
      ##             param.tf = param.tf,
      ##             bwd.tf = bwd.tf,
      ##             safe.param = safe.param,
      ##             lag.vectors = lag.vectors,
      ##             XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
      ##             yy = yy, TT = TT, 
      ##             psi.function = rho.psi.etc[["psi.function"]], 
      ##             dpsi.function = rho.psi.etc[["dpsi.function"]], 
      ##             tuning.psi = tuning.psi,
      ##             tuning.psi.nr = tuning.psi.nr,
      ##             irwls.initial = irwls.initial,
      ##             irwls.maxiter = irwls.maxiter,
      ##             irwls.reltol = irwls.reltol,
      ##             force.gradient = force.gradient,
      ##             expectations = expectations,
      ##             verbose = verbose
      ##           ) 
      ##
      ##           r.converged <- r.root[["convergence"]] == 0
      ##           r.convergence.code <- r.root[["convergence"]] 
      ##           r.counts <- c( nfcnt = r.root[["feval"]], njcnt = NA_integer_ )
      ## 
      ##           r.gradient <- compute.expanded.estimating.equations(
      ##             allpar = r.root[["par"]],
      ##             slv = TRUE,
      ##             envir = envir,        
      ##             variogram.model = variogram.model,
      ##             fixed.param = c( 
      ##               transformed.param[ !fit.param ], 
      ##               transformed.aniso[ !fit.aniso ]
      ##             ),
      ##             param.name = param.name, 
      ##             aniso.name = aniso.name,
      ##             param.tf = param.tf,
      ##             bwd.tf = bwd.tf,
      ##             safe.param = safe.param,
      ##             lag.vectors = lag.vectors,
      ##             XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
      ##             yy = yy, TT = TT,  
      ##             psi.function = rho.psi.etc[["psi.function"]], 
      ##             dpsi.function = rho.psi.etc[["dpsi.function"]], 
      ##             tuning.psi = tuning.psi,
      ##             tuning.psi.nr = tuning.psi.nr,
      ##             irwls.initial = irwls.initial,
      ##             irwls.maxiter = irwls.maxiter,
      ##             irwls.reltol = irwls.reltol,
      ##             force.gradient = force.gradient,
      ##             expectations = expectations,
      ##             verbose = verbose
      ##           )
      ##           
      ##         }
      
      #         } 
      
      ##         else {
      ##         
      ##         ## minimize sum of squared estimating equations
      ##         
      ##         r.opt.eeq.sq <- optim(
      ##           par = c( 
      ##             transformed.param[ fit.param ], 
      ##             transformed.aniso[ fit.aniso ] 
      ##           ),
      ##           fn = compute.estimating.equations,
      ##           method = control.optim[["method"]], 
      ##           lower = control.optim[["lower"]],
      ##           upper = control.optim[["upper"]],
      ##           control = control.optim[["control"]],
      ##           hessian = FALSE,
      ##           slv = slv,
      ##           envir = envir,        
      ##           variogram.model = variogram.model,
      ##           fixed.param = c( 
      ##             transformed.param[ !fit.param ], 
      ##             transformed.aniso[ !fit.aniso ]
      ##           ),
      ##           param.name = param.name, 
      ##           aniso.name = aniso.name,
      ##           param.tf = param.tf,
      ##           bwd.tf = bwd.tf,
      ##           safe.param = safe.param,
      ##           lag.vectors = lag.vectors,
      ##           XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
      ##           yy = yy, TT = TT, xihat = xihat, 
      ##           psi.function = rho.psi.etc[["psi.function"]], 
      ##           dpsi.function = rho.psi.etc[["dpsi.function"]], 
      ##           tuning.psi = tuning.psi,
      ##           tuning.psi.nr = tuning.psi.nr,
      ##           ml.method = ml.method,
      ##           irwls.initial = irwls.initial,
      ##           irwls.maxiter = irwls.maxiter,
      ##           irwls.reltol = irwls.reltol,
      ##           force.gradient = force.gradient,
      ##           expectations = expectations,
      ##           control.parallel = control.parallel,
      ##           verbose = verbose
      ##         )
      ##         
      ##         r.converged <- r.opt.eeq.sq[["convergence"]] == 0
      ##         r.convergence.code <- r.opt.eeq.sq[["convergence"]]      
      ##         r.counts <- r.opt.eeq.sq[["counts"]]
      ##         
      ##         if( verbose > 0 ){
      ##           cat( 
      ##             "\n  sum(EEQ^2)         :",
      ##             format( 
      ##               signif( r.opt.eeq.sq[["value"]], digits = 7 ), 
      ##               scientific = TRUE, width = 14
      ##             ), sep = "" 
      ##           )
      ##           cat( 
      ##             "\n  convergence code   :", 
      ##             format( 
      ##               signif( r.opt.eeq.sq[["convergence"]], digits = 0 ), 
      ##               scientific = FALSE, width = 14
      ##             ), "\n\n", sep = "" 
      ##           )
      ##         }
      ##                
      ##         #         if( hessian ) r.hessian <- r.opt.eeq.sq[["hessian"]]
      ##         
      ##         r.gradient <- compute.estimating.equations(
      ##           adjustable.param = r.opt.eeq.sq[["par"]],
      ##           slv = TRUE,
      ##           envir = envir,        
      ##           variogram.model = variogram.model,
      ##           fixed.param = c( 
      ##             transformed.param[ !fit.param ], 
      ##             transformed.aniso[ !fit.aniso ]
      ##           ),
      ##           param.name = param.name, 
      ##           aniso.name = aniso.name,
      ##           param.tf = param.tf,
      ##           bwd.tf = bwd.tf,
      ##           safe.param = safe.param,
      ##           lag.vectors = lag.vectors,
      ##           XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
      ##           yy = yy, TT = TT, xihat = xihat, 
      ##           psi.function = rho.psi.etc[["psi.function"]], 
      ##           dpsi.function = rho.psi.etc[["dpsi.function"]], 
      ##           tuning.psi = tuning.psi,
      ##           tuning.psi.nr = tuning.psi.nr,
      ##           ml.method = ml.method,
      ##           irwls.initial = irwls.initial,
      ##           irwls.maxiter = irwls.maxiter,
      ##           irwls.reltol = irwls.reltol,
      ##           force.gradient = force.gradient,
      ##           expectations = expectations,
      ##           control.parallel = control.parallel,
      ##           verbose = verbose
      ##         )
      ## 
      ##       }
      
    } else {
      
      ##  all variogram parameters are fixed
      
      ##  evaluate estimating equations
      
      r.gradient <- compute.estimating.equations(
        adjustable.param = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        #         slv = TRUE,
        envir = envir,        
        variogram.model = variogram.model,
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = param.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
        yy = yy, TT = TT, xihat = xihat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        dpsi.function = rho.psi.etc[["dpsi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        force.gradient = force.gradient,
        expectations = expectations,
        control.parallel = control.parallel,
        verbose = verbose
      )
      
      r.converged <- NA
      r.convergence.code <- NA_integer_
      r.counts <- c( nfcnt = NA_integer_, njcnt = NA_integer_ )
      
    }
    
    r.opt.neg.loglik <- NA_real_
    
  } else {
    
    if( any( c( fit.param, fit.aniso ) ) ){
      
      ##  Gaussian REML estimation
      
      r.opt.neg.restricted.loglik <- optim(
        par = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        fn = negative.restr.loglikelihood,
        gr = gradient.negative.restricted.loglikelihood,
        method = control.optim[["method"]], 
        lower = control.optim[["lower"]],
        upper = control.optim[["upper"]],
        control = control.optim[["control"]],
        hessian = hessian,
        envir = envir,        
        variogram.model = variogram.model,
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = param.tf,
        deriv.fwd.tf = deriv.fwd.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
        yy = yy, TT = TT, xihat = xihat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        dpsi.function = rho.psi.etc[["dpsi.function"]], 
        d2psi.function = rho.psi.etc[["d2psi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        control.parallel = control.parallel,
        verbose = verbose,
        force.gradient = force.gradient
      )    
      
      r.opt.neg.loglik <- r.opt.neg.restricted.loglik[["value"]]     
      r.converged <- r.opt.neg.restricted.loglik[["convergence"]] == 0
      r.convergence.code <- r.opt.neg.restricted.loglik[["convergence"]]      
      r.counts <- r.opt.neg.restricted.loglik[["counts"]]
      
      if( hessian ) r.hessian <- r.opt.neg.restricted.loglik[["hessian"]]
      
      r.gradient <- gradient.negative.restricted.loglikelihood(
        adjustable.param = r.opt.neg.restricted.loglik[["par"]],
        envir = envir,
        variogram.model = variogram.model, 
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = param.tf,
        deriv.fwd.tf = deriv.fwd.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
        yy = yy, TT = TT, xihat = xihat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        dpsi.function = rho.psi.etc[["dpsi.function"]], 
        d2psi.function = rho.psi.etc[["d2psi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        force.gradient = force.gradient,
        control.parallel = control.parallel,
        verbose = verbose
      )
      
    } else {
      
      ##  all variogram parameters are fixed
      
      ##  compute negative restricted loglikelihood and gradient
      
      r.opt.neg.loglik <- negative.restr.loglikelihood(
        adjustable.param = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        envir = envir,
        variogram.model = variogram.model, 
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),,
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = param.tf,
        deriv.fwd.tf = deriv.fwd.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
        yy = yy, TT = TT, xihat = xihat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        dpsi.function = rho.psi.etc[["dpsi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        control.parallel = control.parallel,
        verbose = verbose
      )
      
      r.gradient <- gradient.negative.restricted.loglikelihood(
        adjustable.param = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        envir = envir,
        variogram.model = variogram.model, 
        fixed.param = c( 
          transformed.param[ !fit.param ], 
          transformed.aniso[ !fit.aniso ]
        ),
        param.name = param.name, 
        aniso.name = aniso.name,
        param.tf = param.tf,
        deriv.fwd.tf = deriv.fwd.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, rankdef.x = rankdef.x,
        yy = yy, TT = TT, xihat = xihat, 
        psi.function = rho.psi.etc[["psi.function"]], 
        dpsi.function = rho.psi.etc[["dpsi.function"]], 
        d2psi.function = rho.psi.etc[["d2psi.function"]], 
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        force.gradient = force.gradient,
        control.parallel = control.parallel,
        verbose = verbose
      )
      
      r.converged <- NA
      r.convergence.code <- NA_integer_
      r.counts <- c( nfcnt = NA_integer_, njcnt = NA_integer_ )
      
    }
    
  }
  
  ##  get the other fitted items
  
  lik.item <- get( "lik.item", pos = as.environment( envir ) )
  
  ##  compute covariance matrices of betahat and bhat etc.
  
  if( any( c( 
        cov.bhat, cov.betahat, cov.bhat.betahat, 
        cov.delta.bhat, cov.delta.bhat.betahat, 
        cov.ehat, cov.ehat.p.bhat,
        aux.cov.pred.target
      ) 
    ) 
  ){
    
    ##  compute the covariances
    
    r.cov <- compute.covariances(
      Valpha.objects = lik.item[["Valpha"]],
      Aalpha = lik.item[["effects"]][["Aalpha"]],
      Palpha = lik.item[["effects"]][["Palpha"]],
      rweights = lik.item[["effects"]][["rweights"]],
      XX = XX, TT = TT, names.yy = names( yy ),
      nugget = lik.item[["param"]]["nugget"],
      eta = lik.item[["eta"]],
      expectations = expectations,
      cov.bhat = cov.bhat, full.cov.bhat = full.cov.bhat,
      cov.betahat = cov.betahat, 
      cov.bhat.betahat = cov.bhat.betahat,
      cov.delta.bhat = cov.delta.bhat, full.cov.delta.bhat = full.cov.delta.bhat,
      cov.delta.bhat.betahat = cov.delta.bhat.betahat,
      cov.ehat = cov.ehat, full.cov.ehat = full.cov.ehat, 
      cov.ehat.p.bhat = cov.ehat.p.bhat, full.cov.ehat.p.bhat = full.cov.ehat.p.bhat,
      aux.cov.pred.target = aux.cov.pred.target,
      control.parallel = control.parallel,
      verbose = verbose
    )[-1]
    
  }
  
  ## stop SNOW and snowfall clusters
  
  #   f.stop.cluster()
  
  if( length( lik.item[["defaultCluster"]] ) > 0 ){
    cl <- lik.item[["defaultCluster"]]
    
    junk <- parLapply( cl, 1:length(cl), function( i ) sfStop() )
    junk <- stopCluster( cl )
    sfStop()
    options( error = NULL )  
  } 
  
  if( sfIsRunning() ){
    sfStop()
    options( error = NULL )  
  }
  
  if( file.exists( "SOCKcluster.RData" ) ){
    file.remove( "SOCKcluster.RData" )
  } 
  
  
  attr( r.gradient, "eeq.emp" )    <- lik.item[["eeq"]][["eeq.emp"]]
  attr( r.gradient, "eeq.exp" )    <- lik.item[["eeq"]][["eeq.exp"]]
  
  ##      ##  compute residual degrees of freedom 
  ##      
  ##      r.df <- f.compute.df( 
  ##          Valpha = lik.item[["Valpha"]][["Valpha"]],
  ##          XX = XX, 
  ##          param = lik.item[["param"]]
  ##      )
  
  ##  collect output
  
  result.list <- list(
    loglik = -r.opt.neg.loglik,
    variogram.model = variogram.model,
    param = lik.item[["param"]],
    aniso = lik.item[["aniso"]],
    gradient = r.gradient,
    psi.func = psi.func,
    tuning.psi = tuning.psi,
    coefficients = lik.item[["effects"]][["betahat"]],
    fitted.values = drop( XX %*% lik.item[["effects"]][["betahat"]] )[TT],
    bhat = lik.item[["effects"]][["bhat"]],
    residuals = lik.item[["effects"]][["residuals"]],
    rweights = lik.item[["effects"]][["rweights"]],
    converged = r.converged,
    convergence.code = r.convergence.code,
    iter = r.counts,
    Tmat = TT
  )
  names( result.list[["fitted.values"]] ) <- names( result.list[["residuals"]] )
  
  if( any( c( 
        cov.bhat, cov.betahat, cov.bhat.betahat, 
        cov.delta.bhat, cov.delta.bhat.betahat, 
        cov.ehat, cov.ehat.p.bhat, aux.cov.pred.target
      ) 
    ) 
  ){
    
    result.list[["cov"]] <- compress( r.cov )
    
  }
  
  ## map angles to halfcircle
  
  if( !result.list[["aniso"]][["isotropic"]] ){
    
    if( result.list[["aniso"]][["aniso"]]["omega"] < 0. ){
      result.list[["aniso"]][["aniso"]]["omega"] <- 
      result.list[["aniso"]][["aniso"]]["omega"] + 180.
    }
    if( result.list[["aniso"]][["aniso"]]["omega"] > 180. ){
      result.list[["aniso"]][["aniso"]]["omega"] <- 
      result.list[["aniso"]][["aniso"]]["omega"] - 180.
    }
    if( result.list[["aniso"]][["aniso"]]["phi"] < 0. ){
      result.list[["aniso"]][["aniso"]]["phi"] <- 
      result.list[["aniso"]][["aniso"]]["phi"] + 180.
    }
    if( result.list[["aniso"]][["aniso"]]["phi"] > 180. ){
      result.list[["aniso"]][["aniso"]]["phi"] <- 
      result.list[["aniso"]][["aniso"]]["phi"] - 180.
    }
    if( result.list[["aniso"]][["aniso"]]["zeta"] < 90. ){
      result.list[["aniso"]][["aniso"]]["zeta"] <- 
      result.list[["aniso"]][["aniso"]]["zeta"] + 180.
    }
    if( result.list[["aniso"]][["aniso"]]["zeta"] > 90. ){
      result.list[["aniso"]][["aniso"]]["zeta"] <- 
      result.list[["aniso"]][["aniso"]]["zeta"] - 180.
    }
  
  }
  
  ##      result.list[["df.model"]] <- r.df
  
  if( full.output ){
    
    result.list[["param.tf"]] <- param.tf
    result.list[["fwd.tf"]] <- fwd.tf
    result.list[["bwd.tf"]] <- bwd.tf
    if( !is.null( r.hessian ) ){
      result.list[["hessian"]] <- r.hessian
    }
    result.list[["expectations"]]     <- expectations
    result.list[["Valpha.objects"]]   <- compress( lik.item[["Valpha"]] )
    result.list[["Aalpha"]] <- lik.item[["effects"]][["Aalpha"]]
    result.list[["Palpha"]] <- lik.item[["effects"]][["Palpha"]]

    result.list[["locations.objects"]] <- initial.objects[["locations.objects"]]
    result.list[["locations.objects"]][["lag.vectors"]] <- lag.vectors
    
    result.list[["initial.objects"]] <- list(
      coefficients = initial.objects[["betahat"]],
      bhat = initial.objects[["bhat"]],
      param = param,
      fit.param = fit.param,
      aniso = aniso,
      fit.aniso = fit.aniso
    )
    
    
    
  }
  
  return(result.list)
  
}

#  ##############################################################################

getCall.georob <- 
  function( object )
{
  
  ## Function replaces the name of a formula object in the call component
  ## of a georob object by the formula itself (needed for update.default to
  ## work correctly)

  ## 2013-06-12 AP substituting [["x"]] for $x in all lists

  if( is.null( call <- getElement( object, "call" ) ) )  stop(
    "need an object with call component"
  )
  call[["formula"]] <- update.formula( formula(object), formula( object ) )
  
  return( call )
  
}


################################################################################

compute.semivariance <- 
  function( 
    lag.vectors, variogram.model, param, aniso
  )
{

  ## auxiliary function to compute semivariances for an anisotropic model
  
  ## arguments:
  
  ## param                                        vector with variogram parameters in standard order
  ## aniso                                        list with component rotmat and sclmat for coordinate
  ##                                                                  transformation in 3d
  ## lag.vectors    
  
  ## 2012-04-13 A. Papritz
  ## 2012-05-23 ap correction in model.list for models with more than 4 parameters
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-03-05 AP changes for version 3 of RandomFields
  
  ## matrix for coordinate transformation
  
  A <- with(
    aniso,
    sclmat * rotmat / param["scale"]
  )
  
  ## set up variogram model object
  
  model.list <- list( "+",
    list( "$", 
      var = param["variance"], 
      A = A, 
      if( length( param[-(1:4)] ) > 0 ){
         c( list( variogram.model ) , as.list(param[-(1:4)]) )
      } else {
        list( variogram.model )
      }
    ),
    list( "$", 
      var = sum( param[ c("nugget", "snugget") ] ), 
      list( "nugget" )
    )
  )
  
  ##  semivariance 
  
  ## functions of version 3 of RandomFields
  
  RFoptions(newAniso=FALSE)
  
  r.gamma <- try(
    RFvariogram(
      x = lag.vectors, model = model.list, dim = NCOL( lag.vectors ), grid = FALSE
    ),
    silent = TRUE
  )

  ## functions of version 3 of RandomFields

  ##   RFoldstyle()
  ##   r.gamma <- try(
  ##     Variogram( lag.vectors, model = model.list ),
  ##     silent = TRUE
  ##   )
  
  return( r.gamma )
  
}

################################################################################

f.aux.eeq <- function( 
  param, lik.item, TtT, r.cov,
  lag.vectors, variogram.model, control.pmm, verbose
){
  
  ##  auxiliary function to compute robustified estimating equations
  ##  (called by compute.estimating.equations)
  
  ## 2014-07-29 A. Papritz
  
  switch(
    param,
    nugget = {
      
      ## nugget
      
      # eeq.exp <- sum( 
      #   diag( 
      #     lik.item[["Valpha"]][["Valpha.inverse"]] %*%             
      #     ( 1/TtT * lik.item[["Valpha"]][["Valpha.inverse"]] ) %*% 
      #     r.cov[["cov.bhat"]] 
      #   )
      # )
      aux <- pmm( 
        ( 1/TtT * lik.item[["Valpha"]][["Valpha.inverse"]] ),
        r.cov[["cov.bhat"]],
        control.pmm
      )
      eeq.exp <- sum(
        diag( pmm( lik.item[["Valpha"]][["Valpha.inverse"]], aux, control.pmm ) )
      )
      eeq.emp <- sum( 
        ( lik.item[["effects"]][["z.star"]] )^2 / TtT
      )
      
    },
    snugget = {
      
      ## snaugget
      
      # eeq.exp <- sum(
      #   rowSums( 
      #     (lik.item[["Valpha"]][["Valpha.inverse"]] %*% lik.item[["Valpha"]][["Valpha.inverse"]] ) * 
      #     r.cov[["cov.bhat"]]
      #   )
      # )
      eeq.exp <- sum(
        rowSums( 
          pmm( lik.item[["Valpha"]][["Valpha.inverse"]], lik.item[["Valpha"]][["Valpha.inverse"]], control.pmm ) * 
          r.cov[["cov.bhat"]]
        )
      )
      eeq.emp <- sum( lik.item[["effects"]][["z.star"]]^2 )
      
    },
    variance = {
      
      ## variance
      
      # eeq.exp["variance"] <- sum(
      #   rowSums( 
      #     ( lik.item[["Valpha"]][["Valpha.inverse"]] %*% lik.item[["Valpha"]][["Valpha0"]] %*% lik.item[["Valpha"]][["Valpha.inverse"]] ) * 
      #     r.cov[["cov.bhat"]]
      #   )
      # )
      aux <- pmm( lik.item[["Valpha"]][["Valpha0"]], lik.item[["Valpha"]][["Valpha.inverse"]], control.pmm )
      eeq.exp <- sum(
        rowSums( 
          pmm( lik.item[["Valpha"]][["Valpha.inverse"]], aux, control.pmm ) * r.cov[["cov.bhat"]]
        )
      )
      eeq.emp <- sum( 
        lik.item[["effects"]][["z.star"]] * drop( lik.item[["Valpha"]][["Valpha0"]] %*% lik.item[["effects"]][["z.star"]] )
      )
      
    },
    {
      
      ## extra parameters 
      
      dValpha <- dcorr.dparam(
        x = lag.vectors, variogram.model = variogram.model, param = lik.item[["param"]], 
        d.param = param,
        aniso = lik.item[["aniso"]],
        verbose = verbose
      )
      ##       if( identical( class( dValpha ), "try-error" ) ){
      ##         if( verbose > 0 ) cat( "error in dcorr.dparam\n\n" )
      ##         t.result <- rep( Inf, length( adjustable.param ) )
      ##         names( t.result ) <- names( adjustable.param )
      ##         return( t.result )
      ##       }
      
      # eeq.exp <- sum(
      #   rowSums( 
      #     (lik.item[["Valpha"]][["Valpha.inverse"]] %*% dValpha %*% lik.item[["Valpha"]][["Valpha.inverse"]]) * 
      #     r.cov[["cov.bhat"]]
      #   )
      # )
      aux <- pmm( dValpha, lik.item[["Valpha"]][["Valpha.inverse"]], control.pmm )
      eeq.exp <- sum(
        rowSums( 
          pmm( lik.item[["Valpha"]][["Valpha.inverse"]], aux, control.pmm ) * r.cov[["cov.bhat"]]
        )
      )
      eeq.emp <- sum( 
        lik.item[["effects"]][["z.star"]] * drop( dValpha %*% lik.item[["effects"]][["z.star"]] )
      )
    }
  )
  
  c( eeq.exp = eeq.exp, eeq.emp = eeq.emp )
  
}

################################################################################

f.aux.gradient.nll <- function( 
  param, 
  nme.var.snug,
  TT, TtT, std.res, XX, 
  Qi, Vi,
  n, 
  lag.vectors, variogram.model,
  lik.item,
  param.tf, deriv.fwd.tf, 
  ml.method,
  control.pmm, verbose
){
  
  ##  auxiliary function to compute robustified estimating equations
  ##  (called by compute.estimating.equations)
  
  ## 2014-07-29 A. Papritz
  
  switch(
    param,
    nugget = {
      
      ##  derivative of log( det( Q ) )
      
      if( identical( ml.method, "REML" ) ){
        TtTX <- TtT * XX
        dlogdetQ <- -sum( 
          Qi * rbind( 
            cbind( diag( TtT ),                TtTX   ),
            cbind(     t(TtTX), crossprod( XX, TtTX ) )
          )
        ) / lik.item[["param"]]["nugget"]^2
      } else {
        dlogdetQ <- -sum( TtT * diag(Qi) ) / lik.item[["param"]]["nugget"]^2
      }
      
      ##  derivate of U with respect to nugget
      
      Ttstd.res <- as.vector( tapply( std.res, factor( TT ), sum ) )
      dU <- -0.5 * sum( Ttstd.res^2 / TtT ) / lik.item[["param"]]["nugget"]
      
      t.result <- c( nugget = ( -0.5 * ( n / lik.item[["param"]]["nugget"] + dlogdetQ ) - dU ) / 
        deriv.fwd.tf[[param.tf["nugget"]]]( lik.item[["param"]]["nugget"] )
      )
      
    },
    Vparam = {
      
      ##  compute partial derivative of restricted log-likelihood with
      ##  respect to spatial nugget
      
      Vi2 <- pmm( Vi, Vi, control.pmm )

      t.result <- numeric()
      
      ##  compute partial derivative of restricted log-likelihood with
      ##  respect to variance
      
      if( "variance" %in% nme.var.snug ) {
        
        ## derivative of V with respect to variance
        
        dlogdetV <- sum( 1. - lik.item[["param"]]["snugget"] * diag(Vi) ) / 
        lik.item[["param"]]["variance"]
        
        ##  derivative of log( det( Q ) ) with respect to variance
        
        dlogdetQ <- -sum( Qi[1:n, 1:n] * ( Vi + lik.item[["param"]]["snugget"] * Vi2 ) ) / 
        lik.item[["param"]]["variance"]
        
        ##  derivate of U with respect to variance
        
        dU <- -0.5 * sum( 
          lik.item[["effects"]][["z.star"]] * drop( lik.item[["Valpha"]][["Valpha0"]] %*% lik.item[["effects"]][["z.star"]] ) 
        ) / ( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] )^2
        
        t.result <- c(
          variance = ( -0.5 * ( dlogdetV + dlogdetQ ) - dU ) / 
            deriv.fwd.tf[[param.tf["variance"]]]( lik.item[["param"]]["variance"] )
        )
      }
      
      if( "snugget" %in% nme.var.snug ){
        
        ## derivative of V with respect to spatial nugget
        
        dlogdetV <- sum( diag( Vi ) )
        
        ##  derivative of log( det( Q ) ) with respect to spatial nugget
        
        dlogdetQ <- -sum( Qi[1:n, 1:n] * Vi2 )
        
        ##  derivate of U with respect to spatial nugget
        
        dU <- -0.5 * sum( lik.item[["effects"]][["z.star"]]^2 ) /
        ( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] )^2
        
        t.result <- c(
          t.result,
          snugget = ( -0.5 * ( dlogdetV + dlogdetQ ) - dU ) / 
            deriv.fwd.tf[[param.tf["snugget"]]]( lik.item[["param"]]["snugget"] )
        )
        
      }
      
    },
    {
      
      ## extra parameters 

      dValpha <- dcorr.dparam(
        x = lag.vectors, variogram.model = variogram.model, param = lik.item[["param"]], 
        d.param = param,
        aniso = lik.item[["aniso"]],
        verbose = verbose
      )
      
      dValpha.Valphai <- pmm( dValpha, Vi, control.pmm ) 
      
      ## derivatie of V with respect to extra parameter
      
      dlogdetV <- lik.item[["param"]]["variance"] * sum( diag( dValpha.Valphai ) )
      
      ##  derivative of log( det( Q ) ) with respect to extra parameter
      
      dlogdetQ <- -sum( 
        Qi[1:n, 1:n] * pmm( Vi, dValpha.Valphai, control.pmm )
      ) * lik.item[["param"]]["variance"]
      
      ##  derivate of U with respect to extra parameter
      
      dU <- -0.5 * sum( 
        lik.item[["effects"]][["z.star"]] * drop( dValpha %*% lik.item[["effects"]][["z.star"]] )
      ) * lik.item[["param"]]["variance"] / ( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] )^2
      
      t.result <- ( -0.5 * ( dlogdetV + dlogdetQ ) - dU ) / 
        deriv.fwd.tf[[param.tf[param]]]( 
          c( lik.item[["param"]], lik.item[["aniso"]][["aniso"]] )[param] 
        )
      names( t.result ) <- param
            
    }
  )
  
  return( t.result )
  
}

################################################################################

f.aux.Q <- function( TT, XX,  rankdef.x, lik.item, min.condnum, ml.method, control.pmm ){
  
  ## auxiliary function to compute matrix Q used for Gaussian
  ## log-likelihood (called by prepare.likelihood.calculations)
  
  ## 2014-07-29 A. Papritz
  
  result <- list( error = TRUE, log.det.Q = NULL, Q.inverse = NULL )
  
  TtT <- as.vector( table( TT ) )
  TtTX <- TtT * XX
  
  ##  compute matrix Q
  
  Q <-  lik.item[["Valpha"]][["Valpha.inverse"]] / lik.item[["eta"]]
  diag( Q ) <- diag( Q ) + TtT
  
  if( identical( ml.method, "REML" ) ){
    Q <- rbind( 
      cbind( Q,       TtTX                 ),
      cbind( t(TtTX), crossprod( XX, TtTX) )
    )
  }
  
  Q <- Q / lik.item[["param"]]["nugget"]
  
  if( rankdef.x && ml.method == "REML" ){
    
    ## compute log(pseudo.det(Q)) and (Moore-Penrose) pseudo inverse of Q by svd
    
    result[["error"]] <- FALSE
    s <- svd( Q )
    result[["log.det.Q"]] <- sum( log( s[["d"]][s[["d"]] / max( s[["d"]] ) > min.condnum] ) )
    s[["d"]] <- ifelse( s[["d"]] / max( s[["d"]] ) <= min.condnum, 0., 1. / s[["d"]] )
    #         result[["Q.inverse"]] <- s[["v"]] %*% ( s[["d"]] * t( s[["u"]] ) )
    result[["Q.inverse"]] <- pmm(
      s[["v"]], s[["d"]] * t( s[["u"]] ), control.pmm
    )
    
  } else {
    
    ##  compute log(det(Q)) and inverse of Q by cholesky decomposition
    
    t.chol <- try( chol( Q ), silent = TRUE )
    
    if( !identical( class( t.chol ), "try-error" ) ) {
      
      result[["error"]] <- FALSE
      result[["log.det.Q"]] <- 2 * sum( log( diag( t.chol) ) )
      result[["Q.inverse"]] <- chol2inv( t.chol )
      
    }
    
  }
  
  result
  
}

################################################################################

f.aux.Valpha <- function(
  lag.vectors, variogram.model, param, aniso, control.parallel, verbose
){
  
  ## auxiliary function to compute generalized correlation matrix and
  ## related items (called by prepare.likelihood.calculations)
  
  ## 2014-07-29 A. Papritz
  
  result <- list( 
    error = TRUE, gcr.constant = NULL, Valpha = NULL, Valpha0 = NULL,
    Valpha.ucf = NULL, Valpha.ilcf = NULL, Valpha.inverse = NULL
  )
  
  cormat <- gcr(
    lag.vectors = lag.vectors, variogram.model = variogram.model, param = param, 
    aniso = aniso, control.parallel = control.parallel,
    verbose = verbose
  )
  if( cormat[["error"]] ) return( result )
  
  t.vchol <- try( chol( cormat[["Valpha"]] ), silent = TRUE )
  if( !identical( class( t.vchol ), "try-error" ) ) {
    result[["gcr.constant"]]  <- cormat[["gcr.constant"]]
    result[["Valpha"]]        <- cormat[["Valpha"]]
    result[["Valpha0"]]       <- cormat[["Valpha0"]]
    result[["Valpha.ucf"]]    <- unname( t.vchol )
    result[["Valpha.ilcf"]]   <- try(
      t( 
        backsolve( 
          t.vchol, 
          diag( nrow( result[["Valpha"]] ) ), 
          k = nrow( result[["Valpha"]] ) 
        ) 
      ),
      silent = TRUE
    )
    if( identical( class( result[["Valpha.ilcf"]] ), "try-error" ) ) {
      return( result )
    }
    result[["error"]]         <- FALSE
    result[["Valpha.inverse"]] <- pmm(
      t(result[["Valpha.ilcf"]]),  result[["Valpha.ilcf"]], control.parallel
    )
    attr( result[["Valpha"]], "struc" )         <- "sym"
    attr( result[["Valpha0"]], "struc" )        <- "sym"
    attr( result[["Valpha.inverse"]], "struc" ) <- "sym"
    attr( result[["Valpha.ucf"]], "struc" )     <- "ut"
    attr( result[["Valpha.ilcf"]], "struc" )    <- "lt"
  }
  result
}


################################################################################

f.stop.cluster <- function( cl = NULL ){
  
  ## function to stop snow and snowfall clusters 
  ## 2014-07-31 A. Papritz
  
  if( sfIsRunning() ){
    sfStop()
  }
  
  if( file.exists( "SOCKcluster.RData" ) ){
    if( is.null( cl ) ) load( "SOCKcluster.RData" )
    file.remove( "SOCKcluster.RData" )
  } 
  
  if( !is.null( cl ) ){
    junk <- parLapply( cl, 1:length(cl), function( i ) sfStop() )
    junk <- stopCluster( cl )
  }
  options( error = NULL )  

}

