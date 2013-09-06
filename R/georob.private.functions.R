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
    extended.output = FALSE,
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
  
  
  ## adjust flags for computing covariance matrices
  
  cov.bhat.b       <- FALSE
  cov.bhat.e       <- FALSE
  cov.betahat.b    <- FALSE
  cov.betahat.e    <- FALSE
  
  if( any( c( cov.delta.bhat, aux.cov.pred.target )))                          cov.bhat.b <- TRUE
  if( any( c( cov.ehat, aux.cov.pred.target )))                               cov.bhat.e <- TRUE
  if( any( c( cov.delta.bhat.betahat, cov.ehat.p.bhat, aux.cov.pred.target ))) cov.betahat.b <- TRUE
  if( any( c( cov.ehat, cov.ehat.p.bhat, aux.cov.pred.target )))              cov.betahat.e <- TRUE
  if( any( c( cov.ehat, cov.ehat.p.bhat ) ) )                             cov.betahat <- TRUE
  if( any( c( cov.delta.bhat, cov.delta.bhat.betahat ))){
                                                                          cov.bhat <- TRUE
    if( full.cov.delta.bhat )                                              full.cov.bhat <- TRUE
  }
  if( cov.delta.bhat.betahat )                                             cov.bhat.betahat <- TRUE
  if( cov.ehat ){
                                                                          cov.delta.bhat.betahat <- TRUE
                                                                          cov.delta.bhat <- TRUE
    if( full.cov.ehat )                                                   full.cov.delta.bhat <- TRUE
  }
  
  ## compute required auxiliary items 
    
  result.new <- list( error = FALSE )
  
  a <- expectations["psi2"] 
  b <- expectations["dpsi"]
  
  TtT <- as.vector( table( TT ) )

  V <- eta * nugget * Valpha.objects[["Valpha"]] 
  VTtT <- t( TtT * V )
  
  ## inverse of G_theta (note that Palpha is not symmetric!!)
  
  Gi <- try( 
    solve( t(TtT * Valpha.objects[["Valpha"]]) + Palpha / ( b * eta ) ),
    silent = TRUE
  )
  
  if( identical( class( Gi ), "try-error" ) ){
    result.new[["error"]] <- TRUE
    return( result.new )            
  }
  
  ## factors to compute bhat and betahat from xihat
  
  if( any( c( 
        cov.bhat.b, cov.bhat.e,
        cov.bhat, cov.bhat.betahat ) ) 
  ){
    PaGtiVa <- ( Palpha %*% Gi %*% Valpha.objects[["Valpha"]] )
  }
  
  if( any( c( 
        cov.betahat.b, cov.betahat.e,
        cov.betahat, cov.bhat.betahat
      ) ) 
  ){
    AaGtiVa <- ( Aalpha %*% Gi %*% Valpha.objects[["Valpha"]] )
  }
  
  ## covariance of huberized observations
  
  if( any( c( cov.bhat, cov.betahat, cov.bhat.betahat ) )
  ){
    cov.b.psi <- TtT * VTtT
    diag( cov.b.psi ) <- diag( cov.b.psi ) + (a * nugget / b^2) * TtT
  }
  
  ## covariance of bhat and betahat with B and epsilon
  
  if( cov.bhat.b )    cov.bhat.b      <- PaGtiVa %*% t( VTtT )
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
    aux <- tcrossprod( PaGtiVa, cov.b.psi )
    result.new[["cov.bhat"]] <- if( full.cov.bhat )
    {
      aux <- tcrossprod( aux, PaGtiVa )
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
    result.new[["cov.bhat.betahat"]] <- tcrossprod( tcrossprod( PaGtiVa, cov.b.psi ), AaGtiVa )
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
    result.new[["cov.pred.target"]] <- rbind( cov.bhat.b, cov.betahat.b ) %*%
      Valpha.objects[["Valpha.inverse"]] / eta / nugget
  }
  
  ##   result <- list( error = FALSE )
  ##   
  ##   n <- nrow( XX )
  ##   sel <- 1:n
  ##   
  ##   TtWiT <- b
  ##   TtDT  <- a
  ##   
  ##   ##  ... aggregate elements of Wi and D for replicated observations
  ##   
  ##   if( sum( duplicated( TT ) ) > 0 ){
  ##     TtWiT  <- TtWiT * TtT
  ##     TtDT   <- TtDT * TtT
  ##   }
  ##   
  ##   ##  construct matrix M
  ##   
  ##   TtWiTXX <- TtWiT * XX
  ##   aux <- Valpha.objects[["Valpha.inverse"]] / eta
  ##   diag( aux ) <- diag( aux ) + TtWiT
  ##   
  ##   M <- rbind(
  ##     cbind( aux,        TtWiTXX ),
  ##     cbind( t(TtWiTXX), crossprod( XX , TtWiTXX ) )
  ##   )
  ##   
  ##   ## ... and invert it
  ##   
  ##   M.inverse <- try( solve( M ), silent = TRUE )
  ##   
  ##   if( identical( class( M.inverse ), "try-error" ) ){
  ##     result[["error"]] <- TRUE
  ##     return( result )            
  ##   }
  ##   
  ##   ##  compute auxiliary matrices
  ##   
  ##   PpXQt <- M.inverse[ sel, sel] + XX %*% M.inverse[-sel, sel]
  ##   QpXS <-  M.inverse[ sel,-sel] + XX %*% M.inverse[-sel,-sel]
  ##   sqrtD <- sqrt( TtDT )
  ##   
  ##   ValphaiP <- Valpha.objects[["Valpha.inverse"]] %*% M.inverse[sel, sel] 
  ##   
  ##   if( cov.ehat ){
  ##     Atildet <- ( -Valpha.objects[["Valpha.ilcf"]] %*% t(PpXQt) )[, TT] / eta
  ##     B <- ( PpXQt + QpXS %*% t ( XX ) )[TT, TT]
  ##     Bii <- diag( B )
  ##   }
  ##   
  ##   if( cov.ehat.p.bhat ){
  ##     Ftildet <- ( Valpha.objects[["Valpha.ucf"]] + 
  ##       ( Valpha.objects[["Valpha.ilcf"]] %*% M.inverse[ sel,-sel] %*% t(XX) ) / eta )[, TT]
  ##     Gt <- ( QpXS %*% t(XX) )[TT, TT]
  ##     Gii <- diag( Gt )
  ##   }
  ##   
  ##   if( extended.output ){
  ##     
  ##     sill <- nugget * eta
  ##     E.Dbar <- b
  ##     
  ##     ##  ... aggregate elements of D  for replicated observations
  ##     
  ##     if( sum( duplicated( TT ) ) > 0 ){
  ##       TtT <- as.vector( table( TT ) )
  ##       E.Dbar  <- E.Dbar * TtT
  ##     }
  ##     
  ##     ## compute matrices Mbetastar, Mbstar
  ##     ## debugging status: ok
  ##     
  ##     Astar <- E.Dbar * Valpha.objects[["Valpha"]]
  ##     diag( Astar ) <- diag( Astar ) + 1. / eta
  ##     Astar <- try( solve( Astar ), silent = TRUE )
  ##     if( identical( class( Astar ), "try-error" ) ){
  ##       result[["error"]] <- TRUE
  ##       return( result )            
  ##     }
  ##     
  ##     aux1 <- t(XX) %*% Astar
  ##     aux2 <- E.Dbar * XX
  ##     Mbetastar <- try( solve( aux1 %*% aux2 ), silent = TRUE )
  ##     if( identical( class( Mbetastar ), "try-error" ) ){
  ##       result[["error"]] <- TRUE
  ##       return( result )            
  ##     }
  ##     Mbetastar <- Mbetastar %*% aux1
  ##     
  ##     Mbstar <- -aux2 %*% Mbetastar
  ##     diag( Mbstar ) <- diag( Mbstar ) + 1.
  ##     Mbstar <- Valpha.objects[["Valpha"]] %*% Astar %*% Mbstar
  ##     
  ##     Mestar.t <- (t( XX %*% Mbetastar + Mbstar ))[TT,TT]
  ##     
  ##   }
  ##   
  ##   
  ##   
  ##   ##  compute covariance matrix ...
  ##   
  ##   ##   ##  zur Kontrolle: Gleichung 37, Vorschlag HRK vom 18. Januar 2011
  ##   ##   
  ##   ##   t.bla <- diag( rweights ) %*% Valpha.objects[["Valpha"]] %*% diag( rweights ) * eta + diag( rweights )
  ##   ##   
  ##   ##   t.cov <- nugget * M.inverse %*% rbind(
  ##   ##       cbind( t.bla, t.bla %*% XX ),
  ##   ##       cbind( t(XX) %*% t.bla, t(XX) %*% t.bla %*% XX )
  ##   ##   ) %*% M.inverse
  ##   ##   
  ##   ##   ##  zur Kontrolle: Gleichungen 19, 21 & 22, Entwurf Paper vom 27. Februar 2012
  ##   ##   
  ##   ##   t.A <- b * Valpha.objects[["Valpha"]]
  ##   ##   diag( t.A ) <- diag( t.A ) + 1 / eta
  ##   ##   t.A <- solve( t.A )
  ##   ##   
  ##   ##   t.B <- solve( b* t( XX ) %*% t.A %*% XX )
  ##   ##   
  ##   ##   M.beta <- t.B %*% t(XX) %*% t.A
  ##   ##   
  ##   ##   M.b <- -b * XX %*% M.beta
  ##   ##   diag( M.b ) <- diag( M.b ) + 1.
  ##   ##   M.b <- Valpha.objects[["Valpha"]] %*% t.A %*% M.b
  ##   ##   
  ##   ##   cov.Db.sigmapsi <- b^2 * nugget * eta * Valpha.objects[["Valpha"]]
  ##   ##   diag( cov.Db.sigmapsi ) <- diag( cov.Db.sigmapsi ) + nugget * a
  ##   ##   
  ##   ##   ## Kontrolle der Gleichungen 19, 21 & 22, Entwurf Paper vom 24. Februar 2012
  ##   ##   
  ##   ##   t.A <- b * Valpha.objects[["Valpha"]] * eta
  ##   ##   diag( t.A ) <- diag( t.A ) + 1.
  ##   ##   t.A <- Valpha.objects[["Valpha"]] %*% solve( t.A ) * eta
  ##   ##   
  ##   ##   t.B <- -b * t.A
  ##   ##   diag( t.B ) <- diag( t.B ) + 1.
  ##   ##   t.C <- t.B
  ##   ##   t.B <- solve( t(XX) %*% t.B %*% ( b * XX ) )
  ##   ##   
  ##   ##   M.beta <- t.B %*% t(XX) %*% t.C
  ##   ##   
  ##   ##   M.b <- -(b * XX) %*% M.beta
  ##   ##   diag( M.b ) <- diag( M.b ) + 1.
  ##   ##   M.b <- t.A %*% M.b
  ##   ##   
  ##   ##   cov.Db.sigmapsi <- b^2 * nugget * eta * Valpha.objects[["Valpha"]]
  ##   ##   diag( cov.Db.sigmapsi ) <- diag( cov.Db.sigmapsi ) + nugget * a
  ##   
  ##   ##  ... of bhat (debugging status: ok)
  ##   
  ##   if( cov.bhat ){
  ##     
  ##     
  ##     aux <- -ValphaiP / eta
  ##     diag( aux ) <- diag( aux ) + 1
  ##     
  ##     if( full.cov.bhat ){
  ##       
  ##       ##  full matrix
  ##       
  ##       result[["cov.bhat"]] <- nugget * (
  ##         crossprod( Valpha.objects[["Valpha.ucf"]] %*% aux ) * eta +
  ##         crossprod( sqrtD * PpXQt )
  ##       )
  ##       attr( result[["cov.bhat"]], "struc" ) <- "sym"
  ##       
  ##       ##       ##  zur Kontrolle: Kovarianzmatrix UK-Vorhersage
  ##       ##       
  ##       ##       t.V <- ( nugget * eta ) * Valpha.objects[["Valpha"]]
  ##       ##       t.Sigma <- t.V + nugget * diag( n )
  ##       ##       t.iSigma <- solve( t.Sigma )
  ##       ##       
  ##       ##       t.cov.bhat <- t.V %*% t.iSigma %*% t.V - t.V %*% t.iSigma %*% XX %*% solve( 
  ##       ##           t( XX ) %*% t.iSigma %*% XX 
  ##       ##       ) %*% t(XX) %*% t.iSigma %*% t.V
  ##       ##       
  ##       ##       print( summary( c( result[["cov.bhat"]] - t.cov.bhat ) ) )
  ##       ##             
  ##       ##       ## zur Kontrolle: Gleichung 21, Entwurf Paper vom 9. Februar
  ##       ##       
  ##       ##       t.cov <- M.b %*% cov.Db.sigmapsi %*% t( M.b )
  ##       ##       
  ##       ##       cat( "bhat\n" )
  ##       ##       print( summary( c( result[["cov.bhat"]] - t.cov ) ) )
  ##       
  ##     } else {
  ##       
  ##       ##  diagonal elements only
  ##       
  ##       result[["cov.bhat"]] <- nugget * (
  ##         colSums( 
  ##           drop( 
  ##             Valpha.objects[["Valpha.ucf"]] %*% aux 
  ##           )^2
  ##         ) * eta + colSums( (sqrtD * PpXQt)^2 )
  ##       )
  ##       
  ##       ##       print( summary( result[["cov.bhat"]] - diag( t.cov.bhat ) ) )
  ##       
  ##     }
  ##     
  ##   }
  ##   
  ##   
  ##   ##  ... of betahat (debugging status: ok)
  ##   
  ##   if( cov.betahat ){
  ##     
  ##     result[["cov.betahat"]] <- nugget * (
  ##       crossprod(
  ##         Valpha.objects[["Valpha.ilcf"]] %*% M.inverse[sel, -sel]
  ##       ) / eta + crossprod( sqrtD * QpXS )
  ##     )
  ##     
  ##     attr( result[["cov.betahat"]], "struc" ) <- "sym"
  ##   
  ##     
  ##     ##     ##  zur Kontrolle: Kovarianzmatrix der GLS Schaetzung
  ##     ##     
  ##     ##     t.V <- nugget * eta * Valpha.objects[["Valpha"]]
  ##     ##     t.Sigma <- t.V + nugget * diag( n )
  ##     ##     t.iSigma <- solve( t.Sigma )
  ##     ##     t.cov.betahat <- solve( t(XX) %*% t.iSigma %*% XX )
  ##     ##     print( summary( c( result[["cov.betahat"]] - t.cov.betahat ) ) )
  ##     ##     
  ##     ##     ##  zur Kontrolle: Gleichung 22, Entwurf Paper Februar 2012
  ##     ##     
  ##     ##     t.cov <- M.beta %*% cov.Db.sigmapsi %*% t( M.beta )
  ##     ##     cat( "betahat\n" )
  ##     ##     print( summary( c( result[["cov.betahat"]] - t.cov ) ) )
  ##     
  ##   }
  ##   
  ##   
  ##   ##  ... of bhat and betahat (debugging status: ok)
  ##   
  ##   if( cov.bhat.betahat ){
  ##     
  ##     aux <- t(ValphaiP) / eta
  ##     diag( aux ) <- diag( aux ) - 1
  ##     result[["cov.bhat.betahat"]] <- nugget * (
  ##       aux %*% M.inverse[sel,-sel] +
  ##       crossprod( PpXQt, TtDT * QpXS )
  ##     )
  ##     
  ##     ##     print( summary( result[["cov.bhat.betahat"]] ) )
  ##     ##     print( summary( t.cov[(1:n), -(1:n)] ) )
  ##     ##     print( summary( c( result[["cov.bhat.betahat"]] - t.cov[(1:n), -(1:n)] ) ) )
  ##     
  ##   }
  ##   
  ##   
  ##   ##  ... of delta.z = (z - bhat)  (debugging status: ok)
  ##   
  ##   if( cov.delta.bhat ){
  ##     
  ##     if( full.cov.delta.bhat ){
  ##       
  ##       ##  full matrix
  ##       
  ##       result[["cov.delta.bhat"]] <- nugget * ( 
  ##         M.inverse[sel, sel] %*% ValphaiP / eta +
  ##         crossprod( sqrtD * PpXQt )
  ##       )
  ##       dimnames( result[["cov.delta.bhat"]] ) <- list(
  ##         rownames( XX ), rownames( XX )
  ##       )
  ##       attr( result[["cov.delta.bhat"]], "struc" ) <- "sym"
  ##       
  ##       ##       ##  zur Kontrolle: Kovarianzmatrix UK-Vorhersagefehler
  ##       ##       
  ##       ##       t.V <- nugget * eta * Valpha.objects[["Valpha"]]
  ##       ##       t.Sigma <- t.V + nugget * diag( n )
  ##       ##       t.iSigma <- solve( t.Sigma )
  ##       ##       
  ##       ##       t.cov.delta.bhat <- t.V - t.V %*% t.iSigma %*% t.V + t.V %*% t.iSigma %*% XX %*% solve( 
  ##       ##         t( XX ) %*% t.iSigma %*% XX 
  ##       ##       ) %*% t(XX) %*% t.iSigma %*% t.V
  ##       ##       
  ##       ##       print( summary( c( result[["cov.delta.bhat"]] - t.cov.delta.bhat ) ) )
  ##       
  ##     } else {
  ##       
  ##       ##  diagonal elements only 
  ##       
  ##       result[["cov.delta.bhat"]] <- nugget * (
  ##         colSums( 
  ##           drop( 
  ##             Valpha.objects[["Valpha.ilcf"]] %*% M.inverse[sel, sel] 
  ##           )^2
  ##         ) / eta + colSums( (sqrtD * PpXQt)^2 )
  ##       )
  ##       names( result[["cov.delta.bhat"]] ) <- rownames( XX )
  ##       
  ##       ##       print( summary( c( result[["cov.delta.bhat"]] - diag( t.cov.delta.bhat ) ) ) )
  ##   
  ##     }
  ##     
  ##     
  ##   }
  ##   
  ##   ##  ... of delta.z = (z - bhat) and betahat (debugging status: ok)
  ##   
  ##   if( cov.delta.bhat.betahat ){
  ##     
  ##     result[["cov.delta.bhat.betahat"]] <- -nugget * (
  ##       t(ValphaiP) %*% M.inverse[sel,-sel] / eta +
  ##       crossprod( PpXQt, TtDT * QpXS )
  ##     )
  ##     
  ##     ##     ## zur Kontrolle: Kovarianzmatrix UK-Vorhersagefehler und betagls
  ##     ##     
  ##     ##     t.V <- nugget * eta * Valpha.objects[["Valpha"]]
  ##     ##     t.Sigma <- t.V + nugget * diag( n )
  ##     ##     t.iSigma <- solve( t.Sigma )
  ##     ##     
  ##     ##     t.cov.delta.bhat.betahat <- t.V %*% t.iSigma %*% XX %*% solve( 
  ##     ##         t( XX ) %*% t.iSigma %*% XX
  ##     ##     )
  ##     ##     
  ##     ##     print( summary( c( result[["cov.delta.bhat.betahat"]] - t.cov.delta.bhat.betahat ) ) )
  ##     
  ##   }
  ##   
  ##   ##  ... of residuals ehat  (debugging status: ok)
  ##   ##  vgl. Notizen zu Bias-Korrektur h.9 bis h.11
  ##   
  ##   if( cov.ehat ){
  ##     
  ##     if( full.cov.ehat ){
  ##       
  ##       ##  full matrix
  ##       
  ##       aux <- eta * crossprod( Atildet ) + 
  ##         a * crossprod( B ) - 
  ##         2 * b * B
  ##       diag( aux ) <- diag( aux ) + 1
  ##       
  ##       result[["cov.ehat"]] <- nugget * aux
  ##       
  ##       dimnames( result[["cov.ehat"]] ) <- list( names.yy, names.yy )
  ##       attr( result[["cov.ehat"]], "struc" ) <- "sym"
  ##       
  ##       
  ##       ##       ##  zur Kontrolle Kovarianzmatrix von ehat fuer nicht robusten Fall
  ##       ##       
  ##       ##       t.V <- nugget * eta * Valpha.objects[["Valpha"]]
  ##       ##       t.sigma <- t.V + nugget * diag( n )
  ##       ##       t.isigma <- solve( t.sigma )
  ##       ##       t.A <- t.sigma - XX %*% solve( t(XX) %*% t.isigma %*% XX ) %*% t(XX )
  ##       ##       t.cov.ehat <- nugget^2 * t.isigma %*% t.A %*% t.isigma
  ##       ##       
  ##       ##       print( summary( c( result[["cov.ehat"]] - t.cov.ehat ) ) )
  ##       
  ##     } else {
  ##       
  ##       ##  diagonal elements only
  ##       
  ##       result[["cov.ehat"]] <- nugget * (
  ##         1 + eta * colSums( Atildet^2 ) + 
  ##         a * colSums( B^2 ) -
  ##         2 * b * Bii
  ##       )
  ##       
  ##       names( result[["cov.ehat"]] ) <- names.yy
  ##       
  ##       ##       print( summary( c( result[["cov.ehat"]] - diag(t.cov.ehat) ) ) )
  ##       
  ##     }
  ##   }
  ##   
  ##   
  ##   ##  ... of residuals ehat + bhat  (debugging status: ok)
  ##   ##  vgl. Notizen zu Bias-Korrektur h.25 bis h.28
  ##   
  ##   if( cov.ehat.p.bhat ){
  ##     
  ##     if( full.cov.ehat.p.bhat ){
  ##       
  ##       ##  full matrix
  ##       
  ##       aux <- eta * crossprod( Ftildet ) + 
  ##         a * crossprod( Gt ) - 
  ##         b * ( Gt + t(Gt) )
  ##       diag( aux ) <- diag( aux ) + 1.
  ##       
  ##       result[["cov.ehat.p.bhat"]] <- nugget * aux
  ##       
  ##       dimnames( result[["cov.ehat.p.bhat"]] ) <- list( names.yy, names.yy )
  ##       attr( result[["cov.ehat.p.bhat"]], "struc" ) <- "sym"
  ##       
  ##       
  ##       
  ##       ##  zur Kontrolle Kovarianzmatrix von ehat.p.bhat fuer nicht robusten Fall
  ##       
  ##       t.V <- nugget * eta * Valpha.objects[["Valpha"]]
  ##       t.sigma <- t.V + nugget * diag( n )
  ##       t.isigma <- solve( t.sigma )
  ##       t.A <- t.sigma - XX %*% solve( t(XX) %*% t.isigma %*% XX ) %*% t(XX )
  ##       ##       print( summary( c( result[["cov.ehat.p.bhat"]] -  t.A ) ) )
  ##       
  ##     } else {
  ##       
  ##       ##  diagonal elements only
  ##       
  ##       result[["cov.ehat.p.bhat"]] <- nugget * (
  ##         1. + eta * colSums( Ftildet^2 ) + 
  ##           a * colSums( Gt^2 ) -
  ##           2 * b * Gii
  ##       )
  ##       
  ##       names( result[["cov.ehat.p.bhat"]] ) <- names.yy
  ##       
  ##       ##      print( summary( c( result[["cov.ehat.p.bhat"]] - diag( t.A ) ) ) )
  ##       
  ##     }
  ##     
  ##     
  ##   }
  ##   
  ##   ##  matrix required for computing covariances between kriging predictions and
  ##   ##  prediction targets y (debugging status: ok)
  ##   ##  vgl. Notizen zu robustem Kriging, S. 3--4
  ##   
  ##   if( cov.pred.target ){
  ##     
  ##     aux <- -t(ValphaiP) / eta
  ##     diag( aux ) <- diag( aux ) + 1.
  ##     result[["cov.pred.target"]] <- rbind(
  ##       aux,
  ##       -M.inverse[ -sel, sel] %*% Valpha.objects[["Valpha.inverse"]] / eta
  ##     )
  ##   }
  ##   
  ##   ##   ##  covariance matrix of xihat
  ##   ##   
  ##   ##   l1 <- 1. / ( eta * expectations["dpsi"] )
  ##   ##   l2 <- expectations["psi2"] / ( eta * expectations["dpsi"]^2 )
  ##   ##   
  ##   ##   aux1 <- l1 * Valpha.inverse.Palpha
  ##   ##   diag( aux1 ) <- diag( aux1 ) + TtT
  ##   ##   aux1 <- try( solve( aux1 ), silent = TRUE )
  ##   ##   
  ##   ##   aux2 <- TtT * t( TtT * Valpha.objects[["Valpha"]] )
  ##   ##   diag( aux2 ) <- diag( aux2 ) + l2 * TtT
  ##   ##   
  ##   ##   t.cov.xihat <- aux1 %*% aux2 %*% aux1 * eta * nugget
  ##   ##   
  ##   ##   t.cov.bhat <- Palpha %*% t.cov.xihat %*% t( Palpha )
  ##   ##   
  ##   ##   print( summary( c( t.cov.bhat - result[["cov.bhat"]] ) ) )
  ##   ##   
  ##   ##   stop()
  ##   
  ##   sapply(
  ##     names( result ),
  ##     function( i, old, new  ){ 
  ##       print(i)
  ##       print( summary( c( old[[i]] - new[[i]] ) ) )
  ##     },
  ##     old = result, new =result.new
  ##   )
  
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
  TtT   <- rep( 1, length( Ttpsi ) )      
  
  if( sum( duplicated( TT ) > 0 ) ){
    Ttpsi <- as.vector( tapply( Ttpsi, factor( TT ), sum ) )
    TtT   <- as.vector( table( TT ) )
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
    verbose
  )
{
  
  ## 2013-02-04 AP solving estimating equations for xi
  ## 2013-06-03 AP handling design matrices with rank < ncol(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
   
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
  
  result[["Valpha.inverse.Palpha"]] <- Valpha.inverse %*% result[["Palpha"]]
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
    x, variogram.model, param, 
    aniso, 
    irf.models = georob.control()[["irf.models"]],
    verbose
  )
{
  
  ##  Function computes the generalized correlation (matrix) for the lag
  ##  distances in x.  The result is a generalized correlation matrix
  ##  that is positive definite.
  
  ##  cf. HRK's notes of 2011-06-17 on "Robust Kriging im intrinsischen
  ##  Fall"
  
  ##  2011-12-27 ap
  ##  2012-02-07 AP modified for geometrically anisotropic variograms
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-03-05 AP changes for version 3 of RandomFields
  
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
  
  Valpha0 <- try(
    -RFvariogram(
      x = x, model = model.list, dim = NCOL( x ), grid = FALSE
    ),
    silent = TRUE
  )

  ## functions of version 2.xx of RandomFields
  
  ##   RFoldstyle()
  ##   Valpha0 <- try(
  ##     -Variogram(
  ##       x, 
  ##       model = model.list
  ##     ),
  ##     silent = TRUE
  ##   )
  
  if( !(identical( class( Valpha0 ), "try-error" ) || any( is.na( Valpha0 ) )) ){
    
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
    result[["Valpha0"]]          <- Valpha0 + gcr.constant.Valpha0  # correlation matrix that does not contain spatial nugget
    result[["Valpha"]]           <- Valpha  + gcr.constant.Valpha   # correlation matrix that includes spatial nugget
    
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
    psi.function, dpsi.function, tuning.psi, tuning.psi.nr, 
    irwls.initial, irwls.maxiter, irwls.reltol,
    compute.xihat = TRUE,
    compute.Q,
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
  
  if( 
    isTRUE( 
      all.equal( c( param, aniso ), c( lik.item[["param"]], lik.item[["aniso"]][["aniso"]] ) ) 
    ) && !( compute.Q && is.null( lik.item[["Q"]] ) )
  ) {
    
    ##  return the result of the previous call if the variogram
    ##  parameters are the same
    
    return( lik.item )
    
  } else {
    
    ##  compute updates of required likelihood items if the
    ##  variogram parameters differ
    
    ## check whether variogram parameters are within reasonalble bounds and
    ## return an error otherwise
    
    lik.item[["Valpha"]][["error"]] <- TRUE
    
    if( length( c( param, aniso ) ) && any( c( param, aniso ) > safe.param ) ){
      if( verbose > 1 ){
        t.param <- param
        if( !lik.item[["aniso"]][["isotropic"]] ) t.param <- c( t.param, aniso )
        if( verbose > 1 ) {
          cat( "\n" )
          print( signif( t.param ) )
        }
      }
      return( lik.item )  
    }
    
    ## check whether extra variogram parameters are within allowed bounds and
    ## return an error otherwise
    
    ep <- param.names( model = variogram.model )
    param.bounds <- param.bounds( variogram.model, NROW( lag.vectors ), param )
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
      
#       lik.item[["aniso"]][["rotmat"]] <- with( 
#         lik.item[["aniso"]][["sincos"]],
#         rbind(
#           c(             cp*so,             cp*co,       sp ),
#           c( -cz*co + sz*sp*so,  co*sz*sp + cz*so,   -cp*sz ),
#           c( -co*sz - cz*sp*so, -cz*co*sp + sz*so,    cz*cp )
#         )[ 1:n, 1:n, drop = FALSE ]
#       )
      
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
    
    t.param <- lik.item[["param"]]
    if( !lik.item[["aniso"]][["isotropic"]] ) t.param <- c( t.param, lik.item[["aniso"]][["aniso"]] )
    if( verbose > 1 ) {
      cat( "\n" )
      print( signif( t.param ) )
    }
    
    ##  calculate generalized correlation matrix, its inverse and its
    ##  inverse cholesky factor
    
    cormat <- gcr(
      x = lag.vectors, variogram.model = variogram.model, param = param, 
      aniso = lik.item[["aniso"]], verbose = verbose
    )
    
    if( cormat[["error"]] ) return( lik.item )
    
    t.vchol <- try( chol( cormat[["Valpha"]] ), silent = TRUE )
    
#     print( identical( class( t.vchol ), "try-error" ) )
#     
    if( !identical( class( t.vchol ), "try-error" ) ) {

      lik.item[["Valpha"]][["error"]]         <- FALSE
      lik.item[["Valpha"]][["gcr.constant"]]  <- cormat[["gcr.constant"]]
      lik.item[["Valpha"]][["Valpha"]]            <- cormat[["Valpha"]]
      lik.item[["Valpha"]][["Valpha0"]]           <- cormat[["Valpha0"]]
      lik.item[["Valpha"]][["Valpha.ucf"]]        <- unname( t.vchol )
      
      lik.item[["Valpha"]][["Valpha.ilcf"]]       <- try(
        t( 
          backsolve( 
            t.vchol, 
            diag( nrow( lik.item[["Valpha"]][["Valpha"]] ) ), 
            k = nrow( lik.item[["Valpha"]][["Valpha"]] ) 
          ) 
        ),
        silent = TRUE
      )
      if( identical( class( lik.item[["Valpha"]][["Valpha.ilcf"]] ), "try-error" ) ) {
        lik.item[["Valpha"]][["error"]] <- TRUE
        return( lik.item )
      }
      
      lik.item[["Valpha"]][["Valpha.inverse"]]    <- t( lik.item[["Valpha"]][["Valpha.ilcf"]] ) %*% lik.item[["Valpha"]][["Valpha.ilcf"]]
      
      attr( lik.item[["Valpha"]][["Valpha"]], "struc" )         <- "sym"
      attr( lik.item[["Valpha"]][["Valpha0"]], "struc" )        <- "sym"
      attr( lik.item[["Valpha"]][["Valpha.inverse"]], "struc" ) <- "sym"
      attr( lik.item[["Valpha"]][["Valpha.ucf"]], "struc" )     <- "ut"
      attr( lik.item[["Valpha"]][["Valpha.ilcf"]], "struc" )    <- "lt"
      
    } else {
      
      return( lik.item )          ##  an error occurred
      
    }
    
    ##  estimate fixed and random effects (xihat, betahat, bhat,
    ##  residuals )
    
    ##  either take initial guess of betahat and bhat for the current
    ##  irwls iteration from initial.object or from previous iteration
    
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
      lik.item[["param"]]["nugget"], lik.item[["eta"]], lik.item[["Valpha"]][["Valpha.inverse"]],
      verbose
    )      
    
    if( lik.item[["effects"]][["error"]] ) return( lik.item )     ##  an error occurred
    
    ##  compute Q matrix and its Cholesky factor (required for
    ##  non-robust REML estimate)
    
    if( compute.Q ) {
      
      ##  diagonal matrix with second derivative of rho function
      
      TtDT <- dpsi.function( 
        lik.item[["effects"]][["residuals"]] / sqrt( lik.item[["param"]]["nugget"] ), 
        tuning.psi = tuning.psi
      )
      
      ##  aggregate elements of D for replicated observations
      
      if( sum( duplicated( TT ) > 0 ) ){
        TtDT <- as.vector( tapply( TtDT, factor( TT ), sum ) )
      }
      
      TtDTX <- TtDT * XX
      aux <-  lik.item[["Valpha"]][["Valpha.inverse"]] / lik.item[["eta"]]
      diag( aux ) <- diag( aux ) + TtDT
      
      ##  compute matrix Q
      
      Q <- rbind( 
        cbind( aux,      TtDTX                 ),
        cbind( t(TtDTX), crossprod( XX, TtDTX) )
      ) / lik.item[["param"]]["nugget"]
      
      lik.item[["Q"]] <- list( error = TRUE )
      
      if( rankdef.x ){
        
        ## compute log(pseudo.det(Q)) and (Moore-Penrose) pseudo inverse of Q by svd
        
        lik.item[["Q"]][["error"]] <- FALSE
        s <- svd( Q )
        lik.item[["Q"]][["log.det.Q"]] <- sum( log( s[["d"]][s[["d"]] / max( s[["d"]] ) > min.condnum] ) )
        s[["d"]] <- ifelse( s[["d"]] / max( s[["d"]] ) <= min.condnum, 0., 1. / s[["d"]] )
        lik.item[["Q"]][["Q.inverse"]] <- s[["v"]] %*% ( s[["d"]] * t( s[["u"]] ) )
        
      } else {
        
        ##  compute log(det(Q)) and inverse of Q by cholesky decomposition
        
        t.chol <- try( chol( Q ), silent = TRUE )
        
        if( !identical( class( t.chol ), "try-error" ) ) {
          
          lik.item[["Q"]][["error"]] <- FALSE
          lik.item[["Q"]][["log.det.Q"]] <- 2 * sum( log( diag( t.chol) ) )
          lik.item[["Q"]][["Q.inverse"]] <- chol2inv( t.chol )
          
        } else {
          
          return( lik.item )
          
        }
        
      }
      
    }
    
    ##  store updated lik.item object
    
    assign( "lik.item", lik.item, pos = as.environment( envir ) )
    
    return( lik.item )
    
  }
  
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
  ##  2012-01-24 ap cauchytbm and lgd1 models added
  ##  2012-01-25 ap extra model parameter with same names as in Variogram{RandomFields}
  ##  2012-02-07 AP modified for geometrically anisotropic variograms
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
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
    
    bessel = {
      
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
    }, ##  end case bessel
    
    cauchy = {
      
      A <- unname( param["beta"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -2 * A * hs * ( 1+hs^2 )^(-1-A)        
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        beta = {
          -( 1 + hs^2 )^(-A) * log( 1 + hs^2 )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case cauchy
    
    cauchytbm = {
      
      A <- unname( param["alpha"] )
      B <- unname( param["beta"] )
      C <- unname( param["gamma"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -( 
        B * hs^(-1+A) * (1+hs^A)^(-2-B/A) * ( A + C - (-B + C) * hs^A ) 
      ) / C
      # dgc.dhs <- ifelse(
      #   hs > 0.,
      #   -( 
      #     B * hs^(-1+A) * (1+hs^A)^(-2-B/A) * ( A + C - B * hs^A + C * hs^A) 
      #   ) / C,
      #   if( A > 1. ){
      #     0.
      #   } else if( identical( A, 1. ) ){
      #     -B * (1+C) / C
      #   } else {
      #     -Inf
      #   }
      # )
      
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        # scale = {
        #   ( B * hs^A * (1+hs^A)^(-2-B/A) * (A + C + (-B+C) * hs^A ) ) / ( C * scale )
        # },
        alpha = {
          ( B * (1+hs^A)^(-2 - B/A) * (
              -( A * hs^A * ( A + C + (-B+C) * hs^A ) * log(hs) ) + 
              ( 1 + hs^A) * (C + (-B+C) * hs^A ) * log( 1+hs^A ) 
            ) 
          ) / (A^2 * C )
        },
        # alpha = {
        #   ifelse(
        #     hs > 0.,
        #     ( B * (1+hs^A)^(-2 - B/A) * (
        #         -( A * hs^A * ( A + C + (-B+C) * hs^A ) * log(hs) ) + 
        #         ( 1 + hs^A) * (C + (-B+C) * hs^A ) * log( 1+hs^A ) 
        #       ) 
        #     ) / (A^2 * C ),
        #     0.
        #   )
        # },
        beta = {
          ( -( A * hs^A) - (C + (-B+C) * hs^A ) * log( 1+hs^A ) ) / 
          ( A*C * (1+hs^A)^( (A+B)/A ) )
        },
        gamma = {
          ( B * hs^A ) / ( C^2 * (1+hs^A)^( (A+B)/A) )
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case cauchytbm
    
    circular = {
      
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
    }, ##  end case circular
    
    cubic = {
      
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
    }, ##  end case cubic
    
    dagum = {
      
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
    }, ##  end case dagum
    
    dampedcosine = {
      
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
    }, ##  end case dampedcosine
    
    DeWijsian = {
      
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
    }, ##  end case DeWijsian
    
    
    exponential = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -exp( -hs )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case exponential
    
    fractalB = {
      
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
    }, ##  end case fractalB
    
    gauss = {
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- -2 * hs / exp( hs^2 )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        dgc.dhs * dhs.daniso
      )
    }, ##  end case gauss
    
    genB = {
      
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
    }, ##  end case genB
    
    gencauchy = {
      
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
    }, ##  end case gencauchy
    
    gengneiting = {
      
      
      A <- unname( param["n"] )
      B <- unname( param["alpha"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- if( identical( A, 1 ) ){
        -( (1+B) * (2+B) * (1-hs[sel])^B * hs[sel] )
      } else if( identical( A, 2 ) ){
        -( (3+B) * (4+B) * (1-hs[sel])^(1+B) * hs[sel] * ( 1 + hs[sel] + B*hs[sel]) ) / 3.
      } else if( identical( A, 3 ) ){
        -( 
          (5+B) * (6+B) * (1-hs[sel])^(2+B) * hs[sel] * ( 3 + 3 * (2+B) * hs[sel] + (1+B) * (3+B) * hs[sel]^2 ) 
        ) / 15.
      } else {
        stop( "gengneiting model undefined for 'n' != 1:3" )
      }
      
      result <- rep( 0., length( hs ) )
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        alpha = {
          result[sel] <- if( identical( A, 1 ) ){
            (1-hs[sel])^(1+B) * ( hs[sel] + (1 + hs[sel] + B*hs[sel]) * log( 1-hs[sel]) )
            
          } else if( identical( A, 2 ) ){
            (
              (1-hs[sel])^(2+B) * ( 
                hs[sel] * ( 3 + 2 * (2+B) *hs[sel] ) + 
                ( 3 + 3 * ( 2+B) * hs[sel] + ( 1+B) * (3+B) * hs[sel]^2 ) * log( 1-hs[sel] )
              )
            ) / 3.
          } else if( identical( A, 3 ) ){
            ( 
              (1-hs[sel])^(3+B) * ( 
                hs[sel] * ( 15 + hs[sel] * ( 36 + 23*hs[sel] + 3 * B * ( 4 + (6+B)*hs[sel] ) ) ) + 
                ( 15 + 15 * (3+B) * hs[sel] + ( 45 + 6 * B * (6+B) ) * hs[sel]^2 + (1+B) * (3+B) * (5+B) * hs[sel]^3 ) * 
                log( 1-hs[sel]) 
              ) 
            ) / 15.
          }
          result
        },
        dgc.dhs * dhs.daniso
      )
      
      
    }, ##  end case Gengneiting
    
    gneiting = {
      
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
    }, ##  end case gneiting
    
    lgd1 = {
      
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
    }, ##  end case lgd1
    
    matern = {
      
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
    }, ##  end case matern
    
    penta = {
      
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
    }, ##  end case penta
    
    power = {
      
      A <- unname( param["a"] )
      
      ## derivative with of generalized covariance with respect to
      ## scaled lag distance
      
      dgc.dhs <- rep( 0., length( hs ) )
      sel <- hs < 1.
      dgc.dhs[sel] <- -(A * (1-hs[sel])^(-1+A))         
      
      switch(
        d.param,
        scale = dgc.dhs * dhs.dscale,
        a = {
          result <- rep( 0., length( hs ) )
          sel <- hs < 1.
          result[sel] <- ( 1 - hs[sel] )^A * log( 1 - hs[sel] )
          result
        },
        dgc.dhs * dhs.daniso
      )
    }, ##  end case power
    
    qexponential = {
      
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
    }, ##  end case qexponential
    
    
    spherical = {
      
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
    }, ##  end case spherical
    
    stable = {
      
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
    }, ##  end case stable
    
    wave = {
      
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
    
    whittle = {
      
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
    slv,
    envir,
    variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, dpsi.function, 
    tuning.psi, tuning.psi.nr,
    irwls.initial, irwls.maxiter, irwls.reltol,
    force.gradient,
    expectations,
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

  ##  get lik.item
  
  lik.item <- prepare.likelihood.calculations(
    envir,
    adjustable.param, variogram.model, fixed.param, param.name, aniso.name,
    param.tf, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, rankdef.x, yy, TT, xihat, 
    psi.function, dpsi.function, tuning.psi, tuning.psi.nr,
    irwls.initial, irwls.maxiter, irwls.reltol,
    compute.xihat = TRUE, compute.Q = FALSE,
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
      extended.output = FALSE,
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
    
    ##  initialize estimating equations
    
    eeq.emp <- rep( NA, length( adjustable.param ) )
    names( eeq.emp ) <- names( adjustable.param )
    
    eeq.exp <- rep( NA, length( adjustable.param ) )
    names( eeq.exp ) <- names( adjustable.param )
    
    ##  estimation equation for nugget
    
    if( "nugget" %in% names( adjustable.param ) ) {
      
      ##  compute trace of Cov[ psi( residuals/sqrt(nugget) ) ]
      
      eeq.exp["nugget"] <- sum( 
        diag( 
          lik.item[["Valpha"]][["Valpha.inverse"]] %*%             
          ( 1/TtT * lik.item[["Valpha"]][["Valpha.inverse"]] ) %*% 
          r.cov[["cov.bhat"]] 
        ) 
      )
      eeq.emp["nugget"] <- sum( 
        ( lik.item[["effects"]][["z.star"]] )^2 / TtT
      )
      
    }
    
    ##  estimation equation for spatial nugget
    
    if( "snugget" %in% names( adjustable.param ) ) {
      
      ##  compute trace( Valpha^-1 Cov[bhat] )
      
      eeq.exp["snugget"] <- sum(
        rowSums( 
          (lik.item[["Valpha"]][["Valpha.inverse"]] %*% lik.item[["Valpha"]][["Valpha.inverse"]] ) * 
          r.cov[["cov.bhat"]]
        )
      )
      eeq.emp["snugget"] <- sum( lik.item[["effects"]][["z.star"]]^2 )
      
    }
    
    ##  estimation equation for variance
    
    if( "variance" %in% names( adjustable.param ) ) {
      
      ##  compute trace( Valpha^-1 Cov[bhat] )
      
      eeq.exp["variance"] <- sum(
        rowSums( 
          ( lik.item[["Valpha"]][["Valpha.inverse"]] %*% lik.item[["Valpha"]][["Valpha0"]] %*% lik.item[["Valpha"]][["Valpha.inverse"]] ) * 
          r.cov[["cov.bhat"]]
        )
      )
      eeq.emp["variance"] <- sum( 
        lik.item[["effects"]][["z.star"]] * drop( lik.item[["Valpha"]][["Valpha0"]] %*% lik.item[["effects"]][["z.star"]] )
      )
      
    }
    
    ##  estimation equations for scale, extra variogram and anisotropy
    ##  parameters
    
    extra.par <- names( adjustable.param )[ !( 
      names( adjustable.param ) %in% c( "variance", "snugget", "nugget" )
    )]
    
    for( t.i in extra.par ){
      
      ##  compute trace( Valpha^-1 * dValpha/dalpha * Valpha^-1 * Cov[bhat] )
      
      dValpha <- dcorr.dparam(
        x = lag.vectors, variogram.model = variogram.model, param = lik.item[["param"]], 
        d.param = t.i,
        aniso = lik.item[["aniso"]],
        verbose = verbose
      )
      ##       if( identical( class( dValpha ), "try-error" ) ){
      ##         if( verbose > 0 ) cat( "error in dcorr.dparam\n\n" )
      ##         t.result <- rep( Inf, length( adjustable.param ) )
      ##         names( t.result ) <- names( adjustable.param )
      ##         return( t.result )
      ##       }
      
      eeq.exp[t.i] <- sum(
        rowSums( 
          (lik.item[["Valpha"]][["Valpha.inverse"]] %*% dValpha %*% lik.item[["Valpha"]][["Valpha.inverse"]]) * 
          r.cov[["cov.bhat"]]
        )
      )
      eeq.emp[t.i] <- sum( 
        lik.item[["effects"]][["z.star"]] * drop( dValpha %*% lik.item[["effects"]][["z.star"]] )
      )
      
    }
    
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
    
    if( slv ){
      return( eeq.emp / eeq.exp - 1. )
    } else {
      res <- sum( (eeq.emp / eeq.exp - 1.)^2 )
      if( verbose > 1 ) cat( 
        "  sum(EEQ^2)         :",
        format( 
          signif( res, digits = 7 ), 
          scientific = TRUE, width = 14
        ), "\n", sep = "" 
      )
      return( res )
    }
    
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
##     psi.function, dpsi.function, tuning.psi, tuning.psi.nr,
##     irwls.initial, irwls.maxiter, irwls.reltol,
##     compute.xihat = FALSE, compute.Q = FALSE,
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
##     extended.output = FALSE,
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
    tuning.psi, tuning.psi.nr, 
    irwls.initial, irwls.maxiter, irwls.reltol,
    verbose,
    ...
  )
{
  
  ## function computes laplace approximation to negative restricted
  ## loglikelihood
  
  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2013-06-03 AP changes for estimating xihat
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  
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
    psi.function, dpsi.function, tuning.psi, tuning.psi.nr,
    irwls.initial, irwls.maxiter, irwls.reltol,
    compute.xihat = TRUE, compute.Q = TRUE,
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
  
  ##  compute laplace approximation of negative restricted loglikelihood
  
  t.dim <- dim( XX )
  
  term1 <- -0.5 * (
    diff( t.dim ) * log( 2 * pi ) + t.dim[1] * log( 1./lik.item[["eta"]] ) 
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
    tuning.psi, tuning.psi.nr,
    irwls.initial, irwls.maxiter, irwls.reltol,
    force.gradient,
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
    psi.function, dpsi.function, tuning.psi, tuning.psi.nr,
    irwls.initial, irwls.maxiter, irwls.reltol,
    compute.xihat = TRUE, compute.Q = TRUE,
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
    
    n <- nrow( XX ); sel.z <- 1:n
    
    sigma <- sqrt( lik.item[["param"]]["nugget"] )
    
    std.res <- lik.item[["effects"]][["residuals"]] / sigma
    
    psi   <-   psi.function( std.res, tuning.psi )
    dpsi  <-  dpsi.function( std.res, tuning.psi )
    d2psi <- d2psi.function( std.res, tuning.psi )
    
    Qi <- lik.item[["Q"]][["Q.inverse"]]
    Qi.1.Valphai <- Qi[, sel.z] %*% lik.item[["Valpha"]][["Valpha.inverse"]]
    Valphai.Valpha0 <- lik.item[["Valpha"]][["Valpha.inverse"]] %*% lik.item[["Valpha"]][["Valpha0"]]
    
    
    ##  initialize gradient
    
    r.gradient <- numeric( 0 )
    
    
    ##  compute partial derivative of Laplace approximation of restricted
    ##  log-likelihood with respect to nugget
    
    if( "nugget" %in% names( adjustable.param ) ) {
      
      ##  derivative of bhat and betahat
      
      aux <- tapply( 
        sigma * psi + dpsi * lik.item[["effects"]][["residuals"]],
        factor( TT ),
        sum
      )
      dzbeta <- -0.5 * drop(
        ( Qi[, sel.z] + Qi[, -sel.z] %*% t(XX) ) %*% aux 
      ) / lik.item[["param"]]["nugget"]^2
      
      ##              print( dzbeta[501]/
      ##                  deriv.fwd.tf[[param.tf["nugget"]]]( lik.item[["param"]]["nugget"] ) 
      ##              )
      
      ##  derivative of log( det( Q ) )
      
      Dmat <- dpsi + d2psi * (
        0.5 * std.res + 
        sigma * ( dzbeta[sel.z] + drop( XX %*% dzbeta[-sel.z] ) )[TT]
      )
      Dmat <- as.vector( tapply( Dmat, factor( TT ), sum ) )
      
      DX <- Dmat * XX
      
      dlogdetQ <- -sum( 
        Qi * rbind( 
          cbind( diag( Dmat ),                DX   ),
          cbind(        t(DX), crossprod( XX, DX ) )
        )
      ) / lik.item[["param"]]["nugget"]^2
      
      ##              print( dlogdetQ/
      ##                  deriv.fwd.tf[[param.tf["nugget"]]]( lik.item[["param"]]["nugget"] ) 
      ##              )
      
      ##  derivate of U with respect to nugget
      
      Ttpsi                   <- as.vector( tapply( psi,     factor( TT ), sum ) )
      Ttstd.res <- as.vector( tapply( std.res, factor( TT ), sum ) )
      TtT                       <- as.vector( table( TT ) )
      
      dU <- -0.5 * sum( Ttpsi * Ttstd.res / TtT ) / lik.item[["param"]]["nugget"]
      
      ##              print( dU/
      ##                  deriv.fwd.tf[[param.tf["nugget"]]]( lik.item[["param"]]["nugget"] ) 
      ##              )
      
      r.gradient["nugget"] <- ( -0.5 * ( 
          n / lik.item[["param"]]["nugget"] + dlogdetQ
        ) - dU
      ) / deriv.fwd.tf[[param.tf["nugget"]]]( lik.item[["param"]]["nugget"] )
      
    }
    
    ##  compute partial derivative of Laplace approximation of restricted
    ##  log-likelihood with respect to spatial nugget
    
    if( "snugget" %in% names( adjustable.param ) ) {
      
      ##  derivative of bhat and betahat with respect to spatial nugget
      
      dzbeta <- drop( Qi.1.Valphai %*% lik.item[["effects"]][["z.star"]] ) / 
      sum( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] )^2
      
      ##              print( dzbeta[501]/
      ##                  deriv.fwd.tf[[param.tf["variance"]]]( lik.item[["param"]]["variance"] ) 
      ##              )
      
      ##  derivative of log( det( Q ) ) with respect to spatial nugget
      
      Dmat <- d2psi * sigma * (
        ( dzbeta[sel.z] + drop( XX %*% dzbeta[-sel.z] ) )[TT] 
      )
      Dmat <- as.vector( tapply( Dmat, factor( TT ), sum ) )
      
      DX <- Dmat * XX
      
      aux <- lik.item[["Valpha"]][["Valpha.inverse"]] %*% lik.item[["Valpha"]][["Valpha.inverse"]] / lik.item[["eta"]]^2
      diag( aux ) <- diag( aux ) + Dmat
      
      dlogdetQ <- -sum( 
        Qi * rbind( 
          cbind(   aux,                DX   ),
          cbind( t(DX), crossprod( XX, DX ) )
        )
      ) / lik.item[["param"]]["nugget"]^2
      
      ##              print( dlogdetQ/
      ##                  deriv.fwd.tf[[param.tf["variance"]]]( lik.item[["param"]]["variance"] ) 
      ##              )
      
      ##  derivate of U with respect to spatial nugget
      
      dU <- -0.5 * sum( lik.item[["effects"]][["z.star"]] * lik.item[["effects"]][["z.star"]] ) /
      sum( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] )^2
      
      ##              print( dU /
      ##                  deriv.fwd.tf[[param.tf["variance"]]]( lik.item[["param"]]["variance"] ) 
      ##              )
      
      r.gradient["snugget"] <- ( -0.5 * ( 
          sum( diag( lik.item[["Valpha"]][["Valpha.inverse"]] ) ) / 
          sum( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] ) + dlogdetQ
        ) - dU
      ) / deriv.fwd.tf[[param.tf["snugget"]]]( lik.item[["param"]]["snugget"] )
      
    }
    
    ##  compute partial derivative of Laplace approximation of restricted
    ##  log-likelihood with respect to variance
    
    if( "variance" %in% names( adjustable.param ) ) {
      
      ##  derivative of bhat and betahat with respect to variance
      
      dzbeta <- drop( Qi.1.Valphai %*% lik.item[["Valpha"]][["Valpha0"]] %*% lik.item[["effects"]][["z.star"]] ) / 
      sum( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] )^2
      
      ##              print( dzbeta[501]/
      ##                  deriv.fwd.tf[[param.tf["variance"]]]( lik.item[["param"]]["variance"] ) 
      ##              )
      
      ##  derivative of log( det( Q ) ) with respect to variance
      
      Dmat <- d2psi * sigma * (
        ( dzbeta[sel.z] + drop( XX %*% dzbeta[-sel.z] ) )[TT] 
      )
      Dmat <- as.vector( tapply( Dmat, factor( TT ), sum ) )
      
      DX <- Dmat * XX
      
      aux <- Valphai.Valpha0 %*% lik.item[["Valpha"]][["Valpha.inverse"]] / lik.item[["eta"]]^2
      diag( aux ) <- diag( aux ) + Dmat
      
      dlogdetQ <- -sum( 
        Qi * rbind( 
          cbind(   aux,                DX   ),
          cbind( t(DX), crossprod( XX, DX ) )
        )
      ) / lik.item[["param"]]["nugget"]^2
      
      ##              print( dlogdetQ/
      ##                  deriv.fwd.tf[[param.tf["variance"]]]( lik.item[["param"]]["variance"] ) 
      ##              )
      
      ##  derivate of U with respect to variance
      
      dU <- -0.5 * sum( 
        lik.item[["effects"]][["z.star"]] * drop( lik.item[["Valpha"]][["Valpha0"]] %*% lik.item[["effects"]][["z.star"]] ) 
      ) / sum( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] )^2
      
      ##              print( dU /
      ##                  deriv.fwd.tf[[param.tf["variance"]]]( lik.item[["param"]]["variance"] ) 
      ##              )
      
      r.gradient["variance"] <- ( -0.5 * ( 
          sum( diag( Valphai.Valpha0 ) ) / 
          sum( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] ) + dlogdetQ
        ) - dU
      ) / deriv.fwd.tf[[param.tf["variance"]]]( lik.item[["param"]]["variance"] )
      
    }
    
    ##  compute partial derivatives of Laplace approximation of
    ##  restricted log-likelihood with respect to scale, extra variogram
    ##  and anisotropy parameters
    
    t.extra.par <- !( 
      names( adjustable.param ) %in% 
      c( "variance", "nugget", "snugget" )
    )
    
    if( sum( t.extra.par ) > 0 ) {
      
      t.extra.par <- names( adjustable.param )[t.extra.par]
      
      for( t.i in t.extra.par ){
        
        dValpha <- dcorr.dparam(
          x = lag.vectors, variogram.model = variogram.model, param = lik.item[["param"]], d.param = t.i,
          aniso = lik.item[["aniso"]],
          verbose = verbose
        )
        ##         if( identical( class( dValpha ), "try-error" ) ){
        ##           if( verbose > 0 ) cat( "error in dcorr.dparam\n\n" )
        ##           t.result <- rep( Inf, length( adjustable.param ) )
        ##           names( t.result ) <- names( adjustable.param )
        ##           return( t.result )
        ##         }
        
        dValpha.Valphai <- dValpha %*% lik.item[["Valpha"]][["Valpha.inverse"]]
        
        ##  derivative of bhat and betahat
        
        dzbeta <- drop( 
          Qi.1.Valphai %*% dValpha.Valphai %*% lik.item[["effects"]][["bhat"]] 
        ) * lik.item[["param"]]["variance"] / sum( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] )^2
        
        ##              print( dzbeta[501]/
        ##                  deriv.fwd.tf[[param.tf[t.i]]]( lik.item[["param"]][t.i] ) 
        ##              )
        
        ##  derivative of log( det( Q ) )
        
        Dmat <- d2psi * sigma * ( 
          ( dzbeta[sel.z] + drop( XX %*% dzbeta[-sel.z] ) )[TT]
        )
        Dmat <- as.vector( tapply( Dmat, factor( TT ), sum ) )
        
        DX <- Dmat * XX
        
        aux <- lik.item[["param"]]["variance"] * lik.item[["Valpha"]][["Valpha.inverse"]] %*% dValpha.Valphai / lik.item[["eta"]]^2
        diag( aux ) <- diag( aux ) + Dmat
        
        dlogdetQ <- -sum( 
          Qi * rbind( 
            cbind(   aux,                DX   ),
            cbind( t(DX), crossprod( XX, DX ) )
          )
        ) / lik.item[["param"]]["nugget"]^2
        
        ##              print( dlogdetQ/
        ##                  deriv.fwd.tf[[param.tf[t.i]]]( lik.item[["param"]][t.i] ) 
        ##              )
        ##              
        ##  derivate of U with respect to variance
        
        dU <- -0.5 * sum( 
          lik.item[["effects"]][["z.star"]] * drop( dValpha.Valphai %*% lik.item[["effects"]][["bhat"]] )
        ) * lik.item[["param"]]["variance"] / sum( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] )^2
        
        ##              print( dU /
        ##                  deriv.fwd.tf[[param.tf[t.i]]]( lik.item[["param"]][t.i] ) 
        ##              )
        
        r.gradient[t.i] <- ( -0.5 * ( 
            lik.item[["param"]]["variance"] / 
            sum( lik.item[["param"]]["variance"] + lik.item[["param"]]["snugget"] ) *
            sum( diag( dValpha.Valphai ) ) + dlogdetQ
          ) - dU
        ) / deriv.fwd.tf[[param.tf[t.i]]]( 
          c( lik.item[["param"]], lik.item[["aniso"]][["aniso"]] )[t.i] 
        )
      }
    }
    
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
    slv,
    envir,
    initial.objects,
    variogram.model, param, fit.param,
    aniso, fit.aniso,
    param.tf, 
    fwd.tf, 
    deriv.fwd.tf, 
    bwd.tf,
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
    irwls.initial,
    irwls.maxiter, 
    irwls.reltol, 
    force.gradient,
    zero.dist,
    nleqslv.method, nleqslv.control,
    ## bbsolve.method, bbsolve.control,
    optim.method, optim.lower, optim.upper, hessian, optim.control,
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
  betahat   <- coefficients( initial.objects[["initial.fit"]] )
  bhat      <- initial.objects[["bhat"]]
  coordinates <- initial.objects[["locations.objects"]][["coordinates"]]
  
  ##  check for multiple observations for same location and generate
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
  
  #   indices.pairs <- gtools::combinations( 
  #     1:nrow( coordinates ), 2, repeats.allowed = TRUE 
  #   )
  
  indices.pairs <- combn( NROW( coordinates ), 2 )
  lag.vectors <- coordinates[ indices.pairs[2,], ] - coordinates[ indices.pairs[1,], ]
  
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
    param.name[ !param.name %in% names( param ) ], "'"
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
  
  param.bounds <- param.bounds( variogram.model, NCOL( coordinates ), param )
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
    variogram.model %in% (t.models <- c( "fractalB" ) ) && 
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
  
  ##  initialize values of variogram parameters stored in the environment
  
  lik.item <- get( "lik.item", pos = as.environment( envir ) )

  lik.item[["param"]]   <-  rep( -1., length( param.name ) )
  lik.item[["eta"]]     <- NA
  lik.item[["aniso"]]   <- list( 
    isotropic = initial.objects[["isotropic"]], 
    aniso = rep( -1., length( aniso.name ) )
  )
  lik.item[["Valpha"]]  <- list()
  lik.item[["effects"]] <- list()
  lik.item[["eeq"]]     <- list()
  names( lik.item[["param"]] ) <- param.name
  names( lik.item[["aniso"]][["aniso"]] ) <- aniso.name
  
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
      
      if( slv ){
        
        ##         if( identical( root.finding, "nleqslv" ) ){
      
          r.root <- nleqslv(
            x = c( 
              transformed.param[ fit.param ], 
              transformed.aniso[ fit.aniso ] 
            ),
            fn = compute.estimating.equations,
            method = nleqslv.method,
            control = nleqslv.control,
            slv = slv,
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
            irwls.initial = irwls.initial,
            irwls.maxiter = irwls.maxiter,
            irwls.reltol = irwls.reltol,
            force.gradient = force.gradient,
            expectations = expectations,
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
          
      } else {
        
        ## minimize sum of squared estimating equations
        
        r.opt.eeq.sq <- optim(
          par = c( 
            transformed.param[ fit.param ], 
            transformed.aniso[ fit.aniso ] 
          ),
          fn = compute.estimating.equations,
          #         gr = gradient.negative.restricted.loglikelihood,
          method = optim.method, 
          lower = optim.lower,
          upper = optim.upper,
          control = optim.control,
          hessian = FALSE,
          slv = slv,
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
          irwls.initial = irwls.initial,
          irwls.maxiter = irwls.maxiter,
          irwls.reltol = irwls.reltol,
          force.gradient = force.gradient,
          expectations = expectations,
          verbose = verbose
        )
        
        r.converged <- r.opt.eeq.sq[["convergence"]] == 0
        r.convergence.code <- r.opt.eeq.sq[["convergence"]]      
        r.counts <- r.opt.eeq.sq[["counts"]]
        
        if( verbose > 0 ){
          cat( 
            "\n  sum(EEQ^2)         :",
            format( 
              signif( r.opt.eeq.sq[["value"]], digits = 7 ), 
              scientific = TRUE, width = 14
            ), sep = "" 
          )
          cat( 
            "\n  convergence code   :", 
            format( 
              signif( r.opt.eeq.sq[["convergence"]], digits = 0 ), 
              scientific = FALSE, width = 14
            ), "\n\n", sep = "" 
          )
        }
               
        #         if( hessian ) r.hessian <- r.opt.eeq.sq[["hessian"]]
        
        r.gradient <- compute.estimating.equations(
          adjustable.param = r.opt.eeq.sq[["par"]],
          slv = TRUE,
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
          irwls.initial = irwls.initial,
          irwls.maxiter = irwls.maxiter,
          irwls.reltol = irwls.reltol,
          force.gradient = force.gradient,
          expectations = expectations,
          verbose = verbose
        )

      }
      
    } else {
      
      ##  all variogram parameters are fixed
      
      ##  evaluate estimating equations
      
      r.gradient <- compute.estimating.equations(
        adjustable.param = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        slv = TRUE,
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
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        force.gradient = force.gradient,
        expectations = expectations,
        verbose = verbose
      )
      
      r.converged <- NA
      r.convergence.code <- NA_integer_
      r.counts <- c( nfcnt = NA_integer_, njcnt = NA_integer_ )
      
    }
    
    r.opt.neg.loglik <- NA_real_
    
  } else {
    
    if( any( fit.param ) ){
      
      ##  Gaussian REML estimation
      
      r.opt.neg.restricted.loglik <- optim(
        par = c( 
          transformed.param[ fit.param ], 
          transformed.aniso[ fit.aniso ] 
        ),
        fn = negative.restr.loglikelihood,
        gr = gradient.negative.restricted.loglikelihood,
        method = optim.method, 
        lower = optim.lower,
        upper = optim.upper,
        control = optim.control,
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
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        force.gradient = force.gradient,
        verbose = verbose
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
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        force.gradient = force.gradient,
        verbose = verbose
      )
      
    } else {
      
      ##  all variogram parameters are fixed
      
      ##  compute negative restricted loglikelihood and the gradient
      
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
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
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
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.reltol = irwls.reltol,
        force.gradient = force.gradient,
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
      extended.output = full.output,
      verbose = verbose
    )[-1]
    
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
    result.list[["Valpha.objects"]]    <- compress( 
      lik.item[["Valpha"]][c( "Valpha.inverse", "Valpha", "Valpha.ucf", "Valpha.ilcf", "gcr.constant" )]
    )
    result.list[["Aalpha"]] <- lik.item[["effects"]][["Aalpha"]]
    result.list[["Palpha"]] <- lik.item[["effects"]][["Palpha"]]

    result.list[["locations.objects"]] <- initial.objects[["locations.objects"]]
    
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


