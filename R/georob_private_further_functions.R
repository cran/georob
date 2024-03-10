####################################
#                                  #
#   Hilfsfunktionen fuer georob    #
#                                  #
####################################

##  ##############################################################################
### covariances.fixed.random.effects

covariances.fixed.random.effects <-
  function(
    Valphaxi.objects,
    Aalphaxi, Palphaxi, Valphaxi.inverse.Palphaxi,
    rweights, XX, TT, TtT, names.yy,
    nugget, eta,
    expectations, family = c( "gaussian", "long.tailed" ),
    cov.bhat, full.cov.bhat,
    cov.betahat,
    cov.bhat.betahat,
    cov.delta.bhat, full.cov.delta.bhat,
    cov.delta.bhat.betahat,
    cov.ehat, full.cov.ehat,
    cov.ehat.p.bhat, full.cov.ehat.p.bhat,
    aux.cov.pred.target,
    control.pcmp,
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
  ## 2013-02-05 AP covariance matrix of zhat
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-05-06 AP changes for solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-19 AP correcting error when computing covariance of regression residuals
  ## 2015-03-03 AP correcting error in computing result.new[["cov.bhat"]] for full.cov.bhat == FALSE
  ## 2015-06-30 AP changes for improving efficiency
  ## 2015-07-17 AP TtT passed to function, new name for function
  ## 2015-08-19 AP changes for computing covariances under long-tailed distribution of epsilon
  ## 2016-07-20 AP changes for parallel computations
  ## 2020-02-07 AP sanity checks for switch() and if()
  ## 2023-12-20 AP replacement of identical(class(...), ...) by inherits(..., ...)

#### -- preparation

  family = match.arg( family )

  ## store flags for controlling output

  ccov.bhat               <- cov.bhat
  ffull.cov.bhat          <- full.cov.bhat
  ccov.betahat            <- cov.betahat
  ccov.bhat.betahat       <- cov.bhat.betahat
  ccov.delta.bhat         <- cov.delta.bhat
  ffull.cov.delta.bhat    <- full.cov.delta.bhat
  ccov.delta.bhat.betahat <- cov.delta.bhat.betahat

  ## setting flags for computing required items

  cov.ystar        <- FALSE
  cov.bhat.b       <- FALSE
  cov.bhat.e       <- FALSE
  cov.betahat.b    <- FALSE
  cov.betahat.e    <- FALSE

  if( aux.cov.pred.target ){
    cov.bhat.b <- cov.betahat.b <- TRUE
  }

  if( cov.ehat.p.bhat ){
    cov.betahat <- cov.betahat.b <- cov.betahat.e <- TRUE
  }

  if( cov.ehat ){
    cov.delta.bhat.betahat <- cov.delta.bhat <- cov.betahat <- cov.bhat.e <- cov.betahat.e <- TRUE
    if( full.cov.ehat ) full.cov.delta.bhat <- TRUE
  }

  if( cov.delta.bhat.betahat ){
    cov.betahat.b <- cov.bhat.betahat <- TRUE
  }

  if( cov.delta.bhat ){
    cov.bhat <- cov.bhat.b <- TRUE
    if( full.cov.delta.bhat ) full.cov.bhat <- TRUE
  }

  if( any( cov.bhat, cov.betahat, cov.bhat.betahat, cov.delta.bhat, cov.delta.bhat.betahat, cov.ehat ) ){
    cov.ystar <- TRUE
  }

  ## compute required auxiliary items

  switch(
    family,
    gaussian = {
      var.psi     <- expectations["var.gauss.psi"]
      exp.dpsi    <- expectations["exp.gauss.dpsi"]
      var.eps     <- nugget
      cov.psi.eps <- expectations["exp.gauss.dpsi"]
    },
    long.tailed = {
      var.psi     <- expectations["var.f0.psi"]
      exp.dpsi    <- expectations["var.f0.psi"]
      var.eps     <- nugget * expectations["var.f0.eps"]
      cov.psi.eps <- 1.
    }
  )

  V <- eta * nugget * Valphaxi.objects[["Valphaxi"]]
  VTtT <- t( TtT * V )

  result.new <- list( error = FALSE )

  #   ## auxiliary items for checking computation for Gaussian case
  #
  #   cov.bhat = TRUE
  #   full.cov.bhat = TRUE
  #   cov.betahat = TRUE
  #   cov.bhat.betahat = TRUE
  #   cov.delta.bhat = TRUE
  #   full.cov.delta.bhat = TRUE
  #   cov.delta.bhat.betahat = TRUE
  #   cov.ehat = TRUE
  #   full.cov.ehat = TRUE
  #   cov.ehat.p.bhat = TRUE
  #   full.cov.ehat.p.bhat = TRUE
  #   aux.cov.pred.target = TRUE
  #
  #   Gammat <- VTtT
  #   diag( Gammat ) <- diag( Gammat ) + nugget
  #   Gammati <- solve( Gammat )
  #   VGammati <- VTtT %*% Gammati
  #   tmp1 <- solve( crossprod( XX, Gammati ) %*% XX )
  #   tmp2 <- XX[TT,] %*% tmp1 %*% t(XX[TT,])

#### -- compute S_alphaxi

  aux <- Valphaxi.inverse.Palphaxi / ( exp.dpsi * eta )
  diag( aux ) <- diag( aux ) + TtT
  aux <- try( chol( aux ), silent = TRUE )
  if( inherits( aux, "try-error" ) ){
    result.new[["error"]] <- TRUE
    return( result.new )
  }
  Salphaxi <- chol2inv( aux )

#### -- factors to compute bhat and betahat from zhat

  if( any( c( cov.ystar, cov.bhat.b, cov.bhat.e, cov.bhat, cov.bhat.betahat ) ) ){
    PaSa <- pmm( Palphaxi, Salphaxi, control.pcmp )
  }

  if( any( c( cov.betahat.b, cov.betahat.e, cov.betahat, cov.bhat.betahat) ) ){
    AaSa <- Aalphaxi %*% Salphaxi
  }

#### -- covariance of huberized observations

  if( cov.ystar ){
    cov.ystar <- TtT * VTtT
    diag( cov.ystar ) <- diag( cov.ystar ) + (var.psi * nugget / exp.dpsi^2) * TtT
    PaSa.cov.ystar <- pmm( PaSa, cov.ystar, control.pcmp )
  }

#### -- covariance of bhat and betahat with B and epsilon

  if( cov.bhat.b )    cov.bhat.b      <- pmm( PaSa, t( VTtT ), control.pcmp )
  if( cov.bhat.e )    cov.bhat.e      <- (nugget / exp.dpsi * cov.psi.eps * PaSa)[, TT]
  if( cov.betahat.b ){
    cov.betahat.b <- tcrossprod( AaSa, VTtT )
    #     cov.betahat.b <- AaSa %*% t( VTtT )
    TX.cov.betahat.bT <- (XX %*% cov.betahat.b)[TT,TT]
  }
  if( cov.betahat.e ){
    cov.betahat.e <- (nugget / exp.dpsi * cov.psi.eps * AaSa)[, TT]
    TX.cov.betahat.e <- (XX %*% cov.betahat.e)[TT,]
  }

  ## compute now the requested covariances ...

#### -- covariance of bhat (debugging status ok)

  if( cov.bhat ){
    t.cov.bhat <- if( full.cov.bhat )
    {
      aux <- pmm( PaSa.cov.ystar, t(PaSa ), control.pcmp )
      attr( aux, "struc" ) <- "sym"
      aux
    } else {
      aux <- rowSums( PaSa.cov.ystar * PaSa )
      names( aux ) <- rownames( XX )
      aux
    }
    if( ccov.bhat ) result.new[["cov.bhat"]] <- if( ffull.cov.bhat ){
      t.cov.bhat
    } else {
      f.diag( t.cov.bhat )
    }
  }

  #   print(
  #     summary(
  #       c(
  #         t.cov.betahat -
  #         (VGammati %*% VTtT - VGammati %*% tmp2 %*% t( VGammati ))
  #       )
  #     )
  #   )

#### -- covariance of betahat (debugging status ok)

  if( cov.betahat ){
    t.cov.betahat <- tcrossprod( tcrossprod( AaSa, cov.ystar ), AaSa )
    attr( t.cov.betahat, "struc" ) <- "sym"
    if( ccov.betahat ) result.new[["cov.betahat"]] <- t.cov.betahat
  }

  #   print( summary( c(
  #         t.cov.betahat -
  #         tmp1
  #       )
  #     )
  #   )

  ##  ... of bhat and betahat (debugging status ok)

  if( cov.bhat.betahat ){
    t.cov.bhat.betahat <- tcrossprod( PaSa.cov.ystar, AaSa )
    if( ccov.bhat.betahat ) result.new[["cov.bhat.betahat"]] <- t.cov.bhat.betahat
  }

  #   print( summary( c( t.cov.bhat.betahat ) ) )

#### -- covariance of (b - bhat) (debugging status ok)

  if( cov.delta.bhat ){
    t.cov.delta.bhat <- if( full.cov.delta.bhat ){
      aux <- V + t.cov.bhat - cov.bhat.b - t( cov.bhat.b )
      attr( aux, "struc" ) <- "sym"
      dimnames( aux ) <- list( rownames( XX ), rownames( XX ) )
      aux
    } else {
      #       aux <- diag( V ) - 2 * diag( cov.bhat.b ) + (if( full.cov.bhat ){
      #         diag( t.cov.betahat )
      #       } else {
      #         t.cov.betahat
      #       })
      aux <- diag( V ) - 2. * diag( cov.bhat.b ) + f.diag( t.cov.bhat )
      names( aux ) <- rownames( XX )
      aux
    }
    if( ccov.delta.bhat ) result.new[["cov.delta.bhat"]] <- if( ffull.cov.delta.bhat ){
      t.cov.delta.bhat
    } else {
      f.diag( t.cov.delta.bhat )
    }
  }

  #   print(
  #     summary(
  #       c(
  #         t.cov.delta.bhat -
  #         (VTtT - VGammati %*% VTtT + VGammati %*% tmp2 %*% t( VGammati ))
  #       )
  #     )
  #   )

#### -- covariance of (b - bhat) and betahat (debugging status ok)

  if( cov.delta.bhat.betahat ){
    t.cov.delta.bhat.betahat <- t( cov.betahat.b ) - t.cov.bhat.betahat
    dimnames( t.cov.delta.bhat.betahat ) <- dimnames( XX )
    if( ccov.delta.bhat.betahat ){
      result.new[["cov.delta.bhat.betahat"]] <- t.cov.delta.bhat.betahat
    }
  }

  #   print(
  #     summary(
  #       c(
  #         t.cov.delta.bhat.betahat -
  #         (VGammati %*% XX %*% tmp1)
  #       )
  #     )
  #   )

  ## ... of ehat (debugging status ok)

  if( cov.ehat ){
    aux1 <- tcrossprod( t.cov.delta.bhat.betahat, XX )[TT,TT]
    result.new[["cov.ehat"]] <- if( full.cov.ehat )
    {
      aux <- t.cov.delta.bhat[TT,TT] +
        tcrossprod( tcrossprod( XX, t.cov.betahat ), XX )[TT,TT] -
        aux1 - t(aux1) - cov.bhat.e[TT,] - t(cov.bhat.e)[,TT] -
        TX.cov.betahat.e - t(TX.cov.betahat.e)
      diag( aux ) <- diag( aux ) + var.eps
      attr( aux, "struc" ) <- "sym"
      dimnames( aux ) <- list( names.yy, names.yy )
      aux
    } else {
      #       aux <- (if( full.cov.delta.bhat ){
      #         diag( t.cov.delta.bhat )[TT]
      #       } else {
      #         t.cov.delta.bhat[TT]
      #       }) + rowSums( XX * tcrossprod( XX, t.cov.betahat) )[TT] -
      #         2. * diag( aux1 ) - 2. * diag( cov.bhat.e[TT,] ) - 2. * diag( TX.cov.betahat.e ) +
      #         nugget
      aux <- f.diag( t.cov.delta.bhat )[TT] +
        rowSums( XX * tcrossprod( XX, t.cov.betahat) )[TT] -
        2. * diag( aux1 ) - 2. * diag( cov.bhat.e[TT,] ) - 2. * diag( TX.cov.betahat.e ) +
        var.eps
        #         t.cov.delta.bhat[TT]
        #       }) + rowSums( XX * (XX %*% t.cov.betahat) )[TT] -
        #         2. * diag( aux1 ) - 2. * diag( cov.bhat.e[TT,] ) - 2. * diag( TX.cov.betahat.e ) +
        #         nugget
        names( aux ) <- names.yy
      aux
    }
  }

  #   tmp3 <- -VGammati
  #   diag(tmp3) <- diag(tmp3) + 1.
  #   print(
  #     summary(
  #       c(
  #         result.new[["cov.ehat"]] -
  #         (tmp3 %*% (Gammat - tmp2 ) %*% tmp3)
  #       )
  #     )
  #   )


#### -- covariance of ehat + bhat (debugging status ok)

  if( cov.ehat.p.bhat ){
    result.new[["cov.ehat.p.bhat"]] <- if( full.cov.ehat.p.bhat )
    {
      aux <- tcrossprod( tcrossprod( XX, t.cov.betahat ), XX )[TT,TT] -
        TX.cov.betahat.bT - t(TX.cov.betahat.bT) -
        TX.cov.betahat.e - t(TX.cov.betahat.e) + V[TT,TT]
      diag( aux ) <- diag( aux ) + var.eps
      attr( aux, "struc" ) <- "sym"
      dimnames( aux ) <- list( names.yy, names.yy )
      aux
    } else {
      aux <- rowSums( XX * tcrossprod( XX, t.cov.betahat) )[TT] -
        2. * diag( TX.cov.betahat.bT ) -
        2. * diag( TX.cov.betahat.e ) + diag( V )[TT] + var.eps
#       aux <- rowSums( XX * (XX %*% t.cov.betahat) )[TT] -
#         2. * diag( TX.cov.betahat.bT ) -
#         2. * diag( TX.cov.betahat.e ) + diag( V )[TT] + nugget
      names( aux ) <- names.yy
      aux
    }
  }

  #   print(
  #     summary(
  #       c(
  #         result.new[["cov.ehat.p.bhat"]] -
  #         (Gammat - tmp2)
  #       )
  #     )
  #   )

#### -- auxiliary item to compute covariance of kriging predictions
  ## and observations

  if( aux.cov.pred.target ){
    result.new[["cov.pred.target"]] <- pmm(
      rbind( cov.bhat.b, cov.betahat.b ),
      Valphaxi.objects[["Valphaxi.inverse"]] / eta / nugget,
      control.pcmp
    )
  }

  return( result.new )

}

##   ##############################################################################

update_zhat <-
  function(
    XX, yy, res, TT,
    nugget, eta, reparam,
    Valphaxi.inverse.Palphaxi,
    psi.function, tuning.psi,
    verbose
  )
{

  ## 2013-02-04 AP solving estimating equations for xi
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-06-29 AP solving linear system of equation by cholesky decomposition
  ## 2015-12-02 AP correcting error in try(chol(m))
  ## 2023-12-20 AP replacement of identical(class(...), ...) by inherits(..., ...)
  ## 2024-01-25 AP function update.zhat renamed to update_zhat

  ## function computes (1) updated IRWLS estimates zhat of linearized
  ## normal equations, (2) the associated rweights,
  ## (3) the unstandardized residuals (= estimated epsilons); the results
  ## are returned as a list

  ## compute rweights (cf. p. 7, proposal HRK of 26 July 2010)

  ## 2013-04-23 AP new names for robustness weights
  ## 2015-07-17 AP Gaussian (RE)ML estimation for reparametrized variogram

  if( reparam ){

    Wi <- rep( 1., length( res ) )

  } else {

    std.res <- res / sqrt( nugget )

    ##  construct left-hand side matrix M and right-hand side vector of
    ##  linear equation system

    Wi <- ifelse(
      abs( std.res ) < sqrt( .Machine[["double.eps"]] ),
      1.,
      psi.function( std.res, tuning.psi ) / std.res
    )

  }

  ##  aggregate rweights for replicated observations

  if( sum( duplicated( TT ) ) > 0. ){

    TtWiT  <- as.vector( tapply( Wi, factor( TT ), sum ) )
    TtWiyy <- as.vector( tapply( Wi * yy, factor( TT ), sum ) )

  } else {

    TtWiT <- Wi
    TtWiyy <- Wi * yy

  }

  ##  construct left-hand side matrix M and right-hand side vector b of
  ##  linearized system of equations

  M <- Valphaxi.inverse.Palphaxi / eta
  diag( M ) <- diag( M ) + TtWiT

  b <- TtWiyy

  ##  solve linear system

  result <- list( error = TRUE )

#   r.solve <- try( solve( M, b ), silent = TRUE )
#
#   if( !identical( class( r.solve ), "try-error" ) ) {

  t.chol <- try( chol( M ), silent = TRUE )
  if( !inherits( t.chol, "try-error" ) ){

    r.solve <- forwardsolve( t( t.chol ), b )
    r.solve <- backsolve( t.chol, r.solve )


    ##  collect output

    result[["error"]]      <- FALSE
    result[["zhat"]]       <- r.solve
    result[["residuals"]]  <- yy - result[["zhat"]][TT]
    result[["rweights"]]   <- Wi

  }

  return( result )

}


##    ##############################################################################

estimating.equations.z <- function(
  res, TT, zhat,
  nugget, eta, reparam,
  Valphaxi.inverse.Palphaxi,
  psi.function, tuning.psi
){

  ## auxiliary function to compute estimating equations for zhat

  ## 2015-07-17 AP Gaussian (RE)ML estimation for reparametrized variogram

  if( reparam ){

    Ttpsi <- res
    if( sum( duplicated( TT ) > 0. ) ){
      Ttpsi <- as.vector( tapply( Ttpsi, factor( TT ), sum ) )
    }

  } else {

    Ttpsi <- sqrt( nugget ) * psi.function( res / sqrt( nugget ), tuning.psi )
    if( sum( duplicated( TT ) > 0. ) ){
      Ttpsi <- as.vector( tapply( Ttpsi, factor( TT ), sum ) )
    }

  }

  Ttpsi - drop( Valphaxi.inverse.Palphaxi %*% zhat ) / eta

}

##    ##############################################################################
### estimate.zhat

estimate.zhat <-
  function(
    compute.zhat,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat,
    psi.function, tuning.psi, tuning.psi.nr,
    maxit, ftol,
    nugget, eta, reparam,
    Valphaxi.inverse,
    control.pcmp,
    verbose
  )
{

  ## 2013-02-04 AP solving estimating equations for xi
  ## 2013-06-03 AP handling design matrices with rank < ncol(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
  ## 2015-06-30 AP new method to determine convergence
  ## 2015-07-17 AP computing residual sums of squares (UStar)
  ## 2015-07-17 AP Gaussian (RE)ML estimation for reparametrized variogram
  ## 2016-07-20 AP changes for parallel computations
  ## 2023-12-20 AP replacement of identical(class(...), ...) by inherits(..., ...)
  ## 2024-01-25 AP function update.zhat renamed to update_zhat

  ## function computes (1) estimates zhat, bhat, betahat by
  ## solving robustified estimating equations by IRWLS,
  ## (2) the weights of the IRWLS, (3) the unstandardized residuals
  ## (= estimated epsilons); the results are returned as a list

#### -- compute projection matrix Palphaxi and related items

  result <- list( error = FALSE )

  aux1 <- crossprod( XX, Valphaxi.inverse )

  if( col.rank.XX[["deficient"]] ){
    s <- svd( aux1 %*% XX )
    s[["d"]] <- ifelse( s[["d"]] / max( s[["d"]] ) <= min.condnum, 0., 1. / s[["d"]] )
    aux2 <- s[["v"]] %*% ( s[["d"]] * t( s[["u"]] ) )                # Moore-Penrose inverse
  } else {
    t.chol <- try( chol( aux1 %*% XX ), silent = TRUE )
    if( !inherits( t.chol, "try-error" ) ){
      aux2 <- chol2inv( t.chol )
    } else {
      result[["error"]] <- TRUE
      return( result )
    }
  }

  result[["Aalphaxi"]]             <- aux2 %*% aux1
  dimnames( result[["Aalphaxi"]] ) <- dimnames( t(XX) )

  result[["Palphaxi"]]             <- -XX %*% result[["Aalphaxi"]]
  diag( result[["Palphaxi"]] )     <- diag( result[["Palphaxi"]] ) + 1.
  rownames( result[["Palphaxi"]] ) <- rownames( XX )
  colnames( result[["Palphaxi"]] ) <- rownames( XX )

  result[["Valphaxi.inverse.Palphaxi"]] <- pmm(
    Valphaxi.inverse, result[["Palphaxi"]], control.pcmp
  )
  rownames( result[["Valphaxi.inverse.Palphaxi"]] )      <- rownames( XX )
  colnames( result[["Valphaxi.inverse.Palphaxi"]] )      <- rownames( XX )
  attr( result[["Valphaxi.inverse.Palphaxi"]], "struc" ) <- "sym"

#### -- estimate signal zhat by IRWLS

  if( compute.zhat ){

    res <- yy - zhat[TT]

    eeq.old <- estimating.equations.z(
      res, TT, zhat,
      nugget, eta, reparam,
      result[["Valphaxi.inverse.Palphaxi"]],
      psi.function, tuning.psi
    )
    eeq.old.l2 <- sum( eeq.old^2 ) / 2.

    if( !is.finite( eeq.old.l2 ) ) {
      result[["error"]] <- TRUE
      return( result )
    }

    converged <- FALSE

    if( verbose > 2. ) cat(
      "\n  IRWLS\n",
      "      it     Fnorm.old     Fnorm.new   largest |f|\n", sep = ""
    )

    ##  IRWLS

    for( i in 1L:maxit ){

#### --- compute new estimates

      new <- update_zhat(
        XX, yy, res, TT,
        nugget, eta, reparam,
        result[["Valphaxi.inverse.Palphaxi"]],
        psi.function, tuning.psi,
        verbose
      )

      if( new[["error"]] ) {
        result[["error"]] <- TRUE
        return( result )
      }

#### --- evaluate estimating equations for z and compute its l2 norm

      eeq.new <- estimating.equations.z(
        new[["residuals"]], TT, new[["zhat"]],
        nugget, eta, reparam,
        result[["Valphaxi.inverse.Palphaxi"]],
        psi.function, tuning.psi
      )
      eeq.new.l2 <- sum( eeq.new^2 ) / 2.

      if( !is.finite( eeq.new.l2 ) ) {
        result[["error"]] <- TRUE
        return( result )
      }

      if( verbose > 2. ) cat(
        format( i, width = 8L ),
        format(
          signif(
            c( eeq.old.l2, eeq.new.l2, max(abs(eeq.new) ) ), digits = 7L
          ), scientific = TRUE, width = 14L
        ), "\n", sep = ""
      )

#### --- check for convergence

      if( max( abs( eeq.new ) ) <= ftol && i >= 2L ){
        converged <- TRUE
        break
      }

#### --- update zhat, residuals and eeq.old.l2

      eeq.old.l2 <- eeq.new.l2
      zhat      <- new[["zhat"]]
      res        <- new[["residuals"]]

    }

#### -- compute scaled residuals sum of squares

    ##  collect output

    result[["zhat"]]            <- new[["zhat"]]
    names( result[["zhat"]] )   <- rownames( XX )

    result[["residuals"]]        <- new[["residuals"]]
    result[["rweights"]]         <- new[["rweights"]]
    result[["converged"]]        <- converged
    result[["nit"]]              <- i

  } else {

    result[["zhat"]]            <- zhat
    names( result[["zhat"]] )   <- rownames( XX )

    result[["residuals"]]        <- yy - zhat[TT]

    if( reparam ){
      result[["rweights"]]       <- rep( 1., length( result[["residuals"]] ) )
    } else {
      result[["rweights"]]       <- ifelse(
        abs( std.res <- result[["residuals"]] / sqrt( nugget ) ) < sqrt( .Machine[["double.eps"]] ),
        1.,
        psi.function( std.res, tuning.psi ) / std.res
      )
    }

    result[["converged"]]        <- NA
    result[["nit"]]              <- NA_integer_

  }

#### -- prepare output

  result[["bhat"]]             <- drop( result[["Palphaxi"]] %*% result[["zhat"]] )
  names( result[["bhat"]] )    <- rownames( XX )

  result[["betahat"]]          <- drop( result[["Aalphaxi"]] %*% result[["zhat"]] )
  names( result[["betahat"]] ) <- colnames( XX )

  result[["Valphaxi.inverse.bhat"]] <- drop( Valphaxi.inverse %*% result[["bhat"]] )

  result[["RSS"]] <- f.aux.RSS(
    res = result[["residuals"]],
    TT = TT, TtT = TtT,
    bhat = result[["bhat"]],
    Valphaxi.inverse.bhat = result[["Valphaxi.inverse.bhat"]],
    eta = eta
  )

  return( result )

}


##    ##############################################################################
### likelihood.calculations

likelihood.calculations <-
  function(
    envir,
    adjustable.param.aniso, fixed.param.aniso, name.param.aniso, tf.param.aniso,
    bwd.tf, safe.param, reparam,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat,
    psi.function, tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.ftol,
    compute.zhat = TRUE,
    control.pcmp,
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
  ## 2015-03-16 AP extended variogram parameter transformations, elimination of unused variables
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  ## 2015-07-17 AP TtT passed to function
  ## 2015-07-17 AP new name of function, Gaussian (RE)ML estimation for reparametrized variogram
  ## 2015-07-23 AP changes for avoiding computation of Valphaxi object if not needed
  ## 2015-12-02 AP reparametrized variogram parameters renamed
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-03 AP changes for nested variogram models
  ## 2016-11-14 AP correcting error in 3d rotation matrix for geometrically anisotropic variograms

  ##  function transforms (1) the variogram parameters back to their
  ##  original scale; computes (2) the correlation matrix, its inverse
  ##  and its inverse lower cholesky factor; (3) computes betahat,
  ##  bhat and further associates items; and (4) computes the
  ##  matrices A and the cholesky factor of the matrix Q

#### -- preparations

  d2r <- pi / 180.

  ## load lik.item object

  lik.item <- get( "lik.item", pos = as.environment( envir ) )

  #   print(str(lik.item[c("variogram.object", "eta", "xi")]))

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
    function( i, x, lik.item ){
      c(
        lik.item[["variogram.object"]][[i]][c("variogram.model", "isotropic")],
        x[[i]]
      )
    }, x = variogram.object, lik.item = lik.item
  )

  ##  check whether the current variogram parameters and the variogram
  ##  parameters that were used in the previous call to
  ##  likelihood.calculations are the same

  same.param <- isTRUE( all.equal(
      param.aniso,
      unlist( lapply(
          lik.item[["variogram.object"]],
          function(x) c( x[["param"]], x[["aniso"]] )
        ))
    ))

  if( same.param && !is.null( lik.item[["zhat"]] ) ){

    if( verbose > 4. ) cat(
      "\n     likelihood.calculations: exit without computing any objects\n"
    )
    return( lik.item )
  }

  ## check whether variogram parameters are within reasonable bounds and
  ## return an error otherwise

  if( length( param.aniso ) && any( param.aniso > safe.param ) ){

    if( verbose > 1 ) lapply(
      1L:length(variogram.object),
      function( i, x, d2r, reparam ){

        x <- x[[i]]

        t.param <- x[["param"]]

        if( reparam ){
          t.param <- t.param[!names(t.param) %in% "snugget"]
          tmp <- names( t.param )
          tmp <- gsub(
            "nugget", "eta", gsub(
              "variance", "xi", tmp, fixed = TRUE ), fixed = TRUE )
          names( t.param ) <- tmp
        }

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

      }, x = variogram.object, d2r = d2r, reparam = reparam
    )

    return( lik.item )

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
          return( lik.item )
        },
        param = ep.param,
        bounds = param.bounds
      )
    }, x = variogram.object, lag.vectors = lag.vectors
  )

  ##  update variogram and parameters

#   f.rotate <- function( omega=90, phi=90, zeta=0 ){
#     omega <- omega / 180 * pi
#     phi <- phi / 180 * pi
#     zeta <- zeta / 180 * pi
#
#     co = cos(omega)
#     so = sin(omega)
#     cp = cos(phi)
#     sp = sin(phi)
#     cz = cos(zeta)
#     sz = sin(zeta)
#
#     rbind(
#       c(             so*sp,             co*sp,       cp ),
#       c( -co*cz - so*cp*sz,  so*cz - co*cp*sz,    sp*sz ),
#       c(  co*sz - so*cp*cz, -so*sz - co*cp*cz,    sp*cz )
#     )
#   }

#### -- update objects in lik.item

  lik.item[["variogram.object"]] <- lapply(
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
    x = variogram.object, vo = lik.item[["variogram.object"]],
    n = attr( lag.vectors, "ndim.coords" )
  )

  ## update eta, xi, var.signal

  if( !reparam ){
    reparam.variogram.object <- f.reparam.fwd( lik.item[["variogram.object"]] )
    t.var.signal <- attr( reparam.variogram.object, "var.signal" )
  } else {
    reparam.variogram.object <- lik.item[["variogram.object"]]
    t.var.signal <- NA_real_
  }

  lik.item[["eta"]] <- unname( reparam.variogram.object[[1L]][["param"]]["nugget"] )
  lik.item[["xi"]] <-  unname( sapply(
      reparam.variogram.object, function(x) x[["param"]]["variance"]
    ))
  lik.item[["var.signal"]] <- t.var.signal

  ##  print updated variogram parameters

  if( verbose > 1. ) {

    lapply(
      1L:length(variogram.object),
      function( i, x, d2r, reparam ){

        x <- x[[i]]

        t.param <- x[["param"]]

        if( reparam ){
          t.param <- t.param[!names(t.param) %in% "snugget"]
          tmp <- names( t.param )
          tmp <- gsub(
            "nugget", "eta",
            gsub( "variance", "xi", tmp, fixed = TRUE ), fixed = TRUE
          )
          names( t.param ) <- tmp
        }

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

      }, x = variogram.object, d2r = d2r, reparam = reparam
    )

  }

  #   print(str(lik.item))

#### -- compute updates of required likelihood items if the variogram
  ##  parameters differ and are all within allowed bounds

  if( !same.param || is.null( lik.item[["Valphaxi"]] ) ){

    if( verbose > 4. ) cat(
      "\n     likelihood.calculations: computing 'Valphaxi' object\n"
    )

    lik.item[["Valphaxi"]][["error"]] <- TRUE

#### --- calculate generalized correlation matrix, its inverse and its
    ##  inverse cholesky factor

    t.Valphaxi <- f.aux.Valphaxi(
      lag.vectors = lag.vectors,
      variogram.object = lik.item[["variogram.object"]],
      xi = lik.item[["xi"]],
      control.pcmp = control.pcmp,
      verbose = verbose
    )

    if( !t.Valphaxi[["error"]] ){
      lik.item[["Valphaxi"]] <- t.Valphaxi
    } else {
      return( lik.item )
    }

  }

  #     print(str(lik.item))

#### --- estimate fixed and random effects (zhat, betahat, bhat, residuals )
  ##  and estimate of signal variance for Gaussian (RE)ML

  if( verbose > 4. ) cat(
    "\n     likelihood.calculations: computing 'zhat' object\n"
  )

  ##  either take initial guess of betahat and bhat for the current
  ##  irwls iteration from initial.object or from previous iteration

  if(
    !irwls.initial && !is.null( lik.item[["zhat"]][["zhat"]] )
  ){
    zhat <- lik.item[["zhat"]][["zhat"]]
  }

  lik.item[["zhat"]] <- estimate.zhat(
    compute.zhat,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat,
    psi.function, tuning.psi, tuning.psi.nr,
    irwls.maxiter, irwls.ftol,
    lik.item[["variogram.object"]][[1L]][["param"]]["nugget"], lik.item[["eta"]], reparam,
    lik.item[["Valphaxi"]][["Valphaxi.inverse"]],
    control.pcmp,
    verbose
  )

  if( lik.item[["zhat"]][["error"]] ) return( lik.item )     ##  an error occurred

#### --- compute Q matrix and its Cholesky factor (required for Gaussian
  ##  (RE)ML estimation)

  if( tuning.psi >= tuning.psi.nr ) {

    ## compute matrix Q

    if( verbose > 4. ) cat(
      "\n     likelihood.calculations: computing 'Q' object\n"
    )

    t.Q <- f.aux.Qstar(
      TT = TT, TtT = TtT,
      XX = XX, col.rank.XX = col.rank.XX, min.condnum = min.condnum,
      Vi = lik.item[["Valphaxi"]][["Valphaxi.inverse"]],
      eta = lik.item[["eta"]],
      ml.method = ml.method, control.pcmp = control.pcmp
    )

    if( !t.Q[["error"]] ){
      lik.item[["Q"]] <- t.Q
    } else {
      return( lik.item )
    }

  }

#### -- store updated lik.item object

  assign( "lik.item", lik.item, pos = as.environment( envir ) )

#   print( str( lik.item ) ); stop()

  return( lik.item )

}


##   ##############################################################################
### estimating.equations.theta

estimating.equations.theta <-
  function(
    adjustable.param.aniso,
    envir,
    fixed.param.aniso, name.param.aniso, tf.param.aniso, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat,
    psi.function,
    tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.ftol,
    force.gradient,
    expectations,
    error.family.estimation,
    control.pcmp,
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
  ## 2015-03-16 AP elimination of unused variables
  ## 2015-07-17 AP TtT passed to function
  ## 2015-07-17 AP new function interface, improved efficiency
  ## 2015-07-27 AP changes to further improve efficiency
  ## 2015-07-29 AP changes for elimination of parallelized computation of gradient or estimating equations
  ## 2015-08-19 AP control about error families for computing covariances added
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-04 AP changes for nested variogram models

#### -- get new lik.item

  lik.item <- likelihood.calculations(
    envir,
    adjustable.param.aniso, fixed.param.aniso, name.param.aniso, tf.param.aniso,
    bwd.tf, safe.param, reparam = FALSE,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat,
    psi.function, tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.ftol,
    compute.zhat = TRUE,
    control.pcmp = control.pcmp,
    verbose = verbose
  )

  #   print( str( lik.item ) ); stop()

  ##  check whether generalized covariance matrix is positive definite

  if( lik.item[["Valphaxi"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0. ) cat(
      "\n(generalized) correlation matrix Valphaxi is not positive definite\n"
    )
    t.result <- rep( Inf, length( adjustable.param.aniso ) )
    names( t.result ) <- names( adjustable.param.aniso )
    return( t.result )
  }

  ##  check whether computation of betahat and bhat failed

  if( lik.item[["zhat"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0. ) cat(
      "\nan error occurred when estimating the fixed and random effects\n"
    )
    t.result <- rep( Inf, length( adjustable.param.aniso ) )
    names( t.result ) <- names( adjustable.param.aniso )
    return( t.result )
  }

  ##  check whether estimating equations should be computed for fixed parameters

  if( length( adjustable.param.aniso ) == 0L && force.gradient ){
    adjustable.param.aniso <- fixed.param.aniso
  }

#### -- evaluate estimating equations

  if( length( adjustable.param.aniso ) > 0L ){

#### --- compute Cov[bhat]

    r.cov <- covariances.fixed.random.effects(
      Valphaxi.objects = lik.item[["Valphaxi"]][c("Valphaxi", "Valphaxi.inverse")],
      Aalphaxi = lik.item[["zhat"]][["Aalphaxi"]],
      Palphaxi = lik.item[["zhat"]][["Palphaxi"]],
      Valphaxi.inverse.Palphaxi = lik.item[["zhat"]][["Valphaxi.inverse.Palphaxi"]],
      rweights = lik.item[["zhat"]][["rweights"]],
      XX = XX, TT = TT, TtT = TtT, names.yy = names( yy ),
      nugget = lik.item[["variogram.object"]][[1L]][["param"]]["nugget"],
      eta = lik.item[["eta"]],
      expectations = expectations, family = error.family.estimation,
      cov.bhat = TRUE, full.cov.bhat = TRUE,
      cov.betahat = FALSE,
      cov.bhat.betahat = FALSE,
      cov.delta.bhat = FALSE, full.cov.delta.bhat = FALSE,
      cov.delta.bhat.betahat = FALSE,
      cov.ehat = FALSE, full.cov.ehat = FALSE,
      cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
      aux.cov.pred.target = FALSE,
      control.pcmp = control.pcmp,
      verbose = verbose
    )

    if( r.cov[["error"]] ) {
      warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
      if( verbose > 0. ) cat(
        "\nan error occurred when computing the covariances of fixed and random effects\n"
      )
      t.result <- rep( Inf, length( adjustable.param.aniso ) )
      names( t.result ) <- names( adjustable.param.aniso )
      return( t.result )
    }

    ## compute further required items

    ## Valphaxi^-1 Cov[bhat]

    Valphaxii.cov.bhat <- pmm(
      lik.item[["Valphaxi"]][["Valphaxi.inverse"]], r.cov[["cov.bhat"]],
      control.pcmp
    )

#### --- computation of estimating equations for all elements of adjustable.param.aniso

    t.eeq <- sapply(
      names( adjustable.param.aniso ),
      f.aux.eeq,
      variogram.object = lik.item[["variogram.object"]],
      xi = lik.item[["xi"]],
      Valphaxii = lik.item[["Valphaxi"]][["Valphaxi.inverse"]],
      Valphaxii.cov.bhat = Valphaxii.cov.bhat,
      Valpha = lik.item[["Valphaxi"]][["Valpha"]],
      bh = lik.item[["zhat"]][["bhat"]],
      bhVaxi = lik.item[["zhat"]][["Valphaxi.inverse.bhat"]],
      r.cov = r.cov, lik.item = lik.item,
      TtT = TtT,
      lag.vectors = lag.vectors,
      control.pcmp = control.pcmp, verbose = verbose
    )

    eeq.exp <- t.eeq["eeq.exp", ]
    eeq.emp <- t.eeq["eeq.emp", ]

    names( eeq.exp ) <- names( adjustable.param.aniso )
    names( eeq.emp ) <- names( adjustable.param.aniso )

    if( verbose > 1. ) {

      t.eeq.exp <- eeq.exp
      t.eeq.emp <- eeq.emp

      tmp <- strsplit( names(t.eeq.exp), control.georob()[["sepstr"]], fixed = TRUE )
      names( t.eeq.exp ) <- names( t.eeq.emp ) <- sapply( tmp, function(x) x[1L] )
      cmp <- sapply( tmp, function(x) as.integer(x[2L]) )

      t.eeq.exp <- tapply( t.eeq.exp, factor( cmp ), function( x ) x )
      t.eeq.emp <- tapply( t.eeq.emp, factor( cmp ), function( x ) x )

      lapply(
        1L:length(t.eeq.exp),
        function( i, eeq.exp, eeq.emp ){

          eeq.exp <- eeq.exp[[i]]
          eeq.emp <- eeq.emp[[i]]

          cat( "\n                      ",
            format( names( eeq.emp), width = 14L, justify = "right" ),
            "\n", sep =""
          )
          cat( "  EEQ                :",
            format(
              signif( eeq.emp / eeq.exp - 1., digits = 7L ),
              scientific = TRUE, width = 14L
            ), "\n", sep = ""
          )
          if( verbose > 2. ){
            cat( "      empirical terms:",
              format(
                signif( eeq.emp, digits = 7L ),
                scientific = TRUE, width = 14L
              ), "\n", sep = ""
            )
            cat( "      expected  terms:",
              format(
                signif( eeq.exp, digits = 7L ),
                scientific = TRUE, width = 14L
              ), "\n", sep = ""
            )
          }
          cat("\n")


        }, eeq.exp = t.eeq.exp, eeq.emp = t.eeq.emp
      )

    }

#### -- store objects in lik.item object

    lik.item[["eeq"]] <- list(
      eeq.emp = eeq.emp,
      eeq.exp = eeq.exp
    )

    assign( "lik.item", lik.item, pos = as.environment( envir ) )

    return( eeq.emp / eeq.exp - 1. )

  } else {

    ##  all parameters are fixed

    return( NA_real_ )

  }

}


##   ##############################################################################

negative.loglikelihood <-
  function(
    adjustable.param.aniso,
    envir,
    fixed.param.aniso, name.param.aniso, tf.param.aniso, bwd.tf, safe.param,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat,
    psi.function,
    tuning.psi, tuning.psi.nr, ml.method, reparam,
    irwls.initial, irwls.maxiter, irwls.ftol,
    control.pcmp,
    verbose,
    ...
  )
{

  ## function computes to negative (un)restricted loglikelihood

  ## 2012-04-21 AP scaled psi-function
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-11-04 AP unscaled psi-function
  ## 2012-11-27 AP changes in parameter back-transformation
  ## 2013-06-03 AP changes for estimating zhat
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2015-03-16 AP elimination of unused variables
  ## 2015-07-17 AP new name of function, Gaussian (RE)ML estimation for reparametrized variogram
  ## 2015-07-27 AP correcting error in likelihood for original parametrization (reparam == FALSE)
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-04 AP changes for nested variogram models

  #     sel <- !c( param.name, aniso.name ) %in% names( fixed.param )
  #     names( adjustable.param.aniso ) <- c( param.name, aniso.name )[sel]

  ##  compute required items (param, eta, Valphaxi.inverse, Valphaxi.ilcf,
  ##  betahat, bhat, residuals, etc.)

  lik.item <- likelihood.calculations(
    envir,
    adjustable.param.aniso, fixed.param.aniso, name.param.aniso, tf.param.aniso,
    bwd.tf, safe.param, reparam,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat,
    psi.function, tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.ftol,
    compute.zhat = TRUE,
    control.pcmp = control.pcmp,
    verbose
  )

  ##  check whether generalized covariance matrix is positive definite

 if( lik.item[["Valphaxi"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0. ) cat(
      "\n(generalized) correlation matrix Valphaxi is not positive definite\n"
    )
    return( NA )
  }

  ##  check whether computation of betahat and bhat failed

  if( lik.item[["zhat"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0. ) cat(
      "\nan error occurred when estimating the fixed and random effects\n"
    )
    return( NA )
  }

  ##  check whether Q matrix not positive definite

  if( lik.item[["Q"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0. ) cat(
      "\nan error occurred when determinants required for",
      "Gaussian log-likelihood were computed\n"
    )
    return( NA )
  }

  ##  compute negative (restricted || profile) loglikelihood

  t.q <- t.qmp <- NROW(XX)
  if( identical( ml.method, "REML" ) ){
    t.qmp <- t.q - col.rank.XX[["rank"]]
  }

  if( reparam ){

    ## (restricted) profile loglikelihood

    r.neg.loglik <- 0.5 * (
      t.qmp * (
        1. + log(2.*pi) - log(t.qmp ) +
        log( lik.item[["zhat"]][["RSS"]] )
      ) -
      t.q * log( lik.item[["eta"]] ) +
      lik.item[["Q"]][["log.det.Qstar"]] +
      lik.item[["Valphaxi"]][["log.det.Valphaxi"]]
    )

  } else {

    ## (restricted) loglikelihood

    r.neg.loglik <- 0.5 * (
      t.qmp * (
        log(2.*pi) + log( lik.item[["var.signal"]] )
      ) -
      t.q * log( lik.item[["eta"]] ) +
      lik.item[["Valphaxi"]][["log.det.Valphaxi"]] +
      lik.item[["Q"]][["log.det.Qstar"]] +
      lik.item[["zhat"]][["RSS"]] / lik.item[["var.signal"]]
    )

  }

  attributes( r.neg.loglik ) <- NULL

  if( verbose > 1. ) cat(
    "\n  Negative. restrict. loglikelihood:",
    format(
      signif( r.neg.loglik, digits = 7L ),
      scientific = TRUE, width = 14L
    ), "\n", sep = ""
  )

  return( r.neg.loglik )

}


##   ##############################################################################

gradient.negative.loglikelihood <-
  function(
    adjustable.param.aniso,
    envir,
    fixed.param.aniso, name.param.aniso, tf.param.aniso,
    deriv.fwd.tf, bwd.tf, safe.param, reparam,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat,
    psi.function,
    tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.ftol,
    force.gradient,
    control.pcmp,
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
  ## 2015-03-16 AP elimination of unused variables
  ## 2015-07-17 AP new name of function, Gaussian (RE)ML estimation for reparametrized variogram
  ## 2015-07-29 AP changes for elimination of parallelized computation of gradient or estimating equations
  ## 2015-12-02 AP reparametrized variogram parameters renamed
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-04 AP changes for nested variogram models

  ##  get lik.item

  lik.item <- likelihood.calculations(
    envir,
    adjustable.param.aniso, fixed.param.aniso, name.param.aniso, tf.param.aniso,
    bwd.tf, safe.param, reparam,
    lag.vectors,
    XX, min.condnum, col.rank.XX, yy, TT, TtT, zhat,
    psi.function, tuning.psi, tuning.psi.nr, ml.method,
    irwls.initial, irwls.maxiter, irwls.ftol,
    compute.zhat = TRUE,
    control.pcmp = control.pcmp,
    verbose
  )

  ##  check whether generalized covariance matrix is positive definite

  if( lik.item[["Valphaxi"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0. ) cat(
      "\n(generalized) correlation matrix Valphaxi is not positive definite\n"
    )
    return( rep( NA, length( adjustable.param.aniso ) ) )
  }

  ##  check whether computation of betahat and bhat failed

  if( lik.item[["zhat"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0. ) cat(
      "\nan error occurred when estimating the fixed and random effects\n"
    )
    return( rep( NA, length( adjustable.param.aniso ) ) )
  }

  ##  check whether Q matrix not positive definite

  if( lik.item[["Q"]][["error"]] ) {
    warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
    if( verbose > 0. ) cat(
      "\nan error occurred when determinants required for ",
      "Gaussian log-likelihood were computed\n"
    )
    return( rep( NA, length( adjustable.param.aniso ) ) )
  }

  ##  check whether gradient should be computed for fixed parameters

  if( length( adjustable.param.aniso ) == 0L && force.gradient ){
    adjustable.param.aniso <- fixed.param.aniso
  }

  ##  construct param.tf

  param.tf <- tf.param.aniso[names(adjustable.param.aniso)]
  tmp <- strsplit( names(param.tf), control.georob()[["sepstr"]], fixed = TRUE )
  names(param.tf) <- sapply( tmp, function(x) x[1L] )
  param.tf <- param.tf[!duplicated(names(param.tf))]

  ##  evaluate gradient

  if( length( adjustable.param.aniso ) > 0L ){

    ##  compute auxiliary items

    n <- NROW( XX )

    if( reparam ){

      Qsi       <- lik.item[["Q"]][["Qstar.inverse"]]
      Valphaxii <- lik.item[["Valphaxi"]][["Valphaxi.inverse"]]
      bh        <- lik.item[["zhat"]][["bhat"]]
      bhVaxi    <- lik.item[["zhat"]][["Valphaxi.inverse.bhat"]]
      if( !identical(
          gsub( ".__...__.1", "", names(adjustable.param.aniso), fixed = TRUE ),
          "nugget" )
      ){
        Qst11Vai  <- pmm( Qsi[1L:n, 1L:n], Valphaxii, control.pcmp )
      } else {
        Qst11Vai  <- NULL
      }

    } else {

      Qi     <- lik.item[["Q"]][["Qstar.inverse"]] * lik.item[["var.signal"]]
      Vi     <- lik.item[["Valphaxi"]][["Valphaxi.inverse"]] / lik.item[["var.signal"]]
      bh     <- lik.item[["zhat"]][["bhat"]]
      bhVi   <- lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] / lik.item[["var.signal"]]
      if( !identical(
          gsub( ".__...__.1", "", names(adjustable.param.aniso), fixed = TRUE ),
          "nugget" )
      ){
        Qt11Vi <- pmm( Qi[1L:n, 1L:n], Vi, control.pcmp )
      } else {
        Qt11Vi <- NULL
      }

    }

    ##  compute gradient for elements of adjustable.param.aniso

    if( reparam ){

      ## gradient for reparametrized variogram parameters

      r.gradient <- sapply(
        names( adjustable.param.aniso ),
        f.aux.gradient.npll,
        variogram.object = lik.item[["variogram.object"]],
        eta = lik.item[["eta"]], xi = lik.item[["xi"]],
        TT = TT, TtT = TtT, XX = XX, col.rank.XX = col.rank.XX,
        res = lik.item[["zhat"]][["residuals"]],
        Ustar = lik.item[["zhat"]][["RSS"]],
        Qsi = Qsi, Qst11Vai = Qst11Vai,
        Valphaxii = Valphaxii,
        Valpha = lik.item[["Valphaxi"]][["Valpha"]],
        bh = bh, bhVaxi = bhVaxi,
        lag.vectors = lag.vectors,
        param.tf = param.tf, deriv.fwd.tf = deriv.fwd.tf,
        ml.method = ml.method,
        control.pcmp = control.pcmp, verbose = verbose
      )

    } else {

      ## gradient for original variogram parameters

      r.gradient <- sapply(
        names( adjustable.param.aniso ),
        f.aux.gradient.nll,
        variogram.object = lik.item[["variogram.object"]],
        TT = TT, TtT = TtT, XX = XX,
        res = lik.item[["zhat"]][["residuals"]],
        Qi = Qi, Vi = Vi, Qt11Vi = Qt11Vi,
        Valpha = lik.item[["Valphaxi"]][["Valpha"]],
        bh = bh, bhVi = bhVi,
        lag.vectors = lag.vectors,
        param.tf = param.tf, deriv.fwd.tf = deriv.fwd.tf,
        ml.method = ml.method,
        control.pcmp = control.pcmp, verbose = verbose
      )

    }

    names( r.gradient ) <- names( adjustable.param.aniso )

    ##  rearrange elements of gradient and change sign (for negative
    ##  log-likelihood)

    r.gradient <- -r.gradient[names( adjustable.param.aniso )]

    if( verbose > 1. ) f.aux.print.gradient( r.gradient, reparam )

    return( r.gradient )

  } else {

    ##  all parameters are fixed

    return( NA_real_ )

  }

}


##  ##   ##############################################################################
##
##      f.compute.df <- function( Valphaxi, XX, param ){
##
##          ##  function computes three estimates of the degrees of freedom of
##          ##  the smoothing universal kriging predictor, cf.  Hastie &
##          ##  Tibshirani, 1990, Generalized additive models, pp.52
##
##          ##  2011-07-05
##          ##  Andreas Papritz
##
##          sigma <- param["variance"] * Valphaxi
##          diag( sigma ) <- diag( sigma ) + param["nugget"]
##
##          ##  compute inverse lower cholesky factor of covariance matrix of
##          ##  data
##
##          ilcf <- t( backsolve( chol( sigma ), diag( NROW( Valphaxi ) ), k = NROW( Valphaxi ) ) )
##
##          ##  compute hat matrix
##
##          q <- qr.Q( qr( xtilde <- ilcf %*% XX ) )
##          s <- -tcrossprod( q )
##
##          diag( s ) <- diag( s ) + 1.
##          s <- -param["nugget"] * t( ilcf ) %*% s %*% ilcf
##          diag( s ) <- diag( s ) + 1.
##
##          ##  compute degrees of freedom
##
##          df.1 <- sum( diag( s ) )
##          df.3 <- sum( s^2 )
##          df.2 <- 2. * df.1 - df.3
##
##          return(
##              c(
##                  df.SSt    = t.df.2 <- sum( s^2 ),
##                  df.S      = t.df.1 <- sum( diag( s ) ),
##                  df.2SmSSt = 2. * t.df.1 - t.df.2
##              )
##          )
##
##      }


# #  ##############################################################################
#
# getCall.georob <-
#   function( x, ... )
# {
#
#   ## Function replaces the name of a formula object in the call component
#   ## of a georob object by the formula itself (needed for update.default to
#   ## work correctly)
#
#   ## 2013-06-12 AP substituting [["x"]] for $x in all lists
#   ## 2024-01-26 AP renamed object to x
#
#   if( is.null( call <- getElement( x, "call" ) ) )  stop(
#     "need an object with call component"
#   )
#   call[["formula"]] <- update.formula( formula(x), formula( x ) )
#
#   return( call )
#
# }


################################################################################

f.aux.eeq <- function(
  x,
  variogram.object, xi,
  Valphaxii, Valphaxii.cov.bhat,
  Valpha,
  bh, bhVaxi,
  r.cov, lik.item,
  TtT,
  lag.vectors, control.pcmp, verbose
){

  ##  auxiliary function to compute robustified estimating equations
  ##  (called by estimating.equations.theta)

  ## 2014-07-29 A. Papritz
  ## 2015-03-03 AP changes to optimize computing effort
  ## 2015-07-17 AP new function interface and improved efficiency of computation
  ## 2015-07-27 AP changes to further improve efficiency
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-04 AP changes for nested variogram models
  ## 2020-02-07 AP sanity checks for switch() and if()

  ## correct parameters names and determine model component

  tmp <- unlist( strsplit( x, control.georob()[["sepstr"]], fixed = TRUE ) )
  x <- tmp[1L]
  cmp <- as.integer( tmp[2L] )

  ncmp <- length(variogram.object)

  ## select param and aniso of component cmp

  param <- variogram.object[[cmp]][["param"]]
  aniso <- variogram.object[[cmp]][["aniso"]]

  if( ncmp > 1L ){
    Va <- expand(Valpha[[cmp]][["Valpha"]])
  } else {
    Va <- NULL
  }

  switch(
    x[1],
    nugget = {

      ## nugget

      #       eeq.exp <- sum( diag(
      #           ( 1./TtT * lik.item[["Valphaxi"]][["Valphaxi.inverse"]] ) %*% r.cov[["cov.bhat"]] %*% lik.item[["Valphaxi"]][["Valphaxi.inverse"]]
      #         )
      #       )
      #       eeq.emp <- sum(
      #         ( lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] )^2 / TtT
      #       )

      eeq.exp <- sum( 1./TtT * Valphaxii * Valphaxii.cov.bhat )
      eeq.emp <- sum( bhVaxi^2 / TtT )

    },
    snugget = {

      ## snaugget

      #       eeq.exp <- sum(
      #         diag(
      #           lik.item[["Valphaxi"]][["Valphaxi.inverse"]] %*% lik.item[["Valphaxi"]][["Valphaxi.inverse"]] %*% r.cov[["cov.bhat"]]
      #         )
      #       )
      #       eeq.emp <- sum( lik.item[["zhat"]][["Valphaxi.inverse.bhat"]]^2 )

      eeq.exp <- sum( Valphaxii * Valphaxii.cov.bhat )
      eeq.emp <- sum( bhVaxi^2 )

    },
    variance = {

      ## variance

      #       eeq.exp <- sum(
      #         diag(
      #           lik.item[["Valphaxi"]][["Valphaxi.inverse"]] %*% lik.item[["Valphaxi"]][["Valpha"]] %*%
      #           lik.item[["Valphaxi"]][["Valphaxi.inverse"]] %*% r.cov[["cov.bhat"]]
      #         )
      #       )
      #       eeq.emp <- sum(
      #         lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] * drop( lik.item[["Valphaxi"]][["Valpha"]] %*% lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] )
      #       )

      if( identical( ncmp, 1L ) ){

        if( param["snugget"] > 0. ){

          aux <- -param["snugget"] * Valphaxii
          diag( aux ) <- diag( aux ) + 1.

          eeq.exp <- sum( aux * Valphaxii.cov.bhat )
          eeq.emp <- sum( bhVaxi * ( bh - param["snugget"] * bhVaxi ) )

        } else {

          eeq.exp <- sum( diag( Valphaxii.cov.bhat ) )
          eeq.emp <- sum( bhVaxi * bh )

        }

      } else {

        eeq.exp <- sum( pmm( Va, Valphaxii, control.pcmp ) * Valphaxii.cov.bhat )
        eeq.emp <- sum( bhVaxi * drop( Va %*% bhVaxi ) )

      }

    },
    {

      ## scale and extra parameters

      dVa <- partial.derivatives.variogram(
        d.param = x, x = lag.vectors,
        variogram.model = variogram.object[[cmp]][["variogram.model"]],
        param = param, aniso = aniso,
        sincos = variogram.object[[cmp]][["sincos"]],
        sclmat = variogram.object[[cmp]][["sclmat"]],
        rotmat = variogram.object[[cmp]][["rotmat"]],
        verbose = verbose
      )
      #       aux <- pmm( dVa, lik.item[["Valphaxi"]][["Valphaxi.inverse"]], control.pcmp )
      #       eeq.exp <- sum(
      #         pmm( lik.item[["Valphaxi"]][["Valphaxi.inverse"]], aux, control.pcmp ) * r.cov[["cov.bhat"]]
      #       )
      #       eeq.emp <- sum(
      #         lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] * drop( dVa %*% lik.item[["zhat"]][["Valphaxi.inverse.bhat"]] )
      #       )

      aux <- pmm( dVa, Valphaxii, control.pcmp )
      eeq.exp <- sum( aux * Valphaxii.cov.bhat )
      eeq.emp <- sum( bhVaxi * drop( dVa %*% bhVaxi ) )
    }
  )

  c( eeq.exp = eeq.exp, eeq.emp = eeq.emp )

}

################################################################################

f.aux.gradient.nll <- function(
  x,
  variogram.object,
  TT, TtT, XX,
  res,
  Qi, Vi, Qt11Vi,
  Valpha,
  bh, bhVi,
  lag.vectors,
  param.tf, deriv.fwd.tf,
  ml.method,
  control.pcmp, verbose
){

  ##  auxiliary function to compute gradient of (restricted) log-likelihood
  ##  (called by gradient.negative.loglikelihood)
  ##  original parametrization of variogram

  ## 2014-07-29 A. Papritz
  ## 2015-07-17 AP new parametrization of loglikelihood
  ## 2015-07-17 AP new function interface, improve efficiency
  ## 2015-07-27 AP changes to further improve efficiency
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-04 AP changes for nested variogram models
  ## 2020-02-07 AP sanity checks for switch() and if()

  ## correct parameters names and determine model component

  tmp <- unlist( strsplit( x, control.georob()[["sepstr"]], fixed = TRUE ) )
  x <- tmp[1L]
  cmp <- as.integer( tmp[2L] )

  ncmp <- length(variogram.object)

  ## select param and aniso of component cmp

  param <- variogram.object[[cmp]][["param"]]
  aniso <- variogram.object[[cmp]][["aniso"]]

  if( ncmp > 1L ){
    Va <- expand(Valpha[[cmp]][["Valpha"]])
  } else {
    Va <- NULL
  }


  ## compute constants

  t.q <- NROW(Vi)

  switch(
    x[1],
    nugget = {

      ## compute partial derivative of (restricted) log-likelihood with
      ## respect to nugget

      ## partial derivate of U with respect to nugget

      Ttres <- as.vector( tapply( res, factor( TT ), sum ) )
      dU <- -sum( Ttres^2 / TtT ) / param["nugget"]^2

      ## partial derivative of log(det(Q)) with respect to nugget

      if( identical( ml.method, "REML" ) ){
        TtTX <- TtT * XX
        dlogdetQ <- -sum(
          Qi * rbind(
            cbind( diag( TtT ),                TtTX   ),
            cbind(     t(TtTX), crossprod( XX, TtTX ) )
          )
        ) / param["nugget"]^2
      } else {
        dlogdetQ <- -sum( TtT * diag(Qi) ) / param["nugget"]^2
      }

      ## partial derivative of loglik with respect to transformed nugget

      result <- c( ( -0.5 * ( t.q / param["nugget"] + dlogdetQ + dU ) ) /
        deriv.fwd.tf[[param.tf["nugget"]]]( param["nugget"] )
      )

    },
    snugget = {

      ##  compute partial derivative of (restricted) log-likelihood with
      ##  respect to spatial nugget

      ## partial derivate of U with respect to spatial nugget

      dU <- -sum( bhVi^2 )

      ## partial derivative of log(det(V)) with respect to spatial nugget

      dlogdetV <- sum( diag( Vi ) )

      ## partial derivative of log( det( Q ) ) with respect to spatial nugget

      dlogdetQ <- -sum( Qt11Vi * Vi )

      ## partial derivative of loglik with respect to transformed spatial
      ## nugget

      result <- c(
        ( -0.5 * ( dlogdetV + dlogdetQ + dU ) ) /
        deriv.fwd.tf[[param.tf["snugget"]]]( param["snugget"] )
      )

    },
    variance = {

      ##  compute partial derivative of restricted log-likelihood with
      ##  respect to variance

      if( identical( ncmp, 1L ) ){

        if( param["snugget"] > 0. ){

          aux <- -param["snugget"] * Vi
          diag(aux) <- diag( aux ) + 1.

          ## partial derivate of U with respect to variance

          dU <- -sum( bhVi * ( bh - param["snugget"] * bhVi ) ) / param["variance"]

          # partial derivative of log(det(V)) with respect to variance

          dlogdetV <- sum( diag( aux ) ) / param["variance"]

          ## partial derivative of log(det Q)) with respect to variance

          dlogdetQ <- -sum( Qt11Vi * aux ) / param["variance"]

        } else {

          ## partial derivate of U with respect to variance

          dU <- -sum( bhVi * bh ) / param["variance"]

          # partial derivative of log(det(V)) with respect to variance

          dlogdetV <- NROW( Vi ) / param["variance"]

          ## partial derivative of log(det Q)) with respect to variance

          dlogdetQ <- -sum( diag( Qt11Vi ) ) / param["variance"]
        }

      } else {


        ## partial derivate of U with respect to variance

        dU <- -sum( bhVi * drop( Va %*% bhVi ) )

        # partial derivative of log(det(V)) with respect to variance

        dlogdetV <- sum( Vi * Va )

        ## partial derivative of log(det Q)) with respect to variance

        dlogdetQ <- -sum( Qt11Vi * pmm( Va, Vi, control.pcmp ) )

      }


      ## partial derivative of loglik with respect to transformed variance

      result <- c(
        ( -0.5 * ( dlogdetV + dlogdetQ + dU ) ) /
        deriv.fwd.tf[[param.tf["variance"]]]( param["variance"] )
      )

    },
    {

      ## compute partial derivative of (restricted) log-likelihood with
      ## respect to scale and extra parameters


      ## partial derivative of V_0(alpha) with respect to scale and extra
      ## parameters

      dVa <- partial.derivatives.variogram(
        d.param = x, x = lag.vectors,
        variogram.model = variogram.object[[cmp]][["variogram.model"]],
        param = param, aniso = aniso,
        sincos = variogram.object[[cmp]][["sincos"]],
        sclmat = variogram.object[[cmp]][["sclmat"]],
        rotmat = variogram.object[[cmp]][["rotmat"]],
        verbose = verbose
      )

      dVaVi <- pmm( dVa, Vi, control.pcmp )

      ##  derivate of U

      dU <- -sum( bhVi * drop( dVa %*% bhVi ) )

      ## partial derivatie of V

      dlogdetV <- sum( diag( dVaVi ) )

      ##  derivative of log(det(Q))

      dlogdetQ <- -sum( Qt11Vi * t( dVaVi ) )

      ## partial derivative of log(det(V))

      result <- ( -0.5 * ( dlogdetV + dlogdetQ + dU ) * param["variance"] ) /
        deriv.fwd.tf[[param.tf[x]]](
          c( param, aniso )[x]
        )
      names( result ) <- x

    }
  )

  names( result ) <- x

  return( result )

}

################################################################################

f.aux.gradient.npll <- function(
  x,
  variogram.object,
  eta, xi,
  TT, TtT, XX, col.rank.XX,
  res, Ustar,
  Qsi, Qst11Vai,
  Valphaxii, Valpha,
  bh, bhVaxi,
  lag.vectors,
  param.tf, deriv.fwd.tf,
  ml.method,
  control.pcmp, verbose
){

  ##  auxiliary function to compute gradient of (restricted) profile
  ##  log-likelihood (called by gradient.negative.loglikelihood)
  ## reparametrized variogram

  ## 2015-07-17 A. Papritz
  ## 2015-07-27 AP changes to improve efficiency
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-04 AP changes for nested variogram models
  ## 2020-02-07 AP sanity checks for switch() and if()

  ## correct parameters names and determine model component

  tmp <- unlist( strsplit( x, control.georob()[["sepstr"]], fixed = TRUE ) )
  x <- tmp[1L]
  cmp <- as.integer( tmp[2L] )

  ncmp <- length(variogram.object)

  ## select param and aniso of component cmp

  param <- variogram.object[[cmp]][["param"]]
  aniso <- variogram.object[[cmp]][["aniso"]]

  xi <- xi[cmp]

  if( ncmp > 1L ){
    VamI <- expand(Valpha[[cmp]][["Valpha"]])
    diag(VamI) <- diag(VamI) - 1.
  } else {
    VamI <- NULL
  }

  ## compute constants

  t.q <- t.qmp <- NROW(XX)
  if( identical( ml.method, "REML" ) ){
    t.qmp <- t.q - col.rank.XX[["rank"]]
  }

  switch(
    x[1],
    nugget = {

      ## compute partial derivative of (restricted) profile log-likelihood
      ## with respect to eta

      ## partial derivative of log(Ustar) with respect to eta

      Tt.res <- as.vector( tapply( res, factor( TT ), sum ) )
      dU <- sum( Tt.res^2 / TtT ) / Ustar

      ## partial derivative of log(det(Qstar)) with respect to eta

      if( identical( ml.method, "REML" ) ){
        TtTX <- TtT * XX
        dlogdetQ <- sum(
          Qsi * rbind(
            cbind( diag( TtT ),                TtTX   ),
            cbind(     t(TtTX), crossprod( XX, TtTX ) )
          )
        )
      } else {
        dlogdetQ <- sum( TtT * diag(Qsi) )
      }

      ## partial derivative of profile loglik with respect to transformed eta

      result <- c( -0.5 * ( t.qmp * dU - t.q / eta + dlogdetQ ) /
        deriv.fwd.tf[[param.tf["nugget"]]]( param["nugget"] )
      )

    },
    variance = {

      ## compute partial derivative of (restricted) profile log-likelihood
      ## with respect to xi

      if( identical( ncmp, 1L ) ){

        ## partial derivative of log(Ustar) with respect to xi

        dU <- -sum( bhVaxi * ( bh - bhVaxi ) ) / Ustar / xi

        ## partial derivative of log(det(Valphaxi)) with respect to xi

        dlogdetV <- ( NROW( Valphaxii ) - sum( diag(Valphaxii) ) ) / xi

        ## partial derivative of log(det(Qstar)) with respect to xi

        dlogdetQ <- -( sum( diag(Qst11Vai) ) - sum( Qst11Vai * Valphaxii ) ) / xi

      } else {

        ## partial derivative of log(Ustar) with respect to xi

        dU <- -sum( bhVaxi * drop( VamI %*% bhVaxi ) ) / Ustar

        ## partial derivative of log(det(Valphaxi)) with respect to xi

        dlogdetV <- sum( Valphaxii * VamI )

        ## partial derivative of log(det(Qstar)) with respect to xi

        dlogdetQ <- -sum( Qst11Vai * pmm( VamI, Valphaxii, control.pcmp ) )

      }

      ## partial derivative of profile loglik with respect to transformed xi

      result <- c(
         -0.5 * ( t.qmp * dU + dlogdetV + dlogdetQ ) /
        deriv.fwd.tf[[param.tf["variance"]]]( param["variance"] )
      )

    },
    {

      ## compute partial derivative of (restricted) profile log-likelihood
      ## with respect to scale and extra parameters


      ## partial derivative of V_0(alpha) with respect to scale and extra
      ## parameters

      dVa <- partial.derivatives.variogram(
        d.param = x, x = lag.vectors,
        variogram.model = variogram.object[[cmp]][["variogram.model"]],
        param = param, aniso = aniso,
        sincos = variogram.object[[cmp]][["sincos"]],
        sclmat = variogram.object[[cmp]][["sclmat"]],
        rotmat = variogram.object[[cmp]][["rotmat"]],
        verbose = verbose
      )

      dVaVai <- pmm( dVa, Valphaxii, control.pcmp )  ### !!!dVaVai not symmetric!!!!

      ## partial derivative of log(Ustar) (up to factor xi)

      dU <- -sum( bhVaxi * drop( dVa %*% bhVaxi ) ) / Ustar

      ## partial derivative of log(det(Valphaxi)) (up to factor xi)

      dlogdetV <- sum( diag( dVaVai ) )

      ## partial derivative of log(det(Qstar)) (up to factor (1-xi))

      dlogdetQ <- -sum( Qst11Vai * t( dVaVai ) )

      ## partial derivative of profile loglik

      result <- -0.5 * xi * ( t.qmp * dU + dlogdetV + dlogdetQ ) /
        deriv.fwd.tf[[param.tf[x]]](
          c( param, aniso )[x]
        )

    }
  )

  names( result ) <- x

  return( result )

}

################################################################################

f.aux.Qstar <- function(
  TT, TtT,
  XX, col.rank.XX, min.condnum,
  Vi,
  eta,
  ml.method, control.pcmp
){

  ## auxiliary function to compute matrix Qstar used for Gaussian
  ## log-likelihood (called by likelihood.calculations)

  ## 2014-07-29 A. Papritz
  ## 2015-07-17 AP new function interface and new name
  ## 2016-07-20 AP changes for parallel computations
  ## 2023-12-20 AP replacement of identical(class(...), ...) by inherits(..., ...)

  result <- list( error = TRUE, log.det.Qstar = NULL, Qstar.inverse = NULL )

  TtTX <- eta * TtT * XX

  ##  compute matrix Qstar

  Qstar <-  Vi
  diag( Qstar ) <- diag( Qstar ) + eta *TtT

  if( identical( ml.method, "REML" ) ){
    Qstar <- rbind(
      cbind( Qstar,       TtTX                 ),
      cbind( t(TtTX), crossprod( XX, TtTX) )
    )
  }

  if( col.rank.XX[["deficient"]] && ml.method == "REML" ){

    ## compute log(pseudo.det(Qstar)) and (Moore-Penrose) pseudo inverse of Qstar by svd

    result[["error"]] <- FALSE
    s <- svd( Qstar )
    sel <- s[["d"]] / max( s[["d"]] ) > min.condnum
    #     result[["log.det.Qstar"]] <- sum( log( s[["d"]][s[["d"]] / max( s[["d"]] ) > min.condnum] ) )
    #     s[["d"]] <- ifelse( s[["d"]] / max( s[["d"]] ) <= min.condnum, 0., 1. / s[["d"]] )
    #         result[["Qstar.inverse"]] <- s[["v"]] %*% ( s[["d"]] * t( s[["u"]] ) )
    result[["log.det.Qstar"]] <- sum( log( s[["d"]][sel] ) )
    s[["d"]] <- ifelse( sel,  1. / s[["d"]], 0. )
    result[["Qstar.inverse"]] <- pmm(
      s[["v"]], s[["d"]] * t( s[["u"]] ), control.pcmp
    )

  } else {

    ##  compute log(det(Qstar)) and inverse of Qstar by cholesky decomposition

    t.chol <- try( chol( Qstar ), silent = TRUE )

    if( !inherits( t.chol, "try-error" ) ) {

      result[["error"]] <- FALSE
      result[["log.det.Qstar"]] <- 2. * sum( log( diag( t.chol) ) )
      result[["Qstar.inverse"]] <- chol2inv( t.chol )

    }

  }

  result

}


################################################################################

f.stop.cluster <- function( clstr = NULL, fname ){

  ## function to stop snow and snowfall clusters
  ## 2014-07-31 A. Papritz
  ## 2023-12-20 AP options(error = NULL) deleted
  ## 2024-02-01 AP new argument fname

  if( sfIsRunning() ){
    sfStop()
  }

  if( file.exists( fname ) ){
    if( is.null( clstr ) ) load( fname )
    file.remove( fname )
  }

  ## stop cluster started by child processes in recursive paralellized
  ## computations

  if( !is.null( clstr ) ){
    junk <- parLapply( clstr, 1L:length(clstr), function( i ) snowfall::sfStop() )
    junk <- stopCluster( clstr )
  }

}



################################################################################

## functions to (back)transform variogram parameters for alternative
## parametrization of variogram for Gaussian (RE)ML

## 2016-08-03 A. Papritz
## 2016-08-04 AP changes for nested variogram models
## 2020-02-07 AP sanity checks for switch() and if()

f.reparam.fwd <- function( object, set.fit.param = FALSE ){

  ##  variance of signal process

  var.signal <- sum( unlist( lapply(
        object, function(x)
        x[["param"]][names(x[["param"]]) %in% c("variance", "snugget")]
      )))

  ##  extract fitting flags of all signal variance parameters

  if( set.fit.param ){
    fit.param <- lapply(
      object, function(x)
      x[["fit.param"]][names(x[["fit.param"]]) %in% c("variance", "snugget")]
    )
  } else {
    fit.param = NULL
  }

  res <- lapply(
    1L:length( object ),
    function( i, object, fit.param, var.signal ){

      x <- object[[i]]

      ## reparametrize

      x[["param"]]["variance"] <- x[["param"]]["variance"] / var.signal

      if( i == 1L ){

        x[["param"]]["nugget"] <- var.signal / x[["param"]]["nugget"]
        x[["param"]]["snugget"] <- x[["param"]]["snugget"] / var.signal

      }

      ##  set fitting flags for reparametrized parameters

      if( !is.null( fit.param ) ){

        if( i == 1L ){

          if( x[["fit.param"]]["snugget"] ){
            x[["fit.param"]]["snugget"] <- FALSE
          } else {
            x[["fit.param"]]["variance"] <- FALSE
          }

        } else {

          if( x[["fit.param"]]["variance"] && sum( unlist( fit.param[1L:(i-1L)] ) ) <= 0L ){
            x[["fit.param"]]["variance"] <- FALSE
          }

        }
      }

      x
    }, object = object, fit.param = fit.param, var.signal = var.signal
  )
  attr( res, "var.signal" ) <- var.signal
  if( set.fit.param ) attr( res, "fit.param" ) <- fit.param
  res
}

f.reparam.bwd <- function( object, var.signal, restore.fit.param = FALSE ){

  if( missing( var.signal ) ) var.signal <- attr( object, "var.signal" )

  if( restore.fit.param ) fit.param <- attr( object, "fit.param" )

  res <- lapply(
    1L:length( object ),
    function( i, object, fit.param, var.signal ){

      ##  reparametrize

      x <- object[[i]]
      x[["param"]]["variance"] <- x[["param"]]["variance"] * var.signal
      if( i == 1L ){
        x[["param"]]["nugget"] <- var.signal / x[["param"]]["nugget"]
        x[["param"]]["snugget"] <- x[["param"]]["snugget"] * var.signal
      }

      ##  restore fitting flags

      if( restore.fit.param ){
        x[["fit.param"]]["variance"] <- fit.param[[i]]["variance"]
        if( i == 1L ){
          x[["fit.param"]]["snugget"] <- fit.param[[i]]["snugget"]
        }
      }

      x
    }, object = object, fit.param = fit.param, var.signal = var.signal
  )
  res
}


################################################################################

f.aux.RSS <- function( res, TT, TtT, bhat, Valphaxi.inverse.bhat, eta ){

  ## function computes residual sums of squares required for
  ## likelihood computations

  ## 2015-07-17 A. Papritz

  Ttres <- res

  if( sum( duplicated( TT ) > 0. ) ){
    Ttres <- as.vector( tapply( Ttres, factor( TT ), sum ) )
  }

  eta * sum( Ttres^2 / TtT ) + sum( bhat * Valphaxi.inverse.bhat )

}


##  ###########################################################################

## auxiliary function to extract diagonal of square matrix and to return
## x unchanged otherwise

## 2019-12-13 AP correcting use of class() in if() and switch()

f.diag <- function( x ){
  switch(
    class( x )[1],
    matrix = diag( x ),
    x
  )
}

##  ###########################################################################

## auxiliary function to transform variogram parameters in a variogram.object

f.aux.tf.param.fwd <- function( variogram.object, param.tf, fwd.tf ){

  transformed.param.aniso <- lapply(
    1L:length(variogram.object),
    function( i, x, param.tf ){

      x <- x[[i]]

      ## create local copies of objects

      variogram.model <- x[["variogram.model"]]
      param <- x[["param"]]
      param.name <- names( param )
      aniso <- x[["aniso"]]
      aniso.name <- names( aniso )

      ##  preparation for variogram parameter transformations

      all.param.tf <- param.tf

      t.sel <- match( param.name, names( all.param.tf ) )

      if( any( is.na( t.sel ) ) ){
        stop( "transformation undefined for some variogram parameters" )
      } else {
        param.tf <- all.param.tf[t.sel]
      }
      param.tf <- sapply(
        param.tf,
        function( x, vm ) if( length(x) > 1L ) x[vm] else x,
        vm = variogram.model
      )
      names( param.tf ) <- param.name

      ##  transform initial variogram parameters

      transformed.param <- sapply(
        param.name,
        function( x, fwd.tf, param.tf, param ) fwd.tf[[param.tf[x]]]( param[x] ),
        fwd.tf = fwd.tf,
        param.tf = param.tf,
        param = param
      )

      names( transformed.param ) <- param.name

      ##  preparation for anisotropy parameter transformations

      t.sel <- match( aniso.name, names( all.param.tf ) )

      if( any( is.na( t.sel ) ) ){
        stop( "transformation undefined for some anisotropy parameters" )
      } else {
        aniso.tf <- all.param.tf[t.sel]
      }
      aniso.tf <- sapply(
        aniso.tf,
        function( x ) if( length(x) > 1L ) x[variogram.model] else x
      )
      names( aniso.tf ) <- aniso.name

      ##  transform initial anisotropy parameters

      transformed.aniso <- sapply(
        aniso.name,
        function( x, fwd.tf, aniso.tf, aniso ){
          fwd.tf[[aniso.tf[x]]]( aniso[x] )
        },
        fwd.tf = fwd.tf,
        aniso.tf = aniso.tf,
        aniso = aniso
      )
      names( transformed.aniso ) <- aniso.name

      ## return transformed parameters

      res <- c( transformed.param, transformed.aniso )
      attr.res <- c( param.tf, aniso.tf )
      names(res) <- names(attr.res) <- paste( names(res), i, sep=control.georob()[["sepstr"]])

      attr( res, "tf.param.aniso" ) <- attr.res

      res

    }, x = variogram.object, param.tf = param.tf
  )

  tf.param.aniso <- unlist(lapply(
      transformed.param.aniso,
      function(x) attr( x, "tf.param.aniso" )
    ))

  transformed.param.aniso <- unlist( transformed.param.aniso )

  fit.param.aniso <- unlist( lapply(
      1L:length(variogram.object),
      function( i, object ){
        x <- object[[i]]
        res <- c( x[["fit.param"]], x[["fit.aniso"]] )
        names(res) <- paste( names(res), i, sep=control.georob()[["sepstr"]])
        res
      }, object = variogram.object
    ))

  list(
    transformed.param.aniso = transformed.param.aniso,
    tf.param.aniso = tf.param.aniso,
    fit.param.aniso = fit.param.aniso
  )

}


##  ###########################################################################

## auxiliary function to transform variogram parameters in a variogram.object

f.aux.print.gradient <- function( x, reparam = FALSE ){

  tmp <- strsplit( names(x), control.georob()[["sepstr"]], fixed = TRUE )
  names( x ) <- sapply( tmp, function(x) x[1L] )
  cmp <- sapply( tmp, function(x) as.integer(x[2L]) )

  x <- tapply( x, factor( cmp ), function( x ) x, simplify = FALSE )

  lapply(
    x,
    function( x, reparam ){
      if( reparam ){
        tmp <- names( x )
        tmp <- gsub(
          "nugget", "eta", gsub(
            "variance", "xi", gsub(
              "snugget", "1-sum(xi)", tmp, fixed = TRUE ), fixed = TRUE ), fixed = TRUE )
        names( x ) <- tmp
      }
      cat( "\n                      ",
        format( names( x ), width = 14L, justify = "right" ),
        "\n", sep = ""
      )
      cat( "  Gradient           :",
        format(
          signif( x, digits = 7L ),
          scientific = TRUE, width = 14L
        ), "\n" , sep = ""
      )
    }, reparam = reparam
  )
}


##  ###########################################################################

## auxiliary functions to manipulate georob call

## 2017-08-08 A. Papritz

##  ####################

# ## set main x argument to fixed value (e.g verbose = 1 )

f.call.set_x_to_value <- function( cl, x, value ){

  ## arguments:

  ## cl:      a call to function georob
  ## x:       name of argument to be changed
  ## value:   value passed to x

  if( x %in% names(cl) ) cl <- cl[ -match(x, names(cl)) ]
  res <- c( as.list(cl), x =list(value) )
  tmp <- names(res)
  tmp[length(tmp)] <- x
  names(res) <- tmp
  as.call( res )

}

##  ####################

## set a main argument of a main argument function in a georob call to a
## fixed value

f.call.set_x_to_value_in_fun <- function( cl, arg, fun, x, value ){

  ## Arguments

  ## arg:     name of argument to which fun is passed
  ## fun:     name of function passed to arg
  ## x:       name of argument of fun to be changed
  ## value:   value passed to x

  if( arg %in% names(cl) ){

    ## georob called with cl.arg argument

    cl.arg <- as.list( cl[[arg]] )
    cl <- cl[ -match( arg, names(cl) ) ]
    if( x %in% names(cl.arg) ){
      cl.arg[x] <- list( value )
    } else {
      cl.arg <- c( cl.arg, list(value) )
      tmp <- names(cl.arg)
      tmp[length(tmp)] <- x
      names(cl.arg) <- tmp
    }

  } else {

    ## georob called without cl.arg argument

    cl.arg <- list( as.symbol(fun), value )
    names(cl.arg) <- c( "", x )
  }

  res <- as.call( c( as.list(cl), arg = as.call(cl.arg) ) )
  tmp <- names( res )
  tmp[length(tmp)] <- arg
  names(res) <- tmp
  res

}


##  ####################

## prevent recursive parallelization when updating georob objects

f.call.prevent_recursive_parallelization <- function( cl ){

  ## 2024-02-21 A. Papritz

  ## cl:   a call to function georob

  if( !"pcmp" %in% names( cl[["control"]] ) ){

    ## add pcmp argument if not present in call
    tmp1 <- f.call.set_x_to_value_in_fun(
      cl, "control", "control.georob", "pcmp",
      control.pcmp(pmm.ncores = 1L, gcr.ncores = 1L)
    )
    tmp2 <- c( NA, tmp1[["control"]][["pcmp"]] )
    tmp2[[1]] <- as.symbol("control.pcmp")

  } else {

    ## modifiy pcmp argument if present in call
    tmp1 <- cl
    tmp2 <- as.list(tmp1[["control"]][["pcmp"]])
    if( "pmm.ncores" %in% names( tmp2 ) ){
      tmp2[["pmm.ncores"]] <- 1L
    } else {
      tmp2 <- c( tmp2, pmm.names = 1L )
    }
    if( "gcr.ncores" %in% names( tmp2 ) ){
      tmp2[["gcr.ncores"]] <- 1L
    } else {
      tmp2 <- c( tmp2, gcr.ncores = 1L )
    }

  }
  tmp1[["control"]][["pcmp"]] <- as.call( tmp2 )

  tmp1

}


##  ####################

## set all fit.param and fit.aniso equal to FALSE

f.call.set_allfitxxx_to_false <- function( cl ){

  ## 2023-12-20 AP checking class by inherits()
  ## arguments

  ## cl:   a call to function georob

  if( !inherits( cl, "call" ) ) stop(
    "'cl' must be of class 'call'"
  )

  if( !"variogram.object" %in% names(cl) ){

    x <- as.list(cl)

    sel <- c(
      grep( "^fit.p", names(x), fixed = FALSE ),
      grep( "^fit.a", names(x), fixed = FALSE )

    )

    fit.param <- c( list( as.symbol("default.fit.param" ) ),
      as.list( c( variance = FALSE, nugget = FALSE, scale = FALSE ) )
    )
    fit.aniso <- list( as.symbol("default.fit.aniso" ) )

    cl <- as.call(
      c(
        if( length(sel) ) x[-sel] else x,
        list( fit.param = as.call( fit.param ) ),
        list( fit.aniso = as.call( fit.aniso ) )
      )
    )

  } else {

    cl.vo <- cl[["variogram.object"]]
    cl.m.vo <- cl[ -match("variogram.object", names(cl)) ]

    tmp <- lapply(
      as.list( cl.vo[-1L] ),
      function( x ){

        x <- as.list(x)

        sel <- c(
          grep( "^fit.p", names(x), fixed = FALSE ),
          grep( "^fit.a", names(x), fixed = FALSE )

        )
        fit.param <- c( list( as.symbol("default.fit.param" ) ),
          as.list( c( variance = FALSE, nugget = FALSE, scale = FALSE ) )
        )
        fit.aniso <- list( as.symbol("default.fit.aniso" ) )

        as.call(
          c(
            if( length(sel) ) x[-sel] else x,
            list( fit.param = as.call( fit.param ) ),
            list( fit.aniso = as.call( fit.aniso ) )
          )
        )
      }
    )

    cl.vo.new <- as.call( c( as.list(cl.vo[1L]), as.list(tmp) ) )
    cl <- as.call( c( as.list(cl.m.vo), variogram.object = as.call( cl.vo.new) ) )

  }

  cl

}

##  ####################

## set one fit.param or one fit.aniso equal to FALSE

f.call.set_onefitxxx_to_value <- function( cl, nme, value, i = NULL ){

  ## 2023-12-14 AP checking class by inherits()

  ## arguments

  ## cl:     a call to function georob
  ## nme:    name of parameter that should be change
  ## value:  new value for nme
  ## i:      index of variogram component for which parameter should be fixed


  ## auxiliary function

  f.aux <- function( x, nme, value ){

    if( nme %in% c( "f1", "f2", "omega", "phi", "zeta" ) ){

      ## anisotropy parameter

      sel <- grep( "^fit.a", names(x), fixed = FALSE )

      if( length(sel) ){

        ## fit.aniso argument present in cl

        fit.aniso <- as.list(x[[sel]])

        ## match names of fit.aniso

        if( length(names(fit.aniso)) > 1L ){
          tmp <- sapply( names(fit.aniso)[-1L], match.arg, choices = names(default.fit.aniso()) )
          names(fit.aniso) <- c( "", tmp )
        }

        ## set new value for nme

        if( nme %in% names(fit.aniso) ){

          ## nme present in fit.aniso

          fit.aniso[nme] <- value

        } else {

          ## nme not present in fit.aniso

          tmp <- list( value )
          names( tmp ) <- nme
          fit.aniso <- c( fit.aniso, tmp )

        }

      } else {

        ## no fit.aniso argument in cl

        fit.aniso <- c( list( as.symbol("default.fit.aniso" ) ), list( value ) )
        names(fit.aniso) <- c( "", nme )

      }

      cl <- as.call(
        c(
          if( length(sel) ) x[-sel] else x,
          list( fit.aniso = as.call( fit.aniso ) )
        )
      )


    } else {

      ## variogram parameter

      sel <- grep( "^fit.p", names(x), fixed = FALSE )

      if( length(sel) ){

        ## fit.param argument present in cl

        fit.param <- as.list(x[[sel]])

        ## match names of fit.param

        if( length(names(fit.param)) > 1L ){
          tmp <- sapply( names(fit.param)[-1L], match.arg, choices = names(default.fit.param()) )
          names(fit.param) <- c( "", tmp )
        }

        ## set new value for nme

        if( nme %in% names(fit.param) ){

          ## nme present in fit.param

          fit.param[nme] <- value

        } else {

          ## nme not present in fit.param

          tmp <- list( value )
          names( tmp ) <- nme
          fit.param <- c( fit.param, tmp )
        }


      } else {

        ## no fit.param argument in cl

        fit.param <- c( list( as.symbol("default.fit.param" ) ), list( value ) )
        names(fit.param) <- c( "", nme )

      }

      cl <- as.call(
        c(
          if( length(sel) ) x[-sel] else x,
          list( fit.param = as.call( fit.param ) )
        )
      )

    }
  }

  ## start of main body

  if( !inherits( cl, "call" ) ) stop(
    "'cl' must be of class 'call'"
  )

  if( !"variogram.object" %in% names(cl) ){

    x <- as.list(cl)

    cl.new <- f.aux( x, nme, value )

  } else {

    cl.vo <- cl[["variogram.object"]]
    cl.m.vo <- cl[ -match("variogram.object", names(cl)) ]

    if( is.null(i) || !i %in% 1L:length( as.list( cl.vo[-1L] ) ) ) stop( "'i' wrong" )

    tmp <- lapply(
      1L:length( as.list( cl.vo[-1L] ) ),
      function( i, x, ii, nme, value ){

        x <- as.list(x[[i]])

        if( identical( i, ii ) ){

          f.aux( x, nme, value )

        } else {

          as.call( x )

        }

      }, x = as.list( cl.vo[-1L] ), ii = as.integer(i), nme = nme, value = value
    )

    cl.vo.new <- as.call( c( as.list(cl.vo[1L]), as.list(tmp) ) )
    cl.new <- as.call( c( as.list(cl.m.vo), variogram.object = as.call( cl.vo.new) ) )

  }

  cl.new

}



##  ####################

## set all initial values of variogram parameters in call to fitted values
## and possibly update value for variogram.model

f.call.set_allxxx_to_fitted_values <- function( object ){

  ## arguments

  ## object:   a georob object

  ## 2018-01-17 ap  also update of variogram.model
  ## 2023-12-14 AP checking class by inherits()

  if( !inherits( object, "georob" ) ) stop(
    "'object' must be of class 'georob'"
  )

  cl <- object[["call"]]

  if( !"variogram.object" %in% names(cl) ){

    x <- as.list( cl )

    sel <- c(
      grep( "^v", names(x), fixed = FALSE ),
      grep( "^p", names(x), fixed = FALSE ),
      grep( "^a", names(x), fixed = FALSE )

    )

    variogram.model <- c( list(as.symbol("c")), as.list( object[["variogram.object"]][[1L]][["variogram.model"]] ))
    param <- c( list(as.symbol("c")), as.list( object[["variogram.object"]][[1L]][["param"]] ))
    aniso <- c( list(as.symbol("c")), as.list( object[["variogram.object"]][[1L]][["aniso"]] ))

    cl <- as.call(
      c(
        if( length(sel) ) x[-sel] else x,
        list( variogram.model = as.call( variogram.model ) ),
        list( param = as.call( param ) ),
        list( aniso = as.call( aniso ) )
      )
    )

  } else {

    cl.vo <- cl[["variogram.object"]]
    cl.m.vo <- cl[ -match("variogram.object", names(cl)) ]

    tmp <- lapply(
      1L:length( as.list( cl.vo[-1L] ) ),
      function( i, x, vo ){

        x <- as.list( x[[i]] )
        vo <- vo[[i]]

        sel <- c(
          grep( "^v", names(x), fixed = FALSE ),
          grep( "^p", names(x), fixed = FALSE ),
          grep( "^a", names(x), fixed = FALSE )
        )

        variogram.model <- c( list(as.symbol("c")), as.list( vo[["variogram.model"]] ))
        param <- c( list(as.symbol("c")), as.list( vo[["param"]] ))
        aniso <- c( list(as.symbol("c")), as.list( vo[["aniso"]] ))

        as.call(
          c(
            if( length(sel) ) x[-sel] else x,
            list( variogram.model = as.call( variogram.model )),
            list( param = as.call( param )),
            list( aniso = as.call( aniso ))
          )
        )

      }, x = as.list( cl.vo[-1L] ), vo = object[["variogram.object"]]
    )

    cl.vo.new <- as.call( c( as.list(cl.vo[1L]), as.list(tmp) ) )
    cl <- as.call( c( as.list(cl.m.vo), variogram.object = as.call( cl.vo.new) ) )

  }

  cl

}


##  ####################

## set one param or one aniso equal to FALSE

f.call.set_onexxx_to_value <- function( cl, nme, value, i = NULL ){

  ## 2023-12-14 AP checking class by inherits()

  ## arguments

  ## cl:     a call to function georob
  ## nme:    name of parameter that should be changed
  ## value:  new value for nme
  ## i:      index of variogram component for which parameter should be fixed


  ## auxiliary function

  f.aux <- function( x, nme, value ){

    if( nme %in% c( "f1", "f2", "omega", "phi", "zeta" ) ){

      ## anisotropy parameter

      sel <- grep( "^a", names(x), fixed = FALSE )

      if( length(sel) ){

        ## aniso argument present in cl

        aniso <- as.list(x[[sel]])

        ## match names of aniso

        if( length(names(aniso)) > 1L ){
          tmp <- sapply( names(aniso)[-1L], match.arg, choices = names(default.aniso()) )
          names(aniso) <- c( "", tmp )
        }

        ## set new value for nme

        if( nme %in% names(aniso) ){

          ## nme present in aniso

          aniso[nme] <- value

        } else {

          ## nme not present in aniso

          tmp <- list( value )
          names( tmp ) <- nme
          aniso <- c( aniso, tmp )

        }

      } else {

        ## no aniso argument in cl

        aniso <- c( list( as.symbol("c" ) ), list( value ) )
        names(aniso) <- c( "", nme )

      }

      cl <- as.call(
        c(
          if( length(sel) ) x[-sel] else x,
          list( aniso = as.call( aniso ) )
        )
      )


    } else {

      ## variogram parameter

      sel <- grep( "^p", names(x), fixed = FALSE )

      if( length(sel) ){

        ## param argument present in cl

        param <- as.list(x[[sel]])

        ## match names of param

        if( length(names(param)) > 1L ){
          tmp <- sapply( names(param)[-1L], match.arg,
            choices = names(param.transf()) )
          names(param) <- c( "", tmp )
        }

        ## set new value for nme

        if( nme %in% names(param) ){

          ## nme present in param

          param[nme] <- value

        } else {

          ## nme not present in param

          tmp <- list( value )
          names( tmp ) <- nme
          param <- c( param, tmp )
        }


      } else {

        ## no param argument in cl

        param <- c( list( as.symbol("c" ) ), list( value ) )
        names(param) <- c( "", nme )

      }

      cl <- as.call(
        c(
          if( length(sel) ) x[-sel] else x,
          list( param = as.call( param ) )
        )
      )

    }
  }

  ## start of main body

  if( !inherits( cl, "call" ) ) stop(
    "'cl' must be of class 'call'"
  )

  if( !"variogram.object" %in% names(cl) ){

    x <- as.list(cl)

    cl.new <- f.aux( x, nme, value )

  } else {

    cl.vo <- cl[["variogram.object"]]
    cl.m.vo <- cl[ -match("variogram.object", names(cl)) ]

    if( is.null(i) || !i %in% 1L:length( as.list( cl.vo[-1L] ) ) ) stop( "'i' wrong" )

    tmp <- lapply(
      1L:length( as.list( cl.vo[-1L] ) ),
      function( i, x, ii, nme, value ){

        x <- as.list(x[[i]])

        if( identical( i, ii ) ){

          f.aux( x, nme, value )

        } else {

          as.call( x )

        }

      }, x = as.list( cl.vo[-1L] ), ii = as.integer(i), nme = nme, value = value
    )

    cl.vo.new <- as.call( c( as.list(cl.vo[1L]), as.list(tmp) ) )
    cl.new <- as.call( c( as.list(cl.m.vo), variogram.object = as.call( cl.vo.new) ) )

  }

  cl.new

}
