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
  if( identical( class( aux ), "try-error" ) ){
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

update.zhat <-
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
  if( !identical( class( t.chol ), "try-error" ) ){

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
    if( !identical( class( t.chol ), "try-error" ) ){
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

      new <- update.zhat(
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

#### -- check consistency of arguments

  if( !is.null( gcr.constant ) ){
    if( !is.list( gcr.constant ) ||
      !identical( length(gcr.constant), length(variogram.object) )
    ) stop( "lengths of 'gcr.constant' and 'variogram.object' differ" )
  } else {
    gcr.constant <- as.list( rep( NA_real_, length(variogram.object) ) )
  }

#### -- compute generalized covariances

  res <- lapply(
    1L:length(variogram.object),
    function( i, x, gcr.constant, lag.vectors, control.pcmp ){

      variogram.model <- x[[i]][["variogram.model"]]
      param           <- x[[i]][["param"]]
      param           <- param[!names(param) %in% c( "variance", "snugget", "nugget")]
      aniso           <- x[[i]][c("aniso", "sclmat", "rotmat")]
      gcr.constant    <- gcr.constant[[i]]

      result <- list( error = TRUE )

      #       RFoptions(newAniso=FALSE) ## moved to georob.fit

#### --- preparation: anisotropic model

      if( NCOL( lag.vectors ) > 1L ){

        ## matrix for coordinate transformation

        A <- aniso[["sclmat"]] * aniso[["rotmat"]] / param["scale"]

        ## prepare model

        model.list <- list( variogram.model )
        model.list <- c( model.list, as.list( param[-1L] ) )
        model.list <- list( "$", var = 1., A = A, model.list )

        ##  negative semivariance matrix

        ## functions of version 3 of RandomFields

        ## auxiliary function to compute generalized correlations in parallel

        f.aux <- function( i ){

          ## objects s, e, lag.vectors, model.list taken from parent environment

          #           RandomFields Version 3.1
          #           result <- try(
          #             -RFvariogram(
          #               x = lag.vectors[s[i]:e[i], ], model = model.list,
          #               dim = attr( lag.vectors, "ndim.coords" ), grid = FALSE
          #             ),
          #             silent = TRUE
          #           )

          RFoptions( allow_duplicated_locations = TRUE )

          ## note: RFfctn computes covariance for stationary and negative
          ## semivariance for IRF models, required is negative semivariance
          result <- try(
            RFfctn(
              x = lag.vectors[s[i]:e[i], ], model = model.list,
              dim = attr( lag.vectors, "ndim.coords" ), grid = FALSE
            ),
            silent = TRUE
          )

          RFoptions( allow_duplicated_locations = FALSE )

          if( !(identical( class( result ), "try-error" ) || any( is.na( result ) )) ){
            if(!variogram.model %in% irf.models){
              ## subtract variance for stationary models
              result <- result - model.list[["var"]]
            }
            result
          } else {
            "RFvariogram.error"
          }
        }

      } else {

#### --- preparation: isotropic model

        ## matrix for coordinate transformation

        A <- 1. / param["scale"]

        ## prepare model

        model.list <- list( variogram.model )
        model.list <- c( model.list, as.list( param[-1L] ) )
        model.list <- list( "$", var = 1., A = A, model.list )

        ##  negative semivariance matrix

        ## functions of version 3 of RandomFields

        ## auxiliary function to compute generalized correlations in parallel

        f.aux <- function( i ){

          ## objects s, e, lag.vectors, model.list taken from parent environment

          #           RandomFields Version 3.1
          #           result <- try(
          #             -RFvariogram(
          #               x = lag.vectors[s[i]:e[i]], model = model.list,
          #               grid = FALSE
          #             ),
          #             silent = TRUE
          #           )

          RFoptions( allow_duplicated_locations = TRUE )

          ## note: RFfctn computes covariance for stationary and negative
          ## semivariance for IRF models, required is negative semivariance
          result <- try(
            RFfctn(
              x = lag.vectors[s[i]:e[i]], model = model.list,
              grid = FALSE
            ),
            silent = TRUE
          )

          RFoptions( allow_duplicated_locations = FALSE )

          if( !(identical( class( result ), "try-error" ) || any( is.na( result ) )) ){
            if(!variogram.model %in% irf.models){
              ## subtract variance for stationary models
              result <- result - model.list[["var"]]
            }
            result
          } else {
            "RFvariogram.error"
          }
        }

      }

#### --- compute covariances

      ## determine number of cores

      ncores <- control.pcmp[["gcr.ncores"]]
      if( !control.pcmp[["allow.recursive"]] ) ncores <- 1L

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
          options( error = f.stop.cluster )

          junk <- sfInit( parallel = TRUE, cpus = ncores )

          junk <- sfLibrary( georob, verbose = FALSE )
          #         junk <- sfLibrary( RandomFields, verbose = FALSE )

          junk <- sfExport( "s", "e", "lag.vectors", "model.list" )

         }


        Valpha <- sfLapply( 1L:k, f.aux )

        if( control.pcmp[["sfstop"]] ){
          junk <- sfStop()
          options( error = NULL )
        }

      } else {

        Valpha <- mclapply( 1L:k, f.aux, mc.cores = ncores )

      }

      not.ok <- any( sapply(
          Valpha,
          function( x ) identical( x, "RFvariogram.error" ) || any(is.na(x))
        ))

      if( !not.ok ){

        Valpha <- unlist( Valpha )

        ##  compute additive constant for positive definiteness, this
        ##  implements a sufficient condition for positive definiteness of
        ##  Valpha (strong row sum criterion)

        if( is.na( gcr.constant ) ){
          if( variogram.model %in% irf.models ){
            gcr.constant <- -min( Valpha ) * 2.
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

        warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
        if( verbose > 3. ) cat(
          "\n an error occurred when computing the negative semivariance matrix\n"
        )

      }

      return( result )

    }, x = variogram.object, gcr.constant = gcr.constant,
    lag.vectors = lag.vectors, control.pcmp = control.pcmp
  )

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
    #         -( (1.+B) * (2.+B) * (1.-hs[sel])^B * hs[sel] )
    #       } else if( identical( A, 2 ) ){
    #         -( (3.+B) * (4.+B) * (1.-hs[sel])^(1.+B) * hs[sel] * ( 1. + hs[sel] + B*hs[sel]) ) / 3.
    #       } else if( identical( A, 3 ) ){
    #         -(
    #           (5.+B) * (6.+B) * (1.-hs[sel])^(2.+B) * hs[sel] * ( 3. + 3. * (2.+B) * hs[sel] + (1.+B) * (3.+B) * hs[sel]^2 )
    #         ) / 15.
    #       } else {
    #         stop( "gengneiting model undefined for 'n' != 1:3" )
    #       }
    #
    #       result <- rep( 0., length( hs ) )
    #
    #       switch(
    #         d.param[1],
    #         scale = dgc.dhs * dhs.dscale,
    #         alpha = {
    #           result[sel] <- if( identical( A, 1 ) ){
    #             (1.-hs[sel])^(1.+B) * ( hs[sel] + (1. + hs[sel] + B*hs[sel]) * log( 1.-hs[sel]) )
    #
    #           } else if( identical( A, 2 ) ){
    #             (
    #               (1.-hs[sel])^(2.+B) * (
    #                 hs[sel] * ( 3. + 2. * (2.+B) *hs[sel] ) +
    #                 ( 3. + 3. * ( 2.+B) * hs[sel] + ( 1.+B) * (3.+B) * hs[sel]^2 ) * log( 1.-hs[sel] )
    #               )
    #             ) / 3.
    #           } else if( identical( A, 3 ) ){
    #             (
    #               (1.-hs[sel])^(3.+B) * (
    #                 hs[sel] * ( 15. + hs[sel] * ( 36. + 23.*hs[sel] + 3. * B * ( 4. + (6.+B)*hs[sel] ) ) ) +
    #                 ( 15. + 15. * (3.+B) * hs[sel] + ( 45. + 6. * B * (6.+B) ) * hs[sel]^2 + (1.+B) * (3.+B) * (5.+B) * hs[sel]^3 ) *
    #                 log( 1.-hs[sel])
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
        stop( "RMgengneiting model undefined for 'n' != 1:3" )
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

##   ##############################################################################
### georob.fit

georob.fit <-
  function(
    initial.objects,
    variogram.object,
    param.tf,
    fwd.tf,
    deriv.fwd.tf,
    bwd.tf,
    georob.object,
    safe.param,
    tuning.psi,
    error.family.estimation, error.family.cov.effects, error.family.cov.residuals,
    cov.bhat, full.cov.bhat,
    cov.betahat,
    cov.bhat.betahat,
    cov.delta.bhat, full.cov.delta.bhat,
    cov.delta.bhat.betahat,
    cov.ehat, full.cov.ehat,
    cov.ehat.p.bhat, full.cov.ehat.p.bhat,
    aux.cov.pred.target,
    min.condnum, col.rank.XX,
    psi.func,
    tuning.psi.nr,
    ml.method,
    maximizer,
    reparam,
    irwls.initial,
    irwls.maxiter,
    irwls.ftol,
    force.gradient,
    zero.dist,
    control.nleqslv,
    control.optim,
    control.nlminb,
    hessian,
    control.pcmp,
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
  ## 2013-06-12 AP changes in stored items of Valphaxi object
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-03 AP new transformation of rotation angles
  ## 2013-07-09 AP catching errors occuring when fitting anisotropic
  ##               variograms with default anisotropy parameters
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-02-18 AP correcting error when fitting models with offset
  ## 2014-05-28 AP change in check for initial variogram parameter values
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2015-03-10 AP changes for reparametrization of variogram
  ## 2015-03-16 AP elimination of unused variables, own function for psi function
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  ## 2015-07-17 AP Gaussian (RE)ML estimation for reparametrized variogram,
  ##               nlminb added as maximizer of loglikelihood
  ## 2015-07-23 AP changes for avoiding computation of Valphaxi object if not needed,
  ##               rearrangement of output item (Valpha.objects, zhat.objects)
  ## 2015-08-19 AP variances of eps and psi(eps/sigma) for long-tailed error distribution;
  ##               computing covariances of residuals under long-tailed error model,
  ##               control about error families for computing covariances added
  ## 2015-12-02 AP catching error in computation of covariances
  ## 2016-01-26 AP refined check of initial values of variogram parameters
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-07-22 AP corrected names of gradient components
  ## 2016-08-08 AP changes for nested variogram models
  ## 2016-08-10 AP changes for isotropic variogram models
  ## 2016-11-14 AP correcting error in 3d rotation matrix for geometrically anisotropic variograms
  ## 2016-11-28 AP checking ml.method and presence of intercept for intrinsic models
  ## 2017-02-24 AP warning for negative definite hessian
  ## 2017-07-26 AP changes to allow non-zero snugget if there are no replicated observations
  ## 2020-02-19 AP minor change in computation of hessian
  ## 2020-03-27 AP computing Hessian (observed Fisher information) of untransformed parameters

  ##  ToDos:

#### -- preparations

  ##  main body of georob.fit

  RFoptions(newAniso=FALSE)

  d2r <- pi / 180.

  ##  define rho-function and derivatives (suppress temporarily warnings issued by gamma())

  old.op <- options( warn = -1 )
  rho.psi.etc <- f.psi.function( x = psi.func, tp = tuning.psi )
  options( old.op )

  ##  set number of IRWLS iterations for estimating bhat and betahat to
  ##  1 for non-robust REML case

  if( tuning.psi >= tuning.psi.nr ){
    irwls.maxiter <- 1L
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
  first.dist0 <- unname( apply( dist0, 1L, function( x ) ( (1L:length(x))[x])[1L] ) )


  TT <- matrix( 0L, nrow = length( yy ), ncol = length( yy ) )
  TT[ cbind( 1L:NROW(TT), first.dist0 ) ] <- 1L
  rep.obs <- (1L:NCOL(TT))[ apply( TT, 2L, function( x ) all( x == 0L ) ) ]
  if( length( rep.obs ) > 0L )  TT <- TT[, -rep.obs]

  ## check whether explanatory variables are the identical for the replicated
  ## observations and issue an error if not

  apply(
    TT,
    2L,
    function( i, XX ){
      XX <- XX[as.logical(i), , drop = FALSE]
      apply(
        XX,
        2L,
        function( x ){
          if( length(x) > 1L && any( x[-1L] != x[1L] ) ) stop(
            "explanatory variables differ for some replicated observations"
          )
        }
      )
    },
    XX = XX
  )

  ## store row indices of replicated observations only

  TT <- drop( TT %*% 1L:NCOL( TT ) )
  TtT <- as.vector( table( TT ) )

  ##  omit elements corresponding to replicated observations in XX, bhat
  ##  and coordinates

  if( length( rep.obs ) > 0L ) {
    XX          <- XX[ -rep.obs, , drop = FALSE]
    bhat      <- bhat[ -rep.obs ]
    coordinates <- coordinates[ -rep.obs, , drop = FALSE]
    if( verbose > 0. ) cat( "\n", length(rep.obs), "replicated observations at",
      length( unique( TT[rep.obs] ) ), "sampling locations\n"
    )
  }

#### -- process contents of variogram.object

   variogram.object <- lapply(
    variogram.object,
    function( x, TT, d2r, n, has.intercept ){

      ## create local copies of objects

      variogram.model <- x[["variogram.model"]]
      param <- x[["param"]]
      fit.param <- x[["fit.param"]]
      aniso <- x[["aniso"]]
      fit.aniso <- x[["fit.aniso"]]

      ##  check whether fitting of chosen variogram model is implemented and
      ##  return names of extra parameters (if any)

      ep <- param.names( model = variogram.model )

      ## check whether reml method is chosen to fit intrinsic variogram

      if( variogram.model %in% control.georob()[["irf.models"]] ){
        if( ml.method == "ML" && tuning.psi >= tuning.psi.nr ) stop(
          "models with intrinsic variograms must be estimated by REML"
        )
        if( !has.intercept ) stop(
          "models with intrinsic variograms require a drift model with an intercept"
        )
      }

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

      ## rearrange and check flags controlling variogram parameter fitting

      if(
        variogram.model %in% (t.models <- c( "RMfbm" ) ) &&
        sum( duplicated( TT ) == 0L ) && all(
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

      ## rearrange initial anisotropy parameters

      aniso <- aniso[aniso.name]
      fit.aniso <- fit.aniso[aniso.name]

      ## check whether initial values of anisotropy parameters are valid

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


      ## convert angles to radian

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

    }, TT = TT, d2r = d2r, n = NCOL(coordinates),
    has.intercept = attr(
      terms( initial.objects[["initial.fit"]] ), "intercept"
    ) == 1L
  )

  #   print(str(variogram.object))

#### -- set consistent values for nugget and snugget of nested variogram
  ## models

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

  if( sum( duplicated( TT ) ) == 0L ){
    #     variogram.object[[1L]][["param"]]["nugget"] <- sum(  variogram.object[[1L]][["param"]][c("nugget", "snugget") ])
    #     variogram.object[[1L]][["fit.param"]]["nugget"] <- any(  variogram.object[[1L]][["fit.param"]][c("nugget", "snugget") ])
    #     variogram.object[[1L]][["param"]]["snugget"] <- 0.
    variogram.object[[1L]][["fit.param"]]["snugget"] <- FALSE
  }

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

  ## check whether nugget can be robustly estimated

  if( identical( length(unique( TtT )), 1L ) && tuning.psi < tuning.psi.nr &&
    variogram.object[[1L]][["fit.param"]]["snugget"]
  ) stop( "'snugget' cannot be estimated robustly if all sites have the same number of replicated measurements" )

  #   print(str(variogram.object))

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

  #   print(str(variogram.object))

#### -- reparametrize variograms for Gaussian (RE)ML

  ## reparametrize if there is more than one variance parameter to fit

  original.variogram.object <- variogram.object
  original.param.tf <- param.tf

  fit.param <- unlist( lapply(
    variogram.object, function(x)
    x[["fit.param"]][names(x[["fit.param"]]) %in% c("variance", "snugget", "nugget")]
  ))
  reparam <- reparam && tuning.psi >= tuning.psi.nr && sum( fit.param ) > 1L

  reparam.variogram.object <- f.reparam.fwd( variogram.object, set.fit.param = TRUE )

  if( reparam ){
    variogram.object <- reparam.variogram.object
    param.tf[c("variance", "snugget")] <- "logit1"
  }

  #   print(str(variogram.object))

#### -- transform variogram parameters

  tmp <- f.aux.tf.param.fwd( variogram.object, param.tf, fwd.tf )

  transformed.param.aniso <- tmp[["transformed.param.aniso"]]
  tf.param.aniso <- tmp[["tf.param.aniso"]]
  fit.param.aniso <- tmp[["fit.param.aniso"]]


  #   print(transformed.param.aniso)
  #   print(fit.param.aniso)
  #   print(tf.param.aniso)

  ## compute lag vectors for all pairs of coordinates (or restore from
  ## georob.object if available and coordinates are the same)

  if(
    !is.null( georob.object ) &&
    isTRUE( all.equal( georob.object[["locations.objects"]][["coordinates"]], coordinates ) )
  ){
    lag.vectors <- georob.object[["locations.objects"]][["lag.vectors"]]
  } else {
    if( !all( sapply( variogram.object, function(x) x[["isotropic"]] ) ) ){
      indices.pairs <- combn( NROW( coordinates ), 2L )
      lag.vectors <- coordinates[ indices.pairs[2L,], ] - coordinates[ indices.pairs[1L,], ]
    } else {
      lag.vectors <- as.vector( dist( coordinates ) )
    }
    attr( lag.vectors, "ndim.coords" ) <- NCOL(coordinates)
  }

  #   print(str(lag.vectors))

 #### -- create environment to store items required to compute likelihood and
  ##  estimating equations that are provided by
  ##  likelihood.calculations

  envir <- new.env()
  lik.item <- list()

  ##  initialize values of variogram object item stored in the environment

  lik.item[["variogram.object"]] <- lapply(
    variogram.object,
    function(x) x[c("variogram.model", "param", "isotropic",
      "aniso", "sincos", "sclmat", "rotmat")]
  )
  lik.item[["var.signal"]] <- attr( reparam.variogram.object, "var.signal" )
  lik.item[["eta"]] <- unname( reparam.variogram.object[[1L]][["param"]]["nugget"] )
  lik.item[["xi"]] <-  unname( sapply(
      reparam.variogram.object, function(x) x[["param"]]["variance"]
    ))

#### -- restore Valphaxi object from georob.object if available and variogram
  ## parameters are the same

  if( !is.null( georob.object ) ){

    tmp <- georob.object[["variogram.object"]]
    if( reparam ) tmp <- f.reparam.fwd( tmp )

    tmp <- unlist(
      lapply(
        tmp,
        function( x, d2r ) c( x[["param"]], x[["aniso"]] * c( 1., 1., rep( d2r, 3L ) )
        ), d2r = d2r
      )
    )

    if(
      isTRUE(
        all.equal(
          tmp,
          unlist(
            lapply(
              variogram.object,
              function(x) c( x[["param"]], x[["aniso"]] )
            )
          )
        )) &&
      isTRUE(
        all.equal(
          length( georob.object[["Valphaxi.objects"]][["Valphaxi"]][["diag"]] ),
          NROW( XX )
        ))
    ){
      lik.item[["Valphaxi"]]  <- expand( georob.object[["Valphaxi.objects"]] )
    }

  }


  assign( "lik.item", lik.item, pos = as.environment( envir ) )

#### -- compute various expectations of psi, and its derivative etc.

  expectations <- numeric()

  ##  ... E[ psi'(x) ]  with respect to nominal Gaussian model

  if( is.null( rho.psi.etc[["exp.gauss.dpsi"]] ) ){

    t.exp <- integrate(
      function( x, dpsi.function, tuning.psi ) {
        dnorm( x ) * dpsi.function( x, tuning.psi )
      },
      lower = -Inf, upper = Inf,
      dpsi.function = rho.psi.etc[["dpsi.function"]],
      tuning.psi = tuning.psi
    )
    if( !identical( t.exp[["message"]], "OK" ) ) stop( t.exp[["message"]] )
    expectations["exp.gauss.dpsi"] <- t.exp[["value"]]

  } else {

    expectations["exp.gauss.dpsi"] <- rho.psi.etc[["exp.gauss.dpsi"]]( tuning.psi )

  }

  ##  ... E[ psi(x)^2 ]  with respect to nominal Gaussian model

  if( is.null( rho.psi.etc[["var.gauss.psi"]] ) ){

    t.exp <- integrate(
      function( x, psi.function, tuning.psi ) {
        dnorm( x ) * ( psi.function( x, tuning.psi ) )^2
      },
      lower = -Inf, upper = Inf,
      psi.function = rho.psi.etc[["psi.function"]],
      tuning.psi = tuning.psi
    )
    if( !identical( t.exp[["message"]], "OK" ) ) stop( t.exp[["message"]] )
    expectations["var.gauss.psi"] <- t.exp[["value"]]

  } else {

    expectations["var.gauss.psi"] <- rho.psi.etc[["var.gauss.psi"]]( tuning.psi )

  }

  ## ...  E[ x^2 ] with respect to assumed longtailed distribution
  ## f0 \propto 1/sigma exp( - rho(x/sigma) )

  if( is.null( rho.psi.etc[["var.f0.eps"]] ) ){

    t.exp <- integrate(
      function( x, f0, tuning.psi ) {
        f0( x, tuning.psi, sigma = 1. ) * x^2
      },
      lower = -Inf, upper = Inf,
      f0 = rho.psi.etc[["f0"]],
      tuning.psi = tuning.psi
    )
    if( !identical( t.exp[["message"]], "OK" ) ) stop( t.exp[["message"]] )
    expectations["var.f0.eps"] <- t.exp[["value"]]

  } else {

    expectations["var.f0.eps"] <- rho.psi.etc[["var.f0.eps"]](
      tuning.psi, sigma = 1.
    )

  }

  ## ...  E[ psi(x)^2 ] ( = E[ psi'(x) ]) with respect to assumed longtailed distribution
  ## f0 \propto 1/sigma exp( - rho(x/sigma) )

  expectations["var.f0.psi"] <- rho.psi.etc[["var.f0.psi"]]( tuning.psi )

  if( verbose > 2. ) cat(
    "\n expectation of psi'(epsilon/sigma) under nominal Gaussian model  :",
    signif( expectations["exp.gauss.dpsi"] ), "\n",
    "variance of psi(epsilon/sigma) under nominal Gaussian model      :",
    signif( expectations["var.gauss.psi"] ), "\n",
    "variance of epsilon under long-tailed model                      :",
    signif( expectations["var.f0.eps"] ), "\n",
    "variance of psi(epsilon/sigma) under long-tailed model           :",
    signif( expectations["var.f0.psi"] ), "\n"
  )


#### -- compute signal zhat

  sel <- !is.na (betahat )
  zhat <- drop( XX[, sel, drop=FALSE] %*% betahat[sel] + bhat )
  names( zhat ) <- rownames( XX )

  r.hessian <- NULL

  if( tuning.psi < tuning.psi.nr ) {

#### -- robust REML estimation

    hessian <- FALSE

    if( any( fit.param.aniso ) ){

      ##  find roots of estimating equations

      r.root <- nleqslv(
        x = transformed.param.aniso[fit.param.aniso],
        fn = estimating.equations.theta,
        method = control.nleqslv[["method"]],
        global = control.nleqslv[["global"]],
        xscalm = control.nleqslv[["xscalm"]],
        control = control.nleqslv[["control"]],
        envir = envir,
        fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
        name.param.aniso = names(transformed.param.aniso),
        tf.param.aniso = tf.param.aniso,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat,
        psi.function = rho.psi.etc[["psi.function"]],
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        force.gradient = force.gradient,
        expectations = expectations,
        error.family.estimation = error.family.estimation,
        control.pcmp = control.pcmp,
        verbose = verbose
      )

      #       r.param <- r.root[["x"]] names( r.param ) <- names(
      #       transformed.param[ fit.param ] )

      r.gradient <- r.root[["fvec"]]
      names( r.gradient ) <- names( transformed.param.aniso[fit.param.aniso] )

      r.converged <- r.root[["termcd"]] == 1L
      r.convergence.code <- r.root[["termcd"]]

      r.counts <- c( nfcnt = r.root[["nfcnt"]], njcnt = r.root[["njcnt"]] )

    } else {

      ##  all variogram parameters are fixed

      ##  evaluate estimating equations

      r.gradient <- estimating.equations.theta(
        adjustable.param.aniso = transformed.param.aniso[fit.param.aniso],
        envir = envir,
        fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
        name.param.aniso = names(transformed.param.aniso),
        tf.param.aniso = tf.param.aniso,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat,
        psi.function = rho.psi.etc[["psi.function"]],
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        force.gradient = force.gradient,
        expectations = expectations,
        error.family.estimation = error.family.estimation,
        control.pcmp = control.pcmp,
        verbose = verbose
      )

      r.converged <- NA
      r.convergence.code <- NA_integer_
      r.counts <- c( nfcnt = NA_integer_, njcnt = NA_integer_ )

    }

    r.opt.neg.loglik <- NA_real_

  } else {

    if( any( fit.param.aniso ) ){

#### -- Gaussian (RE)ML estimation

      error.family.cov.effects <- error.family.cov.residuals <- "gaussian"

      if( identical( maximizer, "optim" ) ){

#         print("optim")

        r.opt.neg.restricted.loglik <- optim(
          par = transformed.param.aniso[fit.param.aniso],
          fn = negative.loglikelihood,
          gr = gradient.negative.loglikelihood,
          method = control.optim[["method"]],
          lower = control.optim[["lower"]],
          upper = control.optim[["upper"]],
          control = control.optim[["control"]],
          hessian = FALSE,
          envir = envir,
          fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
          name.param.aniso = names(transformed.param.aniso),
          tf.param.aniso = tf.param.aniso,
          deriv.fwd.tf = deriv.fwd.tf,
          bwd.tf = bwd.tf,
          safe.param = safe.param,
          reparam = reparam,
          lag.vectors = lag.vectors,
          XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
          yy = yy, TT = TT, TtT = TtT, zhat = zhat,
          psi.function = rho.psi.etc[["psi.function"]],
          tuning.psi = tuning.psi,
          tuning.psi.nr = tuning.psi.nr,
          ml.method = ml.method,
          irwls.initial = irwls.initial,
          irwls.maxiter = irwls.maxiter,
          irwls.ftol = irwls.ftol,
          control.pcmp = control.pcmp,
          verbose = verbose,
          force.gradient = force.gradient
        )

        r.opt.neg.loglik <- r.opt.neg.restricted.loglik[["value"]]
        r.converged <- r.opt.neg.restricted.loglik[["convergence"]] == 0L
        r.convergence.code <- r.opt.neg.restricted.loglik[["convergence"]]
        r.counts <- r.opt.neg.restricted.loglik[["counts"]]

      } else {

#         print("nlminb")

        r.opt.neg.restricted.loglik <- nlminb(
          start = transformed.param.aniso[fit.param.aniso],
          objective = negative.loglikelihood,
          gradient = gradient.negative.loglikelihood,
          control = control.nlminb[["control"]],
          lower = control.nlminb[["lower"]],
          upper = control.nlminb[["upper"]],
          envir = envir,
          fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
          name.param.aniso = names(transformed.param.aniso),
          tf.param.aniso = tf.param.aniso,
          deriv.fwd.tf = deriv.fwd.tf,
          bwd.tf = bwd.tf,
          safe.param = safe.param,
          reparam = reparam,
          lag.vectors = lag.vectors,
          XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
          yy = yy, TT = TT, TtT = TtT, zhat = zhat,
          psi.function = rho.psi.etc[["psi.function"]],
          tuning.psi = tuning.psi,
          tuning.psi.nr = tuning.psi.nr,
          ml.method = ml.method,
          irwls.initial = irwls.initial,
          irwls.maxiter = irwls.maxiter,
          irwls.ftol = irwls.ftol,
          control.pcmp = control.pcmp,
          verbose = verbose,
          force.gradient = force.gradient
        )

        r.opt.neg.loglik <- r.opt.neg.restricted.loglik[["objective"]]
        r.converged <- r.opt.neg.restricted.loglik[["convergence"]] == 0L
        r.convergence.code <- r.opt.neg.restricted.loglik[["convergence"]]
        r.counts <- r.opt.neg.restricted.loglik[["evaluations"]]

      }

      r.gradient <- gradient.negative.loglikelihood(
        adjustable.param.aniso = r.opt.neg.restricted.loglik[["par"]],
        envir = envir,
        fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
        name.param.aniso = names(transformed.param.aniso),
        tf.param.aniso = tf.param.aniso,
        deriv.fwd.tf = deriv.fwd.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        reparam = reparam,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat,
        psi.function = rho.psi.etc[["psi.function"]],
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        force.gradient = force.gradient,
        control.pcmp = control.pcmp,
        verbose = verbose
      )

    } else {

      ##  all variogram parameters are fixed

      hessian <- FALSE

      ##  compute negative restricted loglikelihood and gradient

      r.opt.neg.loglik <- negative.loglikelihood(
        adjustable.param.aniso = transformed.param.aniso[fit.param.aniso],
        envir = envir,
        fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
        name.param.aniso = names(transformed.param.aniso),
        tf.param.aniso = tf.param.aniso,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        reparam = reparam,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat,
        psi.function = rho.psi.etc[["psi.function"]],
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        control.pcmp = control.pcmp,
        verbose = verbose
      )

      r.gradient <- gradient.negative.loglikelihood(
        adjustable.param.aniso = transformed.param.aniso[fit.param.aniso],
        envir = envir,
        fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
        name.param.aniso = names(transformed.param.aniso),
        tf.param.aniso = tf.param.aniso,
        deriv.fwd.tf = deriv.fwd.tf,
        bwd.tf = bwd.tf,
        safe.param = safe.param,
        reparam = reparam,
        lag.vectors = lag.vectors,
        XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
        yy = yy, TT = TT, TtT = TtT, zhat = zhat,
        psi.function = rho.psi.etc[["psi.function"]],
        tuning.psi = tuning.psi,
        tuning.psi.nr = tuning.psi.nr,
        ml.method = ml.method,
        irwls.initial = irwls.initial,
        irwls.maxiter = irwls.maxiter,
        irwls.ftol = irwls.ftol,
        force.gradient = force.gradient,
        control.pcmp = control.pcmp,
        verbose = verbose
      )

      r.converged <- NA
      r.convergence.code <- NA_integer_
      r.counts <- c( nfcnt = NA_integer_, njcnt = NA_integer_ )

    }

  }

#   ## set the correct parameter names
#
#   if( reparam ){
#     tmp <- names( r.gradient )
#     tmp <- gsub(
#       "nugget", "eta", gsub(
#         "variance", "xi", gsub(
#           "snugget", "1-sum(xi)", tmp, fixed = TRUE ), fixed = TRUE ), fixed = TRUE )
#     names( r.gradient ) <- tmp
#   }

#### -- get the other fitted items

  lik.item <- get( "lik.item", pos = as.environment( envir ) )

  ## revert if necessary to original parametrization of variogram.object
  ## and create fitted.variogram.object with all angles mapped to
  ## half-circle

  if( reparam ){

    ## compute variance of signal

    t.qmp <- NROW(XX)
    if( identical( ml.method, "REML" ) ){
      t.qmp <- t.qmp - col.rank.XX[["rank"]]
    }
    t.var.signal <- lik.item[["zhat"]][["RSS"]] / t.qmp

    ## revert

    lik.item[["variogram.object"]] <- f.reparam.bwd(
      lik.item[["variogram.object"]],
      var.signal = t.var.signal
    )

  }

  fitted.variogram.object <- lapply(
    1L:length(original.variogram.object),
    function( i, ori, fit, d2r ){
      ori <- ori[[i]]
      fit <- fit[[i]]

      ## map angles to halfcircle

      aniso <- fit[["aniso"]] / c( rep( 1., 2L ), rep( d2r, 3L ) )

      if( !ori[["isotropic"]] ){
        if( aniso["omega"] < 0. ){
          aniso["omega"] <- aniso["omega"] + 180.
        }
        if( aniso["omega"] > 180. ){
          aniso["omega"] <- aniso["omega"] - 180.
        }
        if( aniso["phi"] < 0. ){
          aniso["phi"] <- aniso["phi"] + 180.
        }
        if( aniso["phi"] > 180. ){
          aniso["phi"] <- aniso["phi"] - 180.
        }
        if( aniso["zeta"] < 90. ){
          aniso["zeta"] <- aniso["zeta"] + 180.
        }
        if( aniso["zeta"] > 90. ){
          aniso["zeta"] <- aniso["zeta"] - 180.
        }
      }

      ori[["param"]] <- fit[["param"]]
      ori[["aniso"]] <- aniso
      ori[["sincos"]] <- fit[["sincos"]]
      ori[["sclmat"]] <- fit[["sclmat"]]
      ori[["rotmat"]] <- fit[["rotmat"]]
      ori
    }, ori = original.variogram.object, fit = lik.item[["variogram.object"]], d2r = d2r
  )

#### -- extract or recompute Hessian for Gaussian (RE)ML

  if( hessian ){

    #     if( reparam || identical( maximizer, "nlminb" ) ){

    ## transform variogram parameters

    if( reparam ) param.tf[c("variance", "snugget")] <- original.param.tf[c("variance", "snugget")]

    tmp <- f.aux.tf.param.fwd(
      fitted.variogram.object,
      param.tf,
      fwd.tf
    )

    transformed.param.aniso <- tmp[["transformed.param.aniso"]]
    tf.param.aniso <- tmp[["tf.param.aniso"]]
    fit.param.aniso <- tmp[["fit.param.aniso"]]

    ## Hessian of with respect to transformed parameters
    
    r.hessian <- optimHess(
      par = transformed.param.aniso[ fit.param.aniso ],
      fn = negative.loglikelihood,
      gr = gradient.negative.loglikelihood,
      control = control.optim[["control"]],
      envir = envir,
      fixed.param.aniso = transformed.param.aniso[!fit.param.aniso],
      name.param.aniso = names(transformed.param.aniso),
      tf.param.aniso = tf.param.aniso,
      deriv.fwd.tf = deriv.fwd.tf,
      bwd.tf = bwd.tf,
      safe.param = safe.param,
      reparam = FALSE,
      lag.vectors = lag.vectors,
      XX = XX, min.condnum = min.condnum, col.rank.XX = col.rank.XX,
      yy = yy, TT = TT, TtT = TtT, zhat = zhat,
      psi.function = rho.psi.etc[["psi.function"]],
      tuning.psi = tuning.psi,
      tuning.psi.nr = tuning.psi.nr,
      ml.method = ml.method,
      irwls.initial = irwls.initial,
      irwls.maxiter = irwls.maxiter,
      irwls.ftol = irwls.ftol,
      control.pcmp = control.pcmp,
      verbose = 0.,
      force.gradient = force.gradient
    )

    ## check whether Hessian is positive definite

    if( any( eigen(r.hessian)[["values"]] < 0. ) ) warning(
      "hessian not positive definite, check whether local minimum of log-likelihood has been found"
    )
    
    ## Hessian with respect to untransformed parameters
    ## see email exchange with Victor De Oliveira
    
    ## get vector of untransformed fitted variogram parameters
    
    untransformed.param.aniso <- sapply(
      names(transformed.param.aniso)[fit.param.aniso],
      function(x){
        bwd.tf[[tf.param.aniso[[x]]]](transformed.param.aniso[[x]])
      }
    )
    
    ## compute Jacobian matrix of forward-transformation
    
    jacobian.fwd.tf <- sapply(
      names(transformed.param.aniso)[fit.param.aniso],
      function(x){
        deriv.fwd.tf[[tf.param.aniso[[x]]]](untransformed.param.aniso[[x]])
      }
    )
    
    ## compute Hessian of untransformed parameters
    
    r.hessian.untransformed.param.aniso <- jacobian.fwd.tf * t( 
      jacobian.fwd.tf * t(r.hessian) 
    )


    #   } else {
    #
    #       r.hessian <- r.opt.neg.restricted.loglik[["hessian"]]
    #
    #     }
    #
    #     print(r.hessian)

  }

 #### -- compute the covariances of fixed and random effects under nominal
  ##  Gaussian model

  r.cov <- list()

  if( any( c(
        cov.bhat, cov.betahat, cov.bhat.betahat,
        cov.delta.bhat, cov.delta.bhat.betahat,
        aux.cov.pred.target
      )
    )
  ){

    ##  compute the covariances of fixed and random effects under nominal
    ##  Gaussian model

    r.cov <- covariances.fixed.random.effects(
      Valphaxi.objects = lik.item[["Valphaxi"]][c("Valphaxi", "Valphaxi.inverse")],
      Aalphaxi = lik.item[["zhat"]][["Aalphaxi"]],
      Palphaxi = lik.item[["zhat"]][["Palphaxi"]],
      Valphaxi.inverse.Palphaxi = lik.item[["zhat"]][["Valphaxi.inverse.Palphaxi"]],
      rweights = lik.item[["zhat"]][["rweights"]],
      XX = XX, TT = TT, TtT = TtT, names.yy = names( yy ),
      nugget = lik.item[["variogram.object"]][[1L]][["param"]]["nugget"],
      eta = lik.item[["eta"]],
      expectations = expectations, family = error.family.cov.effects,
      cov.bhat = cov.bhat, full.cov.bhat = full.cov.bhat,
      cov.betahat = cov.betahat,
      cov.bhat.betahat = cov.bhat.betahat,
      cov.delta.bhat = cov.delta.bhat, full.cov.delta.bhat = full.cov.delta.bhat,
      cov.delta.bhat.betahat = cov.delta.bhat.betahat,
      cov.ehat = FALSE, full.cov.ehat = FALSE,
      cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
      aux.cov.pred.target = aux.cov.pred.target,
      control.pcmp = control.pcmp,
      verbose = verbose
    )

    if( r.cov[["error"]] ) {
      warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
      if( verbose > 0. ) cat(
        "\nan error occurred when computing the covariances of fixed and random effects\n"
      )
    } else {
      r.cov <- r.cov[!names(r.cov) %in% "error"]
    }

  }

#### -- compute covariances of residuals under assumed long-tailed error distribution

  if( any( c( cov.ehat, cov.ehat.p.bhat ) ) ){

    ##  compute the covariances of fixed and random effects under nominal
    ##  Gaussian model

    r.cov <- c(
      r.cov,
      covariances.fixed.random.effects(
        Valphaxi.objects = lik.item[["Valphaxi"]][c("Valphaxi", "Valphaxi.inverse")],
        Aalphaxi = lik.item[["zhat"]][["Aalphaxi"]],
        Palphaxi = lik.item[["zhat"]][["Palphaxi"]],
        Valphaxi.inverse.Palphaxi = lik.item[["zhat"]][["Valphaxi.inverse.Palphaxi"]],
        rweights = lik.item[["zhat"]][["rweights"]],
        XX = XX, TT = TT, TtT = TtT, names.yy = names( yy ),
        nugget = lik.item[["variogram.object"]][[1L]][["param"]]["nugget"],
        eta = lik.item[["eta"]],
        expectations = expectations, family = error.family.cov.residuals,
        cov.bhat = FALSE, full.cov.bhat = FALSE,
        cov.betahat = FALSE,
        cov.bhat.betahat = FALSE,
        cov.delta.bhat = FALSE, full.cov.delta.bhat = FALSE,
        cov.delta.bhat.betahat = FALSE,
        cov.ehat = cov.ehat, full.cov.ehat = full.cov.ehat,
        cov.ehat.p.bhat = cov.ehat.p.bhat, full.cov.ehat.p.bhat = full.cov.ehat.p.bhat,
        aux.cov.pred.target = FALSE,
        control.pcmp = control.pcmp,
        verbose = verbose
      )
    )

    if( r.cov[["error"]] ) {
      warning( "there were errors: call 'georob' with argument 'verbose' > 1" )
      if( verbose > 0. ) cat(
        "\nan error occurred when computing the covariances of residuals\n"
      )
    } else {
      r.cov <- r.cov[!names(r.cov) %in% "error"]
    }

  }

#   print(str(r.cov))

  ## stop SNOW and snowfall clusters

  #   f.stop.cluster()

  if( length( lik.item[["defaultCluster"]] ) > 0L ){
    cl <- lik.item[["defaultCluster"]]

    junk <- parLapply( cl, 1L:length(cl), function( i ) sfStop() )
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
  ##          Valphaxi = lik.item[["Valphaxi"]][["Valphaxi"]],
  ##          XX = XX,
  ##          param = lik.item[["param"]]
  ##      )

#### -- prepare output

  result.list <- list(
    loglik = -r.opt.neg.loglik,
    variogram.object = fitted.variogram.object,
    gradient = r.gradient,
    tuning.psi = tuning.psi,
    coefficients = lik.item[["zhat"]][["betahat"]],
    fitted.values = drop( XX %*% lik.item[["zhat"]][["betahat"]] )[TT],
    bhat = lik.item[["zhat"]][["bhat"]],
    residuals = lik.item[["zhat"]][["residuals"]],
    rweights = lik.item[["zhat"]][["rweights"]],
    converged = r.converged,
    convergence.code = r.convergence.code,
    iter = r.counts,
    Tmat = TT
  )
  names( result.list[["fitted.values"]] ) <- names( result.list[["residuals"]] )
  names( result.list[["rweights"]] )      <- names( result.list[["residuals"]] )

  if( any( c(
        cov.bhat, cov.betahat, cov.bhat.betahat,
        cov.delta.bhat, cov.delta.bhat.betahat,
        cov.ehat, cov.ehat.p.bhat, aux.cov.pred.target
      )
    )
  ){

    result.list[["cov"]] <- compress( r.cov )

  }

  result.list[["expectations"]]       <- expectations
  #   result.list[["Valphaxi.objects"]]   <- compress(
  #     lik.item[["Valphaxi"]][!names(lik.item[["Valphaxi"]]) %in% "Valpha"]
  #   )
  result.list[["Valphaxi.objects"]]   <- compress( lik.item[["Valphaxi"]][-1] )
  result.list[["Valphaxi.objects"]][["Valpha"]] <- lapply(
    result.list[["Valphaxi.objects"]][["Valpha"]],
    function(x) x[c("gcr.constant", "Valpha")]
  )
  result.list[["zhat.objects"]]       <- compress(
    lik.item[["zhat"]][c( "Aalphaxi", "Palphaxi", "Valphaxi.inverse.Palphaxi" )]
  )

  result.list[["locations.objects"]] <- initial.objects[["locations.objects"]]
  result.list[["locations.objects"]][["lag.vectors"]] <- lag.vectors

  result.list[["initial.objects"]] <- list(
    coefficients = initial.objects[["betahat"]],
    bhat = initial.objects[["bhat"]],
    variogram.object = original.variogram.object
  )
  if( !is.null( r.hessian ) ){
    result.list[["hessian.tfpa"]] <- r.hessian
    result.list[["hessian.ntfpa"]] <- r.hessian.untransformed.param.aniso
  }
  ##      result.list[["df.model"]] <- r.df

  #   print(str(result.list))

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

    if( !identical( class( t.chol ), "try-error" ) ) {

      result[["error"]] <- FALSE
      result[["log.det.Qstar"]] <- 2. * sum( log( diag( t.chol) ) )
      result[["Qstar.inverse"]] <- chol2inv( t.chol )

    }

  }

  result

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

  if( !identical( class( t.vchol ), "try-error" ) ) {

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

  ## stop cluster started by child processes in recursive paralellized
  ## computations

  if( !is.null( cl ) ){
    junk <- parLapply( cl, 1L:length(cl), function( i ) sfStop() )
    junk <- stopCluster( cl )
  }
  options( error = NULL )

}



################################################################################

f.psi.function <- function( x = c( "logistic", "t.dist", "huber" ), tp ){

  ## define psi-function and derivatives
  ## 2015-03-16 A. Papritz
  ## 2015-08-19 AP pdf and variances of eps and psi(eps/sigma) for long-tailed error distribution,
  ##               new parametrisation for t-dist family
  ## 2020-02-07 AP sanity checks for switch() and if()

  x <- match.arg( x )

  switch(
    x[1],
    logistic = list(
      #       rho.function = function( x, tuning.psi ) {
      #         tuning.psi * (-x + tuning.psi *
      #           ( -log(2.) + log( 1. + exp( 2. * x / tuning.psi ) ) )
      #         )
      #       },
      psi.function = function( x, tuning.psi ) {
        tuning.psi * tanh( x / tuning.psi )
      },
      dpsi.function = function( x, tuning.psi ) {
        1. / (cosh( x / tuning.psi ))^2
      },
      f0 = if( is.finite( gamma((1. + tp^2)/2.) ) ){
        function( x, tuning.psi, sigma = 1. ) {              # pdf f0 propto 1/sigma exp(-rho(eps/sigma))
          ifelse(
            is.finite( exp((2.*x)/(tuning.psi*sigma)) ),
            ( exp( (tuning.psi * ( x - tuning.psi * sigma * log(
                      (1. + exp((2.*x)/(tuning.psi*sigma)))/2.
                    )))/sigma ) * gamma((1. + tuning.psi^2)/2.)) /
            ( tuning.psi * sqrt(pi) * sigma * gamma(tuning.psi^2/2.) ),
            0.
          )
        }
      } else {
        function( x, tuning.psi, sigma = 1. ) dnorm( x, sd = sigma )
      },
      var.f0.eps = NULL,                        # variance of eps under f0 (= expectation of psi'(eps/sigma) under f0)
      var.f0.psi = function( tuning.psi){       # variance of psi(eps/sigma) under f0
        tuning.psi^2 / ( (1. + tuning.psi^2) )
      },
      exp.gauss.dpsi = NULL,                    # expectation of psi'(eps/sigma) under N(0, sigma^2)
      var.gauss.psi = NULL                      # variance of psi(eps/sigma) under N(0, sigma^2)
    ),
    t.dist = list(
      #       rho.function = function( x, tuning.psi ){
      #         0.5 * tuning.psi^2 * log( (x^2 + tuning.psi^2) / tuning.psi^2 )
      #       },
      psi.function = function( x, tuning.psi ){
        tuning.psi^2 * x / ( x^2 + tuning.psi^2 )
      },
      dpsi.function = function( x, tuning.psi ) {
        tuning.psi^2 * ( tuning.psi^2 - x^2 ) / ( x^2 + tuning.psi^2 )^2
      },
      f0 = if( is.finite( gamma((-1.+tp^2)/2.) ) ){
        function( x, tuning.psi, sigma = 1. ) {
          ( (tuning.psi^2)^(tuning.psi^2/2.) * gamma( tuning.psi^2/2.) ) /
          ( tuning.psi*sqrt(pi) * (tuning.psi^2 + (x/sigma)^2 )^(tuning.psi^2/2.) * sigma *
            gamma((-1.+tuning.psi^2)/2.) )
        }
      } else {
        function( x, tuning.psi, sigma = 1. ) dnorm( x, sd = sigma )
      },
      var.f0.eps = if( tp > sqrt(3.) ){
        function( tuning.psi, sigma = 1. ){
          ( tuning.psi*sigma )^2 / ( -3. + tuning.psi^2 )
        }
      } else NULL,
      var.f0.psi = function( tuning.psi ){
        1. - 3./(2. + tuning.psi^2)
      },
      exp.gauss.dpsi = if( is.finite( exp(tp^2/2.) ) ){
        function( tuning.psi ){
          tuning.psi^2 - tuning.psi^3 * exp(tuning.psi^2/2.) * sqrt(pi/2.) *
          2. * (1.-pnorm(tuning.psi))
        }
      } else {
        function( tuning.psi ){ 1. }
      },
      var.gauss.psi = if( is.finite( exp(tp^2/2.) ) ){
        function( tuning.psi ){
          ( tuning.psi^3 * (-2.*tuning.psi + (1 + tuning.psi^2) * exp(tuning.psi^2/2.) *
              sqrt(2.*pi) * 2. * (1.-pnorm(tuning.psi)) ) ) / 4.
        }
      } else {
        function( tuning.psi ){ 1. }
      }
      ),
    huber = list(
      #       rho.function = function( x, tuning.psi ) {
      #         ifelse(
      #           abs(x) <= tuning.psi,
      #           0.5 * x^2,
      #           tuning.psi * abs(x) - 0.5 * tuning.psi^2
      #         )
      #       },
      psi.function = function( x, tuning.psi ) {
        ifelse(
          abs(x) <= tuning.psi,
          x,
          sign(x) * tuning.psi
        )
      },
      dpsi.function = function( x, tuning.psi ) {
        ifelse( abs(x) <= tuning.psi, 1., 0. )
      },
      f0 = if( is.finite( (tp*exp(0.5*tp^2) ) ) ){
        function( x, tuning.psi, sigma = 1. ){
          1. / ( exp(
              ifelse(
                abs(x/sigma) <= tuning.psi,
                0.5 * (x/sigma)^2,
                tuning.psi * abs(x/sigma) - 0.5 * tuning.psi^2
              ))
            * ( (2.*sigma) / (tuning.psi*exp(0.5*tuning.psi^2) ) +
              sqrt(2*pi) * sigma * (2.*pnorm(tuning.psi)-1.))
          )
        }
      } else {
        function( x, tuning.psi, sigma = 1. ) dnorm( x, sd = sigma )
      },
      var.f0.eps = function( tuning.psi, sigma = 1. ){
        sigma^2 * ( 1. +
          ( sqrt(2./pi) * ( 2. + tuning.psi^2 ) ) /
          ( tuning.psi^2 * ( sqrt(2./pi) + tuning.psi *exp(0.5*tuning.psi^2) * (2.*pnorm(tuning.psi) - 1.) ) )
        )
      },
      var.f0.psi = function( tuning.psi ){
        1. + 1. / (
          -1. - 0.5*sqrt(2*pi) * tuning.psi * exp( 0.5*tuning.psi^2 ) * (2.*pnorm(tuning.psi) - 1.)
        )
      },
      exp.gauss.dpsi = function( tuning.psi ){
        (2.*pnorm(tuning.psi) - 1.)
      },
      var.gauss.psi = function( tuning.psi ){
        (2.*pnorm(tuning.psi) - 1.) +  tuning.psi*(-(sqrt(2./pi)/exp(tuning.psi^2/2.)) +
          tuning.psi * 2. * (1.-pnorm(tuning.psi)) )
      }
    )
  )

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

## set all fit.param and fit.aniso equal to FALSE

f.call.set_allfitxxx_to_false <- function( cl ){

  ## arguments

  ## cl:   a call to function georob

  if( !identical( class(cl)[1], "call" ) ) stop(
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

  if( !identical( class(cl)[1], "call" ) ) stop(
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

  if( !identical( class(object)[1], "georob" ) ) stop(
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

  if( !identical( class(cl)[1], "call" ) ) stop(
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
