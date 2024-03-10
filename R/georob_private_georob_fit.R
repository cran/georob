
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
  ## 2023-12-20 AP added on.exit(options(...))
  ## 2024-01-21 AP more efficient calculation of lag.vectors for anisotropic variograms

  ##  ToDos:

#### -- preparations

  ##  main body of georob.fit

  d2r <- pi / 180.

  ##  define rho-function and derivatives (suppress temporarily warnings issued by gamma())

  old.op <- options( warn = -1 )
  on.exit( options( old.op ) )
  rho.psi.etc <- f.psi.function( x = psi.func, tp = tuning.psi )

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
      if( aniso["f2"] < 0. ||  aniso["f2"] > 1. ) stop(
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
      lag.vectors <- apply(
        coordinates, 2,
        function( x ){
          nx <- length( x )
          tmp <- matrix( rep( x, nx ), ncol = nx)
          sel <- lower.tri( tmp )
          tmp[sel] - (t( tmp ))[sel]
        }
      )
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

  ## stop PSOCK and snowfall clusters

  #   f.stop.cluster()

  #   if( length( lik.item[["defaultCluster"]] ) > 0L ){
  #     cl <- lik.item[["defaultCluster"]]
  #
  #     junk <- parLapply( cl, 1L:length(cl), function( i ) sfStop() )
  #     junk <- stopCluster( cl )
  #     sfStop()
  #   }

  #   if( sfIsRunning() ){
  #     sfStop()
  #   }
  #
  #   if( file.exists( "SOCKcluster.RData" ) ){
  #     file.remove( "SOCKcluster.RData" )
  #   }

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
