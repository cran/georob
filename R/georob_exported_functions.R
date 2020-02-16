##  ############################################################################
### georob

georob <-
  function(
    formula, data, subset, weights, na.action,
    model = TRUE, x = FALSE, y = FALSE,
    contrasts = NULL, offset,
    locations,
    variogram.model = c( "RMexp", "RMaskey", "RMbessel", "RMcauchy",
      "RMcircular", "RMcubic", "RMdagum", "RMdampedcos", "RMdewijsian", "RMfbm",
      "RMgauss", "RMgencauchy", "RMgenfbm", "RMgengneiting", "RMgneiting", "RMlgd",
      "RMmatern", "RMpenta", "RMqexp", "RMspheric", "RMstable",
      "RMwave", "RMwhittle"
    ),
    param, fit.param = default.fit.param()[names(param)],
    aniso = default.aniso(), fit.aniso = default.fit.aniso(),
    variogram.object = NULL,
    tuning.psi = 2., control = control.georob(), verbose = 0,
    ...
  )
{

  ## wrapper function to georob.fit with "standard" interface for
  ## statistical modelling

  ## 2012-04-21 AP
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-07 AP correction of error for constant trend
  ## 2012-05-28 AP handle missing names of coefficients after calling update
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-05-23 AP correct handling of missing observations and to construct model.frame
  ## 2013-06-03 AP handling design matrices with rank < ncol(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-02 AP new transformation of rotation angles
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2013-09-06 AP exclusive use of nleqslv for solving estimating equations
  ## 2014-02-18 AP correcting error when fitting models with offset
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2014-05-22 AP correcting error when selecting initial.param
  ## 2014-06-02 AP partial matching of names of variogram parameters
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2014-08-26 AP changes to return ml.method if fitted object
  ## 2015-06-30 AP function and arguments renamed
  ## 2015-07-17 AP change of control of computation of hessian for Gaussian (RE)ML,
  ##               changes for singular design matrices, nlminb optimizer added
  ## 2015-08-19 AP control about error families for computing covariances added
  ## 2015-08-28 AP computation of hessian suppressed; correction of error when using georob.object;
  ##               control arguments hessian, initial.param, initial.fixef newly organized
  ## 2015-11-25 AP new way to control which variogram parameters are fitted
  ## 2016-07-15 AP allowing use of lmrob for computing initial fixed effects for rank-deficient model matrix
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-08-08 AP changes for nested variogram models
  ## 2017-05-07 AP correcting error for interface to rq.fit
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

#### -- check arguments

  ## check whether all mandatory arguments have been provided

  if( missing( formula ) || missing( locations ) ||
    ( is.null( variogram.object ) && missing( param )) ) stop(
    "some mandatory arguments are missing"
  )

  ## check validity of arguments

  if(!missing(data))    stopifnot(is.data.frame(as.data.frame(data)))
  if(!missing(subset))  stopifnot(is.null(subset)  || is.logical(subset) || is.character(subset) || is.numeric(subset))
  if(!missing(weights)) stopifnot(is.null(weights) || is.numeric(weights))
  if(!missing(offset))  stopifnot(is.null(offset)  || is.numeric(offset))

  stopifnot(identical(length(model), 1L) && is.logical(model))
  stopifnot(identical(length(x), 1L)     && is.logical(x))
  stopifnot(identical(length(y), 1L)     && is.logical(y))

  stopifnot(is.logical(fit.param))
  stopifnot(is.logical(fit.aniso))

  stopifnot(identical(length(tuning.psi), 1L) && is.numeric(tuning.psi) && tuning.psi > 0)
  stopifnot(identical(length(verbose), 1L)    && is.numeric(verbose)    && verbose >= 0)

  stopifnot(is.numeric(param))
  stopifnot(is.numeric(aniso))

  stopifnot(is.list(control))

  stopifnot(is.null(contrasts)        || is.list(contrasts))
  stopifnot(is.null(variogram.object) || is.list(variogram.object))

  if( identical( control[["psi.func"]], "t.dist" ) && tuning.psi <= 1. )
    stop( "'tuning.psi' must be greater than 1 for t-dist psi-function" )


#### -- check/setup variogram object

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


#### -- model frame etc.

  ## get model frame, response vector, weights, offset and design
  ## matrix (cf.  lm, lmrob)

  ret.x <- x
  ret.y <- y

#   ## vector with row number of included observations
#
#   in.subset <- 1L:NROW( data )
#   if( !missing( subset ) ) in.subset <- in.subset[subset]

  ## build combined formula for fixed effects and locations

  extended.formula <- update(
    formula,
    paste(
      paste( as.character( formula )[c(2L, 1L, 3L)], collapse = " " ),
      as.character( locations )[2L], sep = " + "
    )
  )

  ## setting-up model frame

  cl <- match.call()
  mf <- match.call( expand.dots = FALSE )
  m <- match(
    c( "formula", "data", "subset", "weights", "na.action", "offset"),
    names(mf), 0L
  )
  mf <- mf[c(1L, m)]
  mf[["formula"]] <- extended.formula
  mf[["drop.unused.levels"]] <- TRUE
  mf[[1L]] <- as.name( "model.frame" )

  mf <- eval( mf, parent.frame() )

  ## setting-up terms objects

  mt     <- terms( formula )
  mt.loc <- terms( locations )

  ## eliminate intercept from mt.loc

  attr( mt.loc, "intercept" ) <- 0L

#   ## ... and assign fixed effects terms object as attribute to model.frame
#
#   attr( mf, "terms" ) <- mt

  ## check whether 'empty' models have been entered

  if( is.empty.model( mt ) )
    stop( "an 'empty' fixed effects model has been specified" )
  if( is.empty.model( mt.loc ) )
    stop( "an 'empty' locations model has been specified" )

  ## check whether fixed effects model includes an intercept if an
  ## intrinsic variogram model is used

  if( identical( attr( mt, "intercept" ), 0L ) &&
    variogram.model %in% control[["irf.model"]] )
  stop(
    "the fixed effects model must include an intercept ",
    "if an unbounded variogram model is used"
  )

  ## extract fixed effects response variable, weights, offset and design matrix

  y <- model.response( mf, "numeric" )

  w <- as.vector( model.weights( mf ) )
  if( !is.null(w) )
    stop( "weights are not yet implemented for this estimator" )

  offset <- as.vector( model.offset(mf) )
  if( !is.null(offset) ) {
    if( length( offset ) != NROW(y) )
      stop( gettextf(
        "number of offsets is %d, should equal %d (number of observations)",
        length(offset), NROW(y) ), domain = NA )
  }

  x <- model.matrix( mt, mf, contrasts )

  ## check if optionally provided bhat has correct length

  if( !is.null( control[["bhat"]] ) && length( y ) != length( control[["bhat"]] ) ) stop(
    "lengths of response vector and 'bhat' do not match"
  )

  ## adjust initial.param if all variogram parameters are fixed

  if( all( sapply(
        variogram.object,
        function( x ) !any( c( x[["fit.param"]], x[["fit.aniso"]] ) )
      ))) control[["initial.param"]] <- FALSE

  ## adjust choice for initial.fixef to compute regression coefficients

  if( tuning.psi < control[["tuning.psi.nr"]] ){
    if( control[["initial.param"]] ) control[["initial.fixef"]] <- "lmrob"
  } else {
    control[["initial.param"]] <- FALSE
  }

  ## check whether design matrix has full column rank

  col.rank.XX <- list( deficient = FALSE, rank = NCOL( x ) )

  sv <- svd( crossprod( x ) )[["d"]]
  min.max.sv <- range( sv )
  condnum <- min.max.sv[1L] / min.max.sv[2L]

  if( condnum <= control[["min.condnum"]] ){
    col.rank.XX[["deficient"]] <- TRUE
    col.rank.XX[["col.rank"]] <- sum( sv / min.max.sv[2L] > control[["min.condnum"]] )
    if( control[["initial.fixef"]] == "rq" ){
      cat(
        "design matrix has not full column rank (condition number of X^T X: ",
        signif( condnum, 2L ), ")\ninitial values of fixed effects coefficients are computed by 'lmrob' instead of 'rq'\n\n"
      )
      control[["initial.fixef"]] <- "lmrob"
      warning(
        "design matrix has not full column rank (condition number of X^T X: ",
        signif( condnum, 2L ), ")\ninitial values of fixed effects coefficients are computed by 'lmrob' instead of 'rq'"
      )
    }
  }

  ## subtract offset

  yy <- y
  if( !is.null( offset ) ) yy <- yy - offset


#### -- initial values of fixed effects parameters

  ## compute initial guess of fixed effects parameters (betahat)

  r.initial.fit <- switch(
    control[["initial.fixef"]],
    rq = {

      Rho <- function( u, tau) u * (tau - (u < 0))
      tau <- control[["rq"]][["tau"]]
      process <- (tau < 0 || tau > 1)

      if(tau == 0) tau <- control[["req"]][["eps"]]
      if(tau == 1) tau <- 1. - control[["req"]][["eps"]]

      fit <- switch(
        control[["rq"]][["method"]],
        br = rq.fit(
          x = x, y = yy, tau = tau, method = control[["rq"]][["method"]],
          alpha = control[["rq"]][["alpha"]], ci = control[["rq"]][["ci"]],
          iid = control[["rq"]][["iid"]], interp = control[["rq"]][["interp"]],
          tcrit = control[["rq"]][["tcrit"]]
        ),
        fnb = rq.fit(
          x = x, y = yy, tau = tau, method = control[["rq"]][["method"]],
          beta = control[["rq"]][["beta"]], eps = control[["rq"]][["eps"]]
        ),
        pfn = rq.fit(
          x = x, y = yy, tau = tau, method = control[["rq"]][["method"]],
          Mm.factor = control[["rq"]][["Mm.factor"]],
          max.bad.fixup = control[["rq"]][["max.bad.fixup"]]
        )
      )

      if (process){
        rho <- list(x = fit$sol[1, ], y = fit$sol[3, ])
      } else {
        if (length(dim(fit$residuals)))
        dimnames(fit$residuals) <- list(dimnames(x)[[1]],
          NULL)
        rho <- sum(Rho(fit$residuals, tau))
      }
      if( control[["rq"]][["method"]] == "lasso" ){
        class(fit) <- c("lassorq", "rq")
      } else if( control[["rq"]][["method"]] == "scad"){
        class(fit) <- c("scadrq", "rq")
      } else {
        class(fit) <- ifelse(process, "rq.process", "rq")
      }
      fit[["na.action"]] <- attr( mf, "na.action" )
      fit[["formula"]] <- formula
      fit[["terms"]] <- mt
      fit[["xlevels"]] <- .getXlevels(mt, mf)
      fit[["call"]] <- cl
      fit[["tau"]] <- tau
      fit[["weights"]] <- w
      fit[["residuals"]] <- drop( fit[["residuals"]] )
      fit[["rho"]] <- rho
      fit[["method"]] <- control[["rq"]][["method"]]
      fit[["fitted.values"]] <- drop( fit[["fitted.values"]] )
      attr(fit, "na.message") <- attr( m, "na.message" )
      if( model ) fit[["model"]] <- mf
      fit

    },
    lmrob = {

      fit <- lmrob.fit( x, yy, control = control[["lmrob"]] )
      fit[["na.action"]] <- attr(mf, "na.action")
      fit[["offset"]] <- offset
      fit[["contrasts"]] <- attr(x, "contrasts")
      fit[["xlevels"]] <- .getXlevels(mt, mf)
      fit[["call"]] <- cl
      fit[["terms"]] <- mt
      if( control[["lmrob"]][["compute.rd"]] && !is.null(x) )
      fit[["MD"]] <- robMD( x, attr( mt, "intercept" ) )
      if( !is.null( offset ) ) fit[["fitted.values"]] + offset
      fit

    },
    lm = {

      fit <- if( is.null(w) ){
        lm.fit(x, y, offset = offset, singular.ok = TRUE )
      } else {
        lm.wfit(x, y, w, offset = offset, singular.ok = TRUE )
      }
      class(fit) <- c(if (is.matrix(y)) "mlm", "lm")
      fit[["na.action"]] <- attr(mf, "na.action")
      fit[["offset"]] <- offset
      fit[["contrasts"]] <- attr(x, "contrasts")
      fit[["xlevels"]] <- .getXlevels(mt, mf)
      fit[["call"]] <- cl
      fit[["terms"]] <- mt
      if (model) fit[["model"]] <- mf
      if (ret.x) fit[["x"]] <- x
      if (ret.y) fit[["y"]] <- y
      fit[["qr"]] <- NULL
      fit

    }
  )


#### -- prepare coordinates

  ## compute coordinates of locations and distance object

  locations.coords <- model.matrix( mt.loc, mf )

  if(
    !( missing( aniso ) || missing( fit.aniso ) ) &&
    ( NCOL( locations.coords ) < 2L || NCOL( locations.coords ) > 3L )
  ) stop(
    "anisotropic variogram models are implemented only for 2 or 3 dimensions"
  )

  names( yy ) <- rownames( mf )

  ##  check whether argument "object."  has been provided in call (e.g. by
  ##  update ) and extract'locations' exists in workspace

  extras <- match.call( expand.dots = FALSE )$...
  georob.object <- extras[names(extras) %in% "object."]
  if(
    length( georob.object ) &&
    exists( as.character( georob.object ), envir = parent.frame() ) #&&
    #     !control[["initial.param"]]
  ){
    if( verbose > 4. ) cat(
      "\n    georob: using some components of 'object.'\n"
    )
    georob.object <- eval(
      georob.object[[1L]], parent.frame()
    )[c( "variogram.object", "locations.objects", "Valphaxi.objects" )]
    georob.object[["Valphaxi.objects"]] <- c(
      list(error=FALSE),
      georob.object[["Valphaxi.objects"]]
    )
  } else {
    georob.object <- NULL
  }


#### -- initial values of variogram parameters

  ## compute initial values of variogram and anisotropy parameters

  if( tuning.psi < control[["tuning.psi.nr"]] && control[["initial.param"]] ){

    if( verbose > 0. ) cat( "\ncomputing robust initial parameter estimates ...\n" )

    t.sel <- r.initial.fit[["rweights"]] > control[["min.rweight"]]

    if( verbose > 0. ) cat(
      "\ndiscarding", sum( !t.sel ), "of", length( t.sel ),
      "observations for computing initial estimates of variogram\nand anisotropy parameter by gaussian reml\n"
    )

    ## collect.items for initial object

    initial.objects <- list(
      x = as.matrix( x[t.sel, ] ),
      y = yy[t.sel],
      betahat = coef( r.initial.fit ),
      bhat = if( is.null( control[["bhat"]] ) ){
        rep( 0., length( yy ) )[t.sel]
      } else {
        control[["bhat"]][t.sel]
      },
      initial.fit = r.initial.fit,
      locations.objects = list(
        locations = locations,
        coordinates = locations.coords[t.sel, ]
      )
    )

    ## estimate model parameters with pruned data set

    t.georob <- georob.fit(
      initial.objects = initial.objects,
      variogram.object = variogram.object,
      param.tf = control[["param.tf"]],
      fwd.tf = control[["fwd.tf"]],
      deriv.fwd.tf = control[["deriv.fwd.tf"]],
      bwd.tf = control[["bwd.tf"]],
      georob.object = georob.object,
      safe.param = control[["safe.param"]],
      tuning.psi = control[["tuning.psi.nr"]],
      error.family.estimation = control[["error.family.estimation"]],
      error.family.cov.effects = control[["error.family.cov.effects"]],
      error.family.cov.residuals = control[["error.family.cov.residuals"]],
      cov.bhat = control[["cov.bhat"]], full.cov.bhat = control[["full.cov.bhat"]],
      cov.betahat = control[["cov.betahat"]],
      cov.bhat.betahat = control[["cov.bhat.betahat"]],
      cov.delta.bhat = control[["cov.delta.bhat"]],
      full.cov.delta.bhat = control[["full.cov.delta.bhat"]],
      cov.delta.bhat.betahat = control[["cov.delta.bhat.betahat"]],
      cov.ehat = control[["cov.ehat"]], full.cov.ehat = control[["full.cov.ehat"]],
      cov.ehat.p.bhat = control[["cov.ehat.p.bhat"]], full.cov.ehat.p.bhat = control[["full.cov.ehat.p.bhat"]],
      aux.cov.pred.target = control[["aux.cov.pred.target"]],
      min.condnum = control[["min.condnum"]],
      col.rank.XX = col.rank.XX,
      psi.func = control[["psi.func"]],
      tuning.psi.nr = tuning.psi,
      ml.method = control[["ml.method"]],
      maximizer = control[["maximizer"]],
      reparam = control[["reparam"]],
      irwls.initial = control[["irwls.initial"]],
      irwls.maxiter = control[["irwls.maxiter"]],
      irwls.ftol = control[["irwls.ftol"]],
      force.gradient = control[["force.gradient"]],
      zero.dist = control[["zero.dist"]],
      control.nleqslv =  control[["nleqslv"]],
      control.optim = control[["optim"]],
      control.nlminb = control[["nlminb"]],
      hessian = FALSE,
      control.pcmp = control[["pcmp"]],
      verbose = verbose
    )

    variogram.object <- lapply(
      1:length(variogram.object),
      function( i, x, fvo ){
        x   <- x[[i]]
        fvo <- fvo[[i]]
        x[["param"]] <- fvo[["param"]][names(x[["param"]])]
        x[["aniso"]] <- fvo[["aniso"]][names(x[["aniso"]])]
        x
      }, x = variogram.object, fvo = t.georob[["variogram.object"]]
    )

  }

  ## collect.items for initial object

  initial.objects <- list(
    x = as.matrix( x ),
    y = yy,
    betahat = coef( r.initial.fit ),
    bhat = if( is.null( control[["bhat"]] ) ){
      rep( 0., length( yy ) )
    } else {
      control[["bhat"]]
    },
    initial.fit = r.initial.fit,
    locations.objects = list(
      locations = locations,
      coordinates = locations.coords
    )
  )


#### -- final estimate of model parameters

  ## estimate model parameters

  if( verbose > 0. ) cat( "\ncomputing final parameter estimates ...\n" )

  r.georob <- georob.fit(
    initial.objects = initial.objects,
    variogram.object = variogram.object,
    param.tf = control[["param.tf"]],
    fwd.tf = control[["fwd.tf"]],
    deriv.fwd.tf = control[["deriv.fwd.tf"]],
    bwd.tf = control[["bwd.tf"]],
    georob.object = georob.object,
    safe.param = control[["safe.param"]],
    tuning.psi = tuning.psi,
    error.family.estimation = control[["error.family.estimation"]],
    error.family.cov.effects = control[["error.family.cov.effects"]],
    error.family.cov.residuals = control[["error.family.cov.residuals"]],
    cov.bhat = control[["cov.bhat"]], full.cov.bhat = control[["full.cov.bhat"]],
    cov.betahat = control[["cov.betahat"]],
    cov.bhat.betahat = control[["cov.bhat.betahat"]],
    cov.delta.bhat = control[["cov.delta.bhat"]],
    full.cov.delta.bhat = control[["full.cov.delta.bhat"]],
    cov.delta.bhat.betahat = control[["cov.delta.bhat.betahat"]],
    cov.ehat = control[["cov.ehat"]], full.cov.ehat = control[["full.cov.ehat"]],
    cov.ehat.p.bhat = control[["cov.ehat.p.bhat"]], full.cov.ehat.p.bhat = control[["full.cov.ehat.p.bhat"]],
    aux.cov.pred.target = control[["aux.cov.pred.target"]],
    min.condnum = control[["min.condnum"]],
    col.rank.XX = col.rank.XX,
    psi.func = control[["psi.func"]],
    tuning.psi.nr = control[["tuning.psi.nr"]],
    ml.method = control[["ml.method"]],
    maximizer = control[["maximizer"]],
    reparam = control[["reparam"]],
    irwls.initial = control[["irwls.initial"]],
    irwls.maxiter = control[["irwls.maxiter"]],
    irwls.ftol = control[["irwls.ftol"]],
    force.gradient = control[["force.gradient"]],
    zero.dist = control[["zero.dist"]],
    control.nleqslv =  control[["nleqslv"]],
    control.optim = control[["optim"]],
    control.nlminb = control[["nlminb"]],
    hessian = control[["hessian"]],
    control.pcmp = control[["pcmp"]],
    verbose = verbose
  )


#### -- prepare output

  ## add offset to fitted values

  if( !is.null( offset ) )
    r.georob[["fitted.values"]] <- r.georob[["fitted.values"]] + offset

  ##

  r.georob[["control"]] <- control

  ## add remaining items to output

  if( control[["lmrob"]][["compute.rd"]] && !is.null( x ) )
    r.georob[["MD"]] <- robMD( x, attr(mt, "intercept") )

  if( model ) r.georob[["model"]] <- mf
  if( ret.x ) r.georob[["x"]] <- x
  if( ret.y ) r.georob[["y"]] <- y

  r.georob[["df.residual"]] <- length(yy) - col.rank.XX[["rank"]]
  r.georob[["na.action"]] <- attr(mf, "na.action")
  r.georob[["offset"]] <- offset
  r.georob[["contrasts"]] <- attr(x, "contrasts")
  r.georob[["xlevels"]] <- .getXlevels(mt, mf)
  r.georob[["rank"]] <- col.rank.XX[["rank"]]
  r.georob[["call"]] <- cl
  r.georob[["terms"]] <- mt

  ## set missing names of coefficients (bug of update)

  if( identical(length( r.georob[["coefficients"]] ), 1L) &&
    is.null( names( r.georob[["coefficients"]] ) ) ){
    names( r.georob[["coefficients"]] ) <- "(Intercept)"
  }


  class( r.georob ) <- c( "georob" )

  invisible( r.georob )

}


##  ##############################################################################

pmm <-
  function(
    A, B, control = control.pcmp()
  )
{

  ## function for parallelized matrix multiplication inspired by function
  ## parMM{snow}

  ## 2014-06-25 A. Papritz
  ## 2015-03-13 AP small changes in f.aux
  ## 2016-07-20 AP arguments renamed, new control of recursive parallelization
  ## 2018-01-05 AP improved memory management in parallel computations

  ## auxiliary function

  f.aux <- function( i ){

    ## s, e, A, B are taken from parent environment

    A %*% B[ , s[i]:e[i], drop = FALSE]

  }

  ## determine number of cores

  ncores <- control[["pmm.ncores"]]
  if( !control[["allow.recursive"]] ) ncores <- 1L

  ## determine columns indices of matrix blocks

  k <- control[["f"]] * ncores
  n <- NCOL(B)
  dn <- floor( n / k )
  s <- ( (0L:(k-1L)) * dn ) + 1L
  e <- (1L:k) * dn
  e[k] <- n

  ## set default value for control of forking if missing (required for backward compatibility)

  if( is.null( control[["fork"]] ) ){
    control[["fork"]] <- !identical( .Platform[["OS.type"]], "windows" )
  }

  if( ncores > 1L ){

    if( !control[["fork"]] ){

      ## create a SNOW cluster on windows OS

      if( !sfIsRunning() ){
        options( error = f.stop.cluster )
        junk <- sfInit( parallel = TRUE, cpus = ncores )
      }

      junk <- sfExport( "s", "e", "A", "B" )

      res <- sfLapply( 1L:k, f.aux )

      if( control[["sfstop"]] ){
        junk <- sfStop()
        options( error = NULL )
      }

    } else {

      res <- mclapply( 1L:k, f.aux, mc.cores = ncores )

    }

    junk <- gc()

    matrix( unlist(res), nrow = NROW( A ) )

  } else A %*% B

}

#  ##############################################################################

control.georob <-
  function(
    ml.method = c( "REML", "ML" ),
    reparam = TRUE,
    maximizer = c( "nlminb", "optim" ),
    initial.param = TRUE,
    initial.fixef = c("lmrob", "rq", "lm"),
    bhat = NULL,
    min.rweight = 0.25,
    param.tf = param.transf(),
    fwd.tf = fwd.transf(),
    deriv.fwd.tf = dfwd.transf(),
    bwd.tf = bwd.transf(),
    psi.func = c( "logistic", "t.dist", "huber" ),
    irwls.maxiter = 50, irwls.ftol = 1.e-5,
    force.gradient = FALSE,
    min.condnum = 1.e-12,
    zero.dist = sqrt( .Machine[["double.eps"]] ),
    error.family.estimation    = c( "gaussian", "long.tailed" ),
    error.family.cov.effects   = c( "gaussian", "long.tailed" ),
    error.family.cov.residuals = c( "gaussian", "long.tailed" ),
    cov.bhat = TRUE, full.cov.bhat = FALSE,
    cov.betahat = TRUE,
    cov.delta.bhat = TRUE, full.cov.delta.bhat = TRUE,
    cov.delta.bhat.betahat = TRUE,
    cov.ehat = TRUE, full.cov.ehat = FALSE,
    cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
    hessian = TRUE,
    rq = control.rq(),
    lmrob = lmrob.control(),
    nleqslv = control.nleqslv(),
    optim = control.optim(),
    nlminb = control.nlminb(),
    pcmp = control.pcmp(),
    ...
  )
{

  ## auxiliary function to set meaningful default values for georob

  ## 2012-04-21 A. Papritz
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-06-12 AP changes in stored items of Valphaxi object
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2015-03-10 AP extended variogram parameter transformations
  ## 2015-06-30 AP function and arguments renamed
  ## 2015-07-17 AP nlminb optimizer added
  ## 2015-08-19 AP control about error families for computing covariances added
  ## 2015-08-28 AP control arguments hessian, initial.param, initial.fixef newly organized
  ## 2016-07-15 AP some arguments hard coded
  ## 2016-08-24 AP new default value for error.family.cov.residuals
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## set values of fixed control arguments

  irwls.initial <- TRUE
  tuning.psi.nr <- 1000.
  cov.bhat.betahat <- FALSE
  aux.cov.pred.target <- FALSE
  safe.param <- 1.e12
  sepstr <- ".__...__."

  ## match arguments

  ml.method     <- match.arg( ml.method )
  maximizer     <- match.arg( maximizer )
  initial.fixef <- match.arg( initial.fixef )
  psi.func      <- match.arg( psi.func )

  error.family.estimation    <- match.arg( error.family.estimation )
  error.family.cov.effects   <- match.arg( error.family.cov.effects )
  error.family.cov.residuals <- match.arg( error.family.cov.residuals )

  ## check arguments

  stopifnot(identical(length(reparam), 1L)                && is.logical(reparam))
  stopifnot(identical(length(initial.param), 1L)          && is.logical(initial.param))
  stopifnot(identical(length(force.gradient), 1L)         && is.logical(force.gradient))
  stopifnot(identical(length(cov.bhat), 1L)               && is.logical(cov.bhat))
  stopifnot(identical(length(full.cov.bhat), 1L)          && is.logical(full.cov.bhat))
  stopifnot(identical(length(cov.betahat), 1L)            && is.logical(cov.betahat))
  stopifnot(identical(length(cov.delta.bhat), 1L)         && is.logical(cov.delta.bhat))
  stopifnot(identical(length(full.cov.delta.bhat), 1L)    && is.logical(full.cov.delta.bhat))
  stopifnot(identical(length(cov.delta.bhat.betahat), 1L) && is.logical(cov.delta.bhat.betahat))
  stopifnot(identical(length(cov.ehat), 1L)               && is.logical(cov.ehat))
  stopifnot(identical(length(full.cov.ehat), 1L)          && is.logical(full.cov.ehat))
  stopifnot(identical(length(cov.ehat.p.bhat), 1L)        && is.logical(cov.ehat.p.bhat))
  stopifnot(identical(length(full.cov.ehat.p.bhat), 1L)   && is.logical(full.cov.ehat.p.bhat))
  stopifnot(identical(length(hessian), 1L)                && is.logical(hessian))

  stopifnot(identical(length(min.rweight), 1L)   && is.numeric(min.rweight)   && min.rweight >= 0)
  stopifnot(identical(length(irwls.maxiter), 1L) && is.numeric(irwls.maxiter) && irwls.maxiter >= 1)
  stopifnot(identical(length(irwls.ftol), 1L)    && is.numeric(irwls.ftol)    && irwls.ftol > 0)
  stopifnot(identical(length(min.condnum), 1L)   && is.numeric(min.condnum)   && min.condnum >= 0)
  stopifnot(identical(length(zero.dist), 1L)     && is.numeric(zero.dist)     && zero.dist > 0)

  stopifnot(is.null(bhat) || is.numeric(bhat))

  stopifnot(is.list(param.tf))
  stopifnot(is.list(fwd.tf))
  stopifnot(is.list(deriv.fwd.tf))
  stopifnot(is.list(bwd.tf))

  stopifnot(is.list(rq))
  stopifnot(is.list(lmrob))
  stopifnot(is.list(nleqslv))
  stopifnot(is.list(optim))
  stopifnot(is.list(nlminb))
  stopifnot(is.list(pcmp))

  if(
    !( all( unlist( param.tf ) %in% names( fwd.tf ) ) &&
       all( unlist( param.tf ) %in% names( deriv.fwd.tf ) ) &&
       all( unlist( param.tf ) %in% names( bwd.tf ) )
    )
  ) stop(
    "undefined transformation of variogram parameters; extend respective function definitions"
  )

  if( !irwls.initial && irwls.ftol >= 1.e-6 ) warning(
    "'irwls.initial == FALSE' and large 'ftol' may create problems for root finding"
  )

  list(
    ml.method = ml.method, reparam = reparam,
    maximizer = maximizer,
    initial.param = initial.param,
    initial.fixef = initial.fixef,
    bhat = bhat,
    min.rweight = min.rweight,
    param.tf = param.tf, fwd.tf = fwd.tf, deriv.fwd.tf = deriv.fwd.tf, bwd.tf = bwd.tf,
    safe.param = safe.param, sepstr = sepstr,
    psi.func = psi.func,
    tuning.psi.nr = tuning.psi.nr,
    irwls.initial = irwls.initial, irwls.maxiter = irwls.maxiter, irwls.ftol = irwls.ftol,
    force.gradient = force.gradient,
    min.condnum = min.condnum,
    zero.dist = zero.dist,
    error.family.estimation     = error.family.estimation,
    error.family.cov.effects    = error.family.cov.effects,
    error.family.cov.residuals  = error.family.cov.residuals,
    cov.bhat = cov.bhat, full.cov.bhat = full.cov.bhat,
    cov.betahat = cov.betahat,
    cov.bhat.betahat = cov.bhat.betahat,
    cov.delta.bhat = cov.delta.bhat, full.cov.delta.bhat = full.cov.delta.bhat,
    cov.delta.bhat.betahat = cov.delta.bhat.betahat,
    cov.ehat = cov.ehat, full.cov.ehat = full.cov.ehat,
    cov.ehat.p.bhat = cov.ehat.p.bhat, full.cov.ehat.p.bhat = full.cov.ehat.p.bhat,
    aux.cov.pred.target = aux.cov.pred.target,
    hessian = hessian,
    irf.models = c( "RMdewijsian", "RMfbm", "RMgenfbm" ),
    rq = rq, lmrob = lmrob, nleqslv = nleqslv,
    optim = optim,
    nlminb = nlminb,
    pcmp = pcmp
  )

}

## ======================================================================
param.transf <-
  function(
    variance = "log", snugget = "log", nugget = "log", scale = "log",
    alpha = c(
      RMaskey = "log", RMdewijsian = "logit2", RMfbm = "logit2", RMgencauchy = "logit2",
      RMgenfbm = "logit2", RMlgd = "identity", RMqexp = "logit1", RMstable = "logit2"
    ),
    beta = c( RMdagum = "logit1", RMgencauchy = "log", RMlgd = "log" ),
    delta = "logit1",
    gamma = c( RMcauchy = "log", RMdagum = "logit1" ),
    kappa = "logit3", lambda = "log",
    mu = "log",
    nu = "log",
    f1 = "log", f2  ="log", omega = "identity", phi = "identity", zeta = "identity"
  )
{

  ## function sets meaningful defaults for transformation of variogram
  ## parameters

  ## 2013-07-02 A. Papritz
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2015-03-10 AP extended transformation
  ## 2015-04-07 AP changes for fitting anisotropic variograms

  list(
    variance = variance, snugget = snugget, nugget = nugget, scale = scale,
    alpha = alpha,
    beta = beta,
    delta = delta,
    gamma = gamma,
    kappa = kappa,
    lambda = lambda,
    mu = mu,
    nu = nu,
    f1 = f1, f2 = f2, omega = omega, phi = phi, zeta = zeta
  )

}

## ======================================================================
fwd.transf <-
  function(
    ...
  )
{

  ## definition of forward transformation of variogram parameters

  ## 2013-07-02 A. Papritz
  ## 2015-03-10 AP extended variogram parameter transformations
  ## 2015-04-07 AP changes for fitting anisotropic variograms

  list(
    log = function(x) log(x),
    logit1 = function(x) log( x / (1. - x) ),
    logit2 = function(x) log( x / (2. - x) ),
    logit3 = function(x) log( (x - 1.) / (3. - x) ),
    identity = function(x) x, ...
  )
}

## ======================================================================
dfwd.transf<-
  function(
    ...
  )
{

  ## definition of first derivative of forward transformation of variogram
  ## parameters

  ## 2013-07-02 A. Papritz
  ## 2015-03-10 AP extended variogram parameter transformations
  ## 2015-04-07 AP changes for fitting anisotropic variograms

  list(
    log = function(x) 1./x,
    logit1 = function(x) 1. / (x - x^2),
    logit2 = function(x) 2. / (2.*x - x^2),
    logit3 = function(x) 2. / (4.*x - 3. - x^2),
    identity = function(x) rep(1., length(x)), ...
  )

}

## ======================================================================
bwd.transf <-
function(
  ...
)
{

  ## definition of backward transformation of variogram parameters

  ## 2013-07-02 A. Papritz
  ## 2015-03-10 AP extended variogram parameter transformations
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  ## 2016-08-03 AP corrections for logitx if argument is +/-Inf
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  list(
    log = function(x) exp(x),
    logit1 = function(x){
      stopifnot(identical(length(x), 1L) && is.numeric(x))
      if( !is.finite(x) ){
        if(sign(x) < 0.) 0. else 1.
      } else exp(x) / (1. + exp(x))
    },
    logit2 = function(x){
      stopifnot(identical(length(x), 1L) && is.numeric(x))
      if( !is.finite(x) ){
        if(sign(x) < 0.) 0. else 2.
      } else 2. * exp(x) / (1. + exp(x))
    },
    logit3 = function(x){
      stopifnot(identical(length(x), 1L) && is.numeric(x))
      if( !is.finite(x) ){
        if(sign(x) < 0.) 1. else 3.
      } else (3. * exp(x) + 1.) / (1. + exp(x))
    },
    identity = function(x) x, ...
  )
}

## ======================================================================
control.rq <-
  function(
    ## specific arguments for rq: tau = 0.5, method = "br"
    tau = 0.5, rq.method = c("br", "fnb", "pfn"),
    ## specific arguments for rq.fit.br: tau = 0.5, alpha = 0.1, ci = FALSE, iid = TRUE, interp = TRUE, tcrit = TRUE
    rq.alpha = 0.1, ci = FALSE, iid = TRUE, interp = TRUE, tcrit = TRUE,
    ## specific arguments for rq.fit.fnb: tau = 0.5, rhs = (1 - tau) * apply(x, 2, sum), beta = 0.99995, eps = 1e-06
    rq.beta = 0.99995, eps = 1.e-06,
    ## specific arguments for rq.fit.pfn: tau = 0.5, Mm.factor = 0.8, max.bad.fixup = 3, eps = 1e-06
    Mm.factor = 0.8, max.bad.fixup = 3,
    ...
  )
{

  ## function sets meaningful defaults for selected arguments of function
  ## rq{quantreg}

  ## 2012-12-14 A. Papritz
  ## 2014-07-29 AP
  ## 2015-06-30 AP function and arguments renamed
  ## 2017-05-09 AP small changes
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## match arguments

  rq.method = match.arg(rq.method)

  ## check values of arguments

  stopifnot(identical(length(ci), 1L)     && is.logical(ci))
  stopifnot(identical(length(iid), 1L)    && is.logical(iid))
  stopifnot(identical(length(interp), 1L) && is.logical(interp))
  stopifnot(identical(length(tcrit), 1L)  && is.logical(tcrit))

  stopifnot(identical(length(tau), 1L)           && is.numeric(tau)           && tau > 0 && tau < 1)
  stopifnot(identical(length(rq.alpha), 1L)      && is.numeric(rq.alpha)      && rq.alpha > 0 && rq.alpha < 1)
  stopifnot(identical(length(rq.beta), 1L)       && is.numeric(rq.beta))
  stopifnot(identical(length(eps), 1L)           && is.numeric(eps))
  stopifnot(identical(length(Mm.factor), 1L)     && is.numeric(Mm.factor)     && Mm.factor > 0)
  stopifnot(identical(length(max.bad.fixup), 1L) && is.numeric(max.bad.fixup) && max.bad.fixup > 0)

  list(
    tau = tau,  method = rq.method,
    alpha = rq.alpha, ci = ci, iid = iid, interp = interp, tcrit = tcrit,
    beta = rq.beta, eps = eps,
    Mm.factor = Mm.factor, max.bad.fixup = max.bad.fixup
  )
}


## ======================================================================
control.nleqslv <-
  function(
    method = c( "Broyden", "Newton"),
    global = c( "dbldog", "pwldog", "qline", "gline", "none" ),
    xscalm = c( "fixed", "auto" ),
    control = list( ftol = 1.e-4 ),
    ...
  )
{

  ## function sets meaningful defaults for selected arguments of function
  ## nleqslv{nleqslv}

  ## 2013-07-12 A. Papritz
  ## 2014-07-29 AP
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check arguments

  stopifnot(is.list(control))

  list(
    method = match.arg( method ),
    global = match.arg( global ),
    xscalm = match.arg( xscalm ),
    control = control
  )
}

## ======================================================================
control.optim <-
  function(
    method = c( "BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN", "Brent" ),
    lower = -Inf, upper = Inf,
    control = list(reltol = 1.e-5),
    ...
  )
{

  ## function sets meaningful defaults for selected arguments of function optim

  ## 2012-12-14 A. Papritz
  ## 2014-07-29 AP
  ## 2015-06-30 AP function and arguments renamed
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check arguments

  stopifnot(is.numeric(lower))
  stopifnot(is.numeric(upper))

  stopifnot(is.list(control))

  list(
    method = match.arg( method ),
    lower = lower, upper = upper,
    control = control
  )
}

## ======================================================================
control.nlminb <-
  function(
    control = list( rel.tol = 1.e-5 ),
    lower = -Inf, upper = Inf,
    ...
  )
{

  ## function sets meaningful defaults for selected arguments of function nlminb

  ## 2015-07-17 A. Papritz
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check arguments

  stopifnot(is.numeric(lower))
  stopifnot(is.numeric(upper))

  stopifnot(is.list(control))

  list(
    control = control,
    lower = lower, upper = upper
  )
}

## ======================================================================
control.pcmp <-
  function(
    pmm.ncores = 1, gcr.ncores = 1, max.ncores = detectCores(),
    f = 1, sfstop = FALSE, allow.recursive = TRUE,
    fork = !identical( .Platform[["OS.type"]], "windows" ),
    ...
  )
{

  ## function sets meaningful defaults for parallelized computations

  ## 2014-07-29 AP
  ## 2015-06-30 AP function and arguments renamed
  ## 2015-07-29 AP changes for elimination of parallelized computation of gradient or estimating equations
  ## 2016-07-20 AP renamed function, separate ncores arguments various parallelized computations
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check arguments

  stopifnot(identical(length(sfstop), 1L)          && is.logical(sfstop))
  stopifnot(identical(length(allow.recursive), 1L) && is.logical(allow.recursive))
  stopifnot(identical(length(fork), 1L)            && is.logical(fork))

  stopifnot(identical(length(pmm.ncores), 1L)      && is.numeric(pmm.ncores) && pmm.ncores >= 1)
  stopifnot(identical(length(gcr.ncores), 1L)      && is.numeric(gcr.ncores) && gcr.ncores >= 1)
  stopifnot(identical(length(max.ncores), 1L)      && is.numeric(max.ncores) && max.ncores >= 1)
  stopifnot(identical(length(f), 1L)               && is.numeric(f)          && f >= 1)

  pmm.ncores <- min( pmm.ncores, max.ncores )
  gcr.ncores <- min( gcr.ncores, max.ncores )

  list(
    pmm.ncores = pmm.ncores, gcr.ncores = gcr.ncores, max.ncores = max.ncores,
    f = f, sfstop = sfstop,
    allow.recursive = allow.recursive,
    fork = fork
  )

}

## ======================================================================

compress <-
  function( m )
{

  ## function stores a list of or a single lower, upper triangular or
  ## symmetric matrix compactly

  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  aux <- function( x ){
    struc <- attr( x, "struc" )
    if( !is.null( struc ) ){
      stopifnot(identical(length(struc), 1L))
      switch(
        struc,
        sym = {
          aux <- list( diag = diag( x ), tri = x[lower.tri(x)] )
          attr( aux, "struc" ) <- "sym"
          aux
        },
        lt = {
          aux <- list( diag = diag( x ), tri = x[lower.tri(x)] )
          attr( aux, "struc" ) <- "lt"
          aux
        }
        ,
        ut = {
          aux <- list( diag = diag( x ), tri = x[upper.tri(x)] )
          attr( aux, "struc" ) <- "ut"
          aux
        }
      )
    } else {
      x
    }
  }

  if( is.list( m ) ){
    lapply( m, aux )
  } else {
    aux ( m )
  }



}

## ======================================================================

expand <-
  function( object )
{

  ## function expands a list of or a compactly stored lower, upper
  ## triangular or symmetric matrices

  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  aux <- function( x ){
    struc <- attr( x, "struc" )
    if( !is.null( struc ) ){
      stopifnot(identical(length(struc), 1L))
      switch(
        struc,
        sym = {
          n <- length( x[["diag"]] )
          dn <- names( x[["diag"]] )
          aux <- matrix( 0., n, n )
          aux[lower.tri( aux )] <- x[["tri"]]
          aux <- aux + t( aux )
          diag( aux ) <- x[["diag"]]
          dimnames( aux ) <- list( dn, dn )
          attr( aux, "struc" ) <- "sym"
          aux
        },
        lt = {
          n <- length( x[["diag"]] )
          dn <- names( x[["diag"]] )
          aux <- matrix( 0., n, n )
          aux[lower.tri( aux )] <- x[["tri"]]
          diag( aux ) <- x[["diag"]]
          dimnames( aux ) <- list( dn, dn )
          attr( aux, "struc" ) <- "lt"
          aux
        }
        ,
        ut = {
          n <- length( x[["diag"]] )
          dn <- names( x[["diag"]] )
          aux <- matrix( 0., n, n )
          aux[upper.tri( aux )] <- x[["tri"]]
          diag( aux ) <- x[["diag"]]
          dimnames( aux ) <- list( dn, dn )
          attr( aux, "struc" ) <- "ut"
          aux
        }
      )
    } else {
      x
    }
  }

  ln <- names( object )
  if( is.list( object ) ){
    if( length( ln ) == 2L && all( ln == c( "diag", "tri" ) ) ){
      aux( object )
    } else {
      lapply( object, aux )
    }
  } else {
    object
  }

}

## ======================================================================

param.names <-
  function( model )
{

  ## function returns names of extra parameters of implemented variogram
  ## models (cf. Variogram{RandomFields})

  ## 2012-01-24 A. Papritz
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check arguments

  stopifnot(identical(length(model), 1L) && is.character(model))

  switch(
    model,
    "RMaskey"         = "alpha",
    "RMbessel"        = "nu",
    "RMcauchy"        = "gamma",
    "RMcircular"      = NULL,
    "RMcubic"         = NULL,
    "RMdagum"         = c( "beta", "gamma" ),
    "RMdampedcos"     = "lambda",
    "RMdewijsian"     = "alpha",
    "RMexp"           = NULL,
    "RMfbm"           = "alpha",
    "RMgauss"         = NULL,
    "RMgencauchy"     = c( "alpha", "beta" ),
    "RMgenfbm"        = c( "alpha", "delta" ),
    "RMgengneiting"   = c( "kappa", "mu" ),
    "RMgneiting"      = NULL,
    "RMlgd"           = c( "alpha", "beta" ),
    "RMmatern"        = "nu",
    "RMpenta"         = NULL,
    "RMqexp"          = "alpha",
    "RMspheric"       = NULL,
    "RMstable"        = "alpha",
    "RMwave"          = NULL,
    "RMwhittle"       = "nu",
    stop( model, " variogram not implemented" )
  )
}

##  ##############################################################################

param.bounds <-
function( model, d )
{

  ## function returns range of parameters for which variogram models are
  ## valid (cf.  Variogram{RandomFields})

  ## 2012-03-30 A. Papritz
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check arguments

  stopifnot(identical(length(model), 1L) && is.character(model))
  stopifnot(identical(length(d), 1L)     && is.numeric(d) && d >= 1)

  switch(
    model,
    "RMaskey"         = list( alpha = c( 0.5 * (d + 1L), Inf ) ),
    "RMbessel"        = list( nu = c( 0.5 * (d - 2L), Inf ) ),
    "RMcauchy"        = list( gamma = c( 1.e-18, Inf ) ),
    "RMcircular"      = NULL,
    "RMcubic"         = NULL,
    "RMdagum"         = list( beta = c( 1.e-18, 1.), gamma = c( 1.e-18, 1.-1.e-18) ),
    "RMdampedcos"     = list( lambda = c( if( d > 2L ) sqrt(3.) else 1., Inf ) ),
    "RMdewijsian"     = list( alpha = c( 1.e-18, 2. ) ),
    "RMexp"           = NULL,
    "RMfbm"           = list( alpha = c( 1.e-18, 2.) ),
    "RMgauss"         = NULL,
    "RMgencauchy"     = list( alpha = c(1.e-18, 2.), beta = c(1.e-18, Inf) ),
    "RMgenfbm"        = list( alpha = c(1.e-18, 2.), delta = c(1.e-18, 1.-1.e-18) ),
    "RMgengneiting"   = list( kappa = c(1, 3), mu = c( d/2, Inf ) ),
    "RMgneiting"      = NULL,
    "RMlgd"           = list(
                        alpha = c(
                          1.e-18,
                          if( d <= 3L ) 0.5 * (3L-d) else stop("dimension > 3 not allowed for RMlgd model" )
                        ),
                        beta = c(1.e-18, Inf)
                      ),
    "RMmatern"        = list( nu = c(1.e-18, Inf) ),
    "RMpenta"         = NULL,
    "RMqexp"          = list( alpha = c(0., 1.) ),
    "RMspheric"       = NULL,
    "RMstable"        = list( alpha = c(1.e-18, 2.) ),
    "RMwave"          = NULL,
    "RMwhittle"       = list( nu = c(1.e-18, Inf) ),
    stop( model, " variogram not implemented" )
  )
}


##  ##############################################################################

default.fit.param <-
function(
  variance = TRUE, snugget = FALSE, nugget = TRUE, scale = TRUE,
  alpha = FALSE, beta = FALSE, delta = FALSE, gamma = FALSE,
  kappa = FALSE, lambda = FALSE, mu = FALSE, nu = FALSE )
{

  ## function sets default flags for fitting variogram parameters

  ## 2015-11-27 A. Papritz
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check arguments

  stopifnot(identical(length(variance), 1L) && is.logical(variance))
  stopifnot(identical(length(snugget), 1L)  && is.logical(snugget))
  stopifnot(identical(length(nugget), 1L)   && is.logical(nugget))
  stopifnot(identical(length(scale), 1L)    && is.logical(scale))
  stopifnot(identical(length(alpha), 1L)    && is.logical(alpha))
  stopifnot(identical(length(beta), 1L)     && is.logical(beta))
  stopifnot(identical(length(delta), 1L)    && is.logical(delta))
  stopifnot(identical(length(gamma), 1L)    && is.logical(gamma))
  stopifnot(identical(length(kappa), 1L)    && is.logical(kappa))
  stopifnot(identical(length(lambda), 1L)   && is.logical(lambda))
  stopifnot(identical(length(mu), 1L)       && is.logical(mu))
  stopifnot(identical(length(nu), 1L)       && is.logical(nu))

  c(
    variance = variance, snugget = snugget, nugget = nugget, scale = scale,
    alpha = alpha, beta = beta, delta = delta, gamma = gamma,
    kappa = kappa, lambda = lambda, mu = mu, nu = nu
  )

}


##  ##############################################################################

default.fit.aniso <-
function( f1 = FALSE, f2 = FALSE, omega = FALSE, phi = FALSE, zeta = FALSE )
{

  ## function sets default flags for fitting anisotropy parameters

  ## 2015-11-27 A. Papritz
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check arguments

  stopifnot(identical(length(f1), 1L)    && is.logical(f1))
  stopifnot(identical(length(f2), 1L)    && is.logical(f2))
  stopifnot(identical(length(omega), 1L) && is.logical(omega))
  stopifnot(identical(length(phi), 1L)   && is.logical(phi))
  stopifnot(identical(length(zeta), 1L)  && is.logical(zeta))

  c( f1 = f1, f2 = f2, omega = omega, phi = phi, zeta = zeta )

}


##  ##############################################################################

default.aniso <-
function(
  f1 = 1., f2 = 1., omega = 90., phi = 90., zeta = 0. )
{

  ## function sets default values for anisotropy parameters

  ## 2015-11-26 A. Papritz
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

  ## check arguments

  stopifnot(identical(length(f1), 1L)    && is.numeric(f1) && f1 > 0)
  stopifnot(identical(length(f2), 1L)    && is.numeric(f2) && f2 > 0)
  stopifnot(identical(length(omega), 1L) && is.numeric(omega))
  stopifnot(identical(length(phi), 1L)   && is.numeric(phi))
  stopifnot(identical(length(zeta), 1L)  && is.numeric(zeta))

  c( f1 = f1, f2 = f2, omega = omega, phi = phi, zeta = zeta )

}


##  ##############################################################################
### profilelogLik

profilelogLik <- function( object, values, use.fitted = TRUE, verbose = 0,
  ncores = min( detectCores(), NROW(values) ) ){

  ## function to compute (restricted) likelihood profile for a georob fit

  ## 2015-03-18 A. Papritz
  ## 2015-04-08 AP changes in returned results
  ## 2016-07-14 AP optimization
  ## 2016-07-20 AP changes for parallel computations
  ## 2016-07-28 AP returns gradient in results
  ## 2016-08-08 AP changes for nested variogra
  ## 2016-08-12 AP changes for nested variogram models
  ## 2017-12-22 AP improved memory management in parallel computations
  ## 2020-02-14 AP sanity checks of arguments and for if() and switch()

#### -- auxiliary function

  ## auxiliary function to fit model and return maximized (pseudo) log-likelihood

  f.aux <- function( i ){

    ## values, object, data are taken from parent environment

    ## set fixed initial values

    values <- values[i, ]
    tmp <- strsplit( names(values), control.georob()[["sepstr"]], fixed = TRUE )
    cl <- object[["call"]]
    for( i in 1L:length(values) ){
      cl <- f.call.set_onexxx_to_value( cl, tmp[[i]][1L], values[i], as.integer(tmp[[i]][2L]) )
    }
    object[["call"]] <- cl

    fit <- update( object, data = data )

    ## extract fitted variogram parameters

    param.aniso <- unlist(lapply(
        1L:length(fit[["variogram.object"]]),
        function( i, x ){

          x <- x[[i]]
          res <- c( x[["param"]][x[["fit.param"]]], x[["aniso"]][x[["fit.aniso"]]] )

          if( length(res) ){
            names(res) <- paste( names(res), i, sep = control.georob()[["sepstr"]] )
            res
          } else {
            NULL
          }
        }, x = fit[["variogram.object"]]
      ))

    ## results

    c(
      loglik = logLik(
        fit, warn = FALSE, REML = identical( object[["control"]][["ml.method"]], "REML" )
      ),
      param.aniso,
      coef( fit ),
      gradient = fit[["gradient"]],
      converged = fit[["converged"]]
    )

  }


#### -- check arguments

  ## check whether all mandatory arguments have been provided

  if( missing(object) || missing(values) ) stop(
    "some mandatory arguments are missing"
  )

  stopifnot(identical(class(object)[1], "georob"))

  stopifnot(identical(length(use.fitted), 1L) && is.logical(use.fitted))

  stopifnot(identical(length(verbose), 1L)    && is.numeric(verbose)  && verbose >= 0)
  stopifnot(identical(length(ncores), 1L)     && is.numeric(ncores)   && ncores >= 1)

  ## warning for robust fits

  if( object[["tuning.psi"]] < object[["control"]][["tuning.psi.nr"]] ){
    warning(
      "likelihood approximated for robustly fitted model by likelihood of\n",
      "  equivalent Gaussian model with heteroscedastic nugget"
    )
  }

  if( !(is.matrix(values) || is.data.frame( values ) || is.list( values ) ) ) stop(
    "'values' must be a dataframe or a matrix"
  )


#### -- prepare data

  ## get data.frame with required variables (note that the data.frame passed
  ## as data argument to georob must exist in GlobalEnv)

  data <- cbind(
    get_all_vars(
      formula( object ), data = eval( getCall(object)[["data"]] )
    ),
    get_all_vars(
      object[["locations.objects"]][["locations"]], eval( getCall(object)[["data"]] )
    )
  )

  if( identical( class( object[["na.action"]] ), "omit" ) ) data <- na.omit(data)

  ## select subset if appropriate

  if( !is.null( getCall(object)[["subset"]] ) ){
    data <- data[eval( getCall(object)[["subset"]] ), ]
  }

  ## extract variogram.object

  ## check values

  if( is.data.frame(values) ) values <- list( values )

  if( !identical( length(values), length(object[["variogram.object"]]) ) ) stop(
    "'values' must be a list of the same length as 'variogram.object'"
  )

  ## check names of fixed variogram parameters and fix respective
  ## parameters in variogram.object

  values <- lapply(
    1L:length(values),
    function( i, v, vo ){

      v <- v[[i]]
      vo <- vo[[i]]

      fixed.param.aniso <- colnames( v )

      ## match fixed.param.aniso
      fixed.param.aniso <- sapply(
        fixed.param.aniso, match.arg,
        choices = c( names(default.fit.param()), names(default.fit.aniso()) )
      )
      if(
        any( !fixed.param.aniso %in% c( names( vo[["param"]] ), names( vo[["aniso"]] ) ) )
      ) stop( "column names of 'values' do not match names of variogram parameters" )

      colnames(v) <- fixed.param.aniso
      v

    }, v = values, vo = object[["variogram.object"]]
  )


#### -- manipulate call to compute only required items

  ## manipulate call so that fitted values are used as initial values

  cl <- f.call.set_allxxx_to_fitted_values( object )

  ## fix initial values for parameters present in values

  for( i in 1L:length(values) ){
    for( nme in colnames(values[[i]]) ){
      cl <- f.call.set_onefitxxx_to_value( cl, nme, FALSE, i )
    }
  }

  ## update object call to avoid computation of covariance matrices and
  ## hessian and to set reparam = FALSE if variance parameters are fitted

  cl <- f.call.set_x_to_value( cl,  "verbose", verbose )

  reparam <- !any(
    unique( sapply( values, colnames ) )  %in% c( "variance", "snugget", "nugget" )
  )
  cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "reparam", reparam )

  cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "hessian", FALSE )

  cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "cov.bhat", FALSE )
  cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "cov.betahat", FALSE )
  cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "cov.delta.bhat", FALSE )
  cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "cov.delta.bhat.betahat", FALSE )
  cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "cov.ehat", FALSE )
  cl <- f.call.set_x_to_value_in_fun( cl, "control", "control.georob", "cov.ehat.p.bhat", FALSE )

  object[["call"]] <- cl


#### -- prepare sets of parameter values for which to compute likelihood

  ## rename variogram parameters in values

  values <- lapply(
    1L:length(values),
    function( i, v, vo ){
      v <- v[[i]]
      names( v ) <- paste( names(v), i, sep = control.georob()[["sepstr"]] )
      v
    }, v = values
  )

  ## expand values to data frame

  if( length(values) > 1L ){
    tmp <- values
    values <- as.list( tmp[[1L]] )

    for( i in 2L:length(tmp) ){
      for( j in colnames( tmp[[i]] ) ){
        values <- c( values, as.list( tmp[[i]][, j, drop=FALSE] ) )
      }
    }
    values <- expand.grid( values )
  } else {
    values <- values[[1L]]
  }


#### -- fit model for sets of parameter values for which to compute likelihood

  ## loop over all elements of values

  values <- as.matrix( values )

  ## set default value for control of forking if missing (required for backward compatibility)

  if( is.null( object[["control"]][["pcmp"]][["fork"]] ) ){
    object[["control"]][["pcmp"]][["fork"]] <- !identical( .Platform[["OS.type"]], "windows" )
  }

  if( ncores > 1L && !object[["control"]][["pcmp"]][["fork"]] ){

    ## create a SNOW cluster on windows OS

    clstr <- makeCluster( ncores, type = "SOCK" )
    save( clstr, file = "SOCKcluster.RData" )
    options( error = f.stop.cluster )

    ## export required items to workers

    junk <- clusterEvalQ( clstr, require( georob, quietly = TRUE ) )
    junk <- clusterExport( clstr, c( "values", "object", "data" ), envir = environment() )

    result <- parLapply(
      clstr,
      1L:NROW(values),
      f.aux
    )

    f.stop.cluster( clstr )

  } else {

    ## fork child processes on non-windows OS

    result <- mclapply(
      1L:NROW(values),
      f.aux,
      mc.cores = ncores,
      mc.allow.recursive = object[["control"]][["pcmp"]][["allow.recursive"]]
    )

  }


#### -- prepare output

  ## collect results

  result <- as.data.frame( cbind( values, t( simplify2array( result ) ) ) )

  if( identical( length( object[["variogram.object"]]), 1L ) ){
    tmp <- colnames(result)
    tmp <- gsub(
      paste( control.georob()[["sepstr"]], "1", sep = "" ),
      "", tmp
    )
    colnames(result) <- tmp
  }
  result

}

