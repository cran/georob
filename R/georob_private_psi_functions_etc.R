################################################################################
### f.psi.function

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
