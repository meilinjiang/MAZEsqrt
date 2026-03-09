##' @title DataGen
##'
##' @description Generate data under zero-inflated mediation models and calculate the true effects
##'
##' @param distM distribution of the mediator. One of '`zilonm`', '`zinbm`', and '`zipm`' for zero-inflated log-normal, negative binomial, and Poisson mediators respectively
##' @param theta vector of true parameter values
##' @param K true number of component \eqn{K} in the zero-inflated mixture mediators. Default is \eqn{K=1} for zero-inflated (non-mixture) mediators
##' @param num_Z number of confounder variables
##' @param n number of observations to generate
##' @param B the upper bound value \eqn{B} to be used in the probability mechanism of observing false zeros
##' @param x1 the first value of independent variable of interest
##' @param x2 the second value of independent variable of interest
##' @param zval the value of confounders to be conditional on when calculating true effects
##' @param mval the fixed value of mediator to be conditional on when calculating true CDE
##' @param true_gen using true data generation pf0
##' @return
##' true_eff: a vector containing true effects (NIE1, NIE2, NIE, NDE, and CDE)
##'
##' dat: a data frame containing variables:
##' - `X`: independent variable,
##' - `Mobs`: observed mediator values (with possibly false zeros)
##' - `M`: true mediator values,
##' - `Y`: outcome,
##' - `Z`: confounder variables (if any)
##' @author Meilin Jiang <meilin.jiang@@ufl.edu> and Zhigang Li <zhigang.li@@ufl.edu>
##' @import stats
##' @export
##' @examples
##' betas.tr <- c(2, 0.12, -6.6, 6.3, -3.8, 0)
##' delta.tr <- 1.1
##' alpha0_k.tr <- c(0.4, 1.1)
##' alpha1_k.tr <- c(0.1, 0.5)
##' alphas.tr <- rbind(alpha0_k.tr,alpha1_k.tr)
##' xi0.tr <- -1.5
##' psi_km1.tr <- c(0.6)
##' gammas.tr <- c(-1.8, 0.5)
##' eta.tr <- 1
##' theta <- c(betas.tr, delta.tr, alphas.tr,
##'            xi0.tr, psi_km1.tr, gammas.tr, eta.tr)
##' out <- DataGen(distM = 'zinbm', theta, K = 2, num_Z=0,
##'                n = 200, B = 20, x1 = 0, x2 = 1, zval = NULL, mval = 0)
##' (true_eff <- out$true_eff)
##' dat <- out$dat


DataGen <- function(distM, theta, K, num_Z = 0, n, B, x1, x2, zval = NULL, mval = 0,true_gen=T) {
    dat_placeholder <- data.frame(NULL)
    class(dat_placeholder) <- c(distM, class(dat_placeholder))
    if(true_gen){
      out <- DataGen_call_true(dat_placeholder, theta, K, num_Z, n, B, x1, x2, zval, mval)
    }else{
      out <- DataGen_call(dat_placeholder, theta, K, num_Z, n, B, x1, x2, zval, mval)
    }

    return(out)
}

