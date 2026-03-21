##' @title Mediation Analysis for ZEro-inflated mediators
##'
##' @description A novel mediation modeling approach to address zero-inflated mediators containing true zeros and false zeros.
##'
##' @details For an independent variable \eqn{X}, a zero-inflated mediator \eqn{M} and a continuous outcome variable \eqn{Y}, the following regression equation is used to model the association between \eqn{Y} and \eqn{(X,M)}:
##' \deqn{Y_{xm1_{(m>0)}}=\beta_0+\beta_1m+\beta_2 1_{(m>0)}+\beta_3x+\beta_4x1_{(m>0)}+\beta_5xm+\epsilon}
##' Users can choose to include either one, both, or none of the two exposure-mediator interaction terms between (i) \eqn{X} and \eqn{1_{(M>0)}} and (ii) \eqn{X} and \eqn{M} using the argument `XMint`.
##'
##' For mediators, zero-inflated log-normal, zero-inflated negative binomial, and zero-inflated Poisson distributions are considered and can be specified through the argument `distM`.
##'
##' The indirect and direct effects (NIE1, NIE2, NIE, NDE, and CDE) are estimated for \eqn{X} changing from `x1` to `x2`. When confounders are present, the conditional effects are estimated given the fixed value `zval`.
##'
##' @param data a data frame containing variables: an independent variable \eqn{X}, a mediator \eqn{M}, an outcome \eqn{Y}, and confounder variables \eqn{Z} (if any). See example dataset: `data(zinb10)`
##' @param distM a vector with choices of the distribution of mediator to try with. One or more of '`zilonm`', '`zinbm`', and '`zipm`' for zero-inflated log-normal, negative binomial, and Poisson mediators respectively. Default is `c('zilonm', 'zinbm', 'zipm')` where all three distributions are fitted and the final mediation model is selected by model selection criterion `selection`
##' @param K a vector with choices of the number of component \eqn{K} in the zero-inflated mixture mediators to try with. Default is \eqn{K=1} for zero-inflated (non-mixture) mediators
##' @param selection model selection criterion when more than one model (combination of different values in `distM` and `K`) is fitted. Either '`AIC`' or '`BIC`'. Default is '`AIC`'
##' @param X name of the independent variable. Can be continuous or discrete
##' @param M name of the mediator variable. Non-negative values
##' @param Y name of the outcome variable. Continuous values
##' @param Z name(s) of confounder variables (if any)
##' @param XMint a logical vector of length 2 indicating whether to include the two exposure-mediator interaction terms between (i) \eqn{X} and \eqn{1_{(M>0)}} and (ii) \eqn{X} and \eqn{M}. Default is `c(TRUE, FALSE)`, which only includes the first
##' @param x1 the first value of independent variable of interest
##' @param x2 the second value of independent variable of interest
##' @param zval a vector of value(s) of confounders to be conditional on when estimating effects
##' @param mval the fixed value of mediator to be conditional on when estimating CDE
##' @param B the upper bound value \eqn{B} to be used in the probability mechanism of observing false zeros
##' @param seed an optional seed number to control randomness
##' @param ncore number of cores available for parallel computing
##' @return a list containing:
##' - `results_effects`: a data frame for the results of estimated effects (NIE1, NIE2, NIE, NDE, and CDE). '_cond' for conditional effects at `zval` and '_avg' for average effects
##' - `results_parameters`: a data frame for the results of model parameters
##' - `selected_model_name`: a string for the distribution of \eqn{M} and number of components \eqn{K} selected in the final mediation model
##' - `BIC`: a numeric value for the BIC of the final mediation model
##' - `AIC`: a numeric value for the AIC of the final mediation model
##' - `models`: a list with all fitted models
##' - `analysis2_out`: a list with output from `analysis2()` function (used for internal check)
##' @author Meilin Jiang <meilin.jiang@@ufl.edu> and Zhigang Li <zhigang.li@@ufl.edu>
##' @import stats foreach doParallel MASS
##' @importFrom flexmix flexmix clusters parameters FLXMRglm
##' @importFrom pracma hessian
##' @importFrom numDeriv grad jacobian
##' @export
##' @examples
##' data(zinb10)
##' \donttest{
##' maze_out <- MAZE(data = zinb10,
##'                  distM = c('zilonm', 'zinbm', 'zipm'),  K = 1,
##'                  selection = 'AIC',
##'                  X = 'X', M = 'Mobs', Y = 'Y', Z = NULL,
##'                  XMint = c(TRUE, FALSE),
##'                  x1 = 0, x2 = 1, zval = NULL, mval = 0,
##'                  B = 20, seed = 1)
##' ## results of selected mediation model
##' maze_out$results_effects # indirect and direct effects
##' maze_out$selected_model_name # selected distribution of the mediator and number of components K
##' maze_out$results_parameters # model parameters
##' maze_out$BIC; maze_out$AIC # BIC and AIC of the selected mediation model
##' }

MAZE <- function(data, distM = c("zilonm", "zinbm", "zipm"), K = 1, selection = "AIC",
    X, M, Y, Z = NULL, XMint = c(TRUE, FALSE), x1, x2, zval = NULL, mval = 0, B = 20,
    seed = 1, ncore = 1) {
    distM_sequence <- distM
    K_sequence <- K
    limits <- 0.001
    # set.seed(seed)

    data <- as.data.frame(data)
    dat <- data.frame(X = data[, X], Y = data[, Y], Mobs = data[, M])
    num_Z <- length(Z)
    if (is.null(Z)) {
        Z_names <- NULL
    } else {
        dat <- data.frame(dat, data[, Z])
        Z_names <- paste0("Z", 1:num_Z)
        names(dat)[3 + 1:num_Z] <- Z_names
    }
    # complete data analysis
    dat <- dat[complete.cases(dat), ]

    out <- tryCatch({
        realanalysis(dat, distM_sequence, K_sequence, selection, XMint, x1, x2, zval,
            num_Z, Z_names, mval, limits, B, seed, ncore)
    }, error = function(e) {
        print(paste0("realanalysis()_error = ", e))
        # list(results_effects = NA, results_parameters = NA, BIC = Inf, AIC =
        # Inf, e = e)
    })
    return(out)
}




