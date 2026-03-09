#' @noRd
DataGen_call <- function(dat_placeholder, ...) {
    UseMethod("DataGen_call")
}

#' @noRd
DataGen_call.zinbm <- function(dat_placeholder, theta, K, num_Z, n, B, x1, x2, zval,
    mval) {
    X <- rnorm(n)
    theta_trans <- trans(dat_placeholder, theta, K, xval = X, num_Z, zval = NULL)
    eps <- rnorm(n, 0, theta_trans[["delta"]])

    ind_nb <- rbinom(n, 1, 1 - theta_trans[["Delstar_i"]])
    # ind_nb = 1: M is from Poisson mixture dist 1 (M is from the corresponding
    # Poisson distribution in the mixture)
    L_ik <- t(rmultinom(n, 1, theta_trans[["psi_k"]]))
    M <- rowSums(L_ik * apply(theta_trans[["mu_ik"]], 2, function(mu) rnbinom(n, size = theta_trans[["r"]],
        mu = mu)))
    # M <- rnbinom(n, size=theta_trans[['r']],
    # mu=rowSums(L_ik*theta_trans[['mu_ik']])) NB could also generate 0's ind_nb
    # = 0: M is an excessive zero
    M[ind_nb == 0] <- 0
    # 1(M > 0)
    ind <- (M > 0) * 1

    L <- apply(L_ik, 1, which.max)
    L[ind_nb == 0] <- NA


    Y <- theta_trans[["beta0"]] + theta_trans[["beta1"]] * M + theta_trans[["beta2"]] *
        ind + theta_trans[["beta3"]] * X + theta_trans[["beta4"]] * X * ind + theta_trans[["beta5"]] *
        X * M + eps
    # probability of observing false zeros
    if (theta_trans[["eta"]] == Inf) {
        pf0 <- rep(0, n)
    } else {
        pf0 <- exp(-M * theta_trans[["eta"]]^2)
    }
    pf0[M > B] <- 0
    R <- rbinom(n, 1, 1 - pf0)  # R=0: false zeros
    dat <- data.frame(X, Y, M, R, Mobs = M * R, L, ind_nb)

    true_eff <- effects(dat_placeholder, theta, x1, x2, K, num_Z, zval, mval, XMint = c(T,
        T), calculate_se = F, vcovar = NULL, Group1 = F)$eff
    names(true_eff) <- c("NIE1", "NIE2", "NIE", "NDE", "CDE")
    out <- list(true_eff = true_eff, dat = dat)
    return(out)
}


#' @noRd
DataGen_call_true <- function(dat_placeholder, ...) {
    UseMethod("DataGen_call")
}
#' @noRd
DataGen_call_true.zinbm <- function(dat_placeholder, theta, K, num_Z, n, B, x1, x2, zval,
                               mval) {
    X <- rnorm(n)
    theta_trans <- trans(dat_placeholder, theta, K, xval = X, num_Z, zval = NULL)
    eps <- rnorm(n, 0, theta_trans[["delta"]])
    
    ind_nb <- rbinom(n, 1, 1 - theta_trans[["Delstar_i"]])
    # ind_nb = 1: M is from Poisson mixture dist 1 (M is from the corresponding
    # Poisson distribution in the mixture)
    L_ik <- t(rmultinom(n, 1, theta_trans[["psi_k"]]))
    M <- rowSums(L_ik * apply(theta_trans[["mu_ik"]], 2, function(mu) rnbinom(n, size = theta_trans[["r"]],
                                                                              mu = mu)))
    # M <- rnbinom(n, size=theta_trans[['r']],
    # mu=rowSums(L_ik*theta_trans[['mu_ik']])) NB could also generate 0's ind_nb
    # = 0: M is an excessive zero
    M[ind_nb == 0] <- 0
    # 1(M > 0)
    ind <- (M > 0) * 1
    
    L <- apply(L_ik, 1, which.max)
    L[ind_nb == 0] <- NA
    
    
    Y <- theta_trans[["beta0"]] + theta_trans[["beta1"]] * M + theta_trans[["beta2"]] *
        ind + theta_trans[["beta3"]] * X + theta_trans[["beta4"]] * X * ind + theta_trans[["beta5"]] *
        X * M + eps
    # probability of observing false zeros
    if (theta_trans[["eta"]] == Inf) {
        pf0 <- rep(0, n)
    } else {
        pf0 <- exp(-sqrt(M) * theta_trans[["eta"]]^2)
    }
    pf0[M > B] <- 0
    R <- rbinom(n, 1, 1 - pf0)  # R=0: false zeros
    dat <- data.frame(X, Y, M, R, Mobs = M * R, L, ind_nb)
    
    true_eff <- effects(dat_placeholder, theta, x1, x2, K, num_Z, zval, mval, XMint = c(T,
                                                                                        T), calculate_se = F, vcovar = NULL, Group1 = F)$eff
    names(true_eff) <- c("NIE1", "NIE2", "NIE", "NDE", "CDE")
    out <- list(true_eff = true_eff, dat = dat)
    return(out)
}

