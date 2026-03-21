#' @noRd
trans.zipm <- function(dat, theta, K, xval, num_Z, zval) {
    beta0 <- theta[1]
    beta1 <- theta[2]
    beta2 <- theta[3]
    beta3 <- theta[4]
    beta4 <- theta[5]
    beta5 <- theta[6]
    delta <- theta[7 + num_Z]
    alpha_k <- theta[7 + num_Z + 1:((2 + num_Z) * K)]

    if (K == 1) {
        psi_k <- 1
    } else {
        psi_k <- c(theta[7 + num_Z + (2 + num_Z) * K + 1:(K - 1)], 1 - sum(theta[7 +
            num_Z + (2 + num_Z) * K + 1:(K - 1)]))
    }
    gammas <- theta[7 + num_Z + (2 + num_Z) * K + (K - 1) + 1:(2 + num_Z)]
    eta <- theta[8 + num_Z + (2 + num_Z) * K + (K - 1) + (2 + num_Z)]

    n0 <- length(xval)
    psi_ik <- k_to_ik(psi_k, n0)
    if (num_Z == 0) {
        beta_T_Z <- 0
        # zval <- NULL
    } else {
        if (n0 == 1) {
            zval <- matrix(zval, nrow = 1)
        }
        beta_Z <- theta[6 + 1:num_Z]
        beta_T_Z <- rowSums(k_to_ik(beta_Z, n0) * zval)
    }
    designMat_M <- cbind(1, xval, zval)
    loglambda_ik <- matrix(apply(matrix(alpha_k, ncol = K), 2, function(t) rowSums(k_to_ik(t,
        n0) * designMat_M)), ncol = K)
    lambda_ik <- exp(loglambda_ik)

    Delstar_i <- expit(rowSums(k_to_ik(gammas, n0) * designMat_M))
    Del_i <- Delstar_i + (1 - Delstar_i) * rowSums(psi_ik * exp(-lambda_ik))
    # neg_log_em1 = -log[exp(lambda) - 1], when lambda too small, log(0) causes
    # issues Taylor expansion: when lambda -> 0 exp(lambda) = 1 + lambda +
    # lambda^2/2 + lambda^3/6 + ...  when lambda is small (< 1e-10), take the
    # first term: exp(lambda) = 1 + lambda, so neg_log_em1 = -log(lambda) when
    # lambda is big: exp(lambda) - 1 = exp(lambda)[1 - exp(-lambda)] more stable
    neg_log_em1 <- -loglambda_ik
    neg_log_em1[lambda_ik > 1e-10] <- (-lambda_ik - log(1 - exp(-lambda_ik)))[lambda_ik >
        1e-10]

    theta_trans <- list(beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3,
        beta4 = beta4, beta5 = beta5, delta = delta, eta = eta, beta_T_Z = beta_T_Z,
        alpha_k = alpha_k, loglambda_ik = loglambda_ik, lambda_ik = lambda_ik, psi_k = psi_k,
        psi_ik = psi_ik, Delstar_i = Delstar_i, Del_i = Del_i, neg_log_em1 = neg_log_em1)
    return(theta_trans)
}
#' @noRd
ComputeInit.zipm <- function(dat, K, num_Z, Z_names, XMint) {
    # group 1, flexmix
    dat_g1 <- dat[which(dat$Mobs > 0), ]

    # fit <- summary(lm(Y_group1 ~ M_group1 + X_group1))
    fm_Y_rhs <- c("Mobs", "X", if (XMint[2]) "Mobs:X" else NULL, Z_names)
    fm_Y <- as.formula(paste0("Y~", paste0(fm_Y_rhs, collapse = "+"), collapse = ""))
    fit <- summary(lm(fm_Y, data = dat_g1))
    # beta02, beta1, beta34, (beta5), beta_Z beta02, beta1, beta34, (beta5),
    # beta_Z
    betas <- fit$coefficients[c("(Intercept)", fm_Y_rhs), "Estimate"]
    delta <- fit$sigma
    # poisson reg has loglink so no need to log(Mobs)
    fm_M <- as.formula(paste0("Mobs~", paste0(c("X", Z_names), collapse = "+"), collapse = ""))
    cl <- flexmix(fm_M, k = K, data = dat_g1, model = FLXMRglm(family = "poisson"))
    cl_size <- table(clusters(cl))
    while (length(cl_size) != K) {
        cl <- flexmix(fm_M, k = K, data = dat_g1, model = FLXMRglm(family = "poisson"))
        cl_size <- table(clusters(cl))
    }
    # method 3: distribution 1 as the one with smaller cluster mean since x=age,
    # intercepts are not accurate enough
    mu <- tapply(dat_g1$Mobs, clusters(cl), mean)
    # mu <- parameters(cl)[1,]+parameters(cl)[2,]*mean(X_group1)
    ord <- order(mu)
    mix <- parameters(cl)[, ord, drop = F]
    # alpha_0k, alpha_1k, xi_0, psik
    mix_init <- c(mix[1:(2 + num_Z), ])
    psi_k <- cl_size[ord]/nrow(dat_g1)
    mix_init <- c(mix_init, psi_k[-K])

    init <- c(betas, delta, mix_init)
    return(init)
}

# negative expectation of log-likelihood function with respect to conditional
# distribution of 1(Ci = k) given data and current estimates for group 1: -Q1
#' @noRd
negQ_G1.zipm <- function(dat, theta, K, num_Z, Z_names, B = NULL, tauG1 = NULL, calculate_tau = F,
    calculate_ll = F) {
    # betas, delta, alpha0k, alpha1k, alphaZ_k, w0k, gammas, eta
    dat_g1 <- dat[which(dat$Mobs > 0), ]
    M_group1 <- dat_g1$Mobs
    Y_group1 <- dat_g1$Y
    X_group1 <- dat_g1$X
    if (num_Z == 0) {
        Z_group1 <- NULL
    } else {
        Z_group1 <- dat_g1[, Z_names]
    }
    theta_trans <- trans(dat, theta, K, X_group1, num_Z, Z_group1)
    log_dpois_nz <- M_group1 * theta_trans[["loglambda_ik"]] + theta_trans[["neg_log_em1"]]

    # group 1
    if (!calculate_tau) {
        # if(K == 1){ calculate_ll <- TRUE }
        if (theta_trans[["eta"]] == Inf) {
            pf0 <- rep(0, length(M_group1))
        } else {
            pf0 <- exp(-theta_trans[["eta"]]^2 * sqrt(M_group1))
        }
        pf0[M_group1 > B] <- 0

        l1_ik <- -log(theta_trans[["delta"]]) + log(1 - pf0) - 0.5 * log(2 * pi) +
            log_dpois_nz - sapply(M_group1, function(t) {
            sum(log(1:t))
        }) - (Y_group1 - theta_trans[["beta0"]] - theta_trans[["beta1"]] * M_group1 -
            theta_trans[["beta2"]] - (theta_trans[["beta3"]] + theta_trans[["beta4"]]) *
            X_group1 - theta_trans[["beta5"]] * X_group1 * M_group1 - theta_trans[["beta_T_Z"]])^2/(2 *
            theta_trans[["delta"]]^2)

        bigpsi_ik <- (1 - theta_trans[["Del_i"]]) * theta_trans[["psi_ik"]]
        pow <- log(bigpsi_ik) + l1_ik
        if (!calculate_ll) {
            # calculate -Q1
            Q1_ik <- tauG1 * pow
            out <- -sum(Q1_ik)
        } else {
            # calculate -l1 alternative method for calculating hessian: from log
            # likelihood function parameter estimation cannot use it (has many
            # saddle points so need EM), but hessian can l1_i <- log(
            # bigpsi1*exp(l1_ik1) + bigpsi2*exp(l1_ik2) )
            pow_max <- apply(pow, 1, max)
            l1_i <- log(rowSums(exp(pow - pow_max))) + pow_max
            out <- -sum(l1_i)
        }
    } else {
        # calculate conditional expectation of 1(Ci = k) given data and current
        # estimates
        if (K == 1) {
            out <- matrix(1, nrow = length(M_group1), ncol = 1)
        } else {
            # num <- theta_trans[['psi_ik']]*exp(log_dpois_nz) out <-
            # num/rowSums(num)
            pow <- log(theta_trans[["psi_ik"]]) + log_dpois_nz
            out <- NULL
            for (k in seq_len(K)) {
                set <- (1:K)[-k]
                pow2 <- pow[, set, drop = F] - pow[, k]
                out <- cbind(out, 1/(1 + rowSums(exp(pow2))))
            }
        }
    }
    return(out)
}
#' @noRd
bounds.zipm <- function(dat, K, group, num_Z, XMint) {
    s <- ifelse(group == 1, 3, 4 + XMint[1]) + XMint[2] + num_Z
    ul <- 1000
    if (K == 1) {
        if (group == 1) {
            ui <- diag(rep(1, s + 2 + num_Z + 1))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, 2 + num_Z))
        } else {
            ui <- diag(rep(1, s + 2 + num_Z + 1 + 2 + num_Z + 1))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, 2 + num_Z + 2 + num_Z), 1e-06)
        }
    } else {
        if (group == 1) {
            ui1 <- diag(rep(1, s + 1 + (2 + num_Z) * K + (K - 1)))
            ui2 <- diag(c(rep(0, s + 1 + (2 + num_Z) * K), rep(-1, K - 1)))
            ui <- rbind(ui1, ui2[s + 1 + (2 + num_Z) * K + 1:(K - 1), ], c(rep(0,
                s + 1 + (2 + num_Z) * K), rep(-1, K - 1)))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, (2 + num_Z) * K), rep(1e-06,
                K - 1), rep(-(1 - 1e-06), K - 1), -(1 - 1e-06))
        } else {
            ui1 <- diag(rep(1, s + 2 + (2 + num_Z) * K + (K - 1) + 2 + num_Z))
            ui2 <- diag(c(rep(0, s + 1 + (2 + num_Z) * K), rep(-1, K - 1), rep(0,
                (2 + num_Z) + 1)))
            ui <- rbind(ui1, ui2[s + 1 + (2 + num_Z) * K + 1:(K - 1), ], c(rep(0,
                s + 1 + (2 + num_Z) * K), rep(-1, K - 1), rep(0, (2 + num_Z) + 1)))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, (2 + num_Z) * K), rep(1e-06,
                K - 1), rep(-1000, 2 + num_Z), 1e-06, rep(-(1 - 1e-06), K - 1), -(1 -
                1e-06))
        }
    }
    return(list(ui = ui, ci = ci))
}
#' @noRd
ComputeInit2.zipm <- function(dat, K, num_Z, Z_names, XMint, B, limits, explicit = T) {
    init2 <- ComputeInit(dat, K, num_Z, Z_names, XMint)
    tau2 <- negQ2_G1(dat, init2, K, num_Z, Z_names, XMint, calculate_tau = T)
    init <- rep(0, length(init2))
    countEM <- 0
    X_group1 <- dat$X[which(dat$Mobs > 0)]
    bd <- bounds(dat, K, group = 1, num_Z, XMint)
    # M-step
    while (euclidean_dist(init, init2) > limits) {
        init <- init2
        tau <- tau2

        # m <- nlminb(init,function(x)negQ2_G1(dat,x,K,num_Z,Z_names,
        # XMint,B,tau),lower = c(rep(-1000,9),1e-6))
        m <- constrOptim(init, function(x) negQ2_G1(dat, x, K, num_Z, Z_names, XMint,
            B, tau), grad = NULL, ui = bd$ui, ci = bd$ci, outer.iterations = 500,
            control = list(maxit = 50000))
        if (m$convergence == 0) {
            init2 <- m$par
        } else {
            print(paste0("dist = zilonm, K = ", K, ", ComputeInit2 = ", m$convergence,
                ", countEM = ", countEM))
            init2 <- rep(NA, length(init))
            break
        }
        switch <- compare_mu(dat, init2, group = 1, K, num_Z, Z_names, XMint)
        init2 <- switch$init2
        tau2 <- switch$tau2
        countEM <- countEM + 1
        if (K == 1) {
            break
        }
    }
    # beta02, beta1, beta34, (beta5), beta_Z, delta, alpha_0k, alpha_1k, w_0k
    initials <- G1_init(dat, init2, K, num_Z, Z_names, XMint, initials_for_full = T)
    return(initials)
}

# log(factorial(m)) = sum(log(m)) more stable using log
#' @noRd
loghik_zipm <- function(m, loglambda, beta0, beta1, beta2, beta3, beta4, beta5, beta_T_Z,
    delta, eta, x, y) {
    m * loglambda - ((y - beta0 - (beta1 * m) - beta2 - (beta3 + beta4) * x - beta5 *
        x * m - beta_T_Z)^2)/(2 * delta^2) - eta^2 * sqrt(m) - sum(log(1:m))
}

# negative expectation of log-likelihood function with respect to conditional
# distribution of 1(Ci = k) given data and current estimates for group 2: -Q2
#' @noRd
negQ_G2.zipm <- function(dat, theta, K, num_Z, Z_names, B, tauG2 = NULL, calculate_tau = F,
    calculate_ll = F) {
    # group 2
    Y_group2 <- dat$Y[which(dat$Mobs == 0)]
    X_group2 <- dat$X[which(dat$Mobs == 0)]
    if (num_Z == 0) {
        zval <- NULL
    } else {
        zval <- Z_group2 <- dat[which(dat$Mobs == 0), Z_names]
    }
    theta_trans <- trans(dat, theta, K, X_group2, num_Z, zval)

    # summation: output is the log[sum(hik)]
    loghik_m <- list(NULL)
    for (i in 1:B) {
        loghik_m[[i]] <- loghik_zipm(m = i, loglambda = theta_trans[["loglambda_ik"]],
            theta_trans[["beta0"]], theta_trans[["beta1"]], theta_trans[["beta2"]],
            theta_trans[["beta3"]], theta_trans[["beta4"]], theta_trans[["beta5"]],
            theta_trans[["beta_T_Z"]], theta_trans[["delta"]], theta_trans[["eta"]],
            x = X_group2, y = Y_group2)
    }
    loghik_mmax <- NULL
    for (k in seq_len(K)) {
        loghik_mmax <- cbind(loghik_mmax, apply(sapply(loghik_m, function(t) t[, k]),
            1, max))
    }
    output <- log(Reduce("+", lapply(loghik_m, function(t) exp(t - loghik_mmax)))) +
        loghik_mmax

    # l2_ik <- cbind(-log(theta_trans[['delta']]) -
    # (Y_group2-theta_trans[['beta0']]-theta_trans[['beta3']]*X_group2-theta_trans[['beta_T_Z']])^2/(2*theta_trans[['delta']]^2)
    # - 0.5*log(2*pi), -log(theta_trans[['delta']]) - 0.5*log(2*pi) +
    # theta_trans[['neg_log_em1']] + output )
    l2_ik <- cbind(-(Y_group2 - theta_trans[["beta0"]] - theta_trans[["beta3"]] *
        X_group2 - theta_trans[["beta_T_Z"]])^2/(2 * theta_trans[["delta"]]^2), theta_trans[["neg_log_em1"]] +
        output) - log(theta_trans[["delta"]]) - 0.5 * log(2 * pi)

    bigpsi_ik <- cbind(theta_trans[["Del_i"]], (1 - theta_trans[["Del_i"]]) * theta_trans[["psi_ik"]])
    pow <- log(bigpsi_ik) + l2_ik
    if (!calculate_tau) {
        # if(K == 1){ calculate_ll <- TRUE }
        if (!calculate_ll) {
            # calculate -Q2
            Q2_ik <- tauG2 * pow
            out <- -sum(Q2_ik)
        } else {
            # calculate -l2 alternative method for calculating hessian: from log
            # likelihood function parameter estimation cannot use it (has many
            # saddle points so need EM), but hessian can l2_i <- log(
            # bigpsi0*exp(l2_ik0) + bigpsi1*exp(l2_ik1) + bigpsi2*exp(l2_ik2) )
            pow_max <- apply(pow, 1, max)
            l2_i <- log(rowSums(exp(pow - pow_max))) + pow_max
            out <- -sum(l2_i)
        }
    } else {
        # calculate conditional expectation of 1(Ci = k) given data and current
        # estiamtes
        ll <- exp(pow)
        out <- ll/rowSums(ll)
    }
    return(out)
}
#' @noRd
effects.zipm <- function(dat, theta, x1, x2, K, num_Z, zval, mval, XMint, calculate_se = F,
    vcovar = NULL, Group1 = F) {
    if (Group1) {
        # theta <- c(theta[1:2], 0, theta[3], 0, theta[3 + 1:(num_Z + (2 +
        # num_Z) * K + 1 + K - 1)], -Inf, 0, rep(0, num_Z), Inf)
    }
    theta <- G12_init(theta, XMint)
    x12 <- c(x1, x2)
    if (num_Z == 0) {
        z12 <- NULL
    } else {
        z12 <- rbind(zval, zval)
    }
    theta_trans <- trans(dat, theta, K, xval = x12, num_Z, zval = z12)

    Del_x12 <- theta_trans[["Del_i"]]
    Delstar_x12 <- theta_trans[["Delstar_i"]]
    m_x12 <- rowSums(theta_trans[["psi_ik"]] * theta_trans[["lambda_ik"]])

    NIE1 <- (theta_trans[["beta1"]] + theta_trans[["beta5"]] * x2) * diff((1 - Delstar_x12) *
        m_x12)
    NIE2 <- (theta_trans[["beta2"]] + theta_trans[["beta4"]] * x2) * (-diff(Del_x12))
    NIE <- NIE1 + NIE2
    NDE <- diff(x12) * (theta_trans[["beta3"]] + theta_trans[["beta4"]] * (1 - Del_x12[1]) +
        theta_trans[["beta5"]] * (1 - Delstar_x12[1]) * m_x12[1])
    CDE <- diff(x12) * (theta_trans[["beta3"]] + theta_trans[["beta4"]] * (mval >
        0) + theta_trans[["beta5"]] * mval)
    out <- data.frame(eff = c(NIE1, NIE2, NIE, NDE, CDE))

    if (calculate_se == T) {
        # asymptotic variance by delta method
        desginMat <- cbind(1, x12, z12)
        g_NIE1_alpha_k <- g_NIE2_alpha_k <- g_NDE_alpha_k <- matrix(NA, 2 + num_Z,
            K)
        g_NIE1_gammas <- g_NIE2_gammas <- g_NDE_gammas <- rep(NA, 2 + num_Z)
        for (var in seq_len(2 + num_Z)) {
            g_NIE1_alpha_k[var, ] <- (theta_trans[["beta1"]] + theta_trans[["beta5"]] *
                x2) * theta_trans[["psi_k"]] * diff((1 - Delstar_x12) * theta_trans[["lambda_ik"]] *
                desginMat[, var])
            g_NIE1_gammas[var] <- (theta_trans[["beta1"]] + theta_trans[["beta5"]] *
                x2) * (-1) * diff(Delstar_x12 * (1 - Delstar_x12) * m_x12 * desginMat[,
                var])

            g_NIE2_alpha_k[var, ] <- (theta_trans[["beta2"]] + theta_trans[["beta4"]] *
                x2) * theta_trans[["psi_k"]] * diff((1 - Delstar_x12) * exp(-theta_trans[["lambda_ik"]]) *
                theta_trans[["lambda_ik"]] * desginMat[, var])
            g_NIE2_gammas[var] <- (theta_trans[["beta2"]] + theta_trans[["beta4"]] *
                x2) * (-1) * diff(Delstar_x12 * (1 - Del_x12) * desginMat[, var])

            g_NDE_alpha_k[var, ] <- diff(x12) * (1 - Delstar_x12[1]) * theta_trans[["psi_k"]] *
                (theta_trans[["beta4"]] * exp(-theta_trans[["lambda_ik"]][1, ]) +
                  theta_trans[["beta5"]]) * theta_trans[["lambda_ik"]][1, ] * desginMat[1,
                var]
            g_NDE_gammas[var] <- diff(x12) * (-1) * desginMat[1, var] * Delstar_x12[1] *
                (theta_trans[["beta4"]] * (1 - Del_x12[1]) + theta_trans[["beta5"]] *
                  (1 - Delstar_x12[1]) * m_x12[1])
        }
        # beta1, (beta5), alpha0_k, alpha1_k, psi_k-1, gammas
        g_NIE1 <- c(0, diff((1 - Delstar_x12) * m_x12), rep(0, 2 + XMint[1]), if (XMint[2]) x2 *
            diff((1 - Delstar_x12) * m_x12) else NULL, rep(0, num_Z + 1), g_NIE1_alpha_k,
            if (K == 1) NULL else (theta_trans[["beta1"]] + theta_trans[["beta5"]] *
                x2) * diff((theta_trans[["lambda_ik"]][, 1:(K - 1)] - theta_trans[["lambda_ik"]][,
                K]) * (1 - Delstar_x12)), g_NIE1_gammas, 0)

        # beta2, (beta4), alpha0_k, alpha1_k, psi_k-1, gammas
        g_NIE2 <- c(rep(0, 2), -diff(Del_x12), 0, if (XMint[1]) x2 * (-diff(Del_x12)) else NULL,
            rep(0, XMint[2] + num_Z + 1), g_NIE2_alpha_k, if (K == 1) NULL else (theta_trans[["beta2"]] +
                theta_trans[["beta4"]] * x2) * (-1) * diff((1 - Delstar_x12) * (exp(-theta_trans[["lambda_ik"]])[,
                1:(K - 1)] - exp(-theta_trans[["lambda_ik"]])[, K])), g_NIE2_gammas,
            0)

        g_NIE <- g_NIE1 + g_NIE2
        # beta3, (beta4,beta5), alphas_k, psi_k-1, gammas
        g_NDE <- c(rep(0, 3), diff(x12), if (XMint[1]) diff(x12) * (1 - Del_x12[1]) else NULL,
            if (XMint[2]) diff(x12) * (1 - Delstar_x12[1]) * m_x12[1] else NULL, rep(0,
                num_Z + 1), g_NDE_alpha_k, if (K == 1) NULL else diff(x12) * (1 -
                Delstar_x12[1]) * (-theta_trans[["beta4"]] * (exp(-theta_trans[["lambda_ik"]][1,
                1:(K - 1)]) - exp(-theta_trans[["lambda_ik"]][1, K])) + theta_trans[["beta5"]] *
                (theta_trans[["lambda_ik"]][1, 1:(K - 1)] - theta_trans[["lambda_ik"]][1,
                  K])), g_NDE_gammas, 0)
        # beta3, (beta4,beta5)
        g_CDE <- c(rep(0, 3), diff(x12), if (XMint[1]) diff(x12) * (mval > 0) else NULL,
            if (XMint[2]) diff(x12) * mval else NULL, rep(0, num_Z + 1 + (2 + num_Z) *
                K + (K - 1) + (2 + num_Z) + 1))

        if (Group1) {
            index <- G1_index(dat, K, num_Z, XMint)
            g_NIE1 <- g_NIE1[index]
            g_NIE2 <- g_NIE2[index]
            g_NIE <- g_NIE[index]
            g_NDE <- g_NDE[index]
            g_CDE <- g_CDE[index]
            NIE2_se <- NA
        } else {
            NIE2_se <- sqrt(c(g_NIE2 %*% vcovar %*% g_NIE2))
        }
        NIE1_se <- sqrt(c(g_NIE1 %*% vcovar %*% g_NIE1))
        NIE_se <- sqrt(c(g_NIE %*% vcovar %*% g_NIE))
        NDE_se <- sqrt(c(g_NDE %*% vcovar %*% g_NDE))
        CDE_se <- sqrt(c(g_CDE %*% vcovar %*% g_CDE))

        out$eff_se <- c(NIE1_se, NIE2_se, NIE_se, NDE_se, CDE_se)
    }
    return(out)
}
