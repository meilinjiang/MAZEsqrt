# Mediation analysis for zero-inflated negative binomial mediators with
# confounders
#' @noRd
trans.zinbm <- function(dat, theta, K, xval, num_Z, zval) {
    beta0 <- theta[1]
    beta1 <- theta[2]
    beta2 <- theta[3]
    beta3 <- theta[4]
    beta4 <- theta[5]
    beta5 <- theta[6]
    delta <- theta[7 + num_Z]
    alpha_k <- theta[7 + num_Z + 1:((2 + num_Z) * K)]
    r <- theta[8 + num_Z + (2 + num_Z) * K]
    if (K == 1) {
        psi_k <- 1
    } else {
        psi_k <- c(theta[8 + num_Z + (2 + num_Z) * K + 1:(K - 1)], 1 - sum(theta[8 +
            num_Z + (2 + num_Z) * K + 1:(K - 1)]))
    }
    gammas <- theta[8 + num_Z + (2 + num_Z) * K + (K - 1) + 1:(2 + num_Z)]
    eta <- theta[9 + num_Z + (2 + num_Z) * K + (K - 1) + (2 + num_Z)]

    n0 <- length(xval)
    psi_ik <- matrix(k_to_ik(psi_k, n0), ncol = K)
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

    logmu_ik <- matrix(apply(matrix(alpha_k, ncol = K), 2, function(t) rowSums(k_to_ik(t,
        n0) * designMat_M)), ncol = K)
    mu_ik <- exp(logmu_ik)

    Delstar_i <- expit(rowSums(k_to_ik(gammas, n0) * designMat_M))
    p_ik <- r/(r + mu_ik)
    Del_i <- Delstar_i + (1 - Delstar_i) * rowSums(psi_ik * p_ik^r)

    theta_trans <- list(beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3,
        beta4 = beta4, beta5 = beta5, delta = delta, alpha_k = alpha_k, r = r, psi_k = psi_k,
        gammas = gammas, eta = eta, beta_T_Z = beta_T_Z, mu_ik = mu_ik, Delstar_i = Delstar_i,
        p_ik = p_ik, psi_ik = psi_ik, Del_i = Del_i)
    return(theta_trans)
}
#' @noRd
ComputeInit.zinbm <- function(dat, K, num_Z, Z_names, XMint) {
    # group 1, flexmix
    dat_g1 <- dat[which(dat$Mobs > 0), ]

    # fit <- summary(lm(Y_group1 ~ M_group1 + X_group1))
    fm_Y_rhs <- c("Mobs", "X", if (XMint[2]) "Mobs:X" else NULL, Z_names)
    fm_Y <- as.formula(paste0("Y~", paste0(fm_Y_rhs, collapse = "+"), collapse = ""))
    fit <- summary(lm(fm_Y, data = dat_g1))
    # beta02, beta1, beta34, (beta5), beta_Z
    betas <- fit$coefficients[c("(Intercept)", fm_Y_rhs), "Estimate"]
    delta <- fit$sigma

    # fm_M <- as.formula(paste0('Mobs~', paste0(c('X', Z_names), collapse =
    # '+'), collapse = '')) cl <- flexmix(fm_M, data = dat_g1, k=K, model =
    # FLXMRnegbin())
    cl <- kmeans(dat_g1$Mobs, centers = K)
    para <- NULL
    # for (j in seq_len(K)) { fit2 <- glm.nb(M_group1[cl$cluster==j] ~
    # X_group1[cl$cluster==j], control = glm.control(maxit=1e15)) para <-
    # cbind(para, c(coef(fit2), fit2$theta )) }
    fm_M <- as.formula(paste0("log(Mobs)~", paste0(c("X", Z_names), collapse = "+"),
        collapse = ""))
    for (j in seq_len(K)) {
        dat_g1_j <- dat_g1[cl$cluster == j, ]
        m_j <- dat_g1_j$Mobs
        # r_init <- mean(m_j)^2/(var(m_j)-mean(m_j)) # numerator larger,
        # denominator smaller than true
        r_init <- mean(m_j)/2
        fit2 <- lm(fm_M, dat_g1_j)
        para <- cbind(para, c(coef(fit2), r_init))
    }

    # method 3: distribution 1 as the one with smaller cluster mean since x=age,
    # intercepts are not accurate enough
    mu <- cl$centers
    # mu <- parameters(cl)[1,]+parameters(cl)[2,]*mean(X_group1)
    ord <- order(mu)
    mix <- para[, ord, drop = F]
    # alpha_0k, alpha_1k, r, psik
    mix_init <- c(mix[1:(2 + num_Z), ], mean(mix[3 + num_Z, ]))
    psi_k <- cl$size[ord]/nrow(dat_g1)
    mix_init <- c(mix_init, psi_k[-K])

    init <- c(betas, delta, mix_init)
    return(init)
}
#' @noRd
logbinom <- function(r, M_i) {
    # sapply(M_i, function(m){lgamma(r+m) - lgamma(r) - lgamma(m+1)})
    lgamma(r + M_i) - lgamma(r) - lgamma(M_i + 1)
}

# negative expectation of log-likelihood function with respect to conditional
# distribution of 1(Ci = k) given data and current estimates for group 1: -Q1
#' @noRd
negQ_G1.zinbm <- function(dat, theta, K, num_Z, Z_names, B = NULL, tauG1 = NULL, calculate_tau = F,
    calculate_ll = F) {
    # betas, delta, alpha0k, alpha1k, alphaZ_k, r, w0k, gammas, eta
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
    log_dnb_nz <- logbinom(theta_trans[["r"]], M_group1) + M_group1 * log(1 - theta_trans[["p_ik"]]) -
        log(theta_trans[["p_ik"]]^(-theta_trans[["r"]]) - 1)
    # group 1
    if (!calculate_tau) {
        # if(K == 1){ calculate_ll <- TRUE }
        if (theta_trans[["eta"]] == Inf) {
            pf0 <- rep(0, length(M_group1))
        } else {
            pf0 <- exp(-theta_trans[["eta"]]^2 * sqrt(M_group1))
        }
        pf0[M_group1 > B] <- 0
        l1_ik <- -log(theta_trans[["delta"]]) - 0.5 * log(2 * pi) + log(1 - pf0) -
            (Y_group1 - theta_trans[["beta0"]] - theta_trans[["beta1"]] * M_group1 -
                theta_trans[["beta2"]] - (theta_trans[["beta3"]] + theta_trans[["beta4"]]) *
                X_group1 - theta_trans[["beta5"]] * X_group1 * M_group1 - theta_trans[["beta_T_Z"]])^2/(2 *
                theta_trans[["delta"]]^2) + log_dnb_nz

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
            out <- rep(1, length(M_group1))
        } else {
            num <- theta_trans[["psi_ik"]] * exp(log_dnb_nz)
            out <- num/rowSums(num)
        }
    }
    return(out)
}
#' @noRd
bounds.zinbm <- function(dat, K, group, num_Z, XMint) {
    # r need to be > 0
    s <- ifelse(group == 1, 3, 4 + XMint[1]) + XMint[2] + num_Z
    ul <- 1000
    if (K == 1) {
        if (group == 1) {
            ui <- diag(rep(1, s + 2 + num_Z + 2))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, 2 + num_Z), 1e-06)
        } else {
            ui <- diag(rep(1, s + 2 + num_Z + 2 + 2 + num_Z + 1))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, 2 + num_Z), 1e-06, rep(-1000,
                2 + num_Z), 1e-06)
        }
    } else {
        if (group == 1) {
            ui1 <- diag(rep(1, s + 2 + (2 + num_Z) * K + (K - 1)))
            ui2 <- diag(c(rep(0, s + 2 + (2 + num_Z) * K), rep(-1, K - 1)))
            ui <- rbind(ui1, ui2[s + 2 + (2 + num_Z) * K + 1:(K - 1), ], c(rep(0,
                s + 2 + (2 + num_Z) * K), rep(-1, K - 1)))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, (2 + num_Z) * K), 1e-06, rep(1e-06,
                K - 1), rep(-(1 - 1e-06), K - 1), -(1 - 1e-06))
        } else {
            ui1 <- diag(rep(1, s + 3 + (2 + num_Z) * K + (K - 1) + 2 + num_Z))
            ui2 <- diag(c(rep(0, s + 2 + (2 + num_Z) * K), rep(-1, K - 1), rep(0,
                (2 + num_Z) + 1)))
            ui <- rbind(ui1, ui2[s + 2 + (2 + num_Z) * K + 1:(K - 1), ], c(rep(0,
                s + 2 + (2 + num_Z) * K), rep(-1, K - 1), rep(0, (2 + num_Z) + 1)))
            ci <- c(rep(-1000, s), 1e-06, rep(-1000, (2 + num_Z) * K), 1e-06, rep(1e-06,
                K - 1), rep(-1000, 2 + num_Z), 1e-06, rep(-(1 - 1e-06), K - 1), -(1 -
                1e-06))
        }
    }
    return(list(ui = ui, ci = ci))
}
#' @noRd
ComputeInit2.zinbm <- function(dat, K, num_Z, Z_names, XMint, B, limits, explicit = T) {
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
            B, tau), grad = NULL, ui = bd$ui, ci = bd$ci, control = list(maxit = 5000))
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

    # beta02, beta1, beta34, (beta5), beta_Z, delta, alpha_0k, alpha_1k, r, w_0k
    initials <- G1_init(dat, init2, K, num_Z, Z_names, XMint, initials_for_full = T)
    s <- 4 + XMint[1] + XMint[2] + num_Z
    initials[s + 1 + (2 + num_Z) * K + 1] <- initials[s + 1 + (2 + num_Z) * K + 1]/2
    return(initials)
}

# log(factorial(m)) = sum(log(m)) more stable using log
#' @noRd
loghik_zinbm <- function(m, r, p, beta0, beta1, beta2, beta3, beta4, beta5, beta_T_Z,
    delta, eta, x, y) {
    logbinom(r, m) + m * log(1 - p) - ((y - beta0 - beta1 * m - beta2 - (beta3 + beta4) *
        x - beta5 * x * m - beta_T_Z)^2)/(2 * delta^2) - eta^2 * sqrt(m)
}

# negative expectation of log-likelihood function with respect to conditional
# distribution of 1(Ci = k) given data and current estimates for group 2: -Q2
#' @noRd
negQ_G2.zinbm <- function(dat, theta, K, num_Z, Z_names, B, tauG2 = NULL, calculate_tau = F,
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
        loghik_m[[i]] <- loghik_zinbm(m = i, theta_trans[["r"]], theta_trans[["p_ik"]],
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

    l2_ik <- cbind(-(Y_group2 - theta_trans[["beta0"]] - theta_trans[["beta3"]] *
        X_group2 - theta_trans[["beta_T_Z"]])^2/(2 * theta_trans[["delta"]]^2), -log(theta_trans[["p_ik"]]^(-theta_trans[["r"]]) -
        1) + output) - log(theta_trans[["delta"]]) - 0.5 * log(2 * pi)

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
        # estimates
        ll <- exp(pow)
        out <- ll/rowSums(ll)
    }
    return(out)
}
#' @noRd
f_NIE2_zinbm <- function(dat, theta, K, num_Z, XMint, x12, z12) {
    # theta <- G12_init(theta, XMint)
    theta_trans <- trans.zinbm(dat, theta, K, xval = x12, num_Z, zval = z12)
    NIE2 <- (theta_trans[["beta2"]] + theta_trans[["beta4"]] * x12[2]) * (-diff(theta_trans[["Del_i"]]))
    return(NIE2)
}
#' @noRd
f_NDE_zinbm <- function(dat, theta, K, num_Z, XMint, x12, z12) {
    # theta <- G12_init(theta, XMint)
    theta_trans <- trans.zinbm(dat, theta, K, xval = x12, num_Z, zval = z12)
    NDE <- diff(x12) * (theta_trans[["beta3"]] + theta_trans[["beta4"]] * (1 - theta_trans[["Del_i"]][1]) +
        theta_trans[["beta5"]] * (1 - theta_trans[["Delstar_i"]][1]) * rowSums(theta_trans[["psi_ik"]] *
            theta_trans[["mu_ik"]])[1])
    return(NDE)
}
#' @noRd
effects.zinbm <- function(dat, theta, x1, x2, K, num_Z, zval, mval, XMint, calculate_se = F,
    vcovar = NULL, Group1 = F) {
    if (Group1) {
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
    m_x12 <- rowSums(theta_trans[["psi_ik"]] * theta_trans[["mu_ik"]])

    NIE1 <- (theta_trans[["beta1"]] + theta_trans[["beta5"]] * x2) * diff((1 - Delstar_x12) *
        m_x12)
    NIE2 <- (theta_trans[["beta2"]] + theta_trans[["beta4"]] * x2) * (-diff(Del_x12))
    # f_NIE2_zinbm(dat, theta, K, num_Z, XMint, x12, z12)
    NIE <- NIE1 + NIE2
    NDE <- diff(x12) * (theta_trans[["beta3"]] + theta_trans[["beta4"]] * (1 - Del_x12[1]) +
        theta_trans[["beta5"]] * (1 - Delstar_x12[1]) * m_x12[1])
    # f_NDE_zinbm(dat, theta, K, num_Z, XMint, x12, z12)
    CDE <- diff(x12) * (theta_trans[["beta3"]] + theta_trans[["beta4"]] * (mval >
        0) + theta_trans[["beta5"]] * mval)
    out <- data.frame(eff = c(NIE1, NIE2, NIE, NDE, CDE))

    if (calculate_se == T) {
        # asymptotic variance by delta method
        desginMat <- cbind(1, x12, z12)
        g_NIE1_alpha_k <- matrix(NA, 2 + num_Z, K)
        g_NIE1_gammas <- g_NIE2_gammas <- rep(NA, 2 + num_Z)
        for (var in seq_len(2 + num_Z)) {
            g_NIE1_alpha_k[var, ] <- (theta_trans[["beta1"]] + theta_trans[["beta5"]] *
                x2) * theta_trans[["psi_k"]] * diff((1 - Delstar_x12) * theta_trans[["mu_ik"]] *
                desginMat[, var])
            g_NIE1_gammas[var] <- (theta_trans[["beta1"]] + theta_trans[["beta5"]] *
                x2) * (-1) * diff(Delstar_x12 * (1 - Delstar_x12) * m_x12 * desginMat[,
                var])
        }
        # beta1, (beta5), alpha0_k, alpha1_k, psi_k-1, gammas
        g_NIE1 <- c(0, diff((1 - Delstar_x12) * m_x12), rep(0, 2 + XMint[1]), if (XMint[2]) x2 *
            diff((1 - Delstar_x12) * m_x12) else NULL, rep(0, num_Z + 1), g_NIE1_alpha_k,
            0, if (K == 1) NULL else (theta_trans[["beta1"]] + theta_trans[["beta5"]] *
                x2) * diff((theta_trans[["mu_ik"]][, 1:(K - 1)] - theta_trans[["mu_ik"]][,
                K]) * (1 - Delstar_x12)), g_NIE1_gammas, 0)

        # beta2, (beta4), alpha0_k, alpha1_k, r, psi_k-1, gammas
        g_NIE2 <- numDeriv::grad(function(x) f_NIE2_zinbm(dat, x, K, num_Z, XMint,
            x12, z12), theta)
        # beta3, (beta4,beta5), alphas_k, r, psi_k-1, gammas
        g_NDE <- numDeriv::grad(function(x) f_NDE_zinbm(dat, x, K, num_Z, XMint, x12,
            z12), theta)
        if (sum(XMint) != 2) {
            g_NIE2 <- g_NIE2[-c(5, 6)[!XMint]]
            g_NDE <- g_NDE[-c(5, 6)[!XMint]]
        }

        g_NIE <- g_NIE1 + g_NIE2

        # beta3, (beta4,beta5)
        g_CDE <- c(rep(0, 3), diff(x12), if (XMint[1]) diff(x12) * (mval > 0) else NULL,
            if (XMint[2]) diff(x12) * mval else NULL, rep(0, num_Z + 1 + (2 + num_Z) *
                K + 1 + (K - 1) + (2 + num_Z) + 1))

        NIE1_se <- sqrt(c(g_NIE1 %*% vcovar %*% g_NIE1))
        NIE2_se <- sqrt(c(g_NIE2 %*% vcovar %*% g_NIE2))
        NIE_se <- sqrt(c(g_NIE %*% vcovar %*% g_NIE))
        NDE_se <- sqrt(c(g_NDE %*% vcovar %*% g_NDE))
        CDE_se <- sqrt(c(g_CDE %*% vcovar %*% g_CDE))

        out$eff_se <- c(NIE1_se, NIE2_se, NIE_se, NDE_se, CDE_se)
    }
    return(out)
}
