# source.R common functions among three distributions

gammas_init <- function(dat, num_Z = 0, Z_names = NULL) {
    if (sum(dat$Mobs == 0) == 0) {
        out <- c(-Inf, rep(0, num_Z + 1))
    } else {
        dat$ind_Del <- (dat$Mobs == 0) * 1
        fm <- as.formula(paste0("ind_Del~", paste0(c("X", Z_names), collapse = "+"),
            collapse = ""))
        logitmodel <- summary(glm(fm, family = "binomial", data = dat))
        out <- logitmodel$coefficients[, "Estimate"]
    }
    return(out)
}

trans <- function(dat, theta, K, xval, ...) {
    UseMethod("trans")
}

ComputeInit <- function(dat, K, ...) {
    UseMethod("ComputeInit")
}

# negative expectation of log-likelihood function with respect to conditional
# distribution of 1(Ci = k) given data and current estimates for group 1: -Q1
negQ_G1 <- function(dat, theta, K, ...) {
    UseMethod("negQ_G1")
}

G1_init <- function(dat, init, K, num_Z, Z_names, XMint, initials_for_full = F) {
    gammas <- gammas_init(dat, num_Z, Z_names)
    # from G1 initials to
    if (initials_for_full) {
        # to set initials for parameters in G2 for full analysis
        initials <- c(init[1:2], 0, init[3], rep(0, XMint[1] * 1), init[4:length(init)],
            gammas, 0.01)
    } else {
        # to set parameters not estimating in G1 at fixed values
        initials <- c(init[1:2], 0, init[3], rep(0, ifelse(XMint[2], 1, 2)), init[4:length(init)],
            gammas, 0.01)
    }
    # zilonm, with beta5 initials <- c(init[1:2], 0, init[3], 0, init[4], init[4
    # + 1:(num_Z + (2 + num_Z) * K + 2 + K - 1)], gammas, 0.01)

    # zinbm, no beta5 initials <- c(init[1:2], 0, init[3], 0, init[3 + 1:(num_Z
    # + (2 + num_Z) * K + 2 + K - 1)], gammas, 0.01)

    # zipm, no beta5 initials <- c(init[1:2], 0, init[3], 0, init[3 + 1:(num_Z +
    # (2 + num_Z) * K + 1 + K - 1)], gammas, 0.01)
    return(initials)
}

# fix parameter not estimating in group 1 analysis
negQ2_G1 <- function(dat, init, K, num_Z, Z_names, XMint, B = NULL, tauG1 = NULL,
    calculate_tau = F, calculate_ll = F) {
    initials <- G1_init(dat, init, K, num_Z, Z_names, XMint)
    out <- negQ_G1(dat, theta = initials, K, num_Z, Z_names, B, tauG1, calculate_tau,
        calculate_ll)
    return(out)
}

compare_mu <- function(dat, init2, group, K, num_Z, Z_names, XMint, B = NULL) {
    if (group == 1) {
        initials <- G1_init(dat, init2, K, num_Z, Z_names, XMint)
        s <- 3 + XMint[2] + num_Z
        dat2 <- dat[which(dat$Mobs > 0), ]
    } else {
        initials <- init2
        s <- 4 + XMint[1] + XMint[2] + num_Z
        dat2 <- dat
    }
    if (num_Z == 0) {
        zval <- NULL
    } else {
        zval <- colMeans(dat2[, Z_names, drop = F])
    }
    theta_trans <- trans(dat, initials, K, xval = mean(dat2$X), num_Z, zval)

    if ("zipm" %in% class(dat)) {
        mean_name <- "lambda_ik"
        n_M2nd_param <- 0
    } else {
        mean_name <- "mu_ik"
        n_M2nd_param <- 1
    }
    ord <- order(theta_trans[[mean_name]])
    if (sum(ord != 1:K) > 0) {
        init2[s + 1 + 1:((2 + num_Z) * K)] <- c(matrix(theta_trans[["alpha_k"]], ncol = K)[,
            ord])
        init2[s + 1 + (2 + num_Z) * K + n_M2nd_param + 1:(K - 1)] <- theta_trans[["psi_k"]][ord][-K]
    }
    if (group == 1) {
        tau2 <- negQ2_G1(dat, init2, K, num_Z, Z_names, XMint, calculate_tau = T)
    } else {
        tau2 <- negQ(dat, init2, K, num_Z, Z_names, XMint, B, calculate_tau = T)
    }
    return(list(init2 = init2, tau2 = tau2))
}

bounds <- function(dat, K, group, num_Z, XMint) {
    UseMethod("bounds")
}

ComputeInit2 <- function(dat, K, B, ...) {
    UseMethod("ComputeInit2")
}

# negative expectation of log-likelihood function with respect to conditional
# distribution of 1(Ci = k) given data and current estimates for group 2: -Q2
negQ_G2 <- function(dat, theta, K, ...) {
    UseMethod("negQ_G2")
}

# fix parameter not estimating in full analysis (group12)
G12_init <- function(theta, XMint) {
    theta <- c(theta[1:4], rep(0, !XMint[1]), theta[5:length(theta)])
    theta <- c(theta[1:5], rep(0, !XMint[2]), theta[6:length(theta)])
    return(theta)
}

# Q for all
negQ <- function(dat, theta, K, num_Z, Z_names, XMint, B, tau = NULL, calculate_tau = F,
    calculate_ll = F) {
    if (sum(dat$Mobs == 0) > 0) {
        # fix parameter not estimating in full analysis (group12)
        theta <- G12_init(theta, XMint)
        if (!calculate_tau) {
            if (!calculate_ll) {
                out <- negQ_G1(dat, theta, K, num_Z, Z_names, B, tauG1 = tau$tauG1,
                  calculate_tau, calculate_ll) + negQ_G2(dat, theta, K, num_Z, Z_names,
                  B, tauG2 = tau$tauG2, calculate_tau, calculate_ll)
            } else {
                out <- negQ_G1(dat, theta, K, num_Z, Z_names, B, tauG1 = NULL, calculate_tau,
                  calculate_ll) + negQ_G2(dat, theta, K, num_Z, Z_names, B, tauG2 = NULL,
                  calculate_tau, calculate_ll)
            }
            # print(theta);print(out)
        } else {
            out <- list(tauG1 = negQ_G1(dat, theta, K, num_Z, Z_names, B = NULL, tauG1 = NULL,
                calculate_tau, calculate_ll), tauG2 = negQ_G2(dat, theta, K, num_Z,
                Z_names, B, tauG2 = NULL, calculate_tau, calculate_ll))
        }
        # if (is.nan(out)) { print(theta);print(out) }
    } else {
        if (!calculate_tau) {
            if (!calculate_ll) {
                out <- negQ2_G1(dat, theta, K, num_Z, Z_names, XMint, B, tauG1 = tau$tauG1,
                  calculate_tau, calculate_ll)
            } else {
                out <- negQ2_G1(dat, theta, K, num_Z, Z_names, XMint, B, tauG1 = NULL,
                  calculate_tau, calculate_ll)
            }
            # print(theta);print(out)
        } else {
            out <- list(tauG1 = negQ2_G1(dat, theta, K, num_Z, Z_names, XMint, B = NULL,
                tauG1 = NULL, calculate_tau, calculate_ll))
        }
    }
    return(out)
}

# for calculating observed I (EM approximation)
g <- function(dat, x, y, K, num_Z, Z_names, XMint, B) {
    tau_f <- negQ(dat, y, K, num_Z, Z_names, XMint, B, calculate_tau = T)
    out <- negQ(dat, x, K, num_Z, Z_names, XMint, B, tau_f)
    return(out)
}

effects <- function(dat, theta, x1, x2, K, ...) {
    UseMethod("effects")
}

parameter_names <- function(dat, K, num_Z, Z_names, XMint) {
    variation <- if ("zipm" %in% class(dat))
        NULL else ifelse("zilonm" %in% class(dat), "xi0", "r")
    out <- c(paste0("beta", c(0:5, Z_names)), "delta", paste0("alpha", rep(1:K, each = 2 +
        num_Z), rep(c(0:1, Z_names), K)), variation, paste0("psi", 1:K)[-K], paste0("gamma",
        c(0:1, Z_names)), "eta")
    if (sum(XMint) != 2) {
        out <- out[-c(5, 6)[!XMint]]
    }
    return(out)
}

G1_index <- function(dat, K, num_Z, XMint) {
    c(1:2, 4, (4 + XMint[1] + XMint[2]) * XMint[2], 4 + XMint[1] + XMint[2] + 1:(num_Z +
        1 + (2 + num_Z) * K + ifelse("zipm" %in% class(dat), 0, 1) + (K - 1)))
}

analysis <- function(dat, K, num_Z, Z_names, XMint, B, limits, seed) {
    cinit <- init2 <- ComputeInit2(dat, K, num_Z, Z_names, XMint, B, limits, explicit = T)
    # cinit <- init2 <- true2
    if (sum(dat$Mobs == 0) == 0) {
        index <- G1_index(dat, K, num_Z, XMint)
        cinit <- init2 <- init2[index]
        tau2 <- negQ(dat, init2, K, num_Z, Z_names, XMint, B, calculate_tau = T)
        p <- length(init2)
        countEM <- NA
    } else {
        tau2 <- negQ(dat, init2, K, num_Z, Z_names, XMint, B, calculate_tau = T)
        p <- length(init2)
        init <- rep(0, p)
        countEM <- 0
        bd <- bounds(dat, K, group = 2, num_Z, XMint)
        # M-step
        while (euclidean_dist(init, init2) > limits) {
            init <- init2
            tau <- tau2

            # m2 <- nlminb(init,function(x)negQ(dat,x,K,num_Z,Z_names,
            # XMint,B,tau), lower = c(rep(-1000,6 + 2*K),rep(0,
            # K-1),-1000,-1000,1e-6,1e-6), upper = c(rep(1000,6 + 2*K),rep(1,
            # K-1),rep(1000,4)), control = list(eval.max=5000,iter.max=5000))
            m <- constrOptim(init, function(x) negQ(dat, x, K, num_Z, Z_names, XMint,
                B, tau), grad = NULL, ui = bd$ui, ci = bd$ci, outer.iterations = 500,
                control = list(maxit = 50000))
            # negQ(dat,true,K,num_Z,Z_names,
            # XMint,B,negQ(dat,true,K,num_Z,Z_names, XMint,B,calculate_tau=T))

            if (m$convergence == 0) {
                init2 <- m$par
            } else {
                print(paste0("dist = ", class(dat)[1], ", K = ", K, ", seed = ", seed,
                  ", analysis()_constrOptim() = ", m$convergence, ", countEM = ",
                  countEM))
                print(m)
                init2 <- rep(NA, length(init))
                break
            }
            switch <- compare_mu(dat, init2, group = 2, K, num_Z, Z_names, XMint,
                B)
            init2 <- switch$init2
            tau2 <- switch$tau2
            countEM <- countEM + 1
            # print(countEM)
        }
    }
    est <- init2
    n <- dim(dat)[1]
    negll <- negQ(dat, est, K, num_Z, Z_names, XMint, B, calculate_ll = T)
    AIC <- 2 * negll + 2 * p
    BIC <- 2 * negll + p * log(n)
    out <- list(init = cinit, est = est, countEM = countEM, AIC = AIC, BIC = BIC,
        tau2 = tau2)
    return(out)
}

analysis2 <- function(dat, K, B, est, tau2, x1, x2, num_Z, Z_names, XMint, zval, mval) {
    # method 1: hessian from log-likelihood function, parameter estimation
    # cannot do this
    hess <- pracma::hessian(function(x) negQ(dat, x, K, num_Z, Z_names, XMint, B,
        calculate_ll = T), x = est)
    vcovar <- solve(hess)
    var <- diag(vcovar)
    # if (sum(var < 0) > 0) { hess <- numDeriv::hessian(function(x) negQ(dat, x,
    # K, num_Z, Z_names, XMint, B, calculate_ll = T), x = est) vcovar <-
    # solve(hess) var <- diag(vcovar)
    if (sum(var < 0) > 0) {
        # method 2: hessian matrix from EM approximation
        h.1 <- pracma::hessian(function(x) negQ(dat, x, K, num_Z, Z_names, XMint,
            B, tau2), x = est)
        # grad() function in numDeriv
        g2 <- function(y) {
            numDeriv::grad(function(x) g(dat, x, y, K, num_Z, Z_names, XMint, B),
                x = est)
        }
        # jacobian() function in numDeriv
        h.2 <- numDeriv::jacobian(g2, est)
        hess <- h.1 + h.2
        vcovar <- solve(hess)
        var <- diag(vcovar)
        if (sum(var < 0) > 0) {
            print("All methods var < 0.")
        }
    }
    # }
    se <- sqrt(var)
    # mediation effects
    Group1 <- sum(dat$Mobs == 0) == 0
    # conditional effects: E(NIE|Z=zval)
    MedZIM <- effects(dat, theta = est, x1, x2, K, num_Z, zval, mval, XMint, calculate_se = T,
        vcovar, Group1)
    # average effects: E(NIE) = E[ E(NIE|Z=zval) ] = 1/n \sum_{z \in Data}
    # E(NIE|Z=zval) Var(NIE) = E[ Var(NIE|Z=zval) ] + Var[ E(NIE|Z=zval) ]
    if (num_Z == 0) {
        MedZIM_avg <- MedZIM
    } else {
        n <- dim(dat)[1]
        eff_i <- eff_i_var <- NULL
        for (i in seq_len(n)) {
            MedZIM_i <- effects(dat, theta = est, x1, x2, K, num_Z, zval = dat[i,
                Z_names], mval, XMint, calculate_se = T, vcovar, Group1)
            eff_i <- cbind(eff_i, MedZIM_i$eff)  # E(NIE|Z=zval)
            eff_i_var <- cbind(eff_i_var, MedZIM_i$eff_se^2)  # Var(NIE|Z=zval)
        }
        MedZIM_avg <- data.frame(eff = rowMeans(eff_i), eff_se = sqrt(rowMeans(eff_i_var) +
            apply(eff_i, 1, var)))
    }

    MedZIM <- rbind(MedZIM, MedZIM_avg)

    out <- list(se = se, MedZIM = MedZIM)
    return(out)
}

realanalysis <- function(dat, distM_sequence, K_sequence, selection, XMint, x1, x2,
    zval, num_Z, Z_names, mval, limits, B, seed, ncore) {
    set.seed(seed)
    iter <- data.frame(distM = rep(distM_sequence, each = length(K_sequence)), K = rep(K_sequence,
        length(distM_sequence)))
    niter <- nrow(iter)

    # t0 <- Sys.time()
    cl <- parallel::makeCluster(ncore, outfile = "")
    fun <- c("analysis", "bounds", "compare_mu", "ComputeInit", "ComputeInit2", "gammas_init",
        "euclidean_dist", "expit", "k_to_ik", "logbinom", "loghik_zinbm",
         "negQ", "negQ_G1", "negQ_G2", "negQ2_G1", "trans")
    parallel::clusterExport(cl = cl, varlist = c(fun), envir = parent.env(environment()))
    doParallel::registerDoParallel(cl)
    parallel::clusterSetRNGStream(cl = cl, seed)
    i <- numeric()
    models <- foreach(i = seq_len(niter), .multicombine = T, .packages = c("stats",
        "flexmix", "MASS"), .combine = "list", .errorhandling = "pass") %dopar% {
        # set.seed(seed)
        dat2 <- dat
        if (iter$distM[i] != "zilonm") {
            dat2$Mobs <- round(dat2$Mobs)
        }
        class(dat2) <- c(iter$distM[i], class(dat2))
        model_i <- tryCatch({
            analysis(dat2, K = iter$K[i], num_Z, Z_names, XMint, B, limits, seed)
        }, error = function(e) {
            print(paste0("dist = ", iter$distM[i], ", K = ", iter$K[i], ", analysis()_error = ",
                e))
            list(init = NA, est = NA, countEM = NA, AIC = Inf, BIC = Inf, tau2 = NA,
                e = e)
        })
        # model_i <- analysis(dat2, K = iter$K[i], num_Z, Z_names, XMint, B,
        # limits, seed)
        print(paste0("distM = ", iter$distM[i], ", K = ", iter$K[i], " completed."))
        return(model_i)
    }
    parallel::stopCluster(cl)
    if (niter == 1) {
        models <- list(models)
    }
    names(models) <- apply(iter, 1, function(t) paste0(t, collapse = "_K"))
    print(models)

    AICs <- sapply(models, function(x) x[["AIC"]])
    BICs <- sapply(models, function(x) x[["BIC"]])
    selection_criteria <- sapply(models, function(x) x[[selection]])
    selected_model <- models[[which.min(selection_criteria)]]
    selected_model_name <- names(models)[[which.min(selection_criteria)]]
    print(paste0("selected_model_name: ", selected_model_name))
    selected_dist <- strsplit(selected_model_name, "_K")[[1]][1]
    selected_K <- as.numeric(strsplit(selected_model_name, "_K")[[1]][2])
    KBIC <- K_sequence[which.min(BICs)]
    KAIC <- K_sequence[which.min(AICs)]

    class(dat) <- c(selected_dist, class(dat))
    # t1 <- Sys.time()
    analysis2_out <- tryCatch({
        analysis2(dat, K = selected_K, B, est = selected_model$est, tau2 = selected_model$tau2,
            x1, x2, num_Z, Z_names, XMint, zval, mval)
    }, error = function(e) {
        print(paste0("dist = ", selected_dist, ", K = ", selected_K, ", analysis2()_error = ",
            e))
        list(se = NA, MedZIM = data.frame(eff = rep(NA, 5), eff_se = rep(NA, 5)),
            e = e)
    })
    # t2 <- Sys.time() print(paste0('EM = ', format(difftime(t1, t0)), ',
    # hessian = ', format(difftime(t2, t1))))

    results_parameters <- results(est = selected_model$est, se = analysis2_out$se,
        init = selected_model$init)
    paranames <- parameter_names(dat, K = selected_K, num_Z, Z_names, XMint)
    if (sum(dat$Mobs == 0) == 0) {
        index <- G1_index(dat, K = selected_K, num_Z, XMint)
        paranames <- paranames[index]
    }

    rownames(results_parameters) <- paranames

    # mediation and direct effects
    out_MedZIM <- analysis2_out$MedZIM
    results_effects <- results(out_MedZIM$eff, out_MedZIM$eff_se)
    effects_names <- c("NIE1", "NIE2", "NIE", "NDE", "CDE")
    rownames(results_effects) <- c(paste0(effects_names, "_cond"), paste0(effects_names,
        "_avg"))

    out <- list(results_effects = results_effects, results_parameters = results_parameters,
        selected_model_name = selected_model_name, BIC = BICs[which.min(selection_criteria)],
        AIC = AICs[which.min(selection_criteria)], models = models, analysis2_out = analysis2_out)
    return(out)
}
