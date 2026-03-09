
k_to_ik <- function(para_k, n) {
    sapply(para_k, function(x) x * rep(1, n))
}

expit <- function(t) {
    # out <- exp(t)/(1 + exp(t)) out[t > 100] <- 1/(1 + exp(-t)) return(out)
    1/(1 + exp(-t))
}

euclidean_dist <- function(p, q) {
    sqrt(sum((q - p)^2))
}

results <- function(est, se, init = NA, d = 5, ci_lb = NULL, ci_ub = NULL, pval = NULL) {
    # 95% CI
    if (is.null(ci_lb) & is.null(ci_ub)) {
        ci_lb <- est - qnorm(0.975) * se
        ci_ub <- est + qnorm(0.975) * se
    }
    if (is.null(pval)) {
        pval <- (1 - pnorm(abs(est/se))) * 2
        # pval <- ( 1-pt(abs(est/se),df=dim(dat)[1]-1) )*2
    }

    out <- data.frame(Initials = init, Estimate = est, SE = se, CI_lower = ci_lb,
        CI_upper = ci_ub, Pvalue = pval)
    if (is.na(init)[1]) {
        out <- out[, -1]
    }
    # out <- round(out, d)
    return(out)
}

