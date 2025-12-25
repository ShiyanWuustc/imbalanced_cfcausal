## Split conformal inference for counterfactuals. See ?conformalCf
conformalCf_split <- function(X, Y,
                              estimand,
                              type, side,
                              quantiles,
                              outfun, outparams,
                              psfun, psparams,
                              trainprop,
                              ps_resample_method = "no",
                              ps_resample_seed = NULL,
                              ps_resample_rho = NULL){

    T <- as.numeric(!is.na(Y))
    inds1 <- which(T == 1)
    inds0 <- which(T == 0)
    n1 <- length(inds1)
    n0 <- length(inds0)

    trainid1 <- sample(n1, floor(n1 * trainprop))
    trainid0 <- sample(n0, floor(n0 * trainprop))
    trainid <- c(inds1[trainid1], inds0[trainid0])

    Xtrain <- X[trainid, , drop = FALSE]
    Ttrain <- T[trainid]

    # ===== NEW: resample only for propensity score fitting (psfun) =====
    Xtrain_ps <- Xtrain
    Ttrain_ps <- Ttrain

    # only needed when ps is used (unconditional or missing)
    if (estimand %in% c("unconditional", "missing") && ps_resample_method != "no") {
        idx_ps <- resample_idx_ps(
            method = ps_resample_method,
            X = Xtrain_ps,
            T = Ttrain_ps,
            rho = ps_resample_rho,
            seed = ps_resample_seed
        )
        Xtrain_ps <- Xtrain_ps[idx_ps, , drop = FALSE]
        Ttrain_ps <- Ttrain_ps[idx_ps]
    }
    # ================================================================

    psparams0 <- psparams
    if (estimand == "unconditional"){
        psparams <- c(list(Y = Ttrain_ps, X = Xtrain_ps), psparams0)
        wtfun <- function(X){
            ps <- do.call(psfun, c(list(Xtest = X), psparams))
            1 / ps
        }
    } else if (estimand == "nonmissing"){
        wtfun <- function(X){
            rep(1, nrow(X))
        }
    } else if (estimand == "missing"){
        psparams <- c(list(Y = Ttrain_ps, X = Xtrain_ps), psparams0)
        wtfun <- function(X){
            ps <- do.call(psfun, c(list(Xtest = X), psparams))
            (1 - ps) / ps
        }
    }

    X <- X[inds1, , drop = FALSE]
    Y <- Y[inds1]
    res <- conformalSplit(X, Y,
                          type, side,
                          quantiles,
                          outfun, outparams,
                          wtfun,
                          trainprop, trainid1)
    return(res)
}
