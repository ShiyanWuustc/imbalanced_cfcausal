## CV+ for counterfactuals. See ?conformalCf
conformalCf_CV <- function(X, Y,
                           estimand,
                           type, side,
                           quantiles,
                           outfun, outparams,
                           psfun, psparams,
                           nfolds,
                           ps_resample_method = "no",
                           ps_resample_seed = NULL,
                           ps_resample_rho = NULL){

    T <- as.numeric(!is.na(Y))
    inds1 <- which(T == 1)
    inds0 <- which(T == 0)
    n1 <- length(inds1)
    n0 <- length(inds0)
    if (n1 < nfolds){
        stop("Insufficient non-missing data")
    }
    idlist1 <- gen_cv_ids(n1, nfolds, offset = 0)
    idlist0 <- gen_cv_ids(n0, nfolds, offset = 0)
    idlist <- lapply(1:nfolds, function(k){
        c(inds1[idlist1[[k]]], inds0[idlist0[[k]]])
    })

    psparams0 <- psparams

    if (estimand == "unconditional"){
        wtfun <- lapply(1:nfolds, function(k){
            testid <- idlist[[k]]
            Xtrain <- X[-testid, , drop = FALSE]
            Ttrain <- T[-testid]

            # ===== NEW: resample only for propensity score fitting =====
            Xtrain_ps <- Xtrain
            Ttrain_ps <- Ttrain
            if (ps_resample_method != "no") {
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
            # ==========================================================

            psparams <- c(list(Y = Ttrain_ps, X = Xtrain_ps), psparams0)
            function(X){
                ps <- do.call(psfun, c(psparams, list(Xtest = X)))
                1 / ps
            }
        })

        # ===== NEW: also resample for the global wtfun_test fit =====
        X_ps <- X
        T_ps <- T
        if (ps_resample_method != "no") {
            idx_ps_all <- resample_idx_ps(
                method = ps_resample_method,
                X = X_ps,
                T = T_ps,
                rho = ps_resample_rho,
                seed = ps_resample_seed
            )
            X_ps <- X_ps[idx_ps_all, , drop = FALSE]
            T_ps <- T_ps[idx_ps_all]
        }
        psparams <- c(list(Y = T_ps, X = X_ps), psparams0)
        # ==========================================================

        wtfun_test <- function(X){
            ps <- do.call(psfun, c(psparams, list(Xtest = X)))
            1 / ps
        }

    } else if (estimand == "nonmissing"){
        wtfun_test <- function(X){
            rep(1, nrow(X))
        }
        wtfun <- lapply(1:nfolds, function(k){
            wtfun_test
        })

    } else if (estimand == "missing"){
        wtfun <- lapply(1:nfolds, function(k){
            testid <- idlist[[k]]
            Xtrain <- X[-testid, , drop = FALSE]
            Ttrain <- T[-testid]

            # ===== NEW: resample only for propensity score fitting =====
            Xtrain_ps <- Xtrain
            Ttrain_ps <- Ttrain
            if (ps_resample_method != "no") {
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
            # ==========================================================

            psparams <- c(list(Y = Ttrain_ps, X = Xtrain_ps), psparams0)
            function(X){
                ps <- do.call(psfun, c(psparams, list(Xtest = X)))
                (1 - ps) / ps
            }
        })

        # ===== NEW: also resample for the global wtfun_test fit =====
        X_ps <- X
        T_ps <- T
        if (ps_resample_method != "no") {
            idx_ps_all <- resample_idx_ps(
                method = ps_resample_method,
                X = X_ps,
                T = T_ps,
                rho = ps_resample_rho,
                seed = ps_resample_seed
            )
            X_ps <- X_ps[idx_ps_all, , drop = FALSE]
            T_ps <- T_ps[idx_ps_all]
        }
        psparams <- c(list(Y = T_ps, X = X_ps), psparams0)
        # ==========================================================

        wtfun_test <- function(X){
            ps <- do.call(psfun, c(psparams, list(Xtest = X)))
            (1 - ps) / ps
        }
    }

    X <- X[inds1, , drop = FALSE]
    Y <- Y[inds1]
    res <- conformalCV(X, Y,
                       type, side,
                       quantiles,
                       outfun, outparams,
                       wtfun,
                       nfolds, idlist1)
    res$wtfun <- wtfun_test
    return(res)
}
