#' Resample indices for propensity score fitting
#' @keywords internal
resample_idx_ps <- function(method, X, T, rho = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  T1_index <- which(T == 1)
  T0_index <- which(T == 0)
  n_T1 <- length(T1_index)
  n_T0 <- length(T0_index)

  if (method == "no") {
    return(seq_along(T))
  }

  if (method == "undersample") {
    if (n_T1 > n_T0) {
      keep_t <- sample(T1_index, n_T0, replace = TRUE)
      keep_c <- T0_index
    } else if (n_T0 > n_T1) {
      keep_t <- T1_index
      keep_c <- sample(T0_index, n_T1, replace = TRUE)
    } else {
      keep_t <- T1_index
      keep_c <- T0_index
    }
    return(c(keep_t, keep_c))
  }

  if (method %in% c("prob1", "prob2", "prob3")) {
    # ---- estimate ehat via logistic regression ----
    df_ps <- data.frame(T = T, X)
    fit_ps <- glm(T ~ ., data = df_ps, family = binomial(link = "logit"))
    ehat <- predict(fit_ps, newdata = df_ps, type = "response")

    if (is.null(rho)) {
      rho <- switch(method, "prob1" = 1, "prob2" = 1.5, "prob3" = 2)
    }

    p0 <- ehat[T0_index]
    w0 <- p0 * rho

    w0[!is.finite(w0)] <- 0
    w0 <- pmax(w0, 0)
    if (all(w0 == 0)) {
      pi0 <- rep(1, length(w0))
    } else {
      pi0 <- pmin(w0, 1)
    }

    keep0 <- rbinom(length(T0_index), size = 1, prob = pi0)
    T0_kept_index <- T0_index[keep0 == 1]

    n_T0_kept <- length(T0_kept_index)
    if (n_T0_kept < n_T1) {
      needed <- n_T1 - n_T0_kept
      remaining_idx <- T0_index[keep0 == 0]
      w_remaining <- w0[keep0 == 0]
      if (length(remaining_idx) > 0) {
        if (sum(w_remaining) <= 0) {
          prob_rem <- rep(1 / length(w_remaining), length(w_remaining))
        } else {
          prob_rem <- w_remaining / sum(w_remaining)
        }
        add_size <- min(needed, length(remaining_idx))
        extra <- sample(remaining_idx, size = add_size, replace = FALSE, prob = prob_rem)
        T0_kept_index <- c(T0_kept_index, extra)
      }
    }

    return(c(T1_index, T0_kept_index))
  }

  stop("Unknown sampling method: ", method)
}
