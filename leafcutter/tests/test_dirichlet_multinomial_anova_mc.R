library(leafcutter)

# Objective that rstan::optimizing reports as $value for dm_glm_multi_conc and
# dm_glm_robust, written from the model: the Dirichlet-multinomial log-likelihood
# without the multinomial coefficient, plus the gamma(concShape, concRate) prior
# on conc and (robust only) the beta(a, b) prior on outlier_prob. The priors
# enter through "~" statements, so Stan drops their normalizing constants.
dm_objective <- function(par, X, y, concShape=1.0001, concRate=1e-4, robust=FALSE, a=1.01, b=100) {
  K=ncol(y)
  beta=matrix(0, ncol(X), K)
  for (p in seq_len(ncol(X)))
    for (k in seq_len(K))
      beta[p, k]=par$beta_scale[p] * (par$beta_raw[p, k] - 1/K)
  total=sum((concShape - 1) * log(par$conc) - concRate * par$conc)
  if (robust)
    total=total + (a - 1) * log(par$outlier_prob) + (b - 1) * log(1 - par$outlier_prob)
  for (n in seq_len(nrow(y))) {
    eta=as.vector(X[n, , drop=FALSE] %*% beta)
    s=exp(eta - max(eta))
    alpha=par$conc * s / sum(s)
    ll=lgamma(sum(alpha)) - lgamma(sum(alpha) + sum(y[n, ])) + sum(lgamma(alpha + y[n, ]) - lgamma(alpha))
    if (robust) {
      outlier=lgamma(K) + sum(lgamma(1 + y[n, ])) - lgamma(K + sum(y[n, ]))
      l=c(log(1 - par$outlier_prob) + ll, log(par$outlier_prob) + outlier)
      ll=max(l) + log(sum(exp(l - max(l))))
    }
    total=total + ll
  }
  total
}

close_to <- function(a, b, tol) abs(a - b) <= tol * max(1, abs(b))

# Six samples per group, three introns; group 1 shifts usage from introns 2-3 to intron 1.
planted_y <- function() {
  y=rbind(c(15, 45, 40), c(25, 35, 40), c(20, 38, 42), c(12, 50, 38), c(28, 30, 42), c(20, 40, 40),
          c(55, 25, 20), c(68, 12, 20), c(58, 22, 20), c(50, 28, 22), c(62, 18, 20), c(57, 21, 22))
  colnames(y)=c("i1", "i2", "i3")
  y
}

# Group 1 holds the same rows as group 0 in another order, so the group effect at the optimum is zero.
null_y <- function() {
  g0=rbind(c(15, 45, 40), c(25, 35, 40), c(20, 38, 42), c(12, 50, 38), c(28, 30, 42), c(20, 40, 40))
  y=rbind(g0, g0[c(4, 1, 6, 2, 5, 3), ])
  colnames(y)=c("i1", "i2", "i3")
  y
}

group <- rep(0:1, each=6)
xFull <- cbind(1, group)
xNull <- xFull[, 1, drop=FALSE]

check_result <- function(r, xFull, xNull, y, robust) {
  K=ncol(y)
  stopifnot(close_to(r$fit_full$value, dm_objective(r$fit_full$par, xFull, y, robust=robust), 1e-8))
  stopifnot(close_to(r$fit_null$value, dm_objective(r$fit_null$par, xNull, y, robust=robust), 1e-8))
  stopifnot(identical(r$loglr, r$fit_full$value - r$fit_null$value))
  stopifnot(r$df == (ncol(xFull) - ncol(xNull)) * (K - 1))
  stopifnot(identical(r$lrtp, pchisq(2 * r$loglr, df=r$df, lower.tail=FALSE)))
  stopifnot(is.logical(r$refit_null_flag), length(r$refit_null_flag) == 1)
  stopifnot(identical(dim(r$fit_full$par$beta_raw), c(ncol(xFull), K)))
  stopifnot(identical(dim(r$fit_null$par$beta_raw), c(ncol(xNull), K)))
  stopifnot(identical(colnames(r$fit_full$par$beta_raw), colnames(y)))
  stopifnot(identical(rownames(r$fit_full$par$beta_raw), colnames(xFull)))
  # a refitted null loses its dimnames, so they are checked only without a refit
  if (!r$refit_null_flag) {
    stopifnot(identical(colnames(r$fit_null$par$beta_raw), colnames(y)))
    stopifnot(identical(rownames(r$fit_null$par$beta_raw), colnames(xNull)))
  }
  stopifnot(all(abs(rowSums(r$fit_full$par$beta_raw) - 1) < 1e-8))
  stopifnot(all(r$fit_full$par$conc > 0), length(r$fit_full$par$conc) == K)
  stopifnot(identical("outlier_prob" %in% names(r$fit_full$par), robust))
  if (robust) stopifnot(r$fit_full$par$outlier_prob > 0, r$fit_full$par$outlier_prob < 1)
}

test_planted_effect <- function() {
  y=planted_y()
  r=dirichlet_multinomial_anova_mc(xFull, xNull, y, robust=FALSE)
  check_result(r, xFull, xNull, y, robust=FALSE)
  stopifnot(r$df == 2)
  stopifnot(r$lrtp < 1e-6)
  # the group effect raises intron 1 relative to the others
  b=beta_real(r$fit_full$par)
  stopifnot(b[2, "i1"] > 0, b[2, "i1"] > b[2, "i2"], b[2, "i1"] > b[2, "i3"])
}

test_null_cluster <- function() {
  y=null_y()
  r=dirichlet_multinomial_anova_mc(xFull, xNull, y, robust=FALSE)
  check_result(r, xFull, xNull, y, robust=FALSE)
  stopifnot(abs(r$loglr) < 1e-4)
  stopifnot(r$lrtp > 0.5)
  stopifnot(identical(r$refit_null_flag, FALSE))
}

test_robust_model <- function() {
  y=planted_y()
  r=dirichlet_multinomial_anova_mc(xFull, xNull, y, robust=TRUE)
  check_result(r, xFull, xNull, y, robust=TRUE)
  stopifnot(r$lrtp < 1e-6)
  # robust=TRUE is the default
  r2=dirichlet_multinomial_anova_mc(xFull, xNull, y)
  check_result(r2, xFull, xNull, y, robust=TRUE)
  rn=dirichlet_multinomial_anova_mc(xFull, xNull, null_y(), robust=TRUE)
  check_result(rn, xFull, xNull, null_y(), robust=TRUE)
  stopifnot(rn$lrtp > 0.5)
}

test_prior_arguments <- function() {
  y=planted_y()
  r=dirichlet_multinomial_anova_mc(xFull, xNull, y, robust=TRUE, concShape=2, concRate=0.05, outlier_prior_a=2, outlier_prior_b=20)
  stopifnot(close_to(r$fit_full$value, dm_objective(r$fit_full$par, xFull, y, concShape=2, concRate=0.05, robust=TRUE, a=2, b=20), 1e-8))
  stopifnot(close_to(r$fit_null$value, dm_objective(r$fit_null$par, xNull, y, concShape=2, concRate=0.05, robust=TRUE, a=2, b=20), 1e-8))
}

test_df_with_covariates <- function() {
  y=planted_y()
  covariate=c(0.3, -1.2, 0.8, 0.1, -0.4, 1.5, -0.7, 0.2, 1.1, -1.6, 0.5, -0.1)
  xF=cbind(1, group, covariate)
  xN=xF[, 1, drop=FALSE]
  r=dirichlet_multinomial_anova_mc(xF, xN, y, robust=FALSE)
  check_result(r, xF, xN, y, robust=FALSE)
  stopifnot(r$df == 4)
  # lrtp < 0.001, so the null is refit from the full fit's intercept row; the
  # reported null is at least as good as that refit
  stopifnot(r$lrtp < 0.001)
  b=r$fit_full$par$beta_raw[1, , drop=FALSE]
  b[b < 1e-6]=1e-6
  b[b > 1 - 1e-6]=1 - 1e-6
  init=list(beta_raw=sweep(b, 1, rowSums(b), "/"), beta_scale=as.array(r$fit_full$par$beta_scale[1]), conc=r$fit_full$par$conc)
  refit=rstan::optimizing(leafcutter:::stanmodels$dm_glm_multi_conc, data=list(N=nrow(y), K=ncol(y), P=1, y=y, x=xN, concShape=1.0001, concRate=1e-4),
                          init=init, as_vector=FALSE)
  stopifnot(r$fit_null$value >= refit$value - 1e-6)
  # the covariate in both models: one tested column
  xN2=xF[, c(1, 3)]
  r2=dirichlet_multinomial_anova_mc(xF, xN2, y, robust=FALSE)
  check_result(r2, xF, xN2, y, robust=FALSE)
  stopifnot(r2$df == 2)
  stopifnot(r2$lrtp < 1e-6)
}

test_random_init <- function() {
  y=planted_y()
  r1=dirichlet_multinomial_anova_mc(xFull, xNull, y, robust=FALSE, init="random", seed=11)
  r2=dirichlet_multinomial_anova_mc(xFull, xNull, y, robust=FALSE, init="random", seed=11)
  check_result(r1, xFull, xNull, y, robust=FALSE)
  stopifnot(identical(r1, r2))
  # smart and random initialization reach the same optimum, up to the optimizer's tolerance
  s=dirichlet_multinomial_anova_mc(xFull, xNull, y, robust=FALSE)
  stopifnot(abs(r1$fit_full$value - s$fit_full$value) < 0.05)
  stopifnot(abs(r1$fit_null$value - s$fit_null$value) < 0.05)
  stopifnot(r1$lrtp < 1e-6)
}

test_cached_null_fit <- function() {
  y=null_y()
  r=dirichlet_multinomial_anova_mc(xFull, xNull, y, robust=FALSE)
  cached=r$fit_null
  cached$par$beta_raw=unname(cached$par$beta_raw)
  r2=dirichlet_multinomial_anova_mc(xFull, xNull, y, robust=FALSE, fit_null=cached)
  check_result(r2, xFull, xNull, y, robust=FALSE)
  stopifnot(identical(r2$fit_null$value, r$fit_null$value))
  stopifnot(identical(r2$fit_null$par$beta_raw, r$fit_null$par$beta_raw))
  stopifnot(identical(r2$loglr, r2$fit_full$value - r$fit_null$value))
}

tests=list(test_planted_effect=test_planted_effect,
           test_null_cluster=test_null_cluster,
           test_robust_model=test_robust_model,
           test_prior_arguments=test_prior_arguments,
           test_df_with_covariates=test_df_with_covariates,
           test_random_init=test_random_init,
           test_cached_null_fit=test_cached_null_fit)
for (name in names(tests)) {
  cat(name, "\n")
  tests[[name]]()
}
