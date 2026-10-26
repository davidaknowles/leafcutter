library(leafcutter)
suppressMessages(library(doMC))
registerDoMC(2)

# See test_dirichlet_multinomial_anova_mc.R: the objective rstan::optimizing
# reports for dm_glm_multi_conc / dm_glm_robust, priors without constants.
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

samples <- paste0("s", 1:10)
x <- rep(0:1, each=5)

# rows: named introns; values: counts for s1..s10
cluster <- function(...) {
  rows=list(...)
  m=do.call(rbind, rows)
  rownames(m)=names(rows)
  m
}

planted <- function(chr, clu, start=100) {
  i1=c(25, 32, 28, 35, 30, 68, 75, 72, 65, 70)
  m=rbind(i1, 100 - i1)
  rownames(m)=sprintf("%s:%d:%d:%s", chr, start, start + c(100, 200), clu)
  m
}

null_cluster <- function(chr, clu, start=100) {
  i1=c(40, 50, 45, 55, 60, 55, 45, 60, 40, 50)
  m=rbind(i1, 100 - i1)
  rownames(m)=sprintf("%s:%d:%d:%s", chr, start, start + c(100, 200), clu)
  m
}

status_counts <- function() {
  m=rbind(
    planted("chr1", "clu_1"),
    null_cluster("chr1", "clu_2", 1000),
    cluster("chr1:5000:5100:clu_3"=rep(50, 10)),
    cluster("chr1:6000:6100:clu_4"=c(30, rep(0, 9)),
            "chr1:6000:6200:clu_4"=c(20, rep(0, 9))),
    # s1 has 20 reads, everyone else 19
    cluster("chr1:7000:7100:clu_5"=rep(10, 10),
            "chr1:7000:7200:clu_5"=c(10, rep(9, 9))),
    # introns 2 and 3 are each used in four samples
    cluster("chr1:8000:8100:clu_6"=rep(50, 10),
            "chr1:8000:8200:clu_6"=c(10, 10, 10, 10, 0, 0, 0, 0, 0, 0),
            "chr1:8000:8300:clu_6"=c(0, 0, 0, 0, 0, 10, 10, 10, 10, 0)),
    # group 1 has 10 reads per sample, below min_coverage
    cluster("chr1:9000:9100:clu_7"=c(rep(50, 5), rep(5, 5)),
            "chr1:9000:9200:clu_7"=c(rep(50, 5), rep(5, 5))),
    do.call(rbind, lapply(1:11, function(i) {
      r=matrix(rep(10 + i, 10), 1)
      rownames(r)=sprintf("chr1:20000:%d:clu_8", 20000 + 100 * i)
      r
    })),
    null_cluster("chr2", "clu_1"))
  colnames(m)=samples
  m
}

expected_status <- c(
  "chr1:clu_1"="Success",
  "chr1:clu_2"="Success",
  "chr1:clu_3"="<=1 junction in cluster",
  "chr1:clu_4"="<=1 sample with coverage>0",
  "chr1:clu_5"="<=1 sample with coverage>min_coverage",
  "chr1:clu_6"="<2 introns used in >=min_samples_per_intron samples",
  "chr1:clu_7"="Not enough valid samples",
  "chr1:clu_8"="Too many introns in cluster",
  "chr2:clu_1"="Success")

result_status <- function(res) sapply(res, function(r) if (is.character(r)) r else "Success")

# The DM objective identity at both fits, evaluated on the samples and introns the cluster test should use.
check_fit <- function(r, y, x_used, robust=FALSE, conf=NULL) {
  stopifnot(is.list(r), !is.character(r))
  xFull=cbind(1, x_used, conf)
  xNull=cbind(rep(1, length(x_used)), conf)
  stopifnot(identical(colnames(r$fit_full$par$beta_raw), colnames(y)))
  stopifnot(close_to(r$fit_full$value, dm_objective(r$fit_full$par, xFull, y, robust=robust), 1e-8))
  stopifnot(close_to(r$fit_null$value, dm_objective(r$fit_null$par, xNull, y, robust=robust), 1e-8))
  stopifnot(identical(r$loglr, r$fit_full$value - r$fit_null$value))
  stopifnot(r$df == ncol(y) - 1)
  stopifnot(identical(r$lrtp, pchisq(2 * r$loglr, df=r$df, lower.tail=FALSE)))
}

test_statuses <- function() {
  counts=status_counts()
  res=differential_splicing(counts, x)
  stopifnot(is.list(res))
  stopifnot(identical(names(res), sort(names(expected_status))))
  stopifnot(identical(result_status(res)[names(expected_status)], expected_status))
  # defaults: max_cluster_size=10, min_samples_per_intron=5, min_samples_per_group=4, min_coverage=20
  y=t(counts[grepl(":clu_1$", rownames(counts)) & startsWith(rownames(counts), "chr1:"), ])
  check_fit(res[["chr1:clu_1"]], y, x)
  stopifnot(res[["chr1:clu_1"]]$lrtp < 1e-6)
  y=t(counts[grepl(":clu_2$", rownames(counts)), ])
  check_fit(res[["chr1:clu_2"]], y, x)
  stopifnot(res[["chr1:clu_2"]]$lrtp > 0.5)
  # clusters are keyed by chromosome and cluster id
  y=t(counts[startsWith(rownames(counts), "chr2:"), ])
  check_fit(res[["chr2:clu_1"]], y, x)
  stopifnot(res[["chr2:clu_1"]]$lrtp > 0.5)
}

test_robust <- function() {
  counts=status_counts()
  res=differential_splicing(counts, x, robust=TRUE)
  stopifnot(identical(result_status(res)[names(expected_status)], expected_status))
  y=t(counts[grepl(":clu_1$", rownames(counts)) & startsWith(rownames(counts), "chr1:"), ])
  check_fit(res[["chr1:clu_1"]], y, x, robust=TRUE)
  stopifnot(res[["chr1:clu_1"]]$lrtp < 1e-6)
  y=t(counts[grepl(":clu_2$", rownames(counts)), ])
  check_fit(res[["chr1:clu_2"]], y, x, robust=TRUE)
}

test_max_cluster_size <- function() {
  counts=rbind(planted("chr1", "clu_1"),
               cluster("chr1:100:400:clu_1"=c(10, 12, 9, 11, 10, 10, 12, 9, 11, 10)))
  colnames(counts)=samples
  stopifnot(identical(result_status(differential_splicing(counts, x, max_cluster_size=3)), c("chr1:clu_1"="Success")))
  stopifnot(identical(result_status(differential_splicing(counts, x, max_cluster_size=2)), c("chr1:clu_1"="Too many introns in cluster")))
}

test_min_samples_per_intron <- function() {
  # intron 3 is used in exactly five samples (three in group 0, two in group 1)
  counts=rbind(planted("chr1", "clu_1"),
               cluster("chr1:100:400:clu_1"=c(5, 5, 5, 0, 0, 5, 5, 0, 0, 0)))
  colnames(counts)=samples
  res=differential_splicing(counts, x, min_samples_per_intron=5)
  check_fit(res[["chr1:clu_1"]], t(counts), x)
  stopifnot(res[["chr1:clu_1"]]$df == 2)
  res=differential_splicing(counts, x, min_samples_per_intron=6)
  check_fit(res[["chr1:clu_1"]], t(counts[1:2, ]), x)
  stopifnot(res[["chr1:clu_1"]]$df == 1)
  stopifnot(identical(colnames(res[["chr1:clu_1"]]$fit_full$par$beta_raw), rownames(counts)[1:2]))
}

test_zero_coverage_samples_dropped <- function() {
  # s10 has no reads in the cluster: it is left out of the fit, and s1..s9 remain
  counts=planted("chr1", "clu_1")
  counts[, 10]=0
  colnames(counts)=samples
  res=differential_splicing(counts, x, min_samples_per_group=4)
  check_fit(res[["chr1:clu_1"]], t(counts[, 1:9]), x[1:9])
  # a zero row adds nothing to the likelihood, so compare with the run on s1..s9
  stopifnot(identical(res[["chr1:clu_1"]], differential_splicing(counts[, 1:9], x[1:9], min_samples_per_group=4)[["chr1:clu_1"]]))
}

test_min_coverage <- function() {
  # every sample has 20 reads: covered at min_coverage=20, not at 21
  counts=cluster("chr1:100:200:clu_1"=c(5, 6, 4, 7, 5, 15, 14, 16, 13, 15),
                 "chr1:100:300:clu_1"=c(15, 14, 16, 13, 15, 5, 6, 4, 7, 5))
  colnames(counts)=samples
  res=differential_splicing(counts, x, min_coverage=20)
  check_fit(res[["chr1:clu_1"]], t(counts), x)
  stopifnot(identical(result_status(differential_splicing(counts, x, min_coverage=21)),
                      c("chr1:clu_1"="<=1 sample with coverage>min_coverage")))
}

test_min_samples_per_group <- function() {
  # group 1 has four samples at or above min_coverage (s10 has 10 reads)
  counts=planted("chr1", "clu_1")
  counts[, 10]=c(6, 4)
  colnames(counts)=samples
  res=differential_splicing(counts, x, min_samples_per_group=4)
  check_fit(res[["chr1:clu_1"]], t(counts), x)
  stopifnot(identical(result_status(differential_splicing(counts, x, min_samples_per_group=5)),
                      c("chr1:clu_1"="Not enough valid samples")))
  # min_coverage decides which samples count towards min_samples_per_group
  res=differential_splicing(counts, x, min_samples_per_group=5, min_coverage=10)
  check_fit(res[["chr1:clu_1"]], t(counts), x)
}

test_confounders <- function() {
  counts=rbind(planted("chr1", "clu_1"), null_cluster("chr1", "clu_2", 1000))
  counts["chr1:1000:1100:clu_2", 10]=0
  counts["chr1:1000:1200:clu_2", 10]=0
  colnames(counts)=samples
  age=c(0.5, -1.2, 0.3, 1.1, -0.7, 0.9, -0.4, 1.4, -1.0, 2.0)
  constant=rep(1, 10)
  # varies only through s10, which has no reads in clu_2
  only_s10=c(rep(0, 9), 1)
  conf=cbind(age=age, constant=constant, only_s10=only_s10)
  res=differential_splicing(counts, x, confounders=conf)
  r1=res[["chr1:clu_1"]]
  stopifnot(identical(rownames(r1$fit_full$par$beta_raw), c("", "x_subset", "age", "only_s10")))
  stopifnot(identical(rownames(r1$fit_null$par$beta_raw), c("", "age", "only_s10")) || r1$refit_null_flag)
  check_fit(r1, t(counts[1:2, ]), x, conf=conf[, c("age", "only_s10")])
  stopifnot(r1$lrtp < 1e-6)
  r2=res[["chr1:clu_2"]]
  stopifnot(identical(rownames(r2$fit_full$par$beta_raw), c("", "x_subset", "age")))
  check_fit(r2, t(counts[3:4, 1:9]), x[1:9], conf=conf[1:9, "age", drop=FALSE])
}

test_bad_x <- function() {
  counts=planted("chr1", "clu_1")
  colnames(counts)=samples
  stopifnot(inherits(try(differential_splicing(counts, x[1:9]), silent=TRUE), "try-error"))
  stopifnot(inherits(try(differential_splicing(counts, c(NA, x[-1])), silent=TRUE), "try-error"))
}

test_error_in_one_cluster <- function() {
  counts=rbind(planted("chr1", "clu_1"), null_cluster("chr1", "clu_2", 1000))
  counts[3, 4]=NA
  colnames(counts)=samples
  res=differential_splicing(counts, x)
  stopifnot(inherits(res[["chr1:clu_2"]], "error"))
  check_fit(res[["chr1:clu_1"]], t(counts[1:2, ]), x)
  tab=cluster_results_table(res)
  stopifnot(identical(tab$cluster, c("chr1:clu_1", "chr1:clu_2")))
  stopifnot(identical(tab$status[1], "Success"))
  stopifnot(startsWith(tab$status[2], "Error"), !grepl("\n", tab$status[2]))
  stopifnot(is.na(tab$loglr[2]), is.na(tab$df[2]), is.na(tab$p[2]), is.na(tab$p.adjust[2]))
  stopifnot(identical(tab$p.adjust[1], tab$p[1]))
  stopifnot(identical(leaf_cutter_effect_sizes(res)$intron, rownames(counts)[1:2]))
}

test_random_init <- function() {
  counts=status_counts()
  r1=differential_splicing(counts, x, init="random", seed=5)
  r2=differential_splicing(counts, x, init="random", seed=5)
  stopifnot(identical(r1, r2))
  stopifnot(identical(result_status(r1)[names(expected_status)], expected_status))
  # smart and random initialization reach the same optimum, up to the optimizer's tolerance
  s=differential_splicing(counts, x)
  for (clu in c("chr1:clu_1", "chr1:clu_2", "chr2:clu_1")) {
    stopifnot(abs(r1[[clu]]$fit_full$value - s[[clu]]$fit_full$value) < 0.05)
    stopifnot(abs(r1[[clu]]$fit_null$value - s[[clu]]$fit_null$value) < 0.05)
  }
  y=t(counts[grepl(":clu_1$", rownames(counts)) & startsWith(rownames(counts), "chr1:"), ])
  check_fit(r1[["chr1:clu_1"]], y, x)
  stopifnot(r1[["chr1:clu_1"]]$lrtp < 1e-6)
}

test_debug_same_results <- function() {
  counts=rbind(planted("chr1", "clu_1"), null_cluster("chr1", "clu_2", 1000))
  colnames(counts)=samples
  stopifnot(identical(differential_splicing(counts, x, debug=TRUE), differential_splicing(counts, x)))
}

test_parallel_matches_serial <- function() {
  counts=status_counts()
  par=differential_splicing(counts, x)
  registerDoMC(1)
  on.exit(registerDoMC(2))
  ser=differential_splicing(counts, x)
  stopifnot(identical(par, ser))
}

tests=list(test_statuses=test_statuses,
           test_robust=test_robust,
           test_max_cluster_size=test_max_cluster_size,
           test_min_samples_per_intron=test_min_samples_per_intron,
           test_zero_coverage_samples_dropped=test_zero_coverage_samples_dropped,
           test_min_coverage=test_min_coverage,
           test_min_samples_per_group=test_min_samples_per_group,
           test_confounders=test_confounders,
           test_bad_x=test_bad_x,
           test_error_in_one_cluster=test_error_in_one_cluster,
           test_random_init=test_random_init,
           test_debug_same_results=test_debug_same_results,
           test_parallel_matches_serial=test_parallel_matches_serial)
for (name in names(tests)) {
  cat(name, "\n")
  tests[[name]]()
}
