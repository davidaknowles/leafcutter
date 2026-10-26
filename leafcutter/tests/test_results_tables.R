library(leafcutter)
suppressMessages(library(doMC))
registerDoMC(2)

fake_fit <- function(loglr, df, lrtp, par=NULL) list(loglr=loglr, df=df, lrtp=lrtp, fit_full=list(par=par))

test_cluster_results_table <- function() {
  results=list("chr1:clu_1"=fake_fit(10, 1, 0.01),
               "chr1:clu_2"="Too many introns in cluster",
               "chr1:clu_3"=fake_fit(2, 2, 0.04),
               "chr2:clu_4"=simpleError("first line\nsecond line"),
               "chr2:clu_5"=fake_fit(0.5, 1, 0.03),
               "chr3:clu_6"="timeout")
  tab=cluster_results_table(results)
  stopifnot(is.data.frame(tab))
  stopifnot(identical(colnames(tab), c("status", "loglr", "df", "p", "cluster", "p.adjust")))
  stopifnot(identical(tab$cluster, names(results)))
  stopifnot(identical(tab$status, c("Success", "Too many introns in cluster", "Success",
                                    "Error: first linesecond line", "Success", "timeout")))
  stopifnot(identical(tab$loglr, c(10, NA, 2, NA, 0.5, NA)))
  stopifnot(identical(tab$df, c(1, NA, 2, NA, 1, NA)))
  stopifnot(identical(tab$p, c(0.01, NA, 0.04, NA, 0.03, NA)))
  # BH over the three tested clusters: 0.01*3/1=0.03, 0.03*3/2=0.045, 0.04*3/3=0.04, cumulative minimum from the top
  stopifnot(isTRUE(all.equal(tab$p.adjust, c(0.03, NA, 0.04, NA, 0.04, NA), tolerance=1e-14)))
}

test_cluster_results_table_no_success <- function() {
  tab=cluster_results_table(list(a="<=1 junction in cluster", b="Not enough valid samples"))
  stopifnot(identical(tab$status, c("<=1 junction in cluster", "Not enough valid samples")))
  stopifnot(identical(tab$cluster, c("a", "b")))
  stopifnot(all(is.na(tab$loglr)), all(is.na(tab$p)), all(is.na(tab$p.adjust)))
}

test_leaf_cutter_effect_sizes <- function() {
  # K=2: beta = 2*(raw - 1/2), so the intercept row is (0.5, -0.5) and the group row (-0.5, 0.5)
  par1=list(beta_raw=matrix(c(0.75, 0.25,
                              0.25, 0.75), 2, byrow=TRUE, dimnames=list(c("", "x"), c("chr1:1:5:clu_1", "chr1:1:9:clu_1"))),
            beta_scale=c(2, 2), conc=c(3, 1))
  # K=3, group row zero: no change between the conditions
  par2=list(beta_raw=matrix(c(0.5, 0.25, 0.25,
                              1/3, 1/3, 1/3), 2, byrow=TRUE, dimnames=list(c("", "x"), c("a", "b", "c"))),
            beta_scale=c(log(4) * 4/3, 5), conc=c(1, 2, 2))
  results=list(c1=fake_fit(1, 1, 0.1, par1), c2="Too many introns in cluster",
               c3=simpleError("boom"), c4=fake_fit(1, 2, 0.1, par2))
  es=leaf_cutter_effect_sizes(results)
  stopifnot(is.data.frame(es))
  stopifnot(identical(colnames(es), c("intron", "logef", "baseline", "perturbed", "deltapsi")))
  stopifnot(identical(es$intron, c("chr1:1:5:clu_1", "chr1:1:9:clu_1", "a", "b", "c")))
  e=exp(1)
  # c1 baseline: softmax(0.5, -0.5) * (3, 1), normalized = (3e, 1)/(3e + 1); perturbed: softmax(0, 0) * (3, 1) -> (0.75, 0.25)
  # c4: intercept row log(4)*4/3*(1/6, -1/12, -1/12) = (2/9, -1/9, -1/9)*log(4); softmax ∝ (4^(2/9), 4^(-1/9), 4^(-1/9))
  w=c(4^(2/9), 2 * 4^(-1/9), 2 * 4^(-1/9))
  expected=data.frame(
    logef=c(-0.5, 0.5, 0, 0, 0),
    baseline=c(3 * e/(3 * e + 1), 1/(3 * e + 1), w/sum(w)),
    perturbed=c(0.75, 0.25, w/sum(w)))
  stopifnot(isTRUE(all.equal(es$logef, expected$logef, tolerance=1e-14)))
  stopifnot(isTRUE(all.equal(es$baseline, expected$baseline, tolerance=1e-14)))
  stopifnot(isTRUE(all.equal(es$perturbed, expected$perturbed, tolerance=1e-14)))
  stopifnot(identical(es$deltapsi, es$perturbed - es$baseline))
}

# Effect sizes recomputed from the fitted parameters of a real run.
test_effect_sizes_from_fit <- function() {
  i1=c(25, 32, 28, 35, 30, 68, 75, 72, 65, 70)
  i2=c(40, 38, 45, 30, 41, 12, 10, 15, 18, 11)
  counts=rbind("chr1:100:200:clu_1"=i1, "chr1:100:300:clu_1"=100 - i1,
               "chr1:500:600:clu_2"=i2, "chr1:500:700:clu_2"=50, "chr1:500:800:clu_2"=100 - i2,
               "chr1:900:1000:clu_3"=10)
  colnames(counts)=paste0("s", 1:10)
  x=rep(0:1, each=5)
  res=differential_splicing(counts, x)
  es=leaf_cutter_effect_sizes(res)
  stopifnot(identical(es$intron, rownames(counts)[1:5]))
  for (clu in c("chr1:clu_1", "chr1:clu_2")) {
    par=res[[clu]]$fit_full$par
    K=ncol(par$beta_raw)
    b0=par$beta_scale[1] * (par$beta_raw[1, ] - 1/K)
    b1=par$beta_scale[2] * (par$beta_raw[2, ] - 1/K)
    psi=function(eta) { a=par$conc * exp(eta)/sum(exp(eta)); a/sum(a) }
    rows=match(colnames(par$beta_raw), es$intron)
    stopifnot(isTRUE(all.equal(es$logef[rows], as.vector(b1), tolerance=1e-12)))
    stopifnot(isTRUE(all.equal(es$baseline[rows], as.vector(psi(b0)), tolerance=1e-12)))
    stopifnot(isTRUE(all.equal(es$perturbed[rows], as.vector(psi(b0 + b1)), tolerance=1e-12)))
    stopifnot(abs(sum(es$baseline[rows]) - 1) < 1e-12, abs(sum(es$perturbed[rows]) - 1) < 1e-12)
  }
  # mean usage of clu_1's first intron goes from 0.30 to 0.70 in group 1, of clu_2's first from 0.26 to 0.09
  stopifnot(es$deltapsi[1] > 0.3, es$deltapsi[1] < 0.5, es$logef[1] > 0)
  stopifnot(es$deltapsi[3] < -0.1, es$deltapsi[3] > -0.25, es$logef[3] < 0)
  stopifnot(abs(es$baseline[1] - 0.30) < 0.05, abs(es$perturbed[1] - 0.70) < 0.05)
  tab=cluster_results_table(res)
  stopifnot(identical(tab$cluster, c("chr1:clu_1", "chr1:clu_2", "chr1:clu_3")))
  stopifnot(identical(tab$status, c("Success", "Success", "<=1 junction in cluster")))
  stopifnot(identical(tab$df, c(1, 2, NA)))
  stopifnot(identical(tab$loglr[1:2], c(res[["chr1:clu_1"]]$loglr, res[["chr1:clu_2"]]$loglr)))
  stopifnot(identical(tab$p[1:2], pchisq(2 * tab$loglr[1:2], df=tab$df[1:2], lower.tail=FALSE)))
  stopifnot(identical(tab$p.adjust, p.adjust(tab$p, method="fdr")))
  stopifnot(all(tab$p[1:2] < 1e-6))
}

tests=list(test_cluster_results_table=test_cluster_results_table,
           test_cluster_results_table_no_success=test_cluster_results_table_no_success,
           test_leaf_cutter_effect_sizes=test_leaf_cutter_effect_sizes,
           test_effect_sizes_from_fit=test_effect_sizes_from_fit)
for (name in names(tests)) {
  cat(name, "\n")
  tests[[name]]()
}
