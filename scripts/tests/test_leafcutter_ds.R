# Command-line tests for scripts/leafcutter_ds.R. The script lives outside the
# package, so R CMD check does not run these; run them with
#   Rscript scripts/tests/test_leafcutter_ds.R
# with the leafcutter package installed on the library path.
library(leafcutter)
suppressMessages(library(doMC))
registerDoMC(2)

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE))
if (length(this_file) != 1) stop("run this file with Rscript")
script <- normalizePath(file.path(dirname(this_file), "..", "leafcutter_ds.R"), mustWork=TRUE)
rscript <- file.path(R.home("bin"), "Rscript")

root <- tempfile("leafcutter_ds_tests")
dir.create(root)

run_ds <- function(args, wd=root) {
  out=tempfile(tmpdir=root)
  err=tempfile(tmpdir=root)
  old=setwd(wd)
  on.exit(setwd(old))
  status=system2(rscript, c(shQuote(script), args), stdout=out, stderr=err)
  list(status=status, stdout=readLines(out), stderr=readLines(err))
}

new_dir <- function() {
  d=tempfile("run", tmpdir=root)
  dir.create(d)
  d
}

write_counts <- function(m, path) {
  con=if (grepl("\\.gz$", path)) gzfile(path, "w") else file(path, "w")
  writeLines(paste(colnames(m), collapse=" "), con)
  write.table(m, con, quote=FALSE, sep=" ", col.names=FALSE)
  close(con)
  path
}

write_groups <- function(df, path) {
  write.table(df, path, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  path
}

read_sig <- function(prefix) read.delim(paste0(prefix, "_cluster_significance.txt"), stringsAsFactors=FALSE, check.names=FALSE)
read_es <- function(prefix) read.delim(paste0(prefix, "_effect_sizes.txt"), stringsAsFactors=FALSE, check.names=FALSE)
file_bytes <- function(path) readBin(path, "raw", file.info(path)$size)
same_outputs <- function(p1, p2) {
  for (suffix in c("_cluster_significance.txt", "_effect_sizes.txt"))
    if (!identical(file_bytes(paste0(p1, suffix)), file_bytes(paste0(p2, suffix)))) return(FALSE)
  TRUE
}
status_of <- function(sig) setNames(sig$status, sig$cluster)
p_of <- function(sig) setNames(sig$p, sig$cluster)
close_rel <- function(a, b, tol) all(abs(a - b) <= tol * pmax(1e-300, abs(b)))

ctrl <- paste0("ctrl", 1:5)
case <- paste0("case", 1:5)
samples <- c(ctrl, case)

two_intron <- function(chr, start, clu, i1, total=100) {
  m=rbind(i1, total - i1)
  rownames(m)=sprintf("%s:%d:%d:%s", chr, start, start + c(100, 200), clu)
  m
}

# Main input, chromosome names without "chr". Columns: ctrl1..5, an unlisted sample, case1..5.
#   clu_1  planted effect: intron 1 at ~30% in ctrl, ~70% in case
#   clu_2  no effect: case rows are the ctrl rows reordered
#   clu_3  three introns, planted effect on intron 1 (for -s)
#   clu_4  intron 2 used in exactly four samples (for -i)
#   clu_5  case samples have 10 reads each (for -c)
#   clu_6  two case samples have 15 reads, three have 100 (for -g)
#   chr2 clu_7  no effect
main_counts <- function() {
  i1=c(25, 32, 28, 35, 30, 68, 75, 72, 65, 70)
  inull=c(40, 50, 45, 55, 60, 55, 45, 60, 40, 50)
  clu3=rbind(i1, c(40, 36, 38, 33, 37, 16, 13, 14, 18, 15), 100 - i1 - c(40, 36, 38, 33, 37, 16, 13, 14, 18, 15))
  rownames(clu3)=sprintf("1:%d:%d:clu_3", 3000, 3000 + c(100, 200, 300))
  clu4=rbind(c(50, 60, 55, 45, 50, 52, 58, 47, 51, 49), c(10, 10, 0, 0, 0, 10, 10, 0, 0, 0))
  rownames(clu4)=sprintf("1:%d:%d:clu_4", 4000, 4000 + c(100, 200))
  m=rbind(two_intron("1", 1000, "clu_1", i1),
          two_intron("1", 2000, "clu_2", inull),
          clu3,
          clu4,
          two_intron("1", 5000, "clu_5", c(40, 50, 45, 55, 60, 5, 4, 6, 5, 3), c(rep(100, 5), rep(10, 5))),
          two_intron("1", 6000, "clu_6", c(25, 32, 28, 35, 30, 68, 75, 72, 10, 11), c(rep(100, 8), 15, 15)),
          two_intron("2", 1000, "clu_7", inull))
  colnames(m)=samples
  # a sample missing from the groups file, with counts that would change every cluster
  cbind(m[, ctrl], unlisted=c(1000, 0)[1 + (seq_len(nrow(m)) %% 2)], m[, case])
}

groups_df <- data.frame(sample=samples, group=rep(c("ctrl", "case"), each=5), stringsAsFactors=FALSE)
counts_file <- write_counts(main_counts(), file.path(root, "main.counts.gz"))
groups_file <- write_groups(groups_df, file.path(root, "groups.txt"))

main_status <- c("chr1:clu_1"="Success", "chr1:clu_2"="Success", "chr1:clu_3"="Success",
                 "chr1:clu_4"="<2 introns used in >=min_samples_per_intron samples",
                 "chr1:clu_5"="Not enough valid samples", "chr1:clu_6"="Success", "chr2:clu_7"="Success")

default_prefix <- file.path(new_dir(), "default")
default_run <- run_ds(c("-p", "2", "-o", shQuote(default_prefix), shQuote(counts_file), shQuote(groups_file)))

check_tables <- function(sig, es, groups) {
  stopifnot(identical(colnames(sig)[1:6], c("status", "loglr", "df", "p", "cluster", "p.adjust")))
  stopifnot(identical(colnames(es), c("intron", "logef", groups, "deltapsi")))
  ok=sig$status == "Success"
  stopifnot(all(is.na(sig$loglr[!ok])), all(is.na(sig$p[!ok])), all(is.na(sig$df[!ok])))
  es_cluster=sub("^(.*):[0-9]+:[0-9]+:(.*)$", "\\1:\\2", es$intron)
  stopifnot(setequal(unique(es_cluster), sig$cluster[ok]))
  n_introns=table(es_cluster)[sig$cluster[ok]]
  stopifnot(all(sig$df[ok] == n_introns - 1))
  stopifnot(close_rel(sig$p[ok], pchisq(2 * sig$loglr[ok], df=sig$df[ok], lower.tail=FALSE), 1e-9))
  stopifnot(close_rel(sig$p.adjust[ok], p.adjust(sig$p, method="fdr")[ok], 1e-9))
  stopifnot(all(abs(es$deltapsi - (es[[groups[2]]] - es[[groups[1]]])) < 1e-12))
  for (clu in sig$cluster[ok]) {
    rows=es_cluster == clu
    stopifnot(abs(sum(es[[groups[1]]][rows]) - 1) < 1e-9, abs(sum(es[[groups[2]]][rows]) - 1) < 1e-9)
    stopifnot(abs(sum(es$logef[rows])) < 1e-9)
  }
}

test_default_run <- function() {
  stopifnot(default_run$status == 0)
  stopifnot(any(default_run$stdout == "Encoding as ctrl =0, case =1"))
  stopifnot(any(default_run$stdout == "No exon_file provided."))
  sig=read_sig(default_prefix)
  es=read_es(default_prefix)
  check_tables(sig, es, c("ctrl", "case"))
  stopifnot(identical(colnames(sig), c("status", "loglr", "df", "p", "cluster", "p.adjust")))
  stopifnot(identical(sig$cluster, names(main_status)))
  stopifnot(identical(status_of(sig), main_status))
  stopifnot(identical(sig$df[sig$status == "Success"], c(1L, 1L, 2L, 1L, 1L)))
  p=p_of(sig)
  stopifnot(p["chr1:clu_1"] < 1e-6, p["chr1:clu_3"] < 1e-6, p["chr1:clu_6"] < 1e-6)
  stopifnot(p["chr1:clu_2"] > 0.5, p["chr2:clu_7"] > 0.5)
  # intron names gain "chr"; introns of untested clusters are absent
  stopifnot(identical(es$intron, c("chr1:1000:1100:clu_1", "chr1:1000:1200:clu_1", "chr1:2000:2100:clu_2", "chr1:2000:2200:clu_2",
                                   "chr1:3000:3100:clu_3", "chr1:3000:3200:clu_3", "chr1:3000:3300:clu_3",
                                   "chr1:6000:6100:clu_6", "chr1:6000:6200:clu_6", "chr2:1000:1100:clu_7", "chr2:1000:1200:clu_7")))
  # planted effect: intron 1 of clu_1 goes from ~0.30 in ctrl to ~0.70 in case
  stopifnot(es$logef[1] > 0, abs(es$ctrl[1] - 0.30) < 0.05, abs(es$case[1] - 0.70) < 0.05)
  stopifnot(es$deltapsi[5] > 0.3)
}

test_unlisted_sample_ignored <- function() {
  m=main_counts()
  counts2=write_counts(m[, samples], file.path(root, "listed_only.counts.txt"))
  prefix=file.path(new_dir(), "listed")
  r=run_ds(c("-p", "2", "-o", shQuote(prefix), shQuote(counts2), shQuote(groups_file)))
  stopifnot(r$status == 0)
  # also checks that a plain-text counts file reads the same as a gzipped one
  stopifnot(same_outputs(prefix, default_prefix))
}

test_default_prefix_and_one_thread <- function() {
  wd=new_dir()
  r=run_ds(c("-p", "1", shQuote(counts_file), shQuote(groups_file)), wd=wd)
  stopifnot(r$status == 0)
  stopifnot(setequal(list.files(wd), c("leafcutter_ds_cluster_significance.txt", "leafcutter_ds_effect_sizes.txt")))
  stopifnot(same_outputs(file.path(wd, "leafcutter_ds"), default_prefix))
}

test_max_cluster_size <- function() {
  prefix=file.path(new_dir(), "s2")
  r=run_ds(c("-p", "2", "-s", "2", "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  sig=read_sig(prefix)
  expected=main_status
  expected["chr1:clu_3"]="Too many introns in cluster"
  stopifnot(identical(status_of(sig), expected))
  stopifnot(!any(grepl(":clu_3$", read_es(prefix)$intron)))
}

test_min_samples_per_intron <- function() {
  prefix=file.path(new_dir(), "i4")
  r=run_ds(c("-p", "2", "-i", "4", "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  sig=read_sig(prefix)
  es=read_es(prefix)
  expected=main_status
  expected["chr1:clu_4"]="Success"
  stopifnot(identical(status_of(sig), expected))
  check_tables(sig, es, c("ctrl", "case"))
  stopifnot(identical(es$intron[grepl(":clu_4$", es$intron)], c("chr1:4000:4100:clu_4", "chr1:4000:4200:clu_4")))
}

test_min_coverage <- function() {
  # case samples of clu_5 have exactly 10 reads
  prefix=file.path(new_dir(), "c10")
  r=run_ds(c("-p", "2", "-c", "10", "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  sig=read_sig(prefix)
  expected=main_status
  expected["chr1:clu_5"]="Success"
  stopifnot(identical(status_of(sig), expected))
  prefix=file.path(new_dir(), "c11")
  r=run_ds(c("-p", "2", "-c", "11", "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  stopifnot(identical(status_of(read_sig(prefix)), main_status))
}

test_min_samples_per_group <- function() {
  # three case samples of clu_6 reach 20 reads
  prefix=file.path(new_dir(), "g4")
  r=run_ds(c("-p", "2", "-g", "4", "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  expected=main_status
  expected["chr1:clu_6"]="Not enough valid samples"
  stopifnot(identical(status_of(read_sig(prefix)), expected))
  prefix=file.path(new_dir(), "g3")
  r=run_ds(c("-p", "2", "-g", "3", "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  stopifnot(same_outputs(prefix, default_prefix))
}

test_timeout_option <- function() {
  prefix=file.path(new_dir(), "t100")
  r=run_ds(c("-p", "2", "-t", "100", "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  i=which(r$stdout == "$timeout")
  stopifnot(length(i) == 1, r$stdout[i + 1] == "[1] 100")
  stopifnot(same_outputs(prefix, default_prefix))
}

exons <- data.frame(chr=c("chr1", "chr1", "chr1", "chr2", "chr2"),
                    start=c(500, 1200, 3300, 1200, 6100),
                    end=c(1000, 1300, 3400, 1300, 6200),
                    strand="+",
                    gene_name=c("GENE_A", "GENE_A", "GENE_B", "GENE_C", "GENE_D"),
                    stringsAsFactors=FALSE)

test_exon_file <- function() {
  exon_file=file.path(root, "exons.txt.gz")
  con=gzfile(exon_file, "w")
  write.table(exons, con, quote=FALSE, sep="\t", row.names=FALSE)
  close(con)
  prefix=file.path(new_dir(), "exons")
  r=run_ds(c("-p", "2", "-e", shQuote(exon_file), "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  sig=read_sig(prefix)
  stopifnot(identical(colnames(sig), c("cluster", "status", "loglr", "df", "p", "p.adjust", "genes")))
  genes=setNames(sig$genes, sig$cluster)
  # clu_1 matches GENE_A at both ends; clu_3 matches GENE_B; chr2's clu_7 matches GENE_C, not the chr1 coordinates
  stopifnot(identical(genes[c("chr1:clu_1", "chr1:clu_3", "chr2:clu_7")], c("chr1:clu_1"="GENE_A", "chr1:clu_3"="GENE_B", "chr2:clu_7"="GENE_C")))
  stopifnot(all(is.na(genes[c("chr1:clu_2", "chr1:clu_4", "chr1:clu_5", "chr1:clu_6")])))
  stopifnot(identical(status_of(sig)[names(main_status)], main_status))
  default_sig=read_sig(default_prefix)
  stopifnot(identical(sig[match(default_sig$cluster, sig$cluster), colnames(default_sig)], default_sig))
  stopifnot(identical(file_bytes(paste0(prefix, "_effect_sizes.txt")), file_bytes(paste0(default_prefix, "_effect_sizes.txt"))))
}

test_exon_file_missing_or_unusable <- function() {
  prefix=file.path(new_dir(), "missing_exons")
  missing=file.path(root, "no_such_exons.txt")
  r=run_ds(c("-p", "2", "-e", shQuote(missing), "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  stopifnot(any(grepl("does not exist", r$stderr)))
  stopifnot(same_outputs(prefix, default_prefix))
  # no gene_name column: the mapping fails, which the script reports as a warning
  bad=file.path(root, "bad_exons.txt")
  write.table(exons[, 1:4], bad, quote=FALSE, sep="\t", row.names=FALSE)
  prefix=file.path(new_dir(), "bad_exons")
  r=run_ds(c("-p", "2", "-e", shQuote(bad), "-o", shQuote(prefix), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r$status == 0)
  stopifnot(any(grepl("Error in", r$stderr, fixed=TRUE)))
  stopifnot(same_outputs(prefix, default_prefix))
}

# The CLI's numbers equal an in-process differential_splicing call with the same settings.
in_process <- function(...) {
  counts=read.table(counts_file, header=TRUE, check.names=FALSE)[, samples]
  res=differential_splicing(counts, rep(0:1, each=5), max_cluster_size=Inf, min_samples_per_intron=5,
                            min_samples_per_group=3, min_coverage=20, timeout=30, ...)
  list(sig=cluster_results_table(res), es=leaf_cutter_effect_sizes(res))
}

matches_in_process <- function(prefix, ref) {
  sig=read_sig(prefix)
  es=read_es(prefix)
  ok=sig$status == "Success"
  identical(sig$status, ref$sig$status) &&
    close_rel(sig$loglr[ok], ref$sig$loglr[ok], 1e-12) &&
    identical(es$intron, add_chr(ref$es$intron)) &&
    close_rel(es$logef, ref$es$logef, 1e-12) &&
    close_rel(es$ctrl, ref$es$baseline, 1e-12) &&
    close_rel(es$case, ref$es$perturbed, 1e-12)
}

test_init_and_seed <- function() {
  p1=file.path(new_dir(), "random7")
  p2=file.path(new_dir(), "random7")
  r1=run_ds(c("-p", "2", "--init", "random", "--seed", "7", "-o", shQuote(p1), shQuote(counts_file), shQuote(groups_file)))
  r2=run_ds(c("-p", "2", "--init", "random", "--seed", "7", "-o", shQuote(p2), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r1$status == 0, r2$status == 0)
  stopifnot(same_outputs(p1, p2))
  sig=read_sig(p1)
  check_tables(sig, read_es(p1), c("ctrl", "case"))
  stopifnot(identical(status_of(sig), main_status))
  p=p_of(sig)
  stopifnot(p["chr1:clu_1"] < 1e-6, p["chr1:clu_2"] > 0.5)
  stopifnot(matches_in_process(p1, in_process(init="random", seed=7)))
  # --seed defaults to 12345
  p3=file.path(new_dir(), "random_default_seed")
  r3=run_ds(c("-p", "2", "--init", "random", "-o", shQuote(p3), shQuote(counts_file), shQuote(groups_file)))
  stopifnot(r3$status == 0)
  stopifnot(matches_in_process(p3, in_process(init="random", seed=12345)))
  # --init defaults to smart
  stopifnot(matches_in_process(default_prefix, in_process(init="smart")))
}

test_group_order <- function() {
  # string labels keep the order of first appearance in the groups file
  g=groups_df[c(6:10, 1:5), ]
  gf=write_groups(g, file.path(root, "groups_case_first.txt"))
  prefix=file.path(new_dir(), "case_first")
  r=run_ds(c("-p", "2", "-o", shQuote(prefix), shQuote(counts_file), shQuote(gf)))
  stopifnot(r$status == 0)
  stopifnot(any(r$stdout == "Encoding as case =0, ctrl =1"))
  sig=read_sig(prefix)
  es=read_es(prefix)
  check_tables(sig, es, c("case", "ctrl"))
  stopifnot(identical(status_of(sig), main_status))
  stopifnot(es$logef[1] < 0, abs(es$case[1] - 0.70) < 0.05, abs(es$ctrl[1] - 0.30) < 0.05)
}

test_numeric_group_labels <- function() {
  # numeric labels are sorted: 1 (case, listed second) is encoded 0
  g=groups_df
  g$group=ifelse(g$group == "ctrl", 2, 1)
  gf=write_groups(g, file.path(root, "groups_numeric.txt"))
  prefix=file.path(new_dir(), "numeric")
  r=run_ds(c("-p", "2", "-o", shQuote(prefix), shQuote(counts_file), shQuote(gf)))
  stopifnot(r$status == 0)
  stopifnot(any(r$stdout == "Encoding as 1 =0, 2 =1"))
  sig=read_sig(prefix)
  es=read_es(prefix)
  check_tables(sig, es, c("1", "2"))
  stopifnot(identical(status_of(sig), main_status))
  stopifnot(es$logef[1] < 0, abs(es[["1"]][1] - 0.70) < 0.05, abs(es[["2"]][1] - 0.30) < 0.05)
}

# Confounder input: clu_1 usage set by batch (b1 20%, b2 50%, b3 80%, plus noise),
# clu_2 usage set by age (-1: 20%, 0: 50%, 1: 80%). Both are correlated with group.
batch <- c("b1", "b1", "b1", "b2", "b2", "b3", "b3", "b3", "b3", "b2")
age <- c(-1, -1, -1, -1, 0, 1, 1, 1, 1, 0)
conf_counts_file <- local({
  i1=c(b1=20, b2=50, b3=80)[batch] + c(-2, 1, 0, 2, -1, 1, -2, 0, 2, -1)
  i2=c(20, 50, 80)[age + 2]
  m=rbind(two_intron("chr1", 100, "clu_1", i1), two_intron("chr1", 900, "clu_2", i2))
  colnames(m)=samples
  write_counts(m, file.path(root, "confounders.counts.gz"))
})

run_conf <- function(extra, name) {
  g=if (is.null(extra)) groups_df else cbind(groups_df, extra)
  gf=write_groups(g, file.path(root, paste0("groups_", name, ".txt")))
  prefix=file.path(new_dir(), name)
  r=run_ds(c("-p", "2", "-o", shQuote(prefix), shQuote(conf_counts_file), shQuote(gf)))
  stopifnot(r$status == 0)
  sig=read_sig(prefix)
  check_tables(sig, read_es(prefix), c("ctrl", "case"))
  stopifnot(all(sig$status == "Success"))
  stopifnot(identical(sig$df, c(1L, 1L)))
  p_of(sig)
}

test_confounders <- function() {
  p=run_conf(NULL, "no_confounders")
  stopifnot(p["chr1:clu_1"] < 1e-3, p["chr1:clu_2"] < 1e-3)
  # categorical: three batches become two indicator columns
  p=run_conf(data.frame(batch=batch), "batch")
  stopifnot(p["chr1:clu_1"] > 0.5)
  # numeric
  p=run_conf(data.frame(age=age), "age")
  stopifnot(p["chr1:clu_2"] > 0.5)
  p=run_conf(data.frame(age=age, batch=batch), "age_batch")
  stopifnot(p["chr1:clu_1"] > 0.5, p["chr1:clu_2"] > 0.5)
}

check_stop <- function(args, pattern) {
  wd=new_dir()
  r=run_ds(c("-o", shQuote(file.path(wd, "out")), args))
  stopifnot(r$status != 0)
  stopifnot(any(grepl(pattern, r$stderr, fixed=TRUE)))
  stopifnot(length(list.files(wd)) == 0)
}

test_stops <- function() {
  missing=file.path(root, "no_such_file.txt")
  check_stop(c(shQuote(missing), shQuote(groups_file)), paste("File", missing, "does not exist"))
  check_stop(c(shQuote(counts_file), shQuote(missing)), paste("File", missing, "does not exist"))
  # smallest group (4) below the default -i 5, and below -g 5 once -i allows it
  g4=write_groups(groups_df[-10, ], file.path(root, "groups_4_case.txt"))
  check_stop(c(shQuote(counts_file), shQuote(g4)), "less than min_samples_per_intron")
  check_stop(c("-i", "4", "-g", "5", shQuote(counts_file), shQuote(g4)), "less than min_samples_per_group")
  g=groups_df
  g$group[c(5, 10)]="other"
  g3=write_groups(g, file.path(root, "groups_three.txt"))
  check_stop(c(shQuote(counts_file), shQuote(g3)), "length(group_names) == 2")
  g=groups_df
  g$group="ctrl"
  g1=write_groups(g, file.path(root, "groups_one.txt"))
  check_stop(c(shQuote(counts_file), shQuote(g1)), "length(group_names) == 2")
}

tests=list(test_default_run=test_default_run,
           test_unlisted_sample_ignored=test_unlisted_sample_ignored,
           test_default_prefix_and_one_thread=test_default_prefix_and_one_thread,
           test_max_cluster_size=test_max_cluster_size,
           test_min_samples_per_intron=test_min_samples_per_intron,
           test_min_coverage=test_min_coverage,
           test_min_samples_per_group=test_min_samples_per_group,
           test_timeout_option=test_timeout_option,
           test_exon_file=test_exon_file,
           test_exon_file_missing_or_unusable=test_exon_file_missing_or_unusable,
           test_init_and_seed=test_init_and_seed,
           test_group_order=test_group_order,
           test_numeric_group_labels=test_numeric_group_labels,
           test_confounders=test_confounders,
           test_stops=test_stops)
for (name in names(tests)) {
  cat(name, "\n")
  tests[[name]]()
}
unlink(root, recursive=TRUE)
