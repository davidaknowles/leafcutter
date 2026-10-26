library(leafcutter)
suppressMessages(library(doMC))
registerDoMC(2)

test_get_intron_meta <- function() {
  meta=get_intron_meta(c("chr1:100:200:clu_1", "chr1:100:350:clu_1", "2:1000:1501:clu_2_+", "chrX:7:9:clu_3_-"))
  stopifnot(is.data.frame(meta))
  stopifnot(identical(colnames(meta), c("chr", "start", "end", "clu", "middle")))
  stopifnot(identical(meta$chr, c("chr1", "chr1", "2", "chrX")))
  stopifnot(identical(meta$start, c(100, 100, 1000, 7)))
  stopifnot(identical(meta$end, c(200, 350, 1501, 9)))
  stopifnot(identical(meta$clu, c("clu_1", "clu_1", "clu_2_+", "clu_3_-")))
  stopifnot(identical(meta$middle, c(150, 225, 1250.5, 8)))
}

test_add_chr <- function() {
  stopifnot(identical(add_chr(c("1", "2", "X")), c("chr1", "chr2", "chrX")))
  stopifnot(identical(add_chr(c("chr1", "chr2", "chrX")), c("chr1", "chr2", "chrX")))
  stopifnot(identical(add_chr(c("1:clu_1", "10:clu_22_+")), c("chr1:clu_1", "chr10:clu_22_+")))
  stopifnot(identical(add_chr(c("chr1:clu_1", "chr10:clu_22_+")), c("chr1:clu_1", "chr10:clu_22_+")))
  stopifnot(identical(add_chr("17:100:200:clu_3"), "chr17:100:200:clu_3"))
}

test_bh <- function() {
  # non-NA p after capping at 1: 0.01, 0.04, 0.03, 1 (n=4)
  # sorted: 0.01*4/1=0.04, 0.03*4/2=0.06, 0.04*4/3=0.0533.., 1*4/4=1, then cumulative minimum from the top
  q=bh(c(0.01, NA, 0.04, 0.03, 2))
  stopifnot(is.na(q[2]))
  stopifnot(isTRUE(all.equal(q[-2], c(0.04, 0.16/3, 0.16/3, 1), tolerance=1e-14)))
  stopifnot(identical(bh(c(NA_real_, NA_real_)), c(NA_real_, NA_real_)))
  stopifnot(isTRUE(all.equal(bh(c(0.5, 0.2)), c(0.5, 0.4), tolerance=1e-14)))
}

test_sanitize_simplex <- function() {
  sanitize_simplex=leafcutter:::sanitize_simplex
  stopifnot(isTRUE(all.equal(sanitize_simplex(c(0, 0.5, 0.5)), c(1e-6, 0.5, 0.5)/(1 + 1e-6), tolerance=1e-14)))
  stopifnot(isTRUE(all.equal(sanitize_simplex(c(1, 0)), c(1 - 1e-6, 1e-6), tolerance=1e-14)))
  stopifnot(isTRUE(all.equal(sanitize_simplex(c(0.2, 0.3, 0.5)), c(0.2, 0.3, 0.5), tolerance=1e-14)))
  stopifnot(isTRUE(all.equal(sanitize_simplex(c(0, 1), eps=0.1), c(0.1, 0.9), tolerance=1e-14)))
  s=sanitize_simplex(c(0, 0, 3))
  stopifnot(abs(sum(s) - 1) < 1e-14)
}

test_beta_real <- function() {
  raw=matrix(c(0.5, 0.25, 0.25,
               0.2, 0.3, 0.5), 2, byrow=TRUE, dimnames=list(c("", "x"), c("i1", "i2", "i3")))
  b=beta_real(list(beta_raw=raw, beta_scale=c(2, -1)))
  # 2 * (raw[1,] - 1/3) and -1 * (raw[2,] - 1/3)
  expected=matrix(c(1/3, -1/6, -1/6,
                    2/15, 1/30, -1/6), 2, byrow=TRUE, dimnames=list(c("", "x"), c("i1", "i2", "i3")))
  stopifnot(isTRUE(all.equal(b, expected, tolerance=1e-14)))
  stopifnot(identical(dimnames(b), dimnames(raw)))
  # rows of beta_real sum to zero for any beta_raw on the simplex
  stopifnot(all(abs(rowSums(b)) < 1e-14))
}

test_map_clusters_to_genes <- function() {
  meta=get_intron_meta(c("chr1:100:200:clu_1", "chr1:100:300:clu_1",
                         "chr1:1000:1200:clu_2", "chr1:1000:1300:clu_2",
                         "chr1:5000:6000:clu_3", "chr1:5000:6500:clu_3",
                         "chr2:100:200:clu_4", "chr2:100:900:clu_4",
                         "chr3:100:200:clu_5", "chr3:150:200:clu_5"))
  exons=data.frame(
    chr=      c("chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr4", "chr3"),
    start=    c(50,     200,    300,    1300,   900,    900,    4900,   200,    6500),
    end=      c(100,    260,    400,    1400,   1000,   950,    5000,   260,    6600),
    strand="+",
    gene_name=c("A",    "A",    "B",    "C",    "C",    "D",    "E",    "F",    "G"),
    stringsAsFactors=FALSE)
  # clu_1: 5' ends at 100 (A), 3' ends at 200 (A) and 300 (B)
  # clu_2: 5' end 1000 matches C's exon end, 3' end 1300 matches C's exon start; C once
  # clu_3: 5000 and 6500 match exons on chr2 (E) and chr3 (G) only, so no gene
  # clu_4: 3' end 900 matches D's start on chr2; A's chr1 coordinates would match too
  # clu_5: 3' end 200 matches exon starts on chr1 (A) and chr4 (F), none on chr3
  r=map_clusters_to_genes(meta, exons)
  stopifnot(identical(class(r), "data.frame"))
  stopifnot(identical(colnames(r), c("clu", "genes")))
  stopifnot(setequal(r$clu, c("chr1:clu_1", "chr1:clu_2", "chr2:clu_4")))
  stopifnot(!any(duplicated(r$clu)))
  genes=setNames(strsplit(r$genes, ","), r$clu)
  stopifnot(identical(sort(genes[["chr1:clu_1"]]), c("A", "B")))
  stopifnot(identical(genes[["chr1:clu_2"]], "C"))
  stopifnot(identical(genes[["chr2:clu_4"]], "D"))
}

tests=list(test_get_intron_meta=test_get_intron_meta,
           test_add_chr=test_add_chr,
           test_bh=test_bh,
           test_sanitize_simplex=test_sanitize_simplex,
           test_beta_real=test_beta_real,
           test_map_clusters_to_genes=test_map_clusters_to_genes)
for (name in names(tests)) {
  cat(name, "\n")
  tests[[name]]()
}
