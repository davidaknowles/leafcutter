#' Parse output of differential_splicing
#'
#' Parse output of \code{\link{differential_splicing}} and make a per cluster results table
#'
#' @param results From \code{\link{differential_splicing}}
#' @return Data.frame with columns status, log likelihood ratio, degrees of freedom, p-value
#' @import dplyr
#' @export
cluster_results_table=function(results) {
  tab=if (length(results) > 0 && all(vapply(results, is_simple_result, logical(1)))) {
    cluster_results_columns(unname(results))
  } else {
    foreach(res=results, .combine=bind_rows) %dopar% {
      if ( is.character(res) | ("error" %in% class(res)) )
        data.frame(status=as.character(res), loglr=NA, df=NA, p=NA, stringsAsFactors = F) else
        data.frame(status="Success", loglr=res$loglr, df=res$df, p=res$lrtp, stringsAsFactors = F)
    }
  }
  tab %>% mutate( cluster=names(results),
                  status   = gsub("\n", "", status),  # strip newlines from any reported errors
                  p.adjust = p.adjust(p, method="fdr") )
}

# A status string, or a fit whose loglr, df and lrtp are plain numbers: what
# differential_splicing() returns for a cluster that did not fail.
is_simple_result=function(res) {
  if (is.character(res)) return(length(res) == 1 && is.null(attributes(res)))
  is.list(res) && !is.object(res) && is_plain_number(res$loglr) && is_plain_number(res$df) && is_plain_number(res$lrtp)
}

is_plain_number=function(x) is.double(x) && length(x) == 1 && is.null(attributes(x))

# The rows cluster_results_table() binds one result at a time, built a column
# at a time. Without a tested cluster the numeric columns stay logical, as
# bind_rows() leaves them.
cluster_results_columns=function(results) {
  tested=!vapply(results, is.character, logical(1))
  fit_values=function(get) {
    if (!any(tested)) return(rep(NA, length(results)))
    vapply(results, function(res) if (is.character(res)) NA_real_ else get(res), numeric(1))
  }
  data.frame(status=vapply(results, function(res) if (is.character(res)) res else "Success", character(1)),
             loglr=fit_values(function(res) res$loglr),
             df=fit_values(function(res) res$df),
             p=fit_values(function(res) res$lrtp),
             stringsAsFactors = F)
}

#' Convert K-1 representation of parameters to real
#'
#' @param r A parameter object from fitting dm_glm or dm_glm_multi_conc. Must have members beta_raw and beta_scale.
#' @return Coefficient matrix.
#' @export
beta_real=function(r)
  sweep(r$beta_raw - 1.0/ncol(r$beta_raw), 1, r$beta_scale, "*")

#' Per intron effect sizes
#'
#' @param results From \code{\link{differential_splicing}}
#' @return data.frame with a row for every tested intron and columns: intron, log effect size, baseline proportions, proportions in the second condition, and resulting deltaPSI
#' @import dplyr
#' @export
leaf_cutter_effect_sizes=function(results) {
  normalize=function(g) { g/sum(g) }
  softmax=function(g) normalize(exp(g))
  to_psi=function(b,conc) { normalize(softmax(b)*conc) }
  cluster_effect_sizes=function(res) {
    if ( is.character(res) | ("error" %in% class(res)) ) return(NULL)
    beta=beta_real( res$fit_full$par )
    columns=list( intron=colnames(beta),
                  logef=beta[2,],
                  baseline=to_psi(beta[1,],res$fit_full$par$conc),
                  perturbed=to_psi(beta[1,]+beta[2,],res$fit_full$par$conc) )
    if (is_plain_effect_sizes(columns)) columns else effect_size_frame(columns)
  }
  tables=tryCatch(lapply(results, cluster_effect_sizes), error=function(e) NULL)
  tab=if (is.null(tables)) {
    # raises the error as the per-cluster fold reports it
    foreach(res=results, .combine = bind_rows) %do% {
      cluster_effect_sizes(res)
    }
  } else bind_effect_sizes(tables)
  tab %>%
  mutate( deltapsi=perturbed-baseline)
}

# TRUE when data.frame() would only drop the names and 1-d dims of a cluster's
# effect sizes, and take its introns as row names without a warning.
is_plain_effect_sizes=function(columns) {
  intron=columns$intron
  plain_proportion=function(p) is.double(p) && length(p) == length(intron) &&
    (identical(attributes(p), list(dim=length(intron))) || identical(attributes(p), list(names=intron)))
  is.character(intron) && length(intron) > 0 && is.null(attributes(intron)) && !anyNA(intron) &&
    all(nzchar(intron)) && !anyDuplicated(intron) &&
    is.double(columns$logef) && identical(attributes(columns$logef), list(names=intron)) &&
    plain_proportion(columns$baseline) && plain_proportion(columns$perturbed)
}

effect_size_frame=function(columns) {
  data.frame( intron=columns$intron,
              logef=columns$logef,
              baseline=columns$baseline,
              perturbed=columns$perturbed,
              stringsAsFactors = F )
}

# The per-cluster effect sizes (NULL for an untested cluster) bound as
# foreach(.combine=bind_rows) folds them. That fold leaves a lone result
# unbound and starts from a 0x0 tibble when the first two results are NULL;
# when the introns are distinct, one bind of all the columns gives the same
# table.
bind_effect_sizes=function(tables) {
  tested=!vapply(tables, is.null, logical(1))
  if (length(tables) >= 2 && any(tested) && !any(vapply(tables[tested], is.data.frame, logical(1)))) {
    column=function(name) unlist(lapply(tables[tested], function(tab) tab[[name]]), use.names=F)
    intron=column("intron")
    if (!anyDuplicated(intron)) {
      tab=structure(list(intron=intron, logef=column("logef"), baseline=column("baseline"), perturbed=column("perturbed")),
                    row.names=intron, class="data.frame")
      tab=tryCatch(if (!tested[1] && !tested[2]) bind_rows(bind_rows(NULL, NULL), tab) else bind_rows(tab),
                   error=function(e) NULL)
      if (!is.null(tab)) return(tab)
    }
  }
  tables=lapply(tables, function(tab) if (is.null(tab) || is.data.frame(tab)) tab else effect_size_frame(tab))
  foreach(tab=tables, .combine = bind_rows) %do% tab
}


#' Perform pairwise differential splicing analysis.
#'
#' Parallelization across tested clusters is achieved using foreach/doMC, so the number of threads that will be used is determined by the cores argument passed to registerDoMC.
#'
#' @param counts An [introns] x [samples] matrix of counts. The rownames must be of the form chr:start:end:cluid. If the counts file comes from the leafcutter clustering code this should be the case already.
#' @param x A [samples] numeric vector, should typically be 0s and 1s, although in principle scaling shouldn't matter.
#' @param confounders A [samples] x [confounders] numeric matrix to be controlled for in the GLM. Factors should already have been converted to a 1-of-(K-1) encoding, e.g. using model.matrix (see scripts/leafcutter_ds.R for how to do this). Can be NULL, implying no covariates are controlled for.
#' @param max_cluster_size Don't test clusters with more introns than this
#' @param min_samples_per_intron Ignore introns used (i.e. at least one supporting read) in fewer than n samples
#' @param min_samples_per_group Require this many samples in each group to have at least min_coverage reads
#' @param min_coverage Require min_samples_per_group samples in each group to have at least this many reads
#' @param timeout Maximum time (in seconds) allowed for a single optimization run
#' @param robust Whether to use the robust model (explicitly models outliers). Generally not required/recommended for differential splicing.
#' @param debug If true writes more output
#' @param init One of 'smart' (default) or 'random'. If 'random' you can pass an additional arg "seed" for reproducibility.
#' @return A per cluster list of results. Clusters that were not tested will be represented by a string saying why.
#' @import foreach
#' @importFrom R.utils withTimeout
#' @export
differential_splicing=function(counts, x, confounders=NULL, max_cluster_size=10, min_samples_per_intron=5, min_samples_per_group=4, min_coverage=20, timeout=10, robust=F, debug=F, init="smart", ...) {

  stopifnot(ncol(counts)==length(x))
  stopifnot(!any(is.na(x)))

  introns=get_intron_meta(rownames(counts))
  cluster_ids=paste(introns$chr,introns$clu,sep = ":")

  cluster_sizes=as.data.frame(table(cluster_ids))
  clu_names=as.character(cluster_sizes$cluster_ids)
  cluster_sizes=cluster_sizes$Freq
  names(cluster_sizes)=clu_names

  cluster_junctions=split(seq_along(cluster_ids), factor(cluster_ids, levels=clu_names))
  count_matrix=as.matrix(counts)

  # R.utils::withTimeout() checks its arguments on every call, which can take
  # longer than a fit. Once it has accepted `timeout` in this process, later
  # fits set the same time limit directly.
  timeout_checked=F
  fit_with_timeout=function(fit) {
    if (timeout_checked) return(with_time_limit(fit, timeout))
    res=R.utils::withTimeout(fit, timeout=timeout, onTimeout="silent")
    timeout_checked <<- is.numeric(timeout) && !is.object(timeout) && length(timeout) == 1 && !is.na(timeout) && timeout >= 0
    res
  }

  test_cluster=function(i) {
    if (cluster_sizes[i] > max_cluster_size)
      return("Too many introns in cluster")
    junctions_in_cluster=cluster_junctions[[i]]
    if (length(junctions_in_cluster) <= 1)
      return("<=1 junction in cluster")
    cluster_counts=t(count_matrix[ junctions_in_cluster, , drop=F ])
    sample_totals=rowSums(cluster_counts)
    samples_to_use=sample_totals>0
    if (sum(samples_to_use)<=1 )
      return("<=1 sample with coverage>0")
    sample_totals=sample_totals[samples_to_use]
    if (sum(sample_totals>=min_coverage)<=1)
      return("<=1 sample with coverage>min_coverage")
    x_subset=x[samples_to_use]
    cluster_counts=cluster_counts[samples_to_use,]
    introns_to_use=colSums(cluster_counts>0)>=min_samples_per_intron # only look at introns used by at least 5 samples
    if ( sum(introns_to_use)<2 )
      return("<2 introns used in >=min_samples_per_intron samples")
    cluster_counts=cluster_counts[,introns_to_use]
    ta=table(x_subset[sample_totals>=min_coverage])
    if (sum(ta >= min_samples_per_group)<2) # at least two groups have this (TODO: continuous x)
      return("Not enough valid samples")
    xFull=cbind(1,x_subset)
    xNull=xFull[,1,drop=F]
    if (!is.null(confounders)) {
        ch=confounders[samples_to_use,,drop=F]
        ch=ch[ , apply(ch,2,sd)>0.0, drop=F ]
        xFull=cbind(xFull,ch)
        xNull=cbind(xNull,ch)
    }
    res <- fit_with_timeout(
      dirichlet_multinomial_anova_mc(xFull,xNull,cluster_counts,robust=robust,debug=debug,init=init,...)
    )
    if (is.null(res)) "timeout" else res
  }

  if (!debug) {
    zz <- file( "/dev/null", open = "wt")
    sink(zz)
    sink(zz, type = "message")
  }

  # Fit times vary a lot between clusters, so clusters go to the workers in
  # chunks that shrink towards the end rather than in one block per worker.
  chunks=cluster_chunks(length(clu_names), max(1, getDoParWorkers()))
  chunk_results=foreach (chunk=chunks, .errorhandling = "pass", .options.multicore=list(preschedule=FALSE)) %dopar% {
    lapply(chunk, function(i) tryCatch(test_cluster(i), error=function(e) e))
  }
  results=vector("list", length(clu_names))
  for (k in seq_along(chunks)) {
    chunk_result=chunk_results[[k]]
    # a chunk whose worker died gives its NULL or error to each of its clusters
    complete=is.list(chunk_result) && !is.object(chunk_result) && length(chunk_result) == length(chunks[[k]])
    results[chunks[[k]]]=if (complete) chunk_result else list(chunk_result)
  }

  if (!debug) {
    sink(type="message")
    sink()
  }

    names(results)=clu_names

    statuses=cluster_results_table(results)$status

  cat("Differential splicing summary:\n")
  print(as.data.frame(table(statuses)))

  results
}

# Contiguous chunks of 1..n, each a 1/(2 * workers) share of what is left but
# at least 32 clusters.
cluster_chunks=function(n, workers) {
  chunks=list()
  start=1
  while (start <= n) {
    size=max(32, ceiling((n - start + 1) / (2 * workers)))
    chunks[[length(chunks) + 1]]=start:min(n, start + size - 1)
    start=start + size
  }
  chunks
}

# The time limit R.utils::withTimeout(expr, timeout=timeout, onTimeout="silent")
# sets, without its argument checks. Errors go through withTimeout()'s own
# handler, which returns NULL for a time limit and re-throws anything else.
with_time_limit=function(expr, timeout) {
  setTimeLimit(cpu=timeout, elapsed=timeout, transient=TRUE)
  on.exit(setTimeLimit(cpu=Inf, elapsed=Inf, transient=FALSE))
  tryCatch(expr, error=function(e) R.utils::withTimeout(stop(e), timeout=Inf, onTimeout="silent"))
}
