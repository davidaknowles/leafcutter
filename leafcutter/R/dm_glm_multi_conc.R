

#' Dirichlet multinomial GLM likelihood ratio test for a single cluster
#'
#' @param xFull [samples] x [covariates] matrix for the alternative model
#' @param xNull [samples] x [covariates] matrix for the null model
#' @param y [samples] x [introns] matrix of intron usage counts
#' @param concShape Gamma shape parameter for concentration parameter
#' @param concRate Gamma rate parameter for concentration parameter
#' @param robust Whether to include an outlier model (i.e. use dm_glm_multi_conc_robust rather than dm_glm_multi_conc)
#' @param outlier_prior_a Only used for robust model. The outlier probability outlier_prob ~ Beta(outlier_prior_a,outlier_prior_b)
#' @param outlier_prior_b Only used for robust model. The outlier probability outlier_prob ~ Beta(outlier_prior_a,outlier_prior_b)
#' @param fit_null Optionally cache the fitted null model to save repeatedly fitting the null for each cis-SNP when sQTL mapping)
#' @param debug Whether to give verbose output from rstan.
#' @param init Can be one of {"smart", "random"}. smart uses an method of moments estimator to get a reasonable initialization. The seed for "random" can be set through the ... arguments passed to rstan::optimizing.
#' @param smart_init_regularizer Used to protect against colinear covariates. 
#' @param ... will be passed on the rstan::optimizing, so can be used for example to set the algorithm used (default is LBFGS) or the random seed if random initialization is requested. 
#' @importFrom rstan optimizing
#' @import foreach
#' @import dplyr
#' @export
dirichlet_multinomial_anova_mc <- function(xFull,xNull,y,concShape=1.0001,concRate=1e-4, robust=T, outlier_prior_a=1.01, outlier_prior_b=100, fit_null=NULL, debug=F, init="smart", smart_init_regularizer=0.001, ...) {
  K=ncol(y)
  
  model_to_use=if (robust) stanmodels$dm_glm_robust else stanmodels$dm_glm_multi_conc
  
  dat_null=list(N=nrow(xNull), K=K, P=ncol(xNull), y=y, x=xNull, concShape=concShape,concRate=concRate, outlier_prior_a=outlier_prior_a, outlier_prior_b=outlier_prior_b)
  
  if (init=="smart") {
    y_norm = sweep(y+1,1,rowSums(y+1),"/") %>% log()
    beta_mm = solve( t(xNull) %*% xNull + smart_init_regularizer * diag(ncol(xNull)), t(xNull) %*% y_norm )
    beta_norm = sweep(beta_mm, 1, rowMeans(beta_mm), "-")
    beta_scale = unlist(lapply(seq_len(nrow(beta_norm)), function(i) {
      beta_row = beta_norm[i,]
      up=beta_row[which.max(abs(beta_row))] / (1-1/K)
      down=beta_row[which.max(-sign(up)*beta_row)] / (1/K)
      up-down
    }))
    beta_raw=sweep(beta_norm,1,beta_scale+1e-20,"/") + 1/K
    beta_raw=sweep(beta_raw,1,rowSums(beta_raw),"/")
    init=list(beta_scale=array(beta_scale), beta_raw=beta_raw, conc=rep(10.0,K))
  }
  
  # fit null model
  if (is.null(fit_null)) 
    fit_null=stan_optimizing(model_to_use, data=dat_null, init=init, as_vector=F, verbose=debug, ...)

  colnames(fit_null$par$beta_raw)=colnames(y)
  rownames(fit_null$par$beta_raw)=colnames(xNull)
  
  dat_full=dat_null
  dat_full$P=ncol(xFull)
  dat_full$x=xFull

  # beta_raw is PxK
  init=list(beta_raw=matrix(1e-4,ncol(xFull),K), beta_scale=rep(1,ncol(xFull)), conc=fit_null$par$conc, outlier_prob=fit_null$par$outlier_prob)

  # beta_raw must live _in_ the simplex
  beta_raw_sanitized=fit_null$par$beta_raw
  beta_raw_sanitized[beta_raw_sanitized<1e-6]=1e-6
  beta_raw_sanitized[beta_raw_sanitized>(1.0-1e-6)]=(1.0-1e-6)
  init$beta_raw[1:ncol(xNull),]=beta_raw_sanitized
  init$beta_raw=sweep(init$beta_raw, 1, rowSums(init$beta_raw), "/") 
  
  init$beta_scale[1:ncol(xNull)]=fit_null$par$beta_scale
  stopifnot(all(is.finite(unlist(init))))
  # fit fit model
  fit_full=stan_optimizing(model_to_use, data=dat_full, init=init, as_vector=F, verbose=debug,...)

  colnames(fit_full$par$beta_raw)=colnames(y)
  rownames(fit_full$par$beta_raw)=colnames(xFull)
  
  loglr=fit_full$value-fit_null$value
  df=( ncol(xFull)-ncol(xNull) )*(K-1)
  
  refit_null_flag=F
  
  lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df )
  if (lrtp < .001) {
    init=fit_full$par
    init$beta_raw=init$beta_raw[seq_len(dat_null$P),,drop=F]
    init$beta_raw[init$beta_raw<1e-6]=1e-6
    init$beta_raw[init$beta_raw>(1.0-1e-6)]=(1.0-1e-6)
    init$beta_raw=sweep(init$beta_raw, 1, rowSums(init$beta_raw), "/") 
    init$beta_scale=as.array(init$beta_scale[seq_len(dat_null$P)])
    refit_null=stan_optimizing(model_to_use, data=dat_null, init=init, as_vector=F, verbose=debug,...)
    if (refit_null$value > fit_null$value) {
      refit_null_flag=T
      fit_null=refit_null
      loglr=fit_full$value-fit_null$value
    }
  }
  
  list( loglr=loglr, df=df, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=df ), fit_null=fit_null, fit_full=fit_full, refit_null_flag=refit_null_flag)
}

# rstan::optimizing(object, data=data, init=init, as_vector=as_vector,
# verbose=verbose, ...) for the models fitted above. rstan re-parses the
# model's C++ code and rebuilds its parameter names and dimensions on every
# call; with rstan 2.21.8, quiet output and nothing but a seed in `...`, this
# takes rstan's steps itself and keeps that work per model and, since the
# parameter dimensions of these models depend only on K and P, per (K, P).
stan_optimizing=function(object, data, init, as_vector, verbose, ...) {
  model=stan_model_info(object)
  dots=as.list(substitute(list(...)))[-1]
  if (!model$fast || !identical(verbose, FALSE) || length(dots) > 1 ||
      (length(dots) == 1 && !identical(names(dots), "seed")) || !has_stan_data(data, model$data_names))
    return(rstan::optimizing(object, data=data, init=init, as_vector=as_vector, verbose=verbose, ...))
  data=try(preprocess_stan_data(data[model$data_names]))
  if (inherits(data, "try-error")) {
    message("failed to preprocess the data; optimization not done")
    return(invisible(list(stanmodel=object)))
  }
  # rstan forces its default seed here
  seed=if (length(dots)) ..1 else sample.int(.Machine$integer.max, 1)
  fit=optimize_with_call(model, data, init, seed)
  if (is.null(fit))
    fit=optimize_with_module(model, data, init, seed)
  if (is.null(fit))
    return(invisible(list(stanmodel=object)))
  optim=fit$optim
  optim$return_code=attr(optim, "return_code")
  # with the call rstan's optimizing() method gives it
  if (optim$return_code != 0)
    warning(simpleWarning("non-zero return code in optimizing", quote(.local(object, ...))))
  attr(optim, "return_code")=NULL
  names(optim$par)=fit$shape$par_names
  theta=rstan_internal("rstan_relist")(optim$par, fit$shape$skeleton)
  theta_tilde=matrix(fit$constrain(theta), nrow=1)
  colnames(theta_tilde)=names(optim$par)
  optim$theta_tilde=theta_tilde
  if (!as_vector)
    optim$par=theta
  optim
}

# optimize_with_module() as one .Call of the model's compiled routine, for a
# model that has one and a (K, P) seen before, with an init and seed that
# rstan's handling leaves without side effects. NULL otherwise, and when the
# C++ code throws: optimize_with_module() then repeats the work and reports
# the error as rstan does.
optimize_with_call=function(model, data, init, seed) {
  if (is.null(model$optimizing))
    return(NULL)
  shape=model$shapes[[paste(data$K, data$P)]]
  simple_init=!is.function(init) && (is.list(init) || is.character(init) || is.numeric(init))
  simple_seed=(is.integer(seed) || is.double(seed)) && !is.object(seed) && length(seed) == 1 &&
    !is.na(seed) && abs(seed) < .Machine$integer.max
  if (is.null(shape) || !simple_init || !simple_seed)
    return(NULL)
  if (is.numeric(init))
    init=as.character(init)
  args=list(init=init, seed=rstan_internal("check_seed")(seed, warn=1), method="optim",
            algorithm="LBFGS", refresh=0L, save_iterations=FALSE)
  out=.Call(model$optimizing, data, as.integer(seed), args, shape$skeleton)
  if (is.null(out))
    return(NULL)
  list(optim=out[[1]], shape=shape, constrain=function(theta) out[[2]])
}

# rstan::optimizing()'s steps from building the C++ sampler to running the
# optimizer. NULL, after rstan's message, where rstan gives up.
optimize_with_module=function(model, data, init, seed) {
  sampler=try(new(model$module, data, as.integer(seed), model$cxxfun))
  if (inherits(sampler, "try-error")) {
    message("failed to create the optimizer; optimization not done")
    return(NULL)
  }
  shape=stan_shape(model, data, sampler)
  if (is.numeric(init))
    init=as.character(init)
  if (is.function(init))
    init=init()
  if (!is.list(init) && !is.character(init)) {
    message("wrong specification of initial values")
    return(NULL)
  }
  seed=rstan_internal("check_seed")(seed, warn=1)
  if (is.null(seed))
    return(NULL)
  args=list(init=init, seed=seed, method="optim", algorithm="LBFGS")
  dotlist=list(refresh=0L, save_iterations=FALSE)
  # The sampler calls read as in rstan: the call of an error they raise ends up
  # in the cluster's status.
  optim=sampler$call_sampler(c(args, dotlist))
  constrain=function(theta) {
    theta=sampler$unconstrain_pars(theta)
    sampler$constrain_pars(theta)
  }
  list(optim=optim, shape=shape, constrain=constrain)
}

# Per Stan model, what rstan::optimizing() would otherwise work out on every
# call; filled on first use in each process.
stan_optimizing_cache=new.env(parent=emptyenv())

stan_model_info=function(object) {
  model=stan_optimizing_cache[[object@model_name]]
  if (!is.null(model))
    return(model)
  model=new.env(parent=emptyenv())
  # the steps below are those of rstan 2.21.8's optimizing()
  model$fast=identical(unname(getNamespaceVersion("rstan")), "2.21.8")
  if (model$fast) {
    model$cxxfun=rstan_internal("grab_cxxfun")(object@dso)
    model$module=object@mk_cppmodule(object)
    model$shapes=new.env(parent=emptyenv())
    cppcode=object@model_cpp$model_cppcode
    model$data_names=stan_data_names(cppcode)
    newer_stan=rstan_internal("test_221")(cppcode) && model$module@constructors[[1]]$nargs == 2L
    model$fast=!newer_stan && all(vapply(model$data_names, rstan_internal("is_legal_stan_vname"), logical(1)))
    if (model$fast && identical(object, stanmodels$dm_glm_multi_conc))
      model$optimizing=`_leafcutter_dm_glm_multi_conc_optimizing`
  }
  assign(object@model_name, model, envir=stan_optimizing_cache)
  model
}

rstan_internal=function(name) get(name, envir=asNamespace("rstan"))

# The data names rstan's parse_data() finds in a model's C++ code.
stan_data_names=function(cppcode) {
  cppcode=scan(what=character(), sep="\n", quiet=TRUE, text=cppcode)
  private=grep("^private:$", cppcode) + 1L
  public=grep("^public:$", cppcode) - 1L
  objects=gsub("^.* ([0-9A-Za-z_]+).*;.*$", "\\1", cppcode[private:public])
  in_data=grep("context__.vals_", cppcode, fixed=TRUE, value=TRUE)
  in_data=gsub("^.*\\(\"(.*)\"\\).*;$", "\\1", in_data)
  intersect(objects, in_data)
}

# TRUE when rstan's parse_data() would return data[data_names]: data is a plain
# list holding each name once. parse_data() skips a value that is not numeric,
# logical or a list and looks for the name in the calling frames instead.
has_stan_data=function(data, data_names) {
  is.list(data) && !is.object(data) && !anyDuplicated(names(data)) && all(data_names %in% names(data)) &&
    all(vapply(data[data_names], typeof, character(1)) %in% c("integer", "double", "logical", "list"))
}

# rstan's data_preprocess() for data whose names are legal Stan names.
preprocess_stan_data=function(data) {
  data_list2array=rstan_internal("data_list2array")
  real_is_integer=rstan_internal("real_is_integer")
  names=names(data)
  data=lapply(names, FUN=function(name) {
    x=data[[name]]
    if (is.data.frame(x)) {
      x=data.matrix(x)
    } else if (is.list(x)) {
      x=data_list2array(x)
    } else if (is.logical(x)) {
      storage.mode(x)="integer"
    }
    if (any(is.na(x)))
      stop("Stan does not support NA (in ", name, ") in data")
    if (!is.numeric(x)) {
      warning("data with name ", name, " is not numeric and not used")
      return(NULL)
    }
    if (is.integer(x))
      return(x)
    if (all(abs(x) < .Machine$integer.max) && real_is_integer(x))
      storage.mode(x)="integer"
    x
  })
  names(data)=names
  data[!sapply(data, is.null)]
}

# The parameter skeleton and flat parameter names rstan::optimizing() gets
# from the sampler, kept per (K, P).
stan_shape=function(model, data, sampler) {
  key=paste(data$K, data$P)
  shape=model$shapes[[key]]
  if (is.null(shape)) {
    pars=sampler$param_names()
    keep=which(pars != "lp__")
    fnames=sampler$param_fnames_oi()
    shape=list(skeleton=rstan_internal("create_skeleton")(pars[keep], sampler$param_dims()[keep]),
               par_names=fnames[-length(fnames)])
    assign(key, shape, envir=model$shapes)
  }
  shape
}
