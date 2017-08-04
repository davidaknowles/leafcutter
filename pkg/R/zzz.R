.onLoad <- function(libname, pkgname) {

  # if (!("methods" %in% .packages())) attachNamespace("methods")

  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)

  # .onLoad <- function(libname, pkgname) { Rcpp::loadRcppModules() }

}

