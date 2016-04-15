myskel=function (name = "anRpackage", list = character(), environment = .GlobalEnv, 
          path = ".", force = FALSE, code_files = character(), stan_files = character()) 
{
  if (length(stan_files) > 0 && !all(grepl("\\.stan$", stan_files))) {
    stop("all files named in 'stan_files' must end with a .stan extension")
  }
  utils::package.skeleton(name=name, code_files = code_files)
  if (R.version$major < 3 || (R.version$major == 3 && R.version$minor < 
                                2.2)) {
    warning("rstan.package.skeleton is only fully operational with R >= 3.2.2", 
            "Follow the package skeleton of the rstanarm package on GitHub")
    return(invisible(NULL))
  }
  DIR <- file.path(path, name)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/cleanup", 
                destfile = file.path(DIR, "cleanup"), quiet = TRUE)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/cleanup.win", 
                destfile = file.path(DIR, "cleanup.win"), quiet = TRUE)
  cat("cleanup*", file = file.path(DIR, ".Rbuildignore"), sep = "\n")
  TOOLS <- file.path(DIR, "tools")
  dir.create(TOOLS)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/tools/make_cpp.R", 
                destfile = file.path(TOOLS, "make_cpp.R"), quiet = TRUE)
  EXEC <- file.path(DIR, "exec")
  dir.create(EXEC)
  file.create(file.path(EXEC, "common_functions.txt"))
  file.copy(stan_files, EXEC)
  SRC <- file.path(DIR, "src")
  dir.create(SRC, showWarnings = FALSE)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/src/Makevars", 
                destfile = file.path(SRC, "Makevars"), quiet = TRUE)
  R <- file.path(DIR, "R")
  dir.create(R, showWarnings = FALSE)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/R/stanmodels.R", 
                destfile = file.path(R, "stanmodels.R"), quiet = TRUE)
  cat(".onLoad <- function(libname, pkgname) { Rcpp::loadRcppModules() }", 
      file = file.path(R, "zzz.R"), sep = "\n", append = TRUE)
  if (length(stan_files) == 0) 
    module_names <- "NAME"
  else module_names <- paste0("stan_fit4", sub("\\.stan$", 
                                               "", basename(stan_files)), "_mod")
  cat("Depends: R (>= 3.0.2), Rcpp (>= 0.11.0)", "Imports: rstan (>= 2.8.1)", 
      "LinkingTo: StanHeaders (>= 2.8.0), rstan (>= 2.8.1), BH (>= 1.58.0), Rcpp (>= 0.11.0), RcppEigen", 
      file = file.path(DIR, "DESCRIPTION"), sep = "\n", append = TRUE)
  cat("RcppModules: ", paste(module_names, collapse = ", "), 
      "\n", file = file.path(DIR, "DESCRIPTION"), append = TRUE)
  cat("\n Stan specific notes:", "If you add any additional .stan files to the exec/ directory, ", 
      "be sure to add an entry in the RcppModules: line of DESCRIPTION.", 
      "You can put into exec/functions.txt any function that is needed by any .stan file, ", 
      "and in that case no .stan file should have its own functions{} block.", 
      "The precompiled stanmodel objects will appear in a named list called 'stanmodels'.", 
      "The 'cleanup' and 'cleanup.win' scripts in the root of the directory must be made executable.", 
      file = file.path(DIR, "Read-and-delete-me"), sep = "\n", 
      append = TRUE)
  return(invisible(NULL))
}

myskel(name="leafcutter", code_files = c("bb_glm.R","differential_splicing.R","dm_glm_multi_conc.R","dm_glm.R","splicing_qtl.R","../junction_plot//junction_plot.R"), stan_files = c("bb_glm.stan","bb_glm_fix_conc.stan","dm_glm_multi_conc.stan","dm_glm.stan"))
