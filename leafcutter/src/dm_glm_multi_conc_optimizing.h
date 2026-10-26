#ifndef LEAFCUTTER_DM_GLM_MULTI_CONC_OPTIMIZING_H
#define LEAFCUTTER_DM_GLM_MULTI_CONC_OPTIMIZING_H

// .Call routine for stan_optimizing() in R/dm_glm_multi_conc.R: the C++ work
// of one rstan::optimizing() call of dm_glm_multi_conc (the stan_fit
// constructor, call_sampler, unconstrain_pars and constrain_pars, with the
// arguments rstan passes them) without building the Rcpp module object.
//
// Returns list(<call_sampler's result>, <constrain_pars' result>), or NULL
// when any step throws; R then takes rstan's own path, which repeats the work
// and reports the error as rstan does. Included once, from
// stanExports_dm_glm_multi_conc.cc, after the model.

#include <Rcpp.h>
#include <stan/math/version.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#if STAN_MATH_MAJOR == 3 && STAN_MATH_MINOR == 0

namespace leafcutter {

// What a stan_fit object holds for one call. It lives behind an external
// pointer: a time limit fired in the optimizer's interrupt check longjmps past
// C++ destructors, and the pointer's finalizer then frees it.
struct dm_optimizer {
    rstan::io::rlist_ref_var_context data;
    stan_model model;
    boost::random::ecuyer1988 base_rng;
    Rcpp::List args_list;
    rstan::stan_args args;
    Rcpp::List holder;
    // call_sampler passes these; the optimizer reads neither
    std::vector<size_t> qoi_idx;
    std::vector<std::string> fnames_oi;

    dm_optimizer(SEXP data_, SEXP seed, SEXP args_)
        : data(data_), model(data, Rcpp::as<boost::uint32_t>(seed), &rstan::io::rcout),
          base_rng(Rcpp::as<boost::uint32_t>(seed)), args_list(args_), args(args_list) {}
};

inline void free_dm_optimizer(SEXP handle) {
    delete static_cast<dm_optimizer*>(R_ExternalPtrAddr(handle));
    R_ClearExternalPtr(handle);
}

// rstan:::rstan_relist(par, skeleton): par's values, in order, in a copy of
// skeleton (a list of double arrays and scalars).
inline Rcpp::List relist(SEXP par, SEXP skeleton) {
    Rcpp::List theta(Rf_duplicate(skeleton));
    const double* from = REAL(par);
    R_xlen_t left = Rf_xlength(par);
    for (R_xlen_t i = 0; i < theta.size(); ++i) {
        SEXP element = theta[i];
        const R_xlen_t n = Rf_xlength(element);
        if (TYPEOF(element) != REALSXP || n > left)
            throw std::invalid_argument("parameter skeleton does not match the model");
        std::copy(from, from + n, REAL(element));
        from += n;
        left -= n;
    }
    return theta;
}

inline Rcpp::List optimize(dm_optimizer& opt, SEXP skeleton) {
    // call_sampler
    int ret = rstan::command(opt.args, opt.model, opt.holder, opt.qoi_idx, opt.fnames_oi,
                             opt.base_rng);
    opt.holder.attr("return_code") = ret;
    // unconstrain_pars(rstan_relist(par, skeleton))
    Rcpp::List theta = relist(opt.holder["par"], skeleton);
    rstan::io::rlist_ref_var_context theta_context(theta);
    std::vector<int> params_i;
    std::vector<double> params_r;
    opt.model.transform_inits(theta_context, params_i, params_r, &rstan::io::rcout);
    // constrain_pars
    if (params_r.size() != opt.model.num_params_r())
        throw std::domain_error("Number of unconstrained parameters does not match that of the model");
    std::vector<int> params_i_constrain(opt.model.num_params_i());
    std::vector<double> constrained;
    opt.model.write_array(opt.base_rng, params_r, params_i_constrain, constrained);
    // Held in an Rcpp object: List::create allocates the list before it stores
    // its arguments, and a bare SEXP from wrap() could be collected meanwhile.
    Rcpp::NumericVector constrained_pars(constrained.begin(), constrained.end());
    return Rcpp::List::create(opt.holder, constrained_pars);
}

}  // namespace leafcutter

#endif

RcppExport SEXP _leafcutter_dm_glm_multi_conc_optimizing(SEXP data, SEXP seed, SEXP args,
                                                         SEXP skeleton) {
#if STAN_MATH_MAJOR == 3 && STAN_MATH_MINOR == 0
    BEGIN_RCPP
    SEXP handle = PROTECT(R_MakeExternalPtr(nullptr, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(handle, leafcutter::free_dm_optimizer);
    Rcpp::RObject result;
    // Rcpp's interrupt and longjump exceptions are not std::exceptions: they
    // go on to END_RCPP, which resumes them.
    try {
        leafcutter::dm_optimizer* opt = new leafcutter::dm_optimizer(data, seed, args);
        R_SetExternalPtrAddr(handle, opt);
        result = leafcutter::optimize(*opt, skeleton);
    } catch (const std::exception&) {
        result = R_NilValue;
    }
    leafcutter::free_dm_optimizer(handle);
    UNPROTECT(1);
    return result;
    END_RCPP
#else
    return R_NilValue;
#endif
}

#endif
