#ifndef LEAFCUTTER_DM_GLM_MULTI_CONC_SAMPLES_H
#define LEAFCUTTER_DM_GLM_MULTI_CONC_SAMPLES_H

// The sample loop of dm_glm_multi_conc's log_prob (statements 24-42 of
// stanExports_dm_glm_multi_conc.h) for reverse-mode autodiff, as one vari per
// sample in place of the 4K + 11 the generated statements put on the tape.
//
// Values and gradients are bit-identical to the generated code: every
// floating-point operation below is one Stan Math 3.0 (StanHeaders 2.21)
// performs, on the same operands, written as the same single expression (so a
// compiler contracting a*b + c treats both alike). With another Stan Math
// version, and for scalar types other than var, add() declines and the
// generated statements run.

#include <stan/math/version.hpp>
#include <stan/math.hpp>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace leafcutter {

template <typename T>
struct dm_samples {
    template <typename... Args>
    static bool add(const Args&...) {
        return false;
    }
};

#if STAN_MATH_MAJOR == 3 && STAN_MATH_MINOR == 0

namespace dm_internal {

using stan::math::ChainableStack;
using stan::math::is_any_nan;
using stan::math::matrix_d;
using stan::math::matrix_vi;
using stan::math::var;
using stan::math::vari;
using stan::math::vector_d;

// Memory that lives until the gradient has been computed. It is never
// destructed, so nothing placed in it may own heap memory.
template <typename T>
inline T* arena_array(size_t n) {
    return ChainableStack::instance_->memalloc_.alloc_array<T>(n);
}

inline std::uint64_t bits(double x) {
    std::uint64_t b;
    std::memcpy(&b, &x, sizeof b);
    return b;
}

// lgamma of one argument, with its derivative computed the first time the
// reverse pass asks for it.
struct lgamma_slot {
    double x;
    double value;
    double derivative;
    bool has_derivative;

    double digamma() {
        if (!has_derivative) {
            derivative = stan::math::digamma(x);
            has_derivative = true;
        }
        return derivative;
    }
};

// The lgamma values of one log_prob call. Most arguments repeat across
// samples; keys are the argument's bits, so +0/-0 and NaN payloads stay apart.
class lgamma_memo {
public:
    explicit lgamma_memo(size_t max_keys) : shift_(64) {
        size_t size = 1;
        while (size < 2 * max_keys) {
            size *= 2;
            --shift_;
        }
        index_.assign(size, nullptr);
    }

    lgamma_slot* find(double x) {
        const std::uint64_t key = bits(x);
        const size_t mask = index_.size() - 1;
        size_t i = (key * 0x9E3779B97F4A7C15ULL) >> shift_;
        for (; index_[i]; i = (i + 1) & mask)
            if (bits(index_[i]->x) == key)
                return index_[i];
        // stan::math::lgamma throws for NaN; a throw must leave no slot.
        const double value = stan::math::lgamma(x);
        lgamma_slot* slot = arena_array<lgamma_slot>(1);
        slot->x = x;
        slot->value = value;
        slot->has_derivative = false;
        index_[i] = slot;
        return slot;
    }

private:
    std::vector<lgamma_slot*> index_;
    int shift_;
};

// The varis of beta (K x P, column-major) and conc, shared by all samples.
struct call_varis {
    int K;
    int P;
    vari** beta;
    vari** conc;
};

// Statement 31 for one x, and what statements 33-40 derive from x alone.
struct x_state {
    double* x;
    double* s;  // softmax(beta * x)
    lgamma_slot** a;  // a[k] = conc[k] * s[k]
    double suma;
    lgamma_slot* suma_slot;
    double sum_lgamma_a;
};

// Scratch for the reverse pass. Eigen splits vectorized assignments by the
// destination's alignment, so the destinations are plain Eigen vectors as in
// Stan, not maps over arena memory.
struct chain_scratch {
    vector_d s_adj;
    vector_d adj_times_jac;
    std::vector<double> product;
};

inline chain_scratch& scratch() {
    static thread_local chain_scratch s;
    return s;
}

// One sample's statements 31-42: target += lgamma(suma) + sum(lgamma(a + y))
// - lgamma(suma + sum(y)) - sum(lgamma(a)), with a = conc .* softmax(beta * x).
// chain() replays the chain() of every vari the generated code makes for the
// sample, in reverse order of creation. Within a sample those varis are
// consecutive on the tape, so this one vari runs at the same moment they did.
class sample_vari : public vari {
public:
    sample_vari(double value, const call_varis* call, x_state* xs, lgamma_slot** apy, double* y,
                lgamma_slot* t_slot, double sumy, double S1, double X1, double X2)
        : vari(value), call_(call), xs_(xs), apy_(apy), y_(y), t_slot_(t_slot), sumy_(sumy), S1_(S1),
          X1_(X1), X2_(X2) {}

    void chain() {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const int K = call_->K;
        const int P = call_->P;
        lgamma_slot* suma_slot = xs_->suma_slot;
        const double suma = xs_->suma;
        const double S2 = xs_->sum_lgamma_a;

        // Statement 42 in reverse: X2 - S2, X1 - lgamma(t), lgamma(suma) + S1,
        // lgamma(suma), sum(lgamma(a + y)), lgamma(t), suma + sumy, sum(lgamma(a)).
        // Stan's add and subtract varis assign NaN when an operand value is NaN.
        double adj_X2, adj_S2;
        if (is_any_nan(X2_, S2)) {
            adj_X2 = nan;
            adj_S2 = nan;
        } else {
            adj_X2 = 0.0 + adj_;
            adj_S2 = 0.0 - adj_;
        }
        double adj_X1, adj_lgamma_t;
        if (is_any_nan(X1_, t_slot_->value)) {
            adj_X1 = nan;
            adj_lgamma_t = nan;
        } else {
            adj_X1 = 0.0 + adj_X2;
            adj_lgamma_t = 0.0 - adj_X2;
        }
        double adj_lgamma_suma, adj_S1;
        if (is_any_nan(suma_slot->value, S1_)) {
            adj_lgamma_suma = nan;
            adj_S1 = nan;
        } else {
            adj_lgamma_suma = 0.0 + adj_X1;
            adj_S1 = 0.0 + adj_X1;
        }
        double adj_suma = 0.0 + adj_lgamma_suma * suma_slot->digamma();
        const double adj_lgamma_apy = 0.0 + adj_S1;
        if (sumy_ == 0.0) {
            // suma + 0.0 is suma itself
            adj_suma += adj_lgamma_t * suma_slot->digamma();
        } else {
            const double adj_t = 0.0 + adj_lgamma_t * t_slot_->digamma();
            if (is_any_nan(suma, sumy_))
                adj_suma = nan;
            else
                adj_suma += adj_t;
        }
        const double adj_lgamma_a = 0.0 + adj_S2;

        // Statements 40, 39, 37, 36 and 33. Each k has its own targets, so
        // taking k by k the contributions the generated code makes pass by
        // pass gives the same sums.
        chain_scratch& w = scratch();
        w.s_adj.resize(K);
        const double* s = xs_->s;
        for (int k = K - 1; k >= 0; --k) {
            lgamma_slot* a = xs_->a[k];
            double adj_a = 0.0 + adj_lgamma_a * a->digamma();
            if (y_[k] == 0.0) {
                // a + 0.0 is a itself
                adj_a += adj_lgamma_apy * a->digamma();
            } else {
                const double adj_apy = 0.0 + adj_lgamma_apy * apy_[k]->digamma();
                if (is_any_nan(a->x, y_[k]))
                    adj_a = nan;
                else
                    adj_a += adj_apy;
            }
            adj_a += adj_suma;
            vari* c = call_->conc[k];
            if (is_any_nan(c->val_, s[k])) {
                c->adj_ = nan;
                w.s_adj(k) = nan;
            } else {
                c->adj_ += s[k] * adj_a;
                w.s_adj(k) = 0.0 + c->val_ * adj_a;
            }
        }

        // Statement 31: softmax, as Stan's softmax_op, then the product.
        w.adj_times_jac.resize(K);
        Eigen::Map<vector_d> s_vec(xs_->s, K);
        w.adj_times_jac = -s_vec * w.s_adj.dot(s_vec) + s_vec.cwiseProduct(w.s_adj);
        // Products first, sums second, so the two cannot be fused into one
        // operation where Stan's matrix product rounds them apart.
        w.product.resize(K * P);
        for (int p = 0; p < P; ++p)
            for (int k = 0; k < K; ++k)
                w.product[k + K * p] = (0.0 + w.adj_times_jac(k)) * xs_->x[p];
        for (int i = 0; i < K * P; ++i)
            call_->beta[i]->adj_ += w.product[i];
    }

private:
    const call_varis* call_;
    x_state* xs_;
    lgamma_slot** apy_;  // lgamma of a[k] + y[k], the a[k] slot where y[k] == 0
    double* y_;
    lgamma_slot* t_slot_;  // lgamma of suma + sum(y)
    double sumy_;
    double S1_;
    double X1_;
    double X2_;
};

// Statement 31 for a new x, in Stan's order of checks and operations.
inline x_state* new_x_state(const Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>& beta,
                            double* beta_val, const vector_d& x, matrix_d& AB, vector_d& v,
                            vector_d& theta, vector_d& yv) {
    const int K = beta.rows();
    const int P = beta.cols();
    stan::math::check_multiplicable("multiply", "A", beta, "B", x);
    stan::math::check_not_nan("multiply", "A", beta);
    stan::math::check_not_nan("multiply", "B", x);
    x_state* xs = arena_array<x_state>(1);
    xs->x = arena_array<double>(P);
    Eigen::Map<matrix_d>(xs->x, P, 1) = x;
    // The expression of Stan's multiply_mat_vari, so Eigen picks the same
    // product kernel for the size.
    AB.noalias() = Eigen::Map<matrix_d>(beta_val, K, P) * Eigen::Map<matrix_d>(xs->x, P, 1);
    stan::math::check_nonzero_size("softmax", "alpha", AB);
    v.resize(K);
    for (int k = 0; k < K; ++k)
        v(k) = AB(k);
    // stan::math::softmax(v), into reused vectors
    stan::math::check_nonzero_size("softmax", "v", v);
    theta.resize(K);
    theta = (v.array() - v.maxCoeff()).exp();
    yv = theta.array() / theta.sum();
    xs->s = arena_array<double>(K);
    for (int k = 0; k < K; ++k)
        xs->s[k] = yv(k);
    return xs;
}

}  // namespace dm_internal

template <>
struct dm_samples<stan::math::var> {
    // Adds each sample's term to lp_accum, as the generated loop does, and
    // sets statement before each step that can throw, as it does.
    static bool add(const Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>& beta,
                    const std::vector<Eigen::VectorXd>& x, const std::vector<Eigen::VectorXd>& y,
                    const std::vector<stan::math::var>& conc,
                    stan::math::accumulator<stan::math::var>& lp_accum, int& statement) {
        using namespace dm_internal;
        const size_t N = x.size();
        if (N == 0)
            return true;
        const int K = beta.rows();
        const int P = beta.cols();

        call_varis* call = arena_array<call_varis>(1);
        call->K = K;
        call->P = P;
        call->beta = arena_array<vari*>(K * P);
        Eigen::Map<matrix_vi>(call->beta, K, P) = beta.vi();
        call->conc = arena_array<vari*>(K);
        for (int k = 0; k < K; ++k)
            call->conc[k] = conc[k].vi_;
        double* beta_val = arena_array<double>(K * P);
        Eigen::Map<matrix_d>(beta_val, K, P) = beta.val();

        lgamma_memo memo(N * (2 * K + 2));
        matrix_d AB;
        vector_d v, theta, yv;
        std::vector<double> a(K);
        // Samples share x (an intercept, or intercept and group); the first
        // few distinct x keep their statement 31-40 results for the others.
        const int max_kept = 8;
        x_state* kept[max_kept];
        int n_kept = 0;

        for (size_t n = 0; n < N; ++n) {
            statement = 31;
            x_state* xs = nullptr;
            for (int i = 0; i < n_kept && !xs; ++i)
                if (std::memcmp(kept[i]->x, x[n].data(), P * sizeof(double)) == 0)
                    xs = kept[i];
            const bool repeated = xs != nullptr;
            if (!repeated)
                xs = new_x_state(beta, beta_val, x[n], AB, v, theta, yv);

            const double sumy = stan::math::sum(y[n]);
            double* yn = arena_array<double>(K);
            for (int k = 0; k < K; ++k)
                yn[k] = y[n](k);
            lgamma_slot** apy = arena_array<lgamma_slot*>(K);
            if (repeated) {
                for (int k = 0; k < K; ++k) {
                    if (yn[k] == 0.0) {
                        apy[k] = xs->a[k];
                    } else {
                        statement = 39;
                        apy[k] = memo.find(xs->a[k]->x + yn[k]);
                    }
                }
            } else {
                for (int k = 0; k < K; ++k)
                    a[k] = call->conc[k]->val_ * xs->s[k];
                double suma = a[0];
                for (int k = 1; k < K; ++k)
                    suma = suma + a[k];
                xs->suma = suma;
                xs->a = arena_array<lgamma_slot*>(K);
                for (int k = 0; k < K; ++k) {
                    statement = 39;
                    if (yn[k] == 0.0) {
                        xs->a[k] = apy[k] = memo.find(a[k]);
                    } else {
                        apy[k] = memo.find(a[k] + yn[k]);
                        statement = 40;
                        xs->a[k] = memo.find(a[k]);
                    }
                }
                double S2 = xs->a[0]->value;
                for (int k = 1; k < K; ++k)
                    S2 = S2 + xs->a[k]->value;
                xs->sum_lgamma_a = S2;
            }
            double S1 = apy[0]->value;
            for (int k = 1; k < K; ++k)
                S1 = S1 + apy[k]->value;

            statement = 42;
            const double t = sumy == 0.0 ? xs->suma : xs->suma + sumy;
            lgamma_slot* t_slot = memo.find(t);
            if (!repeated) {
                xs->suma_slot = memo.find(xs->suma);
                if (n_kept < max_kept)
                    kept[n_kept++] = xs;
            }
            const double X1 = xs->suma_slot->value + S1;
            const double X2 = X1 - t_slot->value;
            const double term = X2 - xs->sum_lgamma_a;
            lp_accum.add(var(new sample_vari(term, call, xs, apy, yn, t_slot, sumy, S1, X1, X2)));
        }
        return true;
    }
};

#endif

}  // namespace leafcutter

#endif
