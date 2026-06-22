#include <sympp/functions/bessel.hpp>

#include <optional>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>  // sqrt
#include <sympp/functions/trigonometric.hpp>

namespace sympp {

namespace {

// √(2/(π·z)) as √2·π^{−1/2}·z^{−1/2} — split form (matches SymPy; keeping z·π
// under one root would not be provably equal for a generic, possibly-negative z).
[[nodiscard]] Expr sqrt_2_over_pi_z(const Expr& z) {
    return mul(sqrt(integer(2)),
               mul(pow(S::Pi(), rational(-1, 2)), pow(z, rational(-1, 2))));
}

// True when nu is exactly ±1/2 (with `sign` receiving +1 / −1).
[[nodiscard]] bool is_half(const Expr& nu, int& sign) {
    if (nu == rational(1, 2)) { sign = 1; return true; }
    if (nu == rational(-1, 2)) { sign = -1; return true; }
    return false;
}

}  // namespace

// ----- BesselJ ---------------------------------------------------------------

BesselJ::BesselJ(Expr nu, Expr z) : Function(std::vector<Expr>{std::move(nu), std::move(z)}) {
    compute_hash(FunctionId::BesselJ);
}
Expr BesselJ::rebuild(std::vector<Expr> new_args) const {
    return besselj(new_args[0], new_args[1]);
}
std::optional<bool> BesselJ::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true && is_real(args_[1]) == true) {
        return true;
    }
    return std::nullopt;
}
Expr BesselJ::diff_arg(std::size_t i) const {
    if (i == 1) {  // ∂/∂z J_ν = (J_{ν−1} − J_{ν+1})/2
        return mul(rational(1, 2),
                   add(besselj(args_[0] - S::One(), args_[1]),
                       mul(S::NegativeOne(), besselj(args_[0] + S::One(), args_[1]))));
    }
    return Function::diff_arg(i);  // ∂/∂ν non-elementary
}

// ----- BesselY ---------------------------------------------------------------

BesselY::BesselY(Expr nu, Expr z) : Function(std::vector<Expr>{std::move(nu), std::move(z)}) {
    compute_hash(FunctionId::BesselY);
}
Expr BesselY::rebuild(std::vector<Expr> new_args) const {
    return bessely(new_args[0], new_args[1]);
}
Expr BesselY::diff_arg(std::size_t i) const {
    if (i == 1) {  // ∂/∂z Y_ν = (Y_{ν−1} − Y_{ν+1})/2
        return mul(rational(1, 2),
                   add(bessely(args_[0] - S::One(), args_[1]),
                       mul(S::NegativeOne(), bessely(args_[0] + S::One(), args_[1]))));
    }
    return Function::diff_arg(i);
}

// ----- BesselI ---------------------------------------------------------------

BesselI::BesselI(Expr nu, Expr z) : Function(std::vector<Expr>{std::move(nu), std::move(z)}) {
    compute_hash(FunctionId::BesselI);
}
Expr BesselI::rebuild(std::vector<Expr> new_args) const {
    return besseli(new_args[0], new_args[1]);
}
std::optional<bool> BesselI::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true && is_real(args_[1]) == true) {
        return true;
    }
    return std::nullopt;
}
Expr BesselI::diff_arg(std::size_t i) const {
    if (i == 1) {  // ∂/∂z I_ν = (I_{ν−1} + I_{ν+1})/2
        return mul(rational(1, 2),
                   add(besseli(args_[0] - S::One(), args_[1]),
                       besseli(args_[0] + S::One(), args_[1])));
    }
    return Function::diff_arg(i);
}

// ----- BesselK ---------------------------------------------------------------

BesselK::BesselK(Expr nu, Expr z) : Function(std::vector<Expr>{std::move(nu), std::move(z)}) {
    compute_hash(FunctionId::BesselK);
}
Expr BesselK::rebuild(std::vector<Expr> new_args) const {
    return besselk(new_args[0], new_args[1]);
}
Expr BesselK::diff_arg(std::size_t i) const {
    if (i == 1) {  // ∂/∂z K_ν = −(K_{ν−1} + K_{ν+1})/2
        return mul(rational(-1, 2),
                   add(besselk(args_[0] - S::One(), args_[1]),
                       besselk(args_[0] + S::One(), args_[1])));
    }
    return Function::diff_arg(i);
}

// ----- factories -------------------------------------------------------------

Expr besselj(const Expr& nu, const Expr& z) {
    if (z == S::Zero()) {
        if (nu == S::Zero()) return S::One();
        if (is_positive(nu) == true) return S::Zero();
    }
    int sgn = 0;
    if (is_half(nu, sgn)) {
        Expr pref = sqrt_2_over_pi_z(z);  // J_{1/2}=√(2/πz)sin z, J_{−1/2}=…cos z
        return mul(pref, sgn > 0 ? sin(z) : cos(z));
    }
    return make<BesselJ>(nu, z);
}

Expr bessely(const Expr& nu, const Expr& z) {
    int sgn = 0;
    if (is_half(nu, sgn)) {
        // Y_{1/2}(z) = −J_{−1/2}(z) = −√(2/πz)cos z; Y_{−1/2}(z) = J_{1/2}(z).
        Expr pref = sqrt_2_over_pi_z(z);
        return sgn > 0 ? mul(S::NegativeOne(), mul(pref, cos(z))) : mul(pref, sin(z));
    }
    return make<BesselY>(nu, z);
}

Expr besseli(const Expr& nu, const Expr& z) {
    if (z == S::Zero()) {
        if (nu == S::Zero()) return S::One();
        if (is_positive(nu) == true) return S::Zero();
    }
    int sgn = 0;
    if (is_half(nu, sgn)) {
        Expr pref = sqrt_2_over_pi_z(z);  // I_{1/2}=√(2/πz)sinh z, I_{−1/2}=…cosh z
        return mul(pref, sgn > 0 ? sinh(z) : cosh(z));
    }
    return make<BesselI>(nu, z);
}

Expr besselk(const Expr& nu, const Expr& z) {
    int sgn = 0;
    if (is_half(nu, sgn)) {  // K even in ν: K_{±1/2}(z) = √2·√π·e^{−z}/(2·√z)
        Expr pref = mul(mul(sqrt(integer(2)), sqrt(S::Pi())),
                        mul(rational(1, 2), pow(z, rational(-1, 2))));
        return mul(pref, exp(mul(S::NegativeOne(), z)));
    }
    return make<BesselK>(nu, z);
}

}  // namespace sympp
