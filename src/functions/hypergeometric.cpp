#include <sympp/functions/hypergeometric.hpp>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/functions/exponential.hpp>

namespace sympp {

// ----- Hyper -----------------------------------------------------------------

Hyper::Hyper(std::vector<Expr> args) : Function(std::move(args)) {
    compute_hash(FunctionId::Hyper);
}

Expr Hyper::rebuild(std::vector<Expr> new_args) const {
    // The args layout is invariant (p, q, ap..., bq..., z); going through
    // the factory re-applies auto-evaluation rules.
    if (new_args.size() < 3) {
        throw std::invalid_argument("Hyper::rebuild: malformed args");
    }
    auto get_int = [](const Expr& e) -> long {
        if (e->type_id() != TypeId::Integer) {
            throw std::invalid_argument(
                "Hyper::rebuild: parameter count must be a literal Integer");
        }
        try {
            return std::stol(e->str());
        } catch (...) {
            throw std::invalid_argument("Hyper::rebuild: count overflow");
        }
    };
    long p = get_int(new_args[0]);
    long q = get_int(new_args[1]);
    auto exp = static_cast<std::size_t>(2 + p + q + 1);
    if (new_args.size() != exp) {
        throw std::invalid_argument("Hyper::rebuild: arg count mismatch");
    }
    std::vector<Expr> ap_v(new_args.begin() + 2, new_args.begin() + 2 + p);
    std::vector<Expr> bq_v(new_args.begin() + 2 + p,
                              new_args.begin() + 2 + p + q);
    Expr z_v = new_args.back();
    return hyper(ap_v, bq_v, z_v);
}

std::string Hyper::str() const {
    // hyper((a₁, a₂, …), (b₁, b₂, …), z)
    auto p = this->p();
    auto q = this->q();
    auto ap_v = ap();
    auto bq_v = bq();
    Expr z_v = z();
    std::string out = "hyper((";
    for (std::size_t i = 0; i < ap_v.size(); ++i) {
        if (i > 0) out.append(", ");
        out.append(ap_v[i]->str());
    }
    if (p == 1) out.push_back(',');  // singleton-tuple syntax
    out.append("), (");
    for (std::size_t i = 0; i < bq_v.size(); ++i) {
        if (i > 0) out.append(", ");
        out.append(bq_v[i]->str());
    }
    if (q == 1) out.push_back(',');
    out.append("), ");
    out.append(z_v->str());
    out.push_back(')');
    return out;
}

namespace {

[[nodiscard]] long parse_count(const Expr& e) {
    if (e->type_id() != TypeId::Integer) return -1;
    try {
        return std::stol(e->str());
    } catch (...) {
        return -1;
    }
}

}  // namespace

std::size_t Hyper::p() const {
    long v = parse_count(args()[0]);
    return v < 0 ? 0 : static_cast<std::size_t>(v);
}

std::size_t Hyper::q() const {
    long v = parse_count(args()[1]);
    return v < 0 ? 0 : static_cast<std::size_t>(v);
}

std::vector<Expr> Hyper::ap() const {
    using Diff = std::span<const Expr>::difference_type;
    auto pcount = p();
    auto a = args();
    return std::vector<Expr>(a.begin() + 2,
                                a.begin() + 2 + static_cast<Diff>(pcount));
}

std::vector<Expr> Hyper::bq() const {
    using Diff = std::span<const Expr>::difference_type;
    auto pcount = p();
    auto qcount = q();
    auto a = args();
    return std::vector<Expr>(
        a.begin() + 2 + static_cast<Diff>(pcount),
        a.begin() + 2 + static_cast<Diff>(pcount + qcount));
}

Expr Hyper::z() const { return args().back(); }

// d/dz pFq(a; b; z) = (a₁ ⋯ aₚ) / (b₁ ⋯ b_q) · pFq(a₁+1, …; b₁+1, …; z).
Expr Hyper::diff_arg(std::size_t i) const {
    auto a_v = ap();
    auto b_v = bq();
    auto z_v = z();
    auto last_index = 2 + a_v.size() + b_v.size();
    if (i != last_index) {
        // Derivative w.r.t. parameters (rare). Leave unevaluated.
        return S::Zero();
    }
    // Build the parameter-shifted hypergeometric.
    Expr ratio = S::One();
    std::vector<Expr> a_shift, b_shift;
    a_shift.reserve(a_v.size());
    b_shift.reserve(b_v.size());
    for (const auto& a : a_v) {
        ratio = ratio * a;
        a_shift.push_back(a + integer(1));
    }
    for (const auto& b : b_v) {
        ratio = ratio / b;
        b_shift.push_back(b + integer(1));
    }
    return ratio * hyper(a_shift, b_shift, z_v);
}

// ----- MeijerG ---------------------------------------------------------------

MeijerG::MeijerG(std::vector<Expr> args) : Function(std::move(args)) {
    compute_hash(FunctionId::MeijerG);
}

Expr MeijerG::rebuild(std::vector<Expr> new_args) const {
    // Re-build through factory to re-validate the layout.
    if (new_args.size() < 5) {
        throw std::invalid_argument("MeijerG::rebuild: malformed args");
    }
    auto get_count = [](const Expr& e) -> long {
        if (e->type_id() != TypeId::Integer) return -1;
        try { return std::stol(e->str()); } catch (...) { return -1; }
    };
    long n = get_count(new_args[0]);
    long pr = get_count(new_args[1]);
    long m = get_count(new_args[2]);
    long qr = get_count(new_args[3]);
    if (n < 0 || pr < 0 || m < 0 || qr < 0) {
        throw std::invalid_argument(
            "MeijerG::rebuild: parameter counts must be non-negative integers");
    }
    auto total = static_cast<std::size_t>(4 + n + pr + m + qr + 1);
    if (new_args.size() != total) {
        throw std::invalid_argument("MeijerG::rebuild: arg count mismatch");
    }
    auto begin = new_args.begin() + 4;
    std::vector<Expr> an(begin, begin + n);
    begin += n;
    std::vector<Expr> ap_rest(begin, begin + pr);
    begin += pr;
    std::vector<Expr> bm(begin, begin + m);
    begin += m;
    std::vector<Expr> bq_rest(begin, begin + qr);
    Expr z_v = new_args.back();
    return meijerg(an, ap_rest, bm, bq_rest, z_v);
}

std::string MeijerG::str() const {
    auto a = args();
    long n = parse_count(a[0]);
    long pr = parse_count(a[1]);
    long m = parse_count(a[2]);
    long qr = parse_count(a[3]);
    auto fmt_seq = [&](std::size_t start, long count) {
        std::string s = "(";
        for (long i = 0; i < count; ++i) {
            if (i > 0) s.append(", ");
            s.append(a[start + static_cast<std::size_t>(i)]->str());
        }
        if (count == 1) s.push_back(',');
        s.push_back(')');
        return s;
    };
    auto sz = [](long v) { return static_cast<std::size_t>(v); };
    std::string out = "meijerg(";
    out.append(fmt_seq(4, n));
    out.append(", ");
    out.append(fmt_seq(4 + sz(n), pr));
    out.append(", ");
    out.append(fmt_seq(4 + sz(n) + sz(pr), m));
    out.append(", ");
    out.append(fmt_seq(4 + sz(n) + sz(pr) + sz(m), qr));
    out.append(", ");
    out.append(a.back()->str());
    out.push_back(')');
    return out;
}

// ----- Factories with auto-eval ---------------------------------------------

namespace {

// Try to remove a single matching pair (aᵢ, bⱼ) where aᵢ == bⱼ. Returns
// true and shrinks the vectors when a match is found.
[[nodiscard]] bool cancel_one(std::vector<Expr>& ap, std::vector<Expr>& bq) {
    using Diff = std::vector<Expr>::difference_type;
    for (std::size_t i = 0; i < ap.size(); ++i) {
        for (std::size_t j = 0; j < bq.size(); ++j) {
            if (ap[i] == bq[j]) {
                ap.erase(ap.begin() + static_cast<Diff>(i));
                bq.erase(bq.begin() + static_cast<Diff>(j));
                return true;
            }
        }
    }
    return false;
}

[[nodiscard]] Expr build_hyper_node(const std::vector<Expr>& ap,
                                       const std::vector<Expr>& bq,
                                       const Expr& z) {
    std::vector<Expr> args;
    args.reserve(2 + ap.size() + bq.size() + 1);
    args.push_back(integer(static_cast<long>(ap.size())));
    args.push_back(integer(static_cast<long>(bq.size())));
    for (const auto& a : ap) args.push_back(a);
    for (const auto& b : bq) args.push_back(b);
    args.push_back(z);
    return make<Hyper>(std::move(args));
}

}  // namespace

Expr hyper(const std::vector<Expr>& ap_in, const std::vector<Expr>& bq_in,
              const Expr& z) {
    // pFq(a; b; 0) = 1.
    if (z == S::Zero()) return S::One();

    // Numerator/denominator parameter cancellation. Loop because each
    // cancellation can expose another match.
    std::vector<Expr> ap = ap_in;
    std::vector<Expr> bq = bq_in;
    while (cancel_one(ap, bq)) {}

    // ₀F₀(z) = exp(z).
    if (ap.empty() && bq.empty()) {
        return exp(z);
    }
    // ₁F₀(a; ; z) = (1 − z)^(−a).
    if (ap.size() == 1 && bq.empty()) {
        return pow(integer(1) - z, -ap[0]);
    }
    return build_hyper_node(ap, bq, z);
}

Expr hyper(const Expr& a, const Expr& b, const Expr& c, const Expr& z) {
    return hyper(std::vector<Expr>{a, b}, std::vector<Expr>{c}, z);
}

Expr meijerg(const std::vector<Expr>& an, const std::vector<Expr>& ap_rest,
                const std::vector<Expr>& bm,
                const std::vector<Expr>& bq_rest, const Expr& z) {
    std::vector<Expr> args;
    args.reserve(4 + an.size() + ap_rest.size() + bm.size() + bq_rest.size() + 1);
    args.push_back(integer(static_cast<long>(an.size())));
    args.push_back(integer(static_cast<long>(ap_rest.size())));
    args.push_back(integer(static_cast<long>(bm.size())));
    args.push_back(integer(static_cast<long>(bq_rest.size())));
    for (const auto& a : an) args.push_back(a);
    for (const auto& a : ap_rest) args.push_back(a);
    for (const auto& b : bm) args.push_back(b);
    for (const auto& b : bq_rest) args.push_back(b);
    args.push_back(z);
    return make<MeijerG>(std::move(args));
}

}  // namespace sympp
