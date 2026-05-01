#pragma once

// Number — abstract base for the SymPP exact/symbolic number tower.
// Hierarchy:
//     Number
//     ├── Integer       (GMP mpz_class)
//     ├── Rational      (GMP mpq_class)        [Phase 1c]
//     ├── Float         (MPFR mpfr_t)          [Phase 1d]
//     └── ComplexNumber (a + b*I)              [Phase 1d]
//
// Numbers are the leaves of the expression tree that carry value. They are
// always atomic (args() returns empty) and answer the value-introspection
// methods (is_zero, is_one, is_negative, sign).

#include <cstddef>
#include <span>
#include <string>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Number : public Basic {
public:
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return {}; }

    [[nodiscard]] virtual bool is_zero() const noexcept = 0;
    [[nodiscard]] virtual bool is_one() const noexcept = 0;
    [[nodiscard]] virtual bool is_negative() const noexcept = 0;
    [[nodiscard]] virtual bool is_positive() const noexcept = 0;
    [[nodiscard]] virtual bool is_integer() const noexcept = 0;
    [[nodiscard]] virtual bool is_rational() const noexcept = 0;
    [[nodiscard]] virtual bool is_real() const noexcept { return true; }

    // -1 if negative, 0 if zero, +1 if positive. NaN returns 0.
    [[nodiscard]] virtual int sign() const noexcept = 0;
};

}  // namespace sympp
