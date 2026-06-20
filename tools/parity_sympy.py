#!/usr/bin/env python3
"""SymPy ground-truth generator for the assumption parity sweep.

Builds the IDENTICAL battery of expressions as tests/core/parity_probe.cpp and
prints one line per (label, predicate): "<label>\t<predicate>\t<T|F|?>".
Diff against the C++ probe output to find divergences.
"""
import sympy as sp
from sympy import Symbol, Integer, Rational, oo, zoo, I, pi, E, Abs

PREDICATES = [
    "real", "rational", "integer", "positive", "negative", "zero", "nonzero",
    "nonnegative", "nonpositive", "finite", "even", "odd", "complex", "imaginary",
    "prime", "composite", "irrational", "algebraic", "transcendental",
    "extended_real", "infinite", "extended_positive", "extended_negative",
    "extended_nonnegative", "extended_nonpositive", "hermitian", "antihermitian",
    "commutative",
]


def battery():
    b = []
    def add(label, e):
        b.append((label, e))

    add("sym_generic", Symbol("g"))
    add("sym_real", Symbol("s", real=True))
    add("sym_positive", Symbol("s", positive=True))
    add("sym_negative", Symbol("s", negative=True))
    add("sym_nonnegative", Symbol("s", nonnegative=True))
    add("sym_nonpositive", Symbol("s", nonpositive=True))
    add("sym_nonzero", Symbol("s", nonzero=True))
    add("sym_integer", Symbol("s", integer=True))
    add("sym_rational", Symbol("s", rational=True))
    add("sym_even", Symbol("s", even=True))
    add("sym_odd", Symbol("s", odd=True))
    add("sym_prime", Symbol("s", prime=True))
    add("sym_composite", Symbol("s", composite=True))
    add("sym_irrational", Symbol("s", irrational=True))
    add("sym_imaginary", Symbol("s", imaginary=True))
    add("sym_complex", Symbol("s", complex=True))
    add("sym_finite", Symbol("s", finite=True))
    add("sym_extended_real", Symbol("s", extended_real=True))
    add("sym_hermitian", Symbol("s", hermitian=True))
    add("sym_antihermitian", Symbol("s", antihermitian=True))
    add("sym_ext_positive", Symbol("s", extended_positive=True))
    add("sym_ext_negative", Symbol("s", extended_negative=True))

    add("sym_real_nonzero", Symbol("s", real=True, nonzero=True))
    add("sym_integer_positive", Symbol("s", integer=True, positive=True))
    add("sym_integer_nonnegative", Symbol("s", integer=True, nonnegative=True))
    add("sym_ext_positive_finite", Symbol("s", extended_positive=True, finite=True))

    add("int_0", Integer(0))
    add("int_1", Integer(1))
    add("int_2", Integer(2))
    add("int_3", Integer(3))
    add("int_4", Integer(4))
    add("int_neg1", Integer(-1))
    add("int_neg2", Integer(-2))
    add("rat_half", Rational(1, 2))
    add("rat_neg3_2", Rational(-3, 2))
    add("oo", oo)
    add("neg_oo", -oo)
    add("zoo", zoo)
    add("I", I)
    add("pi", pi)
    add("E", E)

    g = Symbol("g")
    p = Symbol("p", positive=True)
    q = Symbol("q", positive=True)
    n = Symbol("n", integer=True)
    r = Symbol("r", real=True)
    add("abs_generic", Abs(g))
    add("neg_positive", -p)
    add("p_squared", p**2)
    add("g_squared", g**2)
    add("g_squared_plus1", g**2 + 1)
    add("p_plus_q", p + q)
    add("p_times_q", p * q)
    add("two_n_plus_1", 2 * n + 1)
    add("i_times_real", I * r)
    return b


def tri(v):
    if v is True:
        return "T"
    if v is False:
        return "F"
    return "?"


def main():
    for label, e in battery():
        for pred in PREDICATES:
            v = getattr(e, "is_" + pred)
            print(f"{label}\t{pred}\t{tri(v)}")


if __name__ == "__main__":
    main()
