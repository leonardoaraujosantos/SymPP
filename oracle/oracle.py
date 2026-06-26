#!/usr/bin/env python3
"""SymPP oracle harness — long-lived SymPy subprocess.

Protocol: line-delimited JSON over stdin/stdout.
Each request: {"id": int, "op": str, ...op-specific kwargs}
Each response: {"id": int, "ok": bool, "result": ..., "error"?: str, "detail"?: str}

Always flushes stdout after each response so the C++ side never blocks.
"""
import json
import sys
import traceback

import sympy


def _to_str(expr):
    """Stringify a SymPy result in a format we can compare downstream."""
    return str(expr)


def _sympify(s):
    """Parse using SymPy. Local symbols come from the caller's expression text;
    we let SymPy infer them via the default sympify rules."""
    return sympy.sympify(s)


def handle(req):
    op = req.get("op")

    # --- Health / introspection ---
    if op == "ping":
        return {
            "ok": True,
            "result": "pong",
            "sympy_version": sympy.__version__,
            "python_version": sys.version.split()[0],
        }

    if op == "shutdown":
        return {"ok": True, "result": "bye", "_shutdown": True}

    # --- Parsing / structural ---
    if op == "parse":
        e = _sympify(req["expr"])
        return {"ok": True, "result": sympy.srepr(e), "str": _to_str(e)}

    if op == "srepr":
        e = _sympify(req["expr"])
        return {"ok": True, "result": sympy.srepr(e)}

    # --- Equivalence ---
    if op == "equivalent":
        a = _sympify(req["a"])
        b = _sympify(req["b"])
        diff = sympy.simplify(a - b)
        is_zero = (diff == 0) or (diff.equals(0) is True)
        return {"ok": True, "result": bool(is_zero), "diff": _to_str(diff)}

    # --- Numeric ---
    if op == "evalf":
        e = _sympify(req["expr"])
        prec = req.get("prec", 15)
        return {"ok": True, "result": _to_str(e.evalf(prec))}

    # --- Algebra / calculus (used in later phases) ---
    if op == "simplify":
        e = _sympify(req["expr"])
        return {"ok": True, "result": _to_str(sympy.simplify(e))}

    if op == "expand":
        e = _sympify(req["expr"])
        return {"ok": True, "result": _to_str(sympy.expand(e))}

    if op == "factor":
        e = _sympify(req["expr"])
        return {"ok": True, "result": _to_str(sympy.factor(e))}

    if op == "diff":
        e = _sympify(req["expr"])
        v = sympy.Symbol(req["var"])
        order = req.get("order", 1)
        return {"ok": True, "result": _to_str(sympy.diff(e, v, order))}

    if op == "integrate":
        e = _sympify(req["expr"])
        v = sympy.Symbol(req["var"])
        if "lower" in req and "upper" in req:
            return {
                "ok": True,
                "result": _to_str(
                    sympy.integrate(e, (v, _sympify(req["lower"]), _sympify(req["upper"])))
                ),
            }
        return {"ok": True, "result": _to_str(sympy.integrate(e, v))}

    if op == "limit":
        e = _sympify(req["expr"])
        v = sympy.Symbol(req["var"])
        target = _sympify(req["to"])
        direction = req.get("dir", "+-")
        return {"ok": True, "result": _to_str(sympy.limit(e, v, target, direction))}

    if op == "series":
        e = _sympify(req["expr"])
        v = sympy.Symbol(req["var"])
        x0 = _sympify(req.get("x0", "0"))
        n = req.get("n", 6)
        return {"ok": True, "result": _to_str(e.series(v, x0, n).removeO())}

    if op == "solve":
        e = _sympify(req["expr"])
        v = sympy.Symbol(req["var"])
        sols = sympy.solve(e, v)
        return {"ok": True, "result": [_to_str(s) for s in sols]}

    if op == "latex":
        e = _sympify(req["expr"])
        return {"ok": True, "result": sympy.latex(e)}

    # --- Phase 4 polynomial ops ---
    if op == "div":
        # Polynomial long division: returns (quotient, remainder) over var.
        p = _sympify(req["p"])
        q = _sympify(req["q"])
        v = sympy.Symbol(req["var"])
        quot, rem = sympy.div(p, q, v)
        return {"ok": True, "quotient": _to_str(quot), "remainder": _to_str(rem)}

    if op == "gcd":
        a = _sympify(req["a"])
        b = _sympify(req["b"])
        return {"ok": True, "result": _to_str(sympy.gcd(a, b))}

    if op == "sqf_list":
        # Returns (content, [(factor, multiplicity), ...]) over var.
        e = _sympify(req["expr"])
        v = sympy.Symbol(req["var"])
        content, factors = sympy.sqf_list(e, v)
        return {
            "ok": True,
            "content": _to_str(content),
            "factors": [[_to_str(f), int(m)] for f, m in factors],
        }

    if op == "apart":
        e = _sympify(req["expr"])
        v = sympy.Symbol(req["var"])
        return {"ok": True, "result": _to_str(sympy.apart(e, v))}

    if op == "together":
        e = _sympify(req["expr"])
        return {"ok": True, "result": _to_str(sympy.together(e))}

    if op == "cancel":
        e = _sympify(req["expr"])
        return {"ok": True, "result": _to_str(sympy.cancel(e))}

    if op == "resultant":
        p = _sympify(req["p"])
        q = _sympify(req["q"])
        v = sympy.Symbol(req["var"])
        return {"ok": True, "result": _to_str(sympy.resultant(p, q, v))}

    if op == "discriminant":
        e = _sympify(req["expr"])
        v = sympy.Symbol(req["var"])
        return {"ok": True, "result": _to_str(sympy.discriminant(e, v))}

    # --- Hypergeometric ---
    if op == "hyperexpand":
        # Run sympy.hyperexpand on the input. Returns the rewritten form
        # as a string. Use this to cross-validate SymPP's hyperexpand
        # against SymPy's reference implementation.
        e = _sympify(req["expr"])
        return {"ok": True, "result": _to_str(sympy.hyperexpand(e))}

    if op == "evalf_is_zero":
        # High-precision numeric check that |expr| < 10**(-tol). Useful when
        # symbolic simplification of nested radicals can't prove zero but the
        # expression is in fact zero (Cardano substitutions, etc.).
        e = _sympify(req["expr"])
        prec = req.get("prec", 50)
        tol = req.get("tol", 25)
        try:
            v = e.evalf(prec)
            mag = abs(v)
            return {"ok": True, "result": bool(mag < sympy.Float(f"1e-{tol}"))}
        except Exception as ex:
            return {"ok": False, "error": type(ex).__name__, "detail": str(ex)}

    # --- Discrete transforms (cross-check against sympy.discrete) ---
    if op == "discrete":
        from sympy.discrete import fft, ifft, ntt, intt, convolution
        from sympy.discrete.transforms import (mobius_transform,
                                               inverse_mobius_transform)
        seq = [sympy.sympify(x) for x in req["seq"]]
        fn = req["fn"]
        if fn == "fft": res = fft(seq)
        elif fn == "ifft": res = ifft(seq)
        elif fn == "ntt": res = ntt(seq, int(req["p"]))
        elif fn == "intt": res = intt(seq, int(req["p"]))
        elif fn == "convolution":
            res = convolution(seq, [sympy.sympify(x) for x in req["seq2"]])
        elif fn == "mobius": res = mobius_transform(seq)
        elif fn == "inv_mobius": res = inverse_mobius_transform(seq)
        else:
            return {"ok": False, "error": "BadFn", "detail": f"unknown discrete fn: {fn!r}"}
        return {"ok": True, "result": ";".join(_to_str(x) for x in res)}

    # --- Tensor algebra (cross-check against sympy.tensor.array) ---
    if op == "tensor":
        from sympy import Array, tensorproduct, tensorcontraction

        def flatten(a):
            if hasattr(a, "shape") and a.shape:
                n = 1
                for d in a.shape:
                    n *= d
                return ";".join(_to_str(x) for x in a.reshape(n))
            return _to_str(a)

        def build(spec):
            arr = Array([sympy.sympify(x) for x in spec["data"]])
            return arr.reshape(*[int(d) for d in spec["shape"]])

        fn = req["fn"]
        if fn == "contract":
            res = tensorcontraction(build(req["a"]), tuple(int(x) for x in req["axes"]))
            return {"ok": True, "result": flatten(res)}
        if fn == "tensorproduct":
            res = tensorproduct(build(req["a"]), build(req["b"]))
            return {"ok": True, "result": flatten(res)}
        return {"ok": False, "error": "BadFn", "detail": f"unknown tensor fn: {fn!r}"}

    # --- Polynomial factorization count (cross-check Zassenhaus) ---
    if op == "polyfactor":
        from sympy import factor_list
        e = _sympify(req["expr"])
        x = sympy.Symbol(req["var"])
        _, facs = factor_list(e, x)
        return {"ok": True, "result": str(sum(int(m) for _, m in facs))}

    # --- Cryptography (cross-check RSA against sympy.crypto) ---
    if op == "crypto":
        from sympy.crypto.crypto import rsa_private_key, rsa_public_key, encipher_rsa
        p, q, e = int(req["p"]), int(req["q"]), int(req["e"])
        fn = req["fn"]
        if fn == "rsa_d":
            return {"ok": True, "result": str(int(rsa_private_key(p, q, e)[1]))}
        if fn == "rsa_encrypt":
            c = encipher_rsa(int(req["m"]), rsa_public_key(p, q, e))
            return {"ok": True, "result": str(int(c))}
        return {"ok": False, "error": "BadFn", "detail": f"unknown crypto fn: {fn!r}"}

    # --- Vector calculus (grad/div/curl/laplacian via sympy.diff) ---
    if op == "vectorcalc":
        vs = [sympy.Symbol(v) for v in req["vars"]]
        fn = req["fn"]
        if fn == "gradient":
            f = sympy.sympify(req["f"])
            return {"ok": True, "result": ";".join(_to_str(sympy.diff(f, v)) for v in vs)}
        if fn == "divergence":
            F = [sympy.sympify(c) for c in req["field"]]
            return {"ok": True,
                    "result": _to_str(sum(sympy.diff(F[i], vs[i]) for i in range(len(vs))))}
        if fn == "curl":
            F = [sympy.sympify(c) for c in req["field"]]
            comps = [sympy.diff(F[2], vs[1]) - sympy.diff(F[1], vs[2]),
                     sympy.diff(F[0], vs[2]) - sympy.diff(F[2], vs[0]),
                     sympy.diff(F[1], vs[0]) - sympy.diff(F[0], vs[1])]
            return {"ok": True, "result": ";".join(_to_str(c) for c in comps)}
        if fn == "laplacian":
            f = sympy.sympify(req["f"])
            return {"ok": True,
                    "result": _to_str(sum(sympy.diff(f, v, 2) for v in vs))}
        return {"ok": False, "error": "BadFn", "detail": f"unknown vectorcalc fn: {fn!r}"}

    # --- Differential geometry (Ricci scalar from a metric) ---
    if op == "diffgeom":
        coords = [sympy.Symbol(c) for c in req["coords"]]
        n = len(coords)
        g = sympy.Matrix([[sympy.sympify(req["metric"][i][j]) for j in range(n)]
                          for i in range(n)])
        gi = g.inv()
        G = [[[sympy.simplify(sum(
                  gi[k, l] * (sympy.diff(g[l, j], coords[i]) + sympy.diff(g[l, i], coords[j])
                              - sympy.diff(g[i, j], coords[l])) for l in range(n)) / 2)
               for j in range(n)] for i in range(n)] for k in range(n)]
        R = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                r = 0
                for k in range(n):
                    r += sympy.diff(G[k][i][j], coords[k]) - sympy.diff(G[k][i][k], coords[j])
                    for l in range(n):
                        r += G[k][k][l] * G[l][i][j] - G[k][j][l] * G[l][i][k]
                R[i][j] = sympy.simplify(r)
        if req["fn"] == "ricci_scalar":
            s = sympy.simplify(sum(gi[i, j] * R[i][j] for i in range(n) for j in range(n)))
            return {"ok": True, "result": _to_str(s)}
        return {"ok": False, "error": "BadFn", "detail": f"unknown diffgeom fn: {req['fn']!r}"}

    # --- Matrix (cross-check singular values) ---
    if op == "matrix":
        rows = [[sympy.sympify(c) for c in row] for row in req["rows"]]
        M = sympy.Matrix(rows)
        if req["fn"] == "singular_values":
            return {"ok": True,
                    "result": ";".join(_to_str(s) for s in M.singular_values())}
        return {"ok": False, "error": "BadFn", "detail": f"unknown matrix fn: {req['fn']!r}"}

    # --- Statistics (cross-check moments / density / cdf) ---
    if op == "stats":
        import sympy.stats as st
        params = [sympy.sympify(p) for p in req["params"]]
        builders = {
            "normal": lambda: st.Normal("X", params[0], params[1]),
            "uniform": lambda: st.Uniform("X", params[0], params[1]),
            "exponential": lambda: st.Exponential("X", params[0]),
            "bernoulli": lambda: st.Bernoulli("X", params[0]),
            "binomial": lambda: st.Binomial("X", params[0], params[1]),
            "poisson": lambda: st.Poisson("X", params[0]),
        }
        rv = builders[req["dist"]]()
        fn = req["fn"]
        if fn == "E":
            return {"ok": True, "result": _to_str(sympy.simplify(st.E(rv)))}
        if fn == "variance":
            return {"ok": True, "result": _to_str(sympy.simplify(st.variance(rv)))}
        if fn == "density":
            xv = sympy.sympify(req["x"])
            return {"ok": True, "result": _to_str(st.density(rv)(xv))}
        if fn == "cdf":
            xv = sympy.sympify(req["x"])
            return {"ok": True, "result": _to_str(st.cdf(rv)(xv))}
        return {"ok": False, "error": "BadFn", "detail": f"unknown stats fn: {fn!r}"}

    # --- Plane geometry (cross-check distance / area / collinearity) ---
    if op == "geometry":
        from sympy import Point, Polygon
        fn = req["fn"]
        coords = [sympy.sympify(c) for c in req["coords"]]
        pts = [Point(coords[i], coords[i + 1]) for i in range(0, len(coords), 2)]
        if fn == "distance":
            return {"ok": True, "result": _to_str(pts[0].distance(pts[1]))}
        if fn == "polygon_area":
            return {"ok": True, "result": _to_str(Polygon(*pts).area)}
        if fn == "collinear":
            return {"ok": True, "result": bool(Point.is_collinear(*pts))}
        return {"ok": False, "error": "BadFn", "detail": f"unknown geometry fn: {fn!r}"}

    # --- 2D pretty printing (ASCII) ---
    if op == "pretty":
        e = _sympify(req["expr"])
        return {"ok": True, "result": sympy.pretty(e, use_unicode=False)}

    # --- Propositional logic (cross-check satisfiable / equivalence) ---
    if op == "logic":
        from sympy.logic import satisfiable as _sat
        fn = req["fn"]
        if fn == "satisfiable":
            return {"ok": True, "result": bool(_sat(_sympify(req["expr"])))}
        if fn == "equivalent":
            a = _sympify(req["a"])
            b = _sympify(req["b"])
            valid = _sat(sympy.Not(sympy.Equivalent(a, b))) is False
            return {"ok": True, "result": bool(valid)}
        return {"ok": False, "error": "BadFn", "detail": f"unknown logic fn: {fn!r}"}

    # --- Combinatorics & group theory ---
    if op == "combinatorics":
        fn = req["fn"]
        if fn == "group_order":
            from sympy.combinatorics import Permutation, PermutationGroup
            gens = [Permutation(list(g)) for g in req["gens"]]
            return {"ok": True, "result": str(PermutationGroup(gens).order())}
        if fn == "is_abelian":
            from sympy.combinatorics import Permutation, PermutationGroup
            gens = [Permutation(list(g)) for g in req["gens"]]
            return {"ok": True, "result": bool(PermutationGroup(gens).is_abelian)}
        if fn == "perm_sign":
            from sympy.combinatorics import Permutation
            return {"ok": True, "result": str(Permutation(list(req["p"])).signature())}
        if fn == "perm_order":
            from sympy.combinatorics import Permutation
            return {"ok": True, "result": str(Permutation(list(req["p"])).order())}
        if fn == "partition":
            from sympy.functions.combinatorial.numbers import partition
            return {"ok": True, "result": str(partition(int(req["n"])))}
        if fn == "orbits":
            from sympy.combinatorics import Permutation, PermutationGroup
            gens = [Permutation(list(g)) for g in req["gens"]]
            os = sorted(sorted(int(x) for x in o) for o in PermutationGroup(gens).orbits())
            return {"ok": True, "result": ";".join(",".join(str(x) for x in o) for o in os)}
        if fn == "burnside":
            from sympy.combinatorics import Permutation, PermutationGroup
            gens = [Permutation(list(g)) for g in req["gens"]]
            G = PermutationGroup(gens); k = int(req["k"])
            def ncyc(g):
                moved = sum(len(c) for c in g.cyclic_form)
                return len(g.cyclic_form) + (g.size - moved)
            v = sum(k ** ncyc(g) for g in G.elements) // G.order()
            return {"ok": True, "result": str(int(v))}
        if fn == "necklaces":
            from sympy import totient, divisors
            n = int(req["n"]); k = int(req["k"])
            v = sum(totient(d) * k ** (n // d) for d in divisors(n)) // n
            return {"ok": True, "result": str(int(v))}
        return {"ok": False, "error": "BadFn", "detail": f"unknown combinatorics fn: {fn!r}"}

    # --- Physics (cross-check Wigner symbols, hydrogen, QHO, spin) ---
    if op == "physics":
        fn = req["fn"]
        S_ = sympy.sympify
        if fn == "wigner_3j":
            from sympy.physics.wigner import wigner_3j as _w3
            v = _w3(*[S_(req[k]) for k in ("j1", "j2", "j3", "m1", "m2", "m3")])
            return {"ok": True, "result": str(sympy.sympify(v))}
        if fn == "wigner_6j":
            from sympy.physics.wigner import wigner_6j as _w6
            v = _w6(*[S_(req[k]) for k in ("j1", "j2", "j3", "j4", "j5", "j6")])
            return {"ok": True, "result": str(sympy.sympify(v))}
        if fn == "wigner_9j":
            from sympy.physics.wigner import wigner_9j as _w9
            ks = ("j1", "j2", "j3", "j4", "j5", "j6", "j7", "j8", "j9")
            v = _w9(*[S_(req[k]) for k in ks])
            return {"ok": True, "result": str(sympy.sympify(v))}
        if fn == "racah":
            from sympy.physics.wigner import racah as _ra
            v = _ra(*[S_(req[k]) for k in ("a", "b", "c", "d", "e", "f")])
            return {"ok": True, "result": str(sympy.sympify(v))}
        if fn == "gaunt":
            from sympy.physics.wigner import gaunt as _ga
            v = _ga(*[S_(req[k]) for k in ("l1", "l2", "l3", "m1", "m2", "m3")])
            return {"ok": True, "result": str(sympy.sympify(v))}
        if fn == "clebsch_gordan":
            from sympy.physics.wigner import clebsch_gordan as _cg
            v = _cg(*[S_(req[k]) for k in ("j1", "j2", "j3", "m1", "m2", "m3")])
            return {"ok": True, "result": str(sympy.sympify(v))}
        if fn == "hydrogen_E":
            from sympy.physics.hydrogen import E_nl
            return {"ok": True, "result": str(E_nl(S_(req["n"]), S_(req["Z"])))}
        if fn == "hydrogen_R":
            from sympy.physics.hydrogen import R_nl
            r = sympy.Symbol("r")
            v = R_nl(int(req["n"]), int(req["l"]), r, S_(req["Z"]))
            return {"ok": True, "result": str(sympy.simplify(v))}
        if fn == "qho_E":
            from sympy.physics.qho_1d import E_n
            from sympy.physics.quantum.constants import hbar
            v = E_n(int(req["n"]), sympy.Symbol("omega")).subs(hbar, 1)
            return {"ok": True, "result": str(v)}
        if fn == "qho_psi":
            from sympy.physics.qho_1d import psi_n
            from sympy.physics.quantum.constants import hbar
            v = psi_n(int(req["n"]), sympy.Symbol("x"), 1, 1).subs(hbar, 1)
            return {"ok": True, "result": str(sympy.simplify(v))}
        return {"ok": False, "error": "BadFn", "detail": f"unknown physics fn: {fn!r}"}

    # --- Number theory (cross-check factorint/divisors/igcdex/jacobi) ---
    if op == "ntheory":
        fn = req["fn"]
        if fn == "factorint":
            d = sympy.factorint(int(req["n"]))
            pairs = sorted((int(p), int(e)) for p, e in d.items())
            return {"ok": True, "result": ";".join(f"{p}^{e}" for p, e in pairs)}
        if fn == "divisors":
            ds = sorted(sympy.divisors(int(req["n"])))
            return {"ok": True, "result": ",".join(str(x) for x in ds)}
        if fn == "igcdex":
            try:
                from sympy import igcdex as _igcdex
            except ImportError:
                from sympy.core.intfunc import igcdex as _igcdex
            x, y, g = _igcdex(int(req["a"]), int(req["b"]))
            return {"ok": True, "result": f"{x},{y},{g}"}
        if fn == "jacobi":
            v = int(sympy.jacobi_symbol(int(req["a"]), int(req["n"])))
            return {"ok": True, "result": str(v)}
        if fn == "continued_fraction":
            cf = sympy.continued_fraction(sympy.sympify(req["n"]))
            return {"ok": True, "result": ",".join(str(x) for x in cf)}
        if fn == "n_order":
            from sympy.ntheory import n_order as _no
            return {"ok": True, "result": str(int(_no(int(req["a"]), int(req["n"]))))}
        if fn == "primitive_root":
            from sympy.ntheory import primitive_root as _pr
            v = _pr(int(req["n"]))
            return {"ok": True, "result": "None" if v is None else str(int(v))}
        if fn == "sqrt_mod":
            from sympy.ntheory import sqrt_mod as _sm
            v = _sm(int(req["a"]), int(req["n"]))
            return {"ok": True, "result": "None" if v is None else str(int(v))}
        if fn == "crt":
            from sympy.ntheory.modular import crt as _crt
            m = [int(x) for x in req["moduli"]]
            v = [int(x) for x in req["residues"]]
            res = _crt(m, v)
            return {"ok": True,
                    "result": "None" if res is None else f"{int(res[0])},{int(res[1])}"}
        if fn == "discrete_log":
            from sympy.ntheory.residue_ntheory import discrete_log as _dl
            try:
                v = _dl(int(req["n"]), int(req["a"]), int(req["b"]))
                return {"ok": True, "result": str(int(v))}
            except ValueError:
                return {"ok": True, "result": "None"}
        if fn == "pell":
            from sympy.solvers.diophantine.diophantine import diop_DN
            sols = diop_DN(int(req["D"]), 1)
            return {"ok": True,
                    "result": "None" if not sols else f"{int(sols[0][0])},{int(sols[0][1])}"}
        if fn == "legendre":
            from sympy import legendre_symbol
            return {"ok": True, "result": str(int(legendre_symbol(int(req["a"]), int(req["p"]))))}
        if fn == "quad_residues":
            from sympy.ntheory.residue_ntheory import quadratic_residues
            qr = quadratic_residues(int(req["n"]))
            return {"ok": True, "result": ",".join(str(int(x)) for x in qr)}
        if fn == "reduced_totient":
            from sympy.functions.combinatorial.numbers import reduced_totient
            return {"ok": True, "result": str(int(reduced_totient(int(req["n"]))))}
        if fn == "nextprime":
            return {"ok": True, "result": str(int(sympy.nextprime(int(req["n"]))))}
        if fn == "prevprime":
            return {"ok": True, "result": str(int(sympy.prevprime(int(req["n"]))))}
        if fn == "primorial":
            return {"ok": True, "result": str(int(sympy.primorial(int(req["n"]))))}
        if fn == "multiplicity":
            from sympy.ntheory import multiplicity as _mult
            return {"ok": True, "result": str(int(_mult(int(req["p"]), int(req["n"]))))}
        return {"ok": False, "error": "BadFn", "detail": f"unknown ntheory fn: {fn!r}"}

    return {"ok": False, "error": "UnknownOp", "detail": f"unknown op: {op!r}"}


def main():
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        req = None
        try:
            req = json.loads(line)
            resp = handle(req)
            resp["id"] = req.get("id")
        except Exception as e:
            resp = {
                "id": req.get("id") if isinstance(req, dict) else None,
                "ok": False,
                "error": type(e).__name__,
                "detail": str(e),
                "traceback": traceback.format_exc(),
            }
        # The C++ harness issues "shutdown" then closes the pipe without
        # reading the reply, so stop here rather than writing into a closed
        # read end (which would raise BrokenPipeError).
        if resp.get("_shutdown"):
            break
        try:
            sys.stdout.write(json.dumps(resp) + "\n")
            sys.stdout.flush()
        except BrokenPipeError:
            # Reader went away mid-stream — nothing left to answer.
            break


if __name__ == "__main__":
    main()
