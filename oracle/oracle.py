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
        sys.stdout.write(json.dumps(resp) + "\n")
        sys.stdout.flush()
        if resp.get("_shutdown"):
            break


if __name__ == "__main__":
    main()
