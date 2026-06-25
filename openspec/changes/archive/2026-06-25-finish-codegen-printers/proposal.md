## Why

Phase 13 (code generation) shipped the C/C++/Fortran/LaTeX/Octave/Rust/Julia/
MathML printers and string function emission, but three SymPy code-generation
capabilities were still missing: the GLSL shader printer, the Graphviz-`dot`
expression-tree printer, the `srepr` constructor printer, a structured codegen
AST, and `autowrap` (compile-to-native). Finishing them closes Tier 3 #12 and
brings the code-generation surface to SymPy parity.

## What Changes

- Add `printing::glsl_code` (GLSL shader syntax), `printing::dot` (Graphviz
  digraph of the expression tree), and `printing::srepr` (SymPy constructor
  form, oracle-validated).
- Add a `sympp::codegen` module with a structured AST — `Assignment`,
  `CodeBlock`, `FunctionDefinition` — plus `cse_codeblock` (CSE via the existing
  `cse` pass) and `emit_c`/`emit_cxx` that render a complete function with
  shared temporaries.
- Add `codegen::autowrap` / `autowrap_available`: emit C, compile it with the
  system C compiler into a shared object, `dlopen`/`dlsym`, and return a
  `CompiledExpr` (SymPy's autowrap C backend).
- Mark Phase 13 ✅ in the README and roadmap.

## Capabilities

- **code-generation** (new): GLSL/dot/srepr printers, a codegen AST with CSE,
  and `autowrap` compile-to-native.

## Non-goals

- Cython/f2py-style autowrap backends (the C/dlopen backend is provided).
- A GLSL/dot oracle in SymPy (those printers are validated by exact-string
  unit tests; `srepr` is cross-checked against SymPy).
