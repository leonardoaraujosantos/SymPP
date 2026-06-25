## 1. Printers

- [x] 1.1 `glsl_code` (GLSL float syntax, `pow`/`sqrt`/`inversesqrt`, PI/E).
- [x] 1.2 `dot` (Graphviz digraph of the expression tree).
- [x] 1.3 `srepr` (SymPy constructor form, oracle-cross-checked).

## 2. Codegen AST

- [x] 2.1 `Assignment` / `CodeBlock` / `FunctionDefinition`.
- [x] 2.2 `cse_codeblock` reusing the `cse` pass for shared temporaries.
- [x] 2.3 `emit_c` / `emit_cxx` rendering a complete function.

## 3. Autowrap

- [x] 3.1 `autowrap_available()` (locates $CC / cc / gcc / clang).
- [x] 3.2 `autowrap(vars, body)` → compile C → shared object → dlopen → callable.
- [x] 3.3 Tests gate on `autowrap_available()` and check vs evalf.
