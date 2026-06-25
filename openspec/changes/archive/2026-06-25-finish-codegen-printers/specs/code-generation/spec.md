# code-generation Specification

## ADDED Requirements

### Requirement: GLSL, Graphviz-dot and srepr printers

SymPP SHALL provide `printing::glsl_code(e)` rendering an expression in GLSL
shader syntax (float literals, `pow`/`sqrt`/`inversesqrt`, named `PI`/`E`),
`printing::dot(e)` rendering a Graphviz `digraph` of the expression tree (one
node per subexpression, directed parent→child edges), and `printing::srepr(e)`
rendering SymPy's constructor form. `srepr` SHALL match SymPy's `srepr` for the
common constructors (`Add`/`Mul`/`Pow`/`Symbol`/`Integer`/`Rational`/function),
modulo documented argument-ordering differences.

#### Scenario: GLSL power
- **WHEN** `glsl_code(x**3)` is evaluated
- **THEN** the result is `pow(x, 3.0)`

#### Scenario: dot is a digraph
- **WHEN** `dot(2*x)` is evaluated
- **THEN** the result begins a `digraph` and contains a `Mul` node with edges to
  its `Integer(2)` and `Symbol(x)` children

#### Scenario: srepr matches SymPy
- **WHEN** `srepr(x**2)` is evaluated
- **THEN** the result equals SymPy's `srepr(x**2)` = `Pow(Symbol('x'), Integer(2))`

### Requirement: Structured codegen AST

SymPP SHALL provide a `sympp::codegen` AST — `Assignment` (a `Symbol` lhs and an
`Expr` rhs), `CodeBlock` (a list of assignments plus a returned expression), and
`FunctionDefinition` (name, argument symbols, body) — together with
`cse_codeblock(expr)` that builds a `CodeBlock` whose repeated subexpressions are
reified into shared temporaries via common-subexpression elimination, and
`emit_c`/`emit_cxx` that render a complete function with those temporaries.

#### Scenario: CSE reifies a repeated subexpression
- **WHEN** `emit_c` is run on `f(x,y) = sin(x*y) + (x*y)**2 + x*y`
- **THEN** the emitted function declares one temporary for `x*y` and returns the
  expression in terms of that temporary

### Requirement: Autowrap compile-to-native

SymPP SHALL provide `codegen::autowrap(vars, body)` that emits C for the body,
compiles it with the system C compiler into a shared object, loads it via
`dlopen`/`dlsym`, and returns a `CompiledExpr` callable whose result matches
numeric evaluation of the expression. `codegen::autowrap_available()` SHALL
report whether a usable C compiler is present; `autowrap` SHALL throw when none
is available.

#### Scenario: Compiled callable matches evalf
- **WHEN** `autowrap_available()` is true and `f = autowrap(x, x*x + 1)` is built
- **THEN** `f({3.0})` equals `10.0` (matching `evalf` of `x**2 + 1` at `x = 3`)

#### Scenario: No compiler
- **WHEN** no C compiler is available
- **THEN** `autowrap_available()` is false and `autowrap(...)` throws
