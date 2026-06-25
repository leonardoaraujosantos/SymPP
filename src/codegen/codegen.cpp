#include <sympp/codegen/codegen.hpp>

#include <sstream>
#include <string>

#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/printing/printing.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::codegen {

CodeBlock cse_codeblock(const Expr& expr) {
    auto result = cse(expr);

    CodeBlock block;
    block.body.reserve(result.substitutions.size());

    // simplify::cse names its temporaries "_cse_0", "_cse_1", ...; rename them
    // to friendlier "t0", "t1", ... and propagate the renaming through later
    // definitions and the returned expression.
    Expr ret = result.expr;
    int idx = 0;
    for (const auto& [old_temp, definition] : result.substitutions) {
        Expr fresh = symbol("t" + std::to_string(idx++));
        Expr rhs = subs(definition, old_temp, fresh);
        // Rewrite already-collected definitions too, in case a later temp
        // referenced an earlier one under its old name.
        for (auto& assignment : block.body) {
            assignment.rhs = subs(assignment.rhs, old_temp, fresh);
        }
        ret = subs(ret, old_temp, fresh);
        block.body.push_back(Assignment{fresh, std::move(rhs)});
    }
    block.ret = std::move(ret);
    return block;
}

namespace {

using PrinterFn = std::string (*)(const Expr&);

std::string emit_function(const FunctionDefinition& fn, PrinterFn print) {
    std::ostringstream os;
    os << "double " << fn.name << "(";
    for (std::size_t i = 0; i < fn.args.size(); ++i) {
        if (i > 0) os << ", ";
        os << "double " << fn.args[i]->str();
    }
    os << ") {\n";
    for (const auto& assignment : fn.block.body) {
        os << "    const double " << assignment.lhs->str() << " = "
           << print(assignment.rhs) << ";\n";
    }
    os << "    return " << print(fn.block.ret) << ";\n}";
    return os.str();
}

}  // namespace

std::string emit_c(const FunctionDefinition& fn) {
    return emit_function(fn, printing::ccode);
}

std::string emit_cxx(const FunctionDefinition& fn) {
    return emit_function(fn, printing::cxxcode);
}

}  // namespace sympp::codegen
