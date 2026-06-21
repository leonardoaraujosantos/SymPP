#include <sympp/core/lambdify_llvm.hpp>

#include <stdexcept>

#ifndef SYMPP_HAVE_LLVM_JIT

namespace sympp {

bool lambdify_jit_available() noexcept { return false; }

CompiledExpr lambdify_jit(const std::vector<Expr>&, const Expr&) {
    throw std::runtime_error(
        "lambdify_jit: SymPP was built without the LLVM JIT backend "
        "(use lambdify() for the portable interpreter)");
}
CompiledExpr lambdify_jit(const Expr&, const Expr&) {
    throw std::runtime_error("lambdify_jit: SymPP was built without the LLVM JIT backend");
}

}  // namespace sympp

#else  // SYMPP_HAVE_LLVM_JIT

#include <cmath>
#include <memory>
#include <mutex>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include <mpfr.h>

#include <llvm/ExecutionEngine/Orc/LLJIT.h>
#include <llvm/ExecutionEngine/Orc/ExecutionUtils.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/Intrinsics.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/Verifier.h>
#include <llvm/Support/TargetSelect.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

bool lambdify_jit_available() noexcept { return true; }

namespace {

void init_llvm_once() {
    static std::once_flag flag;
    std::call_once(flag, [] {
        llvm::InitializeNativeTarget();
        llvm::InitializeNativeTargetAsmPrinter();
    });
}

[[noreturn]] void unsupported(const std::string& what) {
    throw std::runtime_error("lambdify_jit: unsupported in real evaluation: " + what);
}

[[nodiscard]] double number_symbol_value(NumberSymbolKind k) {
    switch (k) {
        case NumberSymbolKind::Pi: return M_PI;
        case NumberSymbolKind::E: return M_E;
        case NumberSymbolKind::EulerGamma: return 0.57721566490153286060651209;
        case NumberSymbolKind::Catalan: return 0.91596559417721901505460351;
    }
    unsupported("number symbol");
}

// IR generator: walks an expression into LLVM Values (all doubles).
class Codegen {
public:
    Codegen(llvm::Module& m, llvm::IRBuilder<>& b, llvm::Value* args,
            const std::unordered_map<std::string, std::size_t>& idx)
        : m_(m), b_(b), args_(args), idx_(idx), dbl_(b.getDoubleTy()) {}

    llvm::Value* gen(const Expr& e) {
        switch (e->type_id()) {
            case TypeId::Integer:
                return cst(static_cast<const Integer&>(*e).value().get_d());
            case TypeId::Rational:
                return cst(static_cast<const Rational&>(*e).value().get_d());
            case TypeId::Float:
                return cst(mpfr_get_d(static_cast<const Float&>(*e).value(), MPFR_RNDN));
            case TypeId::NumberSymbol:
                return cst(number_symbol_value(static_cast<const NumberSymbol&>(*e).kind()));
            case TypeId::Symbol: {
                auto it = idx_.find(static_cast<const Symbol&>(*e).name());
                if (it == idx_.end()) unsupported("free symbol");
                llvm::Value* gep = b_.CreateConstInBoundsGEP1_64(dbl_, args_, it->second);
                return b_.CreateLoad(dbl_, gep);
            }
            case TypeId::Add: return fold(e, 0.0, /*mul=*/false);
            case TypeId::Mul: return fold(e, 1.0, /*mul=*/true);
            case TypeId::Pow:
                return intr2(llvm::Intrinsic::pow, gen(e->args()[0]), gen(e->args()[1]));
            case TypeId::Function: return gen_function(e);
            case TypeId::Infinity: return cst(HUGE_VAL);
            case TypeId::NegativeInfinity: return cst(-HUGE_VAL);
            case TypeId::ImaginaryUnit: unsupported("imaginary unit I");
            default: unsupported("node type");
        }
    }

private:
    llvm::Value* cst(double v) { return llvm::ConstantFP::get(dbl_, v); }

    llvm::Value* fold(const Expr& e, double identity, bool mul) {
        llvm::Value* acc = cst(identity);
        bool first = true;
        for (const auto& a : e->args()) {
            llvm::Value* v = gen(a);
            if (first) { acc = v; first = false; continue; }
            acc = mul ? b_.CreateFMul(acc, v) : b_.CreateFAdd(acc, v);
        }
        return acc;
    }

    llvm::Value* intr1(llvm::Intrinsic::ID id, llvm::Value* x) {
        llvm::Function* fn = llvm::Intrinsic::getDeclaration(&m_, id, {dbl_});
        return b_.CreateCall(fn, {x});
    }
    llvm::Value* intr2(llvm::Intrinsic::ID id, llvm::Value* x, llvm::Value* y) {
        llvm::Function* fn = llvm::Intrinsic::getDeclaration(&m_, id, {dbl_});
        return b_.CreateCall(fn, {x, y});
    }
    // Call a libm function declared extern (resolved by the JIT from the process).
    llvm::Value* libm(const char* name, llvm::Value* x) {
        llvm::FunctionType* ft = llvm::FunctionType::get(dbl_, {dbl_}, false);
        llvm::FunctionCallee fc = m_.getOrInsertFunction(name, ft);
        return b_.CreateCall(fc, {x});
    }
    llvm::Value* libm2(const char* name, llvm::Value* x, llvm::Value* y) {
        llvm::FunctionType* ft = llvm::FunctionType::get(dbl_, {dbl_, dbl_}, false);
        llvm::FunctionCallee fc = m_.getOrInsertFunction(name, ft);
        return b_.CreateCall(fc, {x, y});
    }
    llvm::Value* recip(llvm::Value* x) { return b_.CreateFDiv(cst(1.0), x); }

    llvm::Value* gen_function(const Expr& e) {
        const auto& fn = static_cast<const Function&>(*e);
        const auto& as = fn.args();
        llvm::Value* a = as.empty() ? nullptr : gen(as[0]);
        switch (fn.function_id()) {
            case FunctionId::Sin: return intr1(llvm::Intrinsic::sin, a);
            case FunctionId::Cos: return intr1(llvm::Intrinsic::cos, a);
            case FunctionId::Exp: return intr1(llvm::Intrinsic::exp, a);
            case FunctionId::Abs: return intr1(llvm::Intrinsic::fabs, a);
            case FunctionId::Floor: return intr1(llvm::Intrinsic::floor, a);
            case FunctionId::Ceiling: return intr1(llvm::Intrinsic::ceil, a);
            case FunctionId::Tan: return libm("tan", a);
            case FunctionId::Asin: return libm("asin", a);
            case FunctionId::Acos: return libm("acos", a);
            case FunctionId::Atan: return libm("atan", a);
            case FunctionId::Sinh: return libm("sinh", a);
            case FunctionId::Cosh: return libm("cosh", a);
            case FunctionId::Tanh: return libm("tanh", a);
            case FunctionId::Asinh: return libm("asinh", a);
            case FunctionId::Acosh: return libm("acosh", a);
            case FunctionId::Atanh: return libm("atanh", a);
            case FunctionId::Gamma: return libm("tgamma", a);
            case FunctionId::LogGamma: return libm("lgamma", a);
            case FunctionId::Erf: return libm("erf", a);
            case FunctionId::Erfc: return libm("erfc", a);
            case FunctionId::Cot: return recip(libm("tan", a));
            case FunctionId::Sec: return recip(intr1(llvm::Intrinsic::cos, a));
            case FunctionId::Csc: return recip(intr1(llvm::Intrinsic::sin, a));
            case FunctionId::Coth: return recip(libm("tanh", a));
            case FunctionId::Sech: return recip(libm("cosh", a));
            case FunctionId::Csch: return recip(libm("sinh", a));
            case FunctionId::Re:
            case FunctionId::Conjugate: return a;
            case FunctionId::Im: return cst(0.0);
            case FunctionId::Factorial:
                return libm("tgamma", b_.CreateFAdd(a, cst(1.0)));
            case FunctionId::Log: {
                if (as.size() == 1) return intr1(llvm::Intrinsic::log, a);
                return b_.CreateFDiv(intr1(llvm::Intrinsic::log, a),
                                     intr1(llvm::Intrinsic::log, gen(as[1])));
            }
            case FunctionId::Atan2:
                return libm2("atan2", a, gen(as[1]));
            case FunctionId::Min:
                return fold_minmax(as, llvm::Intrinsic::minnum);
            case FunctionId::Max:
                return fold_minmax(as, llvm::Intrinsic::maxnum);
            default:
                unsupported("function '" + std::string(fn.name()) + "'");
        }
    }

    llvm::Value* fold_minmax(const std::span<const Expr>& as, llvm::Intrinsic::ID id) {
        llvm::Value* acc = gen(as[0]);
        for (std::size_t i = 1; i < as.size(); ++i) acc = intr2(id, acc, gen(as[i]));
        return acc;
    }

    llvm::Module& m_;
    llvm::IRBuilder<>& b_;
    llvm::Value* args_;
    const std::unordered_map<std::string, std::size_t>& idx_;
    llvm::Type* dbl_;
};

}  // namespace

CompiledExpr lambdify_jit(const std::vector<Expr>& vars, const Expr& body) {
    init_llvm_once();

    std::unordered_map<std::string, std::size_t> idx;
    for (std::size_t i = 0; i < vars.size(); ++i) {
        if (!vars[i] || vars[i]->type_id() != TypeId::Symbol) {
            throw std::runtime_error("lambdify_jit: variables must be Symbols");
        }
        idx.emplace(static_cast<const Symbol&>(*vars[i]).name(), i);
    }

    auto ctx = std::make_unique<llvm::LLVMContext>();
    auto mod = std::make_unique<llvm::Module>("sympp_lambdify", *ctx);
    llvm::IRBuilder<> builder(*ctx);

    llvm::Type* dbl = builder.getDoubleTy();
    llvm::Type* dblPtr = llvm::PointerType::getUnqual(*ctx);
    llvm::FunctionType* ft = llvm::FunctionType::get(dbl, {dblPtr}, false);
    llvm::Function* f =
        llvm::Function::Create(ft, llvm::Function::ExternalLinkage, "sympp_fn", mod.get());
    llvm::BasicBlock* bb = llvm::BasicBlock::Create(*ctx, "entry", f);
    builder.SetInsertPoint(bb);

    Codegen cg(*mod, builder, f->getArg(0), idx);
    builder.CreateRet(cg.gen(body));

    if (llvm::verifyFunction(*f, &llvm::errs())) {
        throw std::runtime_error("lambdify_jit: generated invalid IR");
    }

    auto jit_or = llvm::orc::LLJITBuilder().create();
    if (!jit_or) throw std::runtime_error("lambdify_jit: failed to create LLJIT");
    std::shared_ptr<llvm::orc::LLJIT> jit = std::move(*jit_or);

    // Resolve libm symbols (tan, tgamma, …) from the host process.
    auto gen = llvm::orc::DynamicLibrarySearchGenerator::GetForCurrentProcess(
        jit->getDataLayout().getGlobalPrefix());
    if (!gen) throw std::runtime_error("lambdify_jit: symbol generator failed");
    jit->getMainJITDylib().addGenerator(std::move(*gen));

    if (auto err = jit->addIRModule(
            llvm::orc::ThreadSafeModule(std::move(mod), std::move(ctx)))) {
        throw std::runtime_error("lambdify_jit: failed to add module");
    }

    auto sym = jit->lookup("sympp_fn");
    if (!sym) throw std::runtime_error("lambdify_jit: symbol lookup failed");
    auto fptr = sym->toPtr<double (*)(const double*)>();

    std::size_t arity = vars.size();
    return [jit, fptr, arity](const std::vector<double>& xs) {
        if (xs.size() != arity) {
            throw std::runtime_error("lambdify_jit: wrong number of arguments");
        }
        return fptr(xs.data());
    };
}

CompiledExpr lambdify_jit(const Expr& var, const Expr& body) {
    return lambdify_jit(std::vector<Expr>{var}, body);
}

}  // namespace sympp

#endif  // SYMPP_HAVE_LLVM_JIT
