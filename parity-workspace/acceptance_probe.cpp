// Parity acceptance probe for Tier 1 + Tier 2. Prints PASS/FAIL per case and a
// final "REMAINING=N / TOTAL=M". The completion loop runs until REMAINING=0.
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <sympp/sympp.hpp>
#include <sympp/core/rewrite.hpp>
#include <sympp/functions/hypergeometric.hpp>
#include <sympp/simplify/hyperexpand.hpp>
#include <sympp/combinatorics/combinatorics.hpp>
#include <sympp/physics/physics.hpp>
#include "oracle/oracle.hpp"
using namespace sympp;
static auto& O(){ return sympp::testing::Oracle::instance(); }
static int total=0, fail=0;
static bool eqv(const std::string&a,const std::string&b){ try{return O().equivalent(a,b);}catch(...){return false;} }
static void want(const std::string&label,bool ok){
  ++total; if(!ok){++fail; std::cout<<"FAIL  "<<label<<"\n";} else std::cout<<"pass  "<<label<<"\n";
}
static bool coeff_eq_zero(const Expr&e){ return e && simplify(e)==S::Zero(); }
// factor must equal expected AND differ from input
static void factor_is(const std::string&label,const Expr&e,const Expr&var,const std::string&expected){
  std::string g; try{ g=factor(e,var)->str(); }catch(...){ g="THROW"; }
  want(label, eqv(g,expected) && g!=e->str());
}
// indefinite integral: no Integral marker and derivative matches
static void integ_closed(const std::string&label,const Expr&e,const Expr&x){
  std::string g; try{ g=integrate(e,x)->str(); }catch(...){ g="THROW"; }
  bool un = g.find("Integral")!=std::string::npos || g=="THROW";
  want(label, !un && eqv("diff("+g+",x)", e->str()));
}
// definite/expr equals expected closed form
static void expr_is(const std::string&label,const Expr&got,const std::string&expected){
  // Direct string match first: oracle.equivalent computes a-b and cannot
  // certify infinite values (oo - oo = nan), so compare those textually.
  want(label, got->str()==expected || eqv(got->str(), expected));
}
#ifdef PROBE_GROUP
static void probe_group_theory(){
  using namespace sympp::combinatorics;
  std::cout<<"== T2#8 group theory (conjugacy/center/derived) ==\n";
  want("S3 #conj=3", num_conjugacy_classes(symmetric_group(3))==3);
  want("S4 #conj=5", num_conjugacy_classes(symmetric_group(4))==5);
  want("A4 #conj=4", num_conjugacy_classes(alternating_group(4))==4);
  want("D4 #conj=5", num_conjugacy_classes(dihedral_group(4))==5);
  want("S3 center=1", group_center(symmetric_group(3)).order()==1);
  want("D4 center=2", group_center(dihedral_group(4)).order()==2);
  want("S3 derived=3", derived_subgroup(symmetric_group(3)).order()==3);
  want("S4 derived=12", derived_subgroup(symmetric_group(4)).order()==12);
  want("A4 derived=4", derived_subgroup(alternating_group(4)).order()==4);
  want("D4 derived=2", derived_subgroup(dihedral_group(4)).order()==2);
  want("S4 solvable", is_solvable(symmetric_group(4)));
  want("A5 not solvable", !is_solvable(alternating_group(5)));
  want("D4 nilpotent", is_nilpotent(dihedral_group(4)));
  want("S3 not nilpotent", !is_nilpotent(symmetric_group(3)));
  want("C6 nilpotent", is_nilpotent(cyclic_group(6)));
  want("A5 simple", is_simple(alternating_group(5)));
  want("A4 not simple", !is_simple(alternating_group(4)));
  want("C6 abelian_invariants={2,3}",
       abelian_invariants(cyclic_group(6)) == std::vector<long>{2, 3});
  want("D4 abelian_invariants={2,2}",
       abelian_invariants(dihedral_group(4)) == std::vector<long>{2, 2});
  want("A5 abelian_invariants={}", abelian_invariants(alternating_group(5)).empty());
}
#endif

#ifdef PROBE_PHYS
static void probe_physics(){
  using namespace sympp::physics;
  std::cout<<"== T2#9 multi-mode second quantization ==\n";
  // [a_i, a_i^dagger] = 1 (same mode).
  MultiModeFockState s{{2,5}, S::One()};
  want("[a0,a0+]=1 same mode", coeff_eq_zero(apply_boson_commutator(s,0,0).coefficient - S::One()));
  // [a_1, a_2^dagger] = 0 (cross mode commute).
  want("[a0,a1+]=0 cross mode", apply_boson_commutator(s,0,1).is_zero());
  // Number eigenvalue: N_1|2 5> = 5|2 5>.
  want("N1|2 5>=5", coeff_eq_zero(apply_number(s,1).coefficient - integer(5)));
  // Vacuum annihilation a_i|0> = 0.
  MultiModeFockState vac{{0,0}, S::One()};
  want("a0|0>=0 vacuum", apply_annihilation(vac,0).is_zero());
}
#endif

#ifdef PROBE_PERF
// Performance regression for the cyclotomic fast path: factor(x^N − 1) is
// decomposed via ∏_{d|N} Φ_d(x) instead of the super-linear Kronecker search,
// so high-divisor N must complete near-instantly (was >30s for N=60).
static void probe_perf(){
  std::cout<<"== PERF cyclotomic factor(x^N-1) ==\n";
  auto x=symbol("x");
  auto P=[&](const Expr&b,int e){return pow(b,integer(e));};
  auto bench=[&](int N){
    Expr e=expand(P(x,N)-S::One());
    auto t0=std::chrono::steady_clock::now();
    Expr g=factor(e,x);
    auto t1=std::chrono::steady_clock::now();
    long ms=std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
    std::string label="factor(x^"+std::to_string(N)+"-1) "+std::to_string(ms)+"ms";
    bool fast = ms < 2000;
    bool correct = eqv(g->str(), "x**"+std::to_string(N)+"-1");
    want(label+" <2000ms", fast);
    want(label+" expands back", correct);
  };
  bench(36);
  bench(60);
}
#endif

int main(){
  auto x=symbol("x"),y=symbol("y"),z=symbol("z"),n=symbol("n");
  auto oo=S::Infinity(); auto I=[&](int i){return integer(i);};
  auto P=[&](const Expr&b,int e){return pow(b,integer(e));};

  std::cout<<"== T1#1 erf ==\n";
  expr_is("erfc(x)^2 0..oo", integrate(P(erfc(x),2),x,I(0),oo), "(2-sqrt(2))/sqrt(pi)");
  expr_is("x e^-x^2 erf(x) 0..oo", integrate(mul(mul(x,exp(mul(I(-1),P(x,2)))),erf(x)),x,I(0),oo), "sqrt(2)/4");
  expr_is("erfc(x) e^-x^2 0..oo", integrate(mul(erfc(x),exp(mul(I(-1),P(x,2)))),x,I(0),oo), "sqrt(pi)/4");

  std::cout<<"== T2#4 Wang ==\n";
  factor_is("x^2-y^2", P(x,2)-P(y,2), x, "(x-y)*(x+y)");
  factor_is("x^3-y^3", P(x,3)-P(y,3), x, "(x-y)*(x**2+x*y+y**2)");
  factor_is("x^4-y^4", P(x,4)-P(y,4), x, "(x-y)*(x+y)*(x**2+y**2)");
  factor_is("x^2-(y+1)x+y", add({P(x,2), mul(mul(I(-1),y+I(1)),x), y}), x, "(x-1)*(x-y)");
  factor_is("xy-x-y+1", add({mul(x,y),mul(I(-1),x),mul(I(-1),y),I(1)}), x, "(x-1)*(y-1)");
  factor_is("x^3+y^3+z^3-3xyz", add({P(x,3),P(y,3),P(z,3),mul({I(-3),x,y,z})}), x, "(x+y+z)*(x**2-x*y-x*z+y**2-y*z+z**2)");

  std::cout<<"== T2#5 Risch (special-fn antiderivatives) ==\n";
  integ_closed("e^x/x -> Ei", mul(exp(x),P(x,-1)), x);
  integ_closed("1/log x -> li", P(log(x),-1), x);
  integ_closed("e^{x^2} -> erfi", exp(P(x,2)), x);
  integ_closed("sin(x)/x -> Si", mul(sin(x),P(x,-1)), x);
  integ_closed("x e^x/(x+1)^2", mul(mul(x,exp(x)),P(x+I(1),-2)), x);

  std::cout<<"== T2#6 hyperexpand ==\n";
  expr_is("2F1(1,1;2;z)", hyperexpand(hyper({I(1),I(1)},{I(2)},x)), "-log(1-x)/x");
  expr_is("1F1(1;2;z)", hyperexpand(hyper({I(1)},{I(2)},x)), "(exp(x)-1)/x");
  expr_is("0F1(;3/2;-z^2/4)", hyperexpand(hyper({},{rational(3,2)},mul(rational(-1,4),P(x,2)))), "sin(x)/x");
  expr_is("2F1(1/2,1;3/2;z^2)", hyperexpand(hyper({rational(1,2),I(1)},{rational(3,2)},P(x,2))), "atanh(x)/x");

  std::cout<<"== T2#7 Gruntz ==\n";
  expr_is("Gamma(x+1)/Gamma(x+1/2)->oo", limit(mul(gamma(x+I(1)),P(gamma(x+rational(1,2)),-1)),x,oo), "oo");
  expr_is("x^(1/x)->1", limit(pow(x,P(x,-1)),x,oo), "1");
  expr_is("(1+1/x)^x->E", limit(pow(I(1)+P(x,-1),x),x,oo), "E");
  expr_is("x^x/(x+1)^x->1/E", limit(mul(pow(x,x),P(pow(x+I(1),x),-1)),x,oo), "exp(-1)");

  std::cout<<"== T2#4b Wang (general multivariate) ==\n";
  factor_is("x^2-y^2+x-y", add({P(x,2),mul(I(-1),P(y,2)),x,mul(I(-1),y)}), x, "(x-y)*(x+y+1)");
  factor_is("x^2 y-x^2-y+1", add({mul(P(x,2),y),mul(I(-1),P(x,2)),mul(I(-1),y),I(1)}), x, "(x-1)*(x+1)*(y-1)");
  factor_is("(x+y)^2-z^2", P(x+y,2)-P(z,2), x, "(x+y-z)*(x+y+z)");
  factor_is("x^3-x y^2", P(x,3)-mul(x,P(y,2)), x, "x*(x-y)*(x+y)");
  factor_is("x^2+2xy+y^2+2x+2y+1", add({P(x,2),mul({I(2),x,y}),P(y,2),mul(I(2),x),mul(I(2),y),I(1)}), x, "(x+y+1)**2");
  factor_is("x^4+x^2 y^2+y^4", add({P(x,4),mul(P(x,2),P(y,2)),P(y,4)}), x, "(x**2-x*y+y**2)*(x**2+x*y+y**2)");

  std::cout<<"== T2#6b hyperexpand ==\n";
  expr_is("1F1(2;1;z)", hyperexpand(hyper({I(2)},{I(1)},x)), "(x+1)*exp(x)");
  expr_is("1F1(-1;1;z)", hyperexpand(hyper({I(-1)},{I(1)},x)), "1-x");
  expr_is("2F1(-2,1;1;z)", hyperexpand(hyper({I(-2),I(1)},{I(1)},x)), "(1-x)**2");
  expr_is("1F1(2;3;z)", hyperexpand(hyper({I(2)},{I(3)},x)), "(2*x-2)*exp(x)/x**2+2/x**2");
  expr_is("0F1(;1/2;z^2/4)", hyperexpand(hyper({},{rational(1,2)},mul(rational(1,4),P(x,2)))), "cosh(x)");

  std::cout<<"== T2#7b Gruntz (towers) ==\n";
  expr_is("log(x+e^x)/x->1", limit(mul(log(x+exp(x)),P(x,-1)),x,oo), "1");
  expr_is("(e^2x+e^x)^(1/x)->e^2", limit(pow(add({exp(mul(I(2),x)),exp(x)}),P(x,-1)),x,oo), "exp(2)");
  expr_is("x(e^(1/x)-1)->1", limit(mul(x,exp(P(x,-1))+I(-1)),x,oo), "1");
  expr_is("log(log x)/log x->0", limit(mul(log(log(x)),P(log(x),-1)),x,oo), "0");
  expr_is("e^(x^2)/(e^x)^x->1", limit(mul(exp(P(x,2)),P(pow(exp(x),x),-1)),x,oo), "1");

  std::cout<<"== T2#5b Risch / integration ==\n";
  integ_closed("log(x)^2/x", mul(P(log(x),2),P(x,-1)), x);
  integ_closed("x log(x)^2", mul(x,P(log(x),2)), x);
  integ_closed("tan(x)", tan(x), x);
  integ_closed("sec(x)", sec(x), x);
  integ_closed("x/sqrt(1-x^2)", mul(x,pow(I(1)-P(x,2),rational(-1,2))), x);
  integ_closed("e^{2x} sin x", mul(exp(mul(I(2),x)),sin(x)), x);
  integ_closed("(x-1)e^x/(x+1)^3", mul(mul(x-I(1),exp(x)),P(x+I(1),-3)), x);
  integ_closed("1/(x^2-1)", P(P(x,2)-I(1),-1), x);
  integ_closed("(x+1)/(x^2+2x+2)", mul(x+I(1),pow(add({P(x,2),mul(I(2),x),I(2)}),I(-1))), x);
  integ_closed("atan(x)", atan(x), x);
  integ_closed("e^{sqrt x}", exp(pow(x,rational(1,2))), x);
  integ_closed("e^x cos x", mul(exp(x),cos(x)), x);
  integ_closed("x^2 e^x", mul(P(x,2),exp(x)), x);
  integ_closed("log(x^2+1)", log(P(x,2)+I(1)), x);

  std::cout<<"== T2#4c Wang ==\n";
  factor_is("x^4-5x^2y^2+4y^4", add({P(x,4),mul({I(-5),P(x,2),P(y,2)}),mul(I(4),P(y,4))}), x, "(x-2*y)*(x-y)*(x+y)*(x+2*y)");
  factor_is("x^2y^2-x^2-y^2+1", add({mul(P(x,2),P(y,2)),mul(I(-1),P(x,2)),mul(I(-1),P(y,2)),I(1)}), x, "(x-1)*(x+1)*(y-1)*(y+1)");
  factor_is("x^3+x^2-xy^2-y^2", add({P(x,3),P(x,2),mul(I(-1),mul(x,P(y,2))),mul(I(-1),P(y,2))}), x, "(x+1)*(x-y)*(x+y)");

  std::cout<<"== T2#6c hyperexpand ==\n";
  expr_is("2F1(2,1;1;z)", hyperexpand(hyper({I(2),I(1)},{I(1)},x)), "(1-x)**(-2)");

  std::cout<<"== T2#7c Gruntz (adversarial) ==\n";
  expr_is("e^e^x/e^(e^x+x)->0", limit(mul(exp(exp(x)),pow(exp(add({exp(x),x})),I(-1))),x,oo), "0");
  expr_is("log(x)^(1/x)->1", limit(pow(log(x),P(x,-1)),x,oo), "1");
  expr_is("e^x(e^(1/x)-1)->oo", limit(mul(exp(x),add({exp(P(x,-1)),I(-1)})),x,oo), "oo");
  expr_is("x^(1/log log x)->oo", limit(pow(x,pow(log(log(x)),I(-1))),x,oo), "oo");

  std::cout<<"== T2#4d Wang (4-variable) ==\n";
  factor_is("ab-ac+bd-cd", add({mul(symbol("a"),symbol("b")),mul(I(-1),mul(symbol("a"),symbol("c"))),mul(symbol("b"),symbol("d")),mul(I(-1),mul(symbol("c"),symbol("d")))}), symbol("a"), "(a+d)*(b-c)");
  factor_is("wx+wy+zx+zy", add({mul(symbol("w"),x),mul(symbol("w"),y),mul(z,x),mul(z,y)}), symbol("w"), "(w+z)*(x+y)");
  factor_is("x^2y^2-z^2w^2", P(mul(x,y),2)-P(mul(z,symbol("w")),2), x, "(x*y-w*z)*(x*y+w*z)");

  std::cout<<"== T2#5c Risch residual ==\n";
  integ_closed("x^2 e^{-x^2}", mul(P(x,2),exp(mul(I(-1),P(x,2)))), x);
  integ_closed("1/(x log^2 x)", P(mul(x,P(log(x),2)),-1), x);
  integ_closed("cos(x)/x -> Ci", mul(cos(x),P(x,-1)), x);
  expr_is("x^2 e^{-x^2} 0..oo", integrate(mul(P(x,2),exp(mul(I(-1),P(x,2)))),x,I(0),oo), "sqrt(pi)/4");

  std::cout<<"== T2#6d hyperexpand (polylog) ==\n";
  expr_is("3F2(1,1,1;2,2;z)", hyperexpand(hyper({I(1),I(1),I(1)},{I(2),I(2)},x)), "polylog(2,x)/x");

#ifdef PROBE_GROUP
  probe_group_theory();
#endif
#ifdef PROBE_PHYS
  probe_physics();
#endif
#ifdef PROBE_PERF
  probe_perf();
#endif

  std::cout<<"\nREMAINING="<<fail<<" / TOTAL="<<total<<"\n";
  return 0;
}
