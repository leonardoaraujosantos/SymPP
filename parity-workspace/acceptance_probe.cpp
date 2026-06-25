// Parity acceptance probe for Tier 1 + Tier 2. Prints PASS/FAIL per case and a
// final "REMAINING=N / TOTAL=M". The completion loop runs until REMAINING=0.
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <sympp/sympp.hpp>
#include <sympp/core/rewrite.hpp>
#include <sympp/functions/hypergeometric.hpp>
#include <sympp/simplify/hyperexpand.hpp>
#include "oracle/oracle.hpp"
using namespace sympp;
static auto& O(){ return sympp::testing::Oracle::instance(); }
static int total=0, fail=0;
static bool eqv(const std::string&a,const std::string&b){ try{return O().equivalent(a,b);}catch(...){return false;} }
static void want(const std::string&label,bool ok){
  ++total; if(!ok){++fail; std::cout<<"FAIL  "<<label<<"\n";} else std::cout<<"pass  "<<label<<"\n";
}
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
  want(label, eqv(got->str(), expected));
}
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

  std::cout<<"\nREMAINING="<<fail<<" / TOTAL="<<total<<"\n";
  return 0;
}
