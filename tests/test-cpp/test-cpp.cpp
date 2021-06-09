#include <Rcpp.h>
#include <iostream>
#include <typeinfo>

using namespace Rcpp;

struct A{};

// [[Rcpp::export]]
void printClass(SEXP x){
  std::cout << typeid(x).name() << std::endl;
}