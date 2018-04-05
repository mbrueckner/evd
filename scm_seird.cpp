#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
imat scmSEIRDcxx(vec beta, vec delta, double gamma, vec f, double cp, int tbegin, int tend,
		 int N, int E0, int I0) {
  
  irowvec state = { N - E0 - I0, E0, I0, 0, 0 };
    
  imat steps = { {-1,  1,  0, 0, 0 },
		 { 0, -1,  1, 0, 0 },
		 { 0,  0, -1, 1, 0 },
		 { 0,  0, -1, 0, 1 } };
  
  beta = beta / N;

  double time = tbegin;
  
  unsigned int i = 0;
  NumericVector a(4);
  
  unsigned int nr = 1.5*N;
  imat out(nr, 7);
  
  while((time < tend) & (state(1) > 0 | state(2) > 0)) {
    // calculate process probabilities for current state
    // S -> E, E -> I, I -> R, I -> D
    if(time < cp) {
      a = { beta(0)*state(0)*state(2),
	    gamma*state(1),
	    (1-f(0))*delta(0)*state(2),
	    f(0)*delta(0)*state(2) };
    } else {    
      a = { beta(1)*state(0)*state(2),
	    gamma*state(1),
	    (1-f(1))*delta(1)*state(2),
	    f(1)*delta(1)*state(2) };
    }
        
    // WHEN does the next process happen?
    double rate = sum(a);      
    time = time + rexp(1, rate)(0);
    
    // WHICH process happens after tau?
    NumericVector p = a / rate;
    unsigned int tt = sample(4, 1, false, p)(0) - 1;
    state = state + steps.row(tt);

    // grow output matrix if needed
    if(i >= nr) {
      nr = nr + 100000;
      out.resize(nr, 7);
    }

    out(i, 0) = round(time);
    out(i, span(1, 5)) = state;
    out(i, 6) = (tt == 1);
    
    i = i+1;
  }

  return resize(out, i, 7);
}
