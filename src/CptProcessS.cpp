#include <Rcpp.h>
using namespace Rcpp;

// function to compute the process S, given a vector of wake and sleep

//' @export
// [[Rcpp::export]]
NumericVector CptSprocess(NumericVector Sleep, NumericVector Wake, NumericVector time,
                          double U, double L, double t_w, double t_s, double init){
  
  
  int Nstep = time.length();
  
  // Vector of process S
  NumericVector ProcS(Nstep+1);
  

  // Fill with initial value
  ProcS[0] = init;
  
  // variable for updated S
  double St;
  
  // variable for time spent in wake
  double TinW;
  
  // variable for time spent in sleep
  double TinS;
  
  for (int i = 0; i < Nstep; i ++ ){
    
    // Time spent in wake during this time step
    TinW = Wake[i];
    // Same for sleep
    TinS = Sleep[i];
    
    // First compute S for time spent in wake
    St = U - (U - ProcS[i]) * exp(-TinW / t_w);
    
    // Then compute S for time spent in sleep
    St = L + (St - L) * exp(-TinS / t_s);
    
    
    ProcS[i+1] = St;
    
    
  }
  
  return(ProcS);
  
}