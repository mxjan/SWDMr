#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// See http://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf for NumericVector
//
// derivative of second-order ODE: Damped, driven harmonic oscillator
// Function defining derivatives dx/dt and dv/dt
// http://matlab-monkey.com/ODE/resonance/resonance.html


// SDDHO_RungeKutta function to solve the second order differential equation using either
// Runge-Kutta 4th order

// In this model, the mass is considered as 1.

// y: is a 2 value vector with initial position and speed
// time: is a vector of time step
// force: is a vector of applied force to the system (required to be the same length as time)
// gamma: damping factor
// k: spring constant (squared natural oscillation with a mass of 1)

// A sin force of 24h period is added to the model

//' @export
// [[Rcpp::export]]
List SDDHO_SinF_RungeKutta(NumericVector y, NumericVector time, NumericVector force, double gamma, double k, double AmpSin, double PhiSin, double PerSin){
  
  // define pi
  const double pi = 3.14159265358979323846;
  
  // number of step
  int Nstep = time.length();
  
  // Vector of position
  NumericVector Vx(Nstep+1);
  
  // Vector of speed
  NumericVector Vv(Nstep+1);
  
  // Compute Iteration step (h)
  double h  = time[2]-time[1];
  double h_half = h/2;
  double h_sq = pow(h,2);
  
  // Vector of time that will be return
  NumericVector TimeReturn(Nstep+1);
  
  
  // position and speed value
  double x = y[0];
  double v = y[1];
  // value for x or v that will be updated during steps
  double xstep;
  double vstep;
  
  // Insert initial value
  Vx[0] = x;
  Vv[0] = v;
  TimeReturn[0] = time[0]-h;
  
  // Omega used in sin-wave force
  double omegaSin;
  omegaSin = 2*pi/PerSin;
  
  //  ---------- Run RK4 ----------
  double k1;
  double k2;
  double k3;
  double k4;
  double F; // force for time t
  double sddho_sinF_derivsecfun(double x,double v, double gamma, double k, double Force); // second derivative
  
  for (int i = 0; i < Nstep; i ++ ){
    
    // force applied for this time step
    F = force[i];
    
    // add sin-wave force
    F = F + AmpSin*sin(time[i]*omegaSin+PhiSin);
    
    // Start RK4
    k1 = sddho_sinF_derivsecfun(x, v, gamma, k, F);
    
    // k2
    xstep = x + (h_half)*v;
    vstep = v + (h_half)*k1;
    k2 = sddho_sinF_derivsecfun (xstep, vstep, gamma, k, F);
    
    // k3
    xstep = x + (h_half)*v + ( (h_sq/4) * k1);
    vstep = v + (h_half)*k2;
    k3 = sddho_sinF_derivsecfun (xstep, vstep, gamma, k, F);
    
    // k4
    xstep = x + h*v + ( (h_sq/2) * k2);
    vstep = v + h*k3;
    k4 = sddho_sinF_derivsecfun (xstep, vstep, gamma, k, F);
    
    // compute new x and v
    x = x + h*v + (h_sq/6) * (k1+k2+k3);
    v = v + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    
    Vx[i+1] = x;
    Vv[i+1] = v;
    TimeReturn[i+1] = time[i];
  }
  
  //  ---------- end ----------
  
  return List::create(
    _["time"] = TimeReturn,
    _["y1"] = Vx,
    _["y2"] = Vv
  );
  
}

// Derivative for RK4
double sddho_sinF_derivsecfun(double x, double v, double gamma, double k, double Force) {
  return Force-gamma*v-k*x;
}


// Version to compute RK4 using two first order deriv.
// NumericVector sddho_derivfun(NumericVector y, double gamma, double omega, double Force) {
//  NumericVector newy(2);
//  newy[0] = y[1]; // dx/dt
//  newy[1] = Force - gamma*y[1] - pow(omega,2)*y[0]; // dv/dt
//  return newy;
//}
//