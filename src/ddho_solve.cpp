#include <Rcpp.h>
#include <math.h>
#include <complex.h>
using namespace Rcpp;

// Authors: Anne Skeldon, Maxime Jan
// Solution from Anne Skeldon <a.skeldon@surray.ac.uk> 
// Adapted by Maxime Jan <Maxime.Jan@unil.ch>

// Original equation have some precision trouble for high time and gamma
// adapted by evaluating transient term between t0 = 0 and t_end = t0 + step. 
// Phase of sinusoidal force is adapted at each step

//' @export
// [[Rcpp::export]]
List Solve_SDDHO_longdouble(NumericVector y, NumericVector time, NumericVector force, long double gamma, long double k,long double AmpSin, long double PhiSin, long double PerSin){
  
  // define pi
  const long double pi = 3.14159265358979323846264338327950288419716939937510L;
  const long double Twopi = 2.0L * pi;
  
  // number of step
  int Nstep = time.length();
  
  // Vector of position
  NumericVector Valpha(Nstep+1);
  
  // Vector of speed
  NumericVector Vbeta(Nstep+1);
  
  // Compute Iteration step (h)
  long double h  = time[1]-time[0];
  
  // Vector of time that will be returned
  NumericVector TimeReturn(Nstep+1);
  
  // inital position and speed value
  long double alpha = y[0];
  long double beta = y[1];
  
  // Insert initial values
  Valpha[0] = alpha;
  Vbeta[0] = beta;
  TimeReturn[0] = time[0]-h;
  
  //transcient reponse frequency
  long double omega0 = sqrtl(k);
  long double gamma_sq = pow(gamma,2);
  std::complex<long double> omega1 = std::sqrt(std::complex<long double>(k - gamma_sq/4,0.0));
  
  long double omega = 2*pi/PerSin;
  long double omega_sq = pow(omega,2);
  long double omega_sq_diff = k - pow(omega,2);
  long double omega_sq_diff_sq = pow(omega_sq_diff,2);
  
  // Solution part
  NumericVector circ_sol(Nstep);
  NumericVector SW_response(Nstep);
  NumericVector trans_sol(Nstep);
  
  
  // always use time as e.g. 0 and 0.1
  long double t0 = 0.0L;
  long double tend = t0 + h;
  // problems with numerical errors for large t - 
  // should look at taking care of gamma_tilde in combination with gamma.
  long double gamma_tilde = expl(-gamma * t0 / 2);
  
  std::complex<long double> C1 = cos(omega1 * t0);
  std::complex<long double> S1 = sin(omega1 * t0);
  long double C = cos(omega*t0);
  long double S = sin(omega*t0);
  std::complex<long double> C1_end = cos(omega1*tend);
  std::complex<long double> S1_end = sin(omega1*tend);
  long double C_end = cos(omega*tend);
  long double S_end = sin(omega*tend);
  std::complex<long double> a21 = ( (- gamma * C1 / 2.0L) - omega1 * S1);
  std::complex<long double> a22 = ( (- gamma * S1 / 2.0L) + omega1 * C1);
  
  for (int i = 0; i < Nstep; i ++ ){
    
    long double teval = TimeReturn[i];
    long double angle = (PhiSin+(Twopi/(PerSin))*(teval)); // new angle given time
    long double PhiSin_new =  (angle - Twopi * floor(angle/(Twopi))); // modulo
    long double Fsw = force[i]; // force at time i
    
    // SOLVE
    long double force_A = AmpSin * cos(PhiSin_new);
    long double force_B = AmpSin * sin(PhiSin_new);
    
    // Response sinusoidal forcing
    long double D = ( - omega * gamma * force_A + omega_sq_diff * force_B ) / ( omega_sq_diff_sq + omega_sq*gamma_sq );
    long double E = (   omega * gamma * force_B + omega_sq_diff * force_A ) / ( omega_sq_diff_sq + omega_sq*gamma_sq );  
    
    // Coefficient of transient term
    std::complex<long double> Atilde = ( 1.0L / (omega1 * gamma_tilde) ) *
      ( a22 * alpha - S1 * beta +
      - a22 * Fsw / (omega0*omega0) +
      D * ( a22 * (-C) - omega * S1 * S ) +
      E * ( a22 * (-S) + omega * S1 * C ));
    
    std::complex<long double> Btilde  = ( 1.0L / (omega1 * gamma_tilde) ) *
      ( -a21 * alpha + C1 * beta +
      + a21 * Fsw / (omega0*omega0) +
      D * ( a21 * C + omega * C1 * S ) +
      E * ( a21 * S - omega * C1 * C ));
    
    std::complex<long double>  alpha_new = expl( - gamma * tend / 2.0L ) * (Atilde * C1_end + Btilde * S1_end) +
      Fsw / (omega0*omega0) + D * C_end + E * S_end ;
    
    std::complex<long double>  beta_new = (-gamma / 2.0L) * expl( - gamma * tend / 2.0L ) * (Atilde * C1_end + Btilde * S1_end) +
      expl( - gamma * tend / 2.0L ) * ( - omega1 * Atilde * S1_end + omega1 * Btilde * C1_end) +
      - omega * D * S_end + omega * E * C_end ;
      
    
      // insert new alpha and beta
      alpha = alpha_new.real();
      beta = beta_new.real();
      Valpha[i+1] = alpha;
      Vbeta[i+1] = beta;
      TimeReturn[i+1] = time[i];
      
      // Sinusoidal response
      circ_sol[i] = D*cos(omega * tend) + E*sin(omega * tend);
      // Sleep-wake response
      SW_response[i] = Fsw/k;
      // transient response
      trans_sol[i]  = (expl(-gamma*(tend)/2.0L) * (Atilde*cos(omega1*tend) + Btilde*sin(omega1*tend))).real();
      
      // total response = circ_sol + SW_response + trans_sol
  }
  
  return List::create(
    _["time"] = TimeReturn,
    _["y1"] = Valpha,
    _["y2"] = Vbeta,
    _["circ_sol"] = circ_sol,
    _["SW_response"] = SW_response,
    _["trans_sol"] = trans_sol
  );
  
}


// Original code translated from matlab 
// with large t and large gamma, we have rapidly troubles caused by numerical errors. 
// changed from double to long double

// List Solve_SDDHO(NumericVector y, NumericVector time, NumericVector force, long double gamma, long double k,long double AmpSin, long double PhiSin, long double PerSin){
//   
//   // define pi
//   const double pi = 3.14159265358979323846264338327950288419716939937510L;
//   long double Twopi = 2.0L * pi;
//   
//   // number of step
//   int Nstep = time.length();
//   
//   // Vector of position
//   NumericVector Valpha(Nstep+1);
//   
//   // Vector of speed
//   NumericVector Vbeta(Nstep+1);
//   
//   // Compute Iteration step (h)
//   long double h  = time[1]-time[0];
//   
//   // Vector of time that will be return
//   NumericVector TimeReturn(Nstep+1);
//   
//   // position and speed value
//   long double alpha = y[0];
//   long double beta = y[1];
//   
//   // Insert initial value
//   Valpha[0] = alpha;
//   Vbeta[0] = beta;
//   TimeReturn[0] = time[0]-h;
//   
//   // SOLVE
//   long double force_A = AmpSin * cos(PhiSin);
//   long double force_B = AmpSin * sin(PhiSin);
//   
//   //transcient reponse frequency
//   long double omega0 = sqrtl(k);
//   long double gamma_sq = pow(gamma,2);
//   std::complex<long double> omega1 = std::sqrt(std::complex<long double>(k - gamma_sq/4,0.0));
//   
//   // Response sinusoidal forcing
//   long double omega = 2*pi/PerSin;
//   long double omega_sq = pow(omega,2);
//   long double omega_sq_diff = k - pow(omega,2);
//   long double omega_sq_diff_sq = pow(omega_sq_diff,2);
//   
//   long double D = ( - omega * gamma * force_A + omega_sq_diff * force_B ) / ( omega_sq_diff_sq + omega_sq*gamma_sq );
//   long double E = (   omega * gamma * force_B + omega_sq_diff * force_A ) / ( omega_sq_diff_sq + omega_sq*gamma_sq );
//   
//   
//   // Solution part
//   NumericVector circ_sol(Nstep);
//   NumericVector SW_response(Nstep);
//   NumericVector trans_sol(Nstep);
//   
//   
//   for (int i = 0; i < Nstep; i ++ ){
//     
//     long double t0 = TimeReturn[i];
//     long double tend = time[i];
//     long double Fsw = force[i];
//     
//     long double gamma_tilde = expl(-gamma * t0 / 2);
//     std::complex<long double> C1 = cos(omega1 * t0);
//     std::complex<long double> S1 = sin(omega1 * t0);
//     long double C = cos(omega*t0);
//     long double S = sin(omega*t0);
//     std::complex<long double> C1_end = cos(omega1*tend);
//     std::complex<long double> S1_end = sin(omega1*tend);
//     long double C_end = cos(omega*tend);
//     long double S_end = sin(omega*tend);
//     
//     
//     std::complex<long double> a21 = ( (- gamma * C1 / 2.0L) - omega1 * S1);
//     std::complex<long double> a22 = ( (- gamma * S1 / 2.0L) + omega1 * C1);
//     
//     std::complex<long double> Atilde = ( 1.0L / (omega1 * gamma_tilde) ) *
//       ( a22 * alpha - S1 * beta +
//       - a22 * Fsw / (omega0*omega0) +
//       D * ( a22 * (-C) - omega * S1 * S ) +
//       E * ( a22 * (-S) + omega * S1 * C ));
//     
//     std::complex<long double> Btilde  = ( 1.0L / (omega1 * gamma_tilde) ) *
//       ( -a21 * alpha + C1 * beta +
//       + a21 * Fsw / (omega0*omega0) +
//       D * ( a21 * C + omega * C1 * S ) +
//       E * ( a21 * S - omega * C1 * C ));
//     
//     std::complex<long double>  alpha_new = expl( - gamma * tend / 2.0L ) * (Atilde * C1_end + Btilde * S1_end) +
//       Fsw / (omega0*omega0) + D * C_end + E * S_end ;
//     
//     std::complex<long double>  beta_new = (-gamma / 2.0L) * expl( - gamma * tend / 2.0L ) * (Atilde * C1_end + Btilde * S1_end) +
//       expl( - gamma * tend / 2.0L ) * ( - omega1 * Atilde * S1_end + omega1 * Btilde * C1_end) +
//       - omega * D * S_end + omega * E * C_end ;
//       
//       // insert new alpha and beta
//       alpha = alpha_new.real();
//       beta = beta_new.real();
//       Valpha[i+1] = alpha;
//       Vbeta[i+1] = beta;
//       TimeReturn[i+1] = time[i];
//       
//       circ_sol[i] = D*cos(omega * tend) + E*sin(omega * tend);
//       SW_response[i] = Fsw/k;
//       trans_sol[i]  = (expl(-gamma*(tend)/2.0L) * (Atilde*cos(omega1*tend) + Btilde*sin(omega1*tend))).real();
//   }
//   
//   return List::create(
//     _["time"] = TimeReturn,
//     _["y1"] = Valpha,
//     _["y2"] = Vbeta,
//     _["circ_sol"] = circ_sol,
//     _["SW_response"] = SW_response,
//     _["trans_sol"] = trans_sol
//   );
//   
// }


