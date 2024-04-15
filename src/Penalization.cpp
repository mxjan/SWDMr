#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
double PenalizationUnstableFit(Rcpp::NumericVector y1,Rcpp::NumericVector Time, double weight, Rcpp::NumericVector TimeInterval, int DayStabilization) {
  
  int Nstep = Time.length();
  IntegerVector idx = seq(0,Nstep-1); 
  
  int h = ceil(24/(Time[3]-Time[2]));
  
  IntegerVector IdxTimeInInterval(h);
  
  int count = 0;
  
  // Get index of element in intervals
  for (int i = 0; i < Nstep; i ++ ){
    if (Time[i] > TimeInterval[0] && Time[i] <= TimeInterval[1]){
      IdxTimeInInterval[count] = i ; 
      count = count + 1;
    }
  }
  
  // Find maxmimum and minimum value in our Time Interval
  NumericVector y1_TimeInterval = y1[IdxTimeInInterval];
  IntegerVector idx_TimeInterval = idx[IdxTimeInInterval];
  
  double y1_max = max(y1_TimeInterval); // maximal value in the interval
  double y1_min = min(y1_TimeInterval); // minimal value in the interval

  int idx_y1_max = idx_TimeInterval[which_max(y1_TimeInterval)]; // index of the maximal value in the interval
  int idx_y1_min = idx_TimeInterval[which_min(y1_TimeInterval)]; // index of the minimal value in the interval
  
  count = 0;
  
  NumericVector y1_Max_Interval(DayStabilization);
  NumericVector y1_Min_Interval(DayStabilization);
  
  for (int i = 1; i < (DayStabilization + 1); i ++ ){
    y1_Max_Interval[count] = y1[idx_y1_max - (h*i)]; // values 24h before maximal value
    y1_Min_Interval[count] = y1[idx_y1_min - (h*i)]; // values 24h before minimal value
    count = count + 1;
  }
  
  double RSS = sum(pow(y1_Min_Interval - y1_min,2) + pow(y1_Max_Interval - y1_max,2)) * weight;
  
  return RSS ; 
  
}