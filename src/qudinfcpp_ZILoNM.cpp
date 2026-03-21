#include <Rcpp.h>
#include <iostream>
#include <math.h>
using namespace Rcpp;

// Sensitivity analysis variant: false-zero probability mechanism exp(-eta^2 * sqrt(M))
// Changes from original (qudinfcpp_ZILoNM.cpp):
//   1. g1_original: -eta*eta*m  ->  -eta*eta*sqrt(m)
//   2. m1 (mode of g1): no closed form, solved via Newton's method
//   3. m3_1, m4_1 (tail bounds from g1): no closed form, solved via Newton's method
//      with bracketing for m4_1 and a safety check for m3_1

double g1_original(double m, double x, double y, double mu, double sig,
                   double ypart, double beta1, double beta5, double delta, double eta){
  // CHANGED: -eta*eta*sqrt(m) instead of -eta*eta*m
  return - pow(y-beta1*m-beta5*x*m-ypart,2)/(2*delta*delta) - eta*eta*sqrt(m);}

double g2_original(double m, double mu, double sig){
  return -log(m) - pow(log(m) - mu,2)/(2*sig*sig);}

double g_original(double m, double x, double y, double mu, double sig,
                  double ypart, double beta1, double beta5, double delta, double eta){
  return  g1_original(m,x,y,mu,sig,ypart,beta1,beta5,delta,eta)+g2_original(m,mu,sig);}

double f_original(double m, double x, double y, double mu, double sig,
                  double ypart, double beta1, double beta5, double delta, double eta, double g_star){
  return exp(g_original(m,x,y,mu,sig,ypart,beta1,beta5,delta,eta)- g_star) ;}

double fquadinf(double yy, double xa, double xb,
                double x, double y, double mu, double sig,
                double ypart, double beta1, double beta5, double delta, double eta, double g_star) {
  double u = -1;
  double v = 1 ;
  double duv = (xb-xa)/(v-u);

  return duv * f_original(xa + duv*(yy-u), x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star);
}

double quadinfcpp(double xa, double xb, List nodes, List weights,
                  double x, double y, double mu, double sig,
                  double ypart, double beta1, double beta5, double delta, double eta, double g_star) {
  double Q;
  double h=0.5;
  NumericVector n=nodes[0];
  NumericVector w=weights[0];
  double s= w[6] * fquadinf(n[6],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star);
  for (int j=0;j<6;j++) {
    s=s + w[j] * (fquadinf(n[j],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star) +
      fquadinf(-n[j],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star) );
  }
  Q=s*h;

  double delta1;
  double newQ;
  double tol=1e-12;
  for (int k=1;k<7;k++) {
    n=nodes[k];
    w=weights[k];
    s=0;
    for (int j=0;j<w.size();j++) {
      s=s+w[j]*(fquadinf(n[j],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star) +
        fquadinf(-n[j],xa, xb,x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star) );
    }
    h = h/2;
    newQ=s*h + Q/2.0;
    delta1 = std::abs(newQ-Q);
    Q=newQ;
    if (delta1<tol) {
      break;
    }
  }

  return Q;}


double loghicpp(double x, double y, double mu, double sig,
                double ypart, double beta1, double beta5, double delta, double eta,
                List nodes, List weights, double C){

  double b = beta1 + beta5*x;

  // CHANGED: mode of g1 for sqrt(M) mechanism via Newton's method
  // dg1/dm  = b*(y-ypart-b*m)/delta^2 - eta^2/(2*sqrt(m))
  // d2g1/dm2 = -b^2/delta^2 + eta^2/(4*m^1.5)
  // Start at m=1.0 (safer than OLS guess for large eta)
  double m1 = 1.0;
  for (int iter=0; iter<100; iter++){
    double f  =  b*(y-ypart-b*m1)/delta/delta - eta*eta/(2*sqrt(m1));
    double fp = -b*b/delta/delta + eta*eta/(4*pow(m1, 1.5));
    double step = f/fp;
    m1 -= step;
    if (m1 <= 0) { m1 = 1e-10; break; }
    if (std::abs(step) < 1e-12) break;
  }
  m1 = std::max(m1, 0.0);

  double m2 = exp(-sig*sig + mu);
  double g2_star = g2_original(m2,mu,sig);
  double start=0;
  double lower=0, upper=0;
  double m3_1=0, m3_2=0, m4_1=0, m4_2=0, g1_star=0,g_star=0,val=0,by=0;
  NumericVector gvalue(1000);

  // situation 1
  if((m1 <= C) and (m2 <= C)){
    if(m1 == 0){start=1e-200;
    }else{start=m1;}
    val=std::min(start,m2);
    by = fabs(m2 - start)/999.0;
    gvalue[0] = g_original(val,x,y,mu,sig,ypart,beta1,beta5,delta,eta) ;
    for (double l=1;l<1000;l++){
      val = val+by;
      gvalue[l] = g_original(val,x,y,mu,sig,ypart,beta1,beta5,delta,eta) ;}
    g_star = max(gvalue);

    double g2_m1=g2_original(m1,mu,sig);
    double g1_m2=g1_original(m2,x,y,mu,sig,ypart,beta1,beta5,delta,eta);

    // CHANGED: m3_1 and m4_1 for sqrt(M) mechanism via Newton's method
    // Solve g1(m) = g1_star - 30, i.e. f(m) = g1(m) - target = 0
    // f'(m) = dg1/dm = b*(y-ypart-b*m)/delta^2 - eta^2/(2*sqrt(m))
    g1_star = g1_original(m1,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
    double target = g1_star - 30;

    // m3_1: left of mode
    // g1(0+) = -(y-ypart)^2/(2*delta^2) [finite], check if left root exists
    double g1_near0 = g1_original(1e-10,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
    if (g1_near0 >= target){
      // g1(0+) > target: function never drops below target on left side
      m3_1 = 0.0;
    } else {
      // genuine left root exists: Newton from 1e-10
      double ml = 1e-10;
      for (int iter=0; iter<100; iter++){
        double f  = g1_original(ml,x,y,mu,sig,ypart,beta1,beta5,delta,eta) - target;
        double fp = b*(y-ypart-b*ml)/delta/delta - eta*eta/(2*sqrt(ml));
        double step = f/fp;
        ml -= step;
        if (ml <= 0) ml = 1e-10;
        if (std::abs(step) < 1e-12) break;
      }
      m3_1 = std::max(ml, 0.0);
    }

    // m4_1: right of mode
    // bracket first: double mr until g1(mr) < target
    double mr = m1*2 + 1.0;
    while (g1_original(mr,x,y,mu,sig,ypart,beta1,beta5,delta,eta) > target) mr *= 2;
    for (int iter=0; iter<100; iter++){
      double f  = g1_original(mr,x,y,mu,sig,ypart,beta1,beta5,delta,eta) - target;
      double fp = b*(y-ypart-b*mr)/delta/delta - eta*eta/(2*sqrt(mr));
      double step = f/fp;
      mr -= step;
      if (mr <= 0) mr = 1e-10;
      if (std::abs(step) < 1e-12) break;
    }
    m4_1 = mr;

    // m3_2, m4_2 unchanged: derived from g2 only
    NumericVector expo30 = NumericVector::create(-sqrt(2*sig*sig *(30-g_star+g1_star-mu+sig*sig/2.0)) -sig*sig + mu ,
                                                 sqrt(2*sig*sig *(30-g_star+g1_star-mu+sig*sig/2.0)) -sig*sig + mu );
    m3_2 = exp(min(expo30));
    m4_2 = exp(max(expo30));

    // m1 < m2
    if (m1 < m2){
      // lower bound
      if (m1==0){
        lower = m3_2 ;
      }else if( g2_m1+g1_star-g_star < -30){ //check m1=m3
        lower = m1;
      }else{
        lower = std::max(m3_1,m3_2);
      }
      // upper bound
      if (g1_m2+g2_star-g_star < -30){//check m2=m4
        upper = m2;
      }else{
        upper = std::min(m4_1,m4_2);
      }
    }else{// m1 > m2
      // lower bound
      if(g1_m2+g2_star-g_star < -30){ //check m2=m4
        lower = m2;
      }else{
        lower = std::max(m3_1,m3_2);
      }
      // upper bound
      if (g2_m1+g1_star-g_star < -30){//check m1=m3
        upper= m1;
      }else{
        upper = std::min(m4_1,m4_2);
      }
    }
    upper = std::min(upper,C);
    // situation 2
  }else if((m1>C) != (m2>C)){
    start = std::min(m1,m2);
    if(start == 0){start=1e-200;}
    val=start;
    by = fabs(C - start)/999.0;
    gvalue[0] = g_original(val,x,y,mu,sig,ypart,beta1,beta5,delta,eta) ;
    for (double l=1;l<1000;l++){
      val = val+by;
      gvalue[l] = g_original(val,x,y,mu,sig,ypart,beta1,beta5,delta,eta) ;}
    g_star = max(gvalue);

    upper = C;
    if(m1<C){g1_star=g1_original(m1,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
    }else{g1_star=g1_original(C,x,y,mu,sig,ypart,beta1,beta5,delta,eta);}
    lower = exp(-sqrt(2*sig*sig *(30-g_star+g1_star-mu+sig*sig/2.0)) -sig*sig + mu);
    // situation 3
  }else{
    g_star = g_original(C,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
    g1_star = g1_original(C,x,y,mu,sig,ypart,beta1,beta5,delta,eta);
    upper = C;
    lower = exp(-sqrt(2*sig*sig *(30-g_star+g1_star-mu+sig*sig/2.0)) -sig*sig + mu);
  }

  double integ =  quadinfcpp(lower,upper,nodes,weights,
                             x,y,mu,sig,ypart,beta1,beta5,delta,eta,g_star);
  double output = log(integ) + g_star;

  return output;
}

// [[Rcpp::export]]
NumericMatrix loghicpp_all(NumericVector X_group2, NumericVector Y_group2,
                           NumericMatrix mu, double sig,
                           NumericVector Ypart, double beta1, double beta5, double delta, double eta,
                           List nodes, List weights, double C){
  NumericMatrix loghik = clone(mu);

  for(int k=0;k < mu.ncol();k++){
    for(int i=0;i < mu.nrow();i++){
      loghik(i,k) = loghicpp(X_group2[i], Y_group2[i], mu(i,k), sig, Ypart[i], beta1, beta5, delta, eta, nodes, weights,C);
    }
  }

  return loghik;
}
