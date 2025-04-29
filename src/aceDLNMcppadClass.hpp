#ifndef ACEDLNMCPPADCLASS_HPP
#define ACEDLNMCPPADCLASS_HPP


// include CppAD
#include <cppad/cppad.hpp>

// autodiff include
// from https://github.com/awstringer1/varcomptest/blob/main/src/reml-ad.cpp
// #include "autodiff/common/meta.hpp"
// #include "autodiff/common/numbertraits.hpp"
// #include "autodiff/common/binomialcoefficient.hpp"
// #include "autodiff/common/vectortraits.hpp"
// #include "autodiff/forward/dual/dual.hpp"
// #include "autodiff/common/eigen.hpp"
// #include "autodiff/common/classtraits.hpp"
// #include "autodiff/forward/utils/derivative.hpp"
// #include "autodiff/forward/real/real.hpp"
// #include "autodiff/forward/utils/taylorseries.hpp"
// #include "autodiff/forward/real.hpp"
// #include "autodiff/forward/real/eigen.hpp"
// #include "autodiff/forward/utils/gradient.hpp"
// using namespace autodiff;



#include "defheader.h"


#include <cmath>

#include <iostream>
using namespace std;



Eigen::VectorXd convertToDouble(const Vec input) {
  int n = input.size();
  Eigen::VectorXd output(n);

  for (int i = 0; i < n; ++i) {
      output(i) = CppAD::Value(input(i));
  }
  return output;
}


Eigen::MatrixXd convertToDouble(const Mat input) {
  int nrow = input.rows();
  int ncol = input.cols();
  Eigen::MatrixXd output(nrow, ncol);

  for (int i = 0; i < nrow; ++i) {
      for (int j = 0; j < ncol; ++j) {
          output(i, j) = CppAD::Value(input(i, j));
      }
  }

  return output;
}


// ************** PART 1: Define some functions **********************


// TODO: the bspline function evaluated at the points outside the boundaries are incorret!
// Bspline(l=0) = 0,0,0,0... It should be a linear function of l, not always equal to 0.
int knotindex(Scalar x,const Vec t) {
  int q = t.size();
  int k=0;
  if (x < t(0)) return -1;
  while(x>=t(k)){
    k++;
    if (k >= q) break;
  }

  return k-1;
}

Scalar weight(Scalar x, const Vec& t,int i,int k) {
  if (t(i+k-1) != t(i-1))
    return((x - t(i-1))/(t(i+k-1)-t(i-1)));
  return 0.;
}


Scalar Bspline(Scalar x, int j, const Vec& t,int p) {
  // Evaluate the jth B-spline
  // B_p(x) of order p (degree p-1) at x
  if (p==1)
    return(x>=t(j-1) && x<t(j+1-1));

  Scalar w1 = weight(x,t,j,p-1);
  Scalar w2 = weight(x,t,j+1,p-1);
  Scalar b1 = Bspline(x,j,t,p-1);
  Scalar b2 = Bspline(x,j+1,t,p-1);

  return w1*b1 + (1.-w2)*b2;
}

Scalar Bspline1st(Scalar x, int j, const Vec& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  Scalar bb1 = Bspline(x,j+1,t,p-1);
  Scalar bb2 = Bspline(x,j,t,p-1);

  Scalar ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  Scalar ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

Scalar Bspline2nd(Scalar x, int j, const Vec& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  Scalar bb1 = Bspline1st(x,j+1,t,p-1);
  Scalar bb2 = Bspline1st(x,j,t,p-1);

  Scalar ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  Scalar ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

Vec Bsplinevec(Scalar x, const Vec& t,int p) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  // int k = knotindex(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return b;
}

Vec BsplinevecCon(Scalar x, const Vec& t, int p, Mat Z) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Vec Bsplinevec1st(Scalar x, const Vec& t,int p) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  // int k = knotindex(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return b;
}

Vec BsplinevecCon1st(Scalar x, const Vec& t, int p, Mat Z) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Vec Bsplinevec2nd(Scalar x, const Vec& t,int p) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return b;
}


Vec BsplinevecCon2nd(Scalar x, const Vec& t,int p, Mat Z) {
  int m = t.size() - p;
  Vec b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}



// Lanczos approximation
// Source code from https://github.com/brianmartens/BetaFunction/blob/master/BetaFunction/bmath.h
Scalar lanczos_lgamma(Scalar z) {
    const Scalar LG_g = 7.0;
    const int LG_N = 9;

    const Scalar ln_sqrt_2_pi = 0.91893853320467274178;

    Vec lct(LG_N+1);
    lct << 0.9999999999998099322768470047347,
    676.520368121885098567009190444019,
   -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
    -176.61502916214059906584551354,
     12.507343278686904814458936853,
    -0.13857109526572011689554707,
    9.984369578019570859563e-6,
    1.50563273514931155834e-7;

    Scalar sum;
    Scalar base;

    // To avoid if condition z < 0.5, we calculate gamma(z+1) which is equal to z*gamma(z).
    // WAS:
    // z = z - 1.0;
    // if (z < 0.5) {
    //   // Use Euler's reflection formula:
    //   // Gamma(z) = Pi / [Sin[Pi*z] * Gamma[1-z]];
    //   out = log(g_pi / sin(g_pi * z)) - lanczos_lgamma(1.0 - z);
    //   return out;
    // }
    // gamma(z) ...

    // New: indeed gamma(z+1)
    base = z + LG_g + 0.5;  // Base of the Lanczos exponential
    sum = 0;
    // We start with the terms that have the smallest coefficients and largest
    // denominator.
    for(int i=LG_N; i>=1; i--) {
      sum += lct[i] / (z + ((double) i));
    }
    sum += lct[0];
    Scalar gammazplus1 = ((ln_sqrt_2_pi + log(sum)) - base) + log(base)*(z+0.5);
    Scalar out = gammazplus1 - log(z);

    return out;
}

// 1st derivative of log gamma function
// https://math.stackexchange.com/questions/481253/differentiate-log-gamma-function
// VERY SLOW!!!
// Scalar lgamma1st (Scalar z) {
//   // Euler's constant https://en.wikipedia.org/wiki/Euler%27s_constant#
//   const Scalar Euler = 0.57721566490153286060651209008240243104215933593992;
//   Scalar out = -1.0 * Euler;
//   int K = 1e6;
//   for (int i = 1; i < K; i++) {
//     out += 1.0/i - 1.0/(i+z-1.0);
//   }
//   std::cout << "z" << CppAD::Value(z) << std::endl;
//   std::cout << "lgamma1st" << CppAD::Value(out) << std::endl;
//   return out;
// }



// https://github.com/tminka/lightspeed/blob/master/digamma.m
Scalar lgamma1st (Scalar x) {
  const Scalar pi = 3.141592653589793238462643383279;
  const Scalar large = 9.5;
  const Scalar d1 = -0.5772156649015328606065121;
  const Scalar d2 = pi*pi/6.0;
  const Scalar small = 1e-6;
  const Scalar s3 = 1.0/12.0;
  const Scalar s4 = 1.0/120.0;
  const Scalar s5 = 1.0/252.0;
  const Scalar s6 = 1.0/240.0;
  const Scalar s7 = 1.0/132.0;
  const Scalar s8 = 691.0/32760.0;
  const Scalar s9 = 1.0/12.0;
  const Scalar s10 = 3617.0/8160.0;

  // Use de Moivre's expansion if x >= large = 9.5
  // calculate lgamma1st(x+10)
  Scalar xplus10 = x + 10.0;
  Scalar y = 0.0;
  Scalar r = 1.0 / xplus10;
  y += log(xplus10) - 0.5 * r;
  r = r * r;
  y = y - r * ( s3 - r * ( s4 - r * (s5 - r * (s6 - r * s7))));

  // lgamma1st(x+10) = (1/x + 1/(x+1) + ... + 1/(x+9)) + lgamma1st(x)
  y = y - 1.0/x - 1.0/(x+1.0) - 1.0/(x+2.0) - 1.0/(x+3.0) - 1.0/(x+4.0) - 1.0/(x+5.0) - 1.0/(x+6.0) - 1.0/(x+7.0) - 1.0/(x+8.0) - 1/(x+9);

  return y;
}



// https://github.com/tminka/lightspeed/blob/master/trigamma.m

// Scalar lgamma2nd (Scalar x) {
//   const Scalar pi = 3.141592653589793238462643383279;
//   const Scalar c = pi*pi/6;
//   const Scalar c1 = -2.404113806319188570799476;
//   const Scalar b2 =  1.0/6.0;
//   const Scalar b4 = -1.0/30.0;
//   const Scalar b6 =  1.0/42.0;
//   const Scalar b8 = -1.0/30.0;
//   const Scalar b10 = 5.0/66.0;

//   // TO DO: % Reduce to trigamma(x+n) where ( X + N ) >= large.

//   Scalar z = 1./(x*x);
//   Scalar y = 0.5*z + (1.0 + z*(b2 + z*(b4 + z*(b6 + z*(b8 + z*b10))))) / x;

//   // std::cout << "x" << CppAD::Value(x) << std::endl;
//   // std::cout << "trigamma" << CppAD::Value(y) << std::endl;
//   return y;
// }

// **************** PART 2: g(mu) = DL term + linear term + smooth term *************************
class Modelcppad {
  // The DLNM model

private:
  // DATA
  const Vec y; // Response
  const Mat Sw_large; // penalty matrix for w(l)
  const std::vector<Mat> Sf_list; // penalty matrix for f(E)
  const std::vector<Mat> B_inner_list;
  const std::vector<Vec> knots_f_list; // knots for f(E) B-spline
  const Mat Dw;  // \int w(l)^2 dl = 1

  const Mat Xfix; // fixed effects
  const Mat Xrand; // random effects
  const Vec r; // rank of each smooth
  const std::vector<Mat> Zf_list;

  const Vec Xoffset; // offset
  
public:
  const Mat K;
  const Vec a; 
  
  int n;
  int kw;
  int kwp; 
  int kwopt;
  int kwpopt;
  int kE;
  int kEp;
  int kbetaR;
  int kbetaF;
  int p; // number of smooth terms in Xrand
  int M;

  // PARAMETERS
  Vec alpha_f;
  Vec phi;
  Vec phiKa;
  Scalar log_theta;
  Vec log_smoothing_f;
  Vec log_smoothing_w;

  Vec betaF; // parameters for fixed effects
  Vec betaR; // parameters for random effects
  Vec logsmoothing; // log smoothing parameters for random effects

  // Components generated
  Scalar theta;
  Vec smoothing_f;
  Vec smoothing_w;
  Vec smoothing;
  Vec phi_long;
  Scalar alpha_w_C_denominator;
  Vec alpha_w_C;
  Mat Bf_matrix;
  Vec Bf;
  Mat E;
  Vec eta;
  Vec eta_remaining; // remaining terms = Xfix * betaF + Xrand * betaR
  Vec mu; // log(mu) = eta + eta_remaining + Xoffset


  // Components for derivatives

  Mat dlogmu_dw_mat;
  Vec dlogdensity_dmu_vec;
  Mat dmu_dw_mat;
  Mat dw_dphi_mat;

  Vec d2logdensity_dmudmu_vec;
  std::vector<Mat> d2w_dphidphi_list;
  Mat he_alpha_w_mat;
  Mat he_alpha_f_alpha_w_mat;


  // gradient and hessian for updating alpha_f, betaR and betaF
  Mat he_inner_mat;



  // full hessian
  Mat he_alpha_f_mat;
  Mat he_betaR_mat;
  Mat he_betaF_mat;
  Mat he_phi_mat;
  Mat he_alpha_f_phi_mat;
  Mat he_alpha_f_betaF_mat;
  Mat he_alpha_f_betaR_mat;
  Mat he_phi_betaF_mat;
  Mat he_phi_betaR_mat;
  Mat he_betaR_betaF_mat;

  Mat he_s_u_mat;


  // AD tape for LAML
  // CppAD::ADFun<double> gr;
  bool ifhastape = false;

  // Constructor
  Modelcppad(const Vec& y_,
             const std::vector<Mat>& B_inner_list_,
             const std::vector<Vec>& knots_f_list_,
             const Mat& Sw_large_,
             const std::vector<Mat>& Sf_list_,
             const Mat& Dw_,
             const Mat& Xrand_,
             const Mat& Xfix_,
             const std::vector<Mat>& Zf_list_,
             const Vec& Xoffset_,
             const Vec& r_,
             const Mat& K_,
             const Vec& a_,
             Vec& alpha_f_,
             Vec& phi_,
             Scalar log_theta_,
             Vec& log_smoothing_f_,
             Vec& log_smoothing_w_,
             Vec& betaR_,
             Vec& betaF_,
             Vec& logsmoothing_) :
    y(y_), B_inner_list(B_inner_list_), knots_f_list(knots_f_list_), Sw_large(Sw_large_), Sf_list(Sf_list_), Dw(Dw_), Xrand(Xrand_), Xfix(Xfix_), Zf_list(Zf_list_), Xoffset(Xoffset_), r(r_),
    K(K_), a(a_),
    alpha_f(alpha_f_), phi(phi_), log_theta(log_theta_), log_smoothing_f(log_smoothing_f_), log_smoothing_w(log_smoothing_w_), betaR(betaR_), betaF(betaF_), logsmoothing(logsmoothing_) {

      n = y.size(); // sample size
      
      M = B_inner_list.size(); 

      kw = a.size() + M;
      kwopt = K.cols(); // conL = TRUE: kwopt = kw-(2*M); conL = FALSE: kwopt = kw-(1*M)

      kwp = Sw_large.cols();
      kwpopt = kwopt/M; 

      kE = alpha_f.size();
      kEp = Sf_list.at(0).cols();


      kbetaR = betaR.size();
      kbetaF = betaF.size();
      p = r.size();

      theta = exp(log_theta);
      
      smoothing_f.resize(M);
      for (int i = 0; i < M; i++) smoothing_f(i) = exp(log_smoothing_f(i));
      
      smoothing_w.resize(M);
      for (int i = 0; i < M; i++) smoothing_w(i) = exp(log_smoothing_w(i));


      smoothing.resize(p);
      for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));

      phiKa = K * phi + a;
      
      phi_long.resize(kw); // phi_long = c(1, phiKa) for each kwp
      for (int i = 0; i < M; i++) {
        phi_long(i*kwp) = 1.0;
        phi_long.segment(i*kwp + 1, kwp-1) = phiKa.segment(i*(kwp-1), (kwp-1));
      }
      

      alpha_w_C.resize(kw);

      for (int i = 0; i < M; i++) {
        alpha_w_C_denominator = sqrt(phi_long.segment(i*kwp, kwp).dot(Dw * phi_long.segment(i*kwp, kwp)));
        alpha_w_C.segment(i*kwp, kwp) = phi_long.segment(i*kwp, kwp) / alpha_w_C_denominator;
      }

      
      E.resize(n, M);
      for (int i = 0; i < M; i++) {
        E.col(i) = B_inner_list.at(i) * alpha_w_C.segment(i*kwp, kwp);
      }
      

      Bf_matrix.resize(n, kE);
      eta.resize(n);
      eta_remaining.resize(n);
      mu.resize(n);
      
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < M; j++) {
          Bf = BsplinevecCon(E(i,j), knots_f_list.at(j), 4, Zf_list.at(j));
          Bf_matrix.block(i, j*kEp, 1, kEp) = Bf.transpose();
        }
        eta(i) = Bf_matrix.row(i).dot(alpha_f);
        eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
        mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
      }

      // Initialize the derivative components and NegativeLogLikelihood
      dw_dphi_mat = dw_dphi(); // d alpha_w / d phi
      d2w_dphidphi_list = d2w_dphidphi(); // d^2 alpha_w / d phi d phi
      
      he_s_u_mat.resize(kwopt+kE+kbetaR+kbetaF, kwopt+kE+kbetaR+kbetaF);
      
      
      derivative_coef(); // comment out for debug
      derivative_he(); // comment out for debug
    }

  // Functions to set parameters
  void setAlphaF(const Vec alpha_f_) {
    alpha_f = alpha_f_;
  }

  void setPhi(const Vec phi_) {
    phi = phi_;
  }
  void setBetaF(const Vec betaF_) {
    betaF = betaF_;
  }
  void setBetaR(const Vec betaR_) {
    betaR = betaR_;
  }
  void setLogTheta(const Scalar log_theta_) {
    log_theta = log_theta_;
    theta = exp(log_theta);
  }

  void setLogSmoothingF(const Vec log_smoothing_f_) {
    log_smoothing_f = log_smoothing_f_;
    for (int i = 0; i < M; i++) smoothing_f(i) = exp(log_smoothing_f(i));
  }

  void setLogSmoothingW(const Vec log_smoothing_w_) {
    log_smoothing_w = log_smoothing_w_;
    for (int i = 0; i < M; i++) smoothing_w(i) = exp(log_smoothing_w(i));
  }
  void setLogsmoothing(const Vec logsmoothing_) { // log smoothing parameters for remaining terms
    logsmoothing = logsmoothing_;
    for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));
  }


  void derivative_coef() {
    // regenerate
    
    phiKa = K * phi + a;

    phi_long.resize(kw); // phi_long = c(1, phiKa) for each kwp
    for (int i = 0; i < M; i++) {
      phi_long(i*kwp) = 1.0;
      phi_long.segment(i*kwp + 1, kwp-1) = phiKa.segment(i*(kwp-1), (kwp-1));
    }
    
    alpha_w_C.resize(kw);

    for (int i = 0; i < M; i++) {
      alpha_w_C_denominator = sqrt(phi_long.segment(i*kwp, kwp).dot(Dw * phi_long.segment(i*kwp, kwp)));
      alpha_w_C.segment(i*kwp, kwp) = phi_long.segment(i*kwp, kwp) / alpha_w_C_denominator;
    }
    
    E.resize(n, M);
    for (int i = 0; i < M; i++) {
      E.col(i) = B_inner_list.at(i) * alpha_w_C.segment(i*kwp, kwp);
    }
    
    Bf_matrix.resize(n, kE);
    eta.resize(n);
    eta_remaining.resize(n);
    mu.resize(n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < M; j++) {
        Bf = BsplinevecCon(E(i,j), knots_f_list.at(j), 4, Zf_list.at(j));
        Bf_matrix.block(i, j*kEp, 1, kEp) = Bf.transpose();
      }
      eta(i) = Bf_matrix.row(i).dot(alpha_f);
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }


    dw_dphi_mat = dw_dphi(); // d alpha_w / d phi
    d2w_dphidphi_list = d2w_dphidphi(); // d^2 alpha_w / d phi d phi


    dlogmu_dw_mat = dlogmu_dw();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_dw_mat = dmu_dw();
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();
   
    he_alpha_w_mat = he_alpha_w();
   
    he_alpha_f_alpha_w_mat = he_alpha_f_alpha_w();


    // obtain hessian
    he_alpha_f_mat = he_alpha_f();
    he_betaR_mat = he_betaR();
    he_betaF_mat = he_betaF();
    he_phi_mat = he_phi();
    he_alpha_f_phi_mat = he_alpha_f_phi();
    he_alpha_f_betaF_mat = he_alpha_f_betaF();
    he_alpha_f_betaR_mat = he_alpha_f_betaR();
    he_phi_betaF_mat = he_phi_betaF();
    he_phi_betaR_mat = he_phi_betaR();
    he_betaR_betaF_mat = he_betaR_betaF();
  }

  // update full gradient and hessian of alpha_f, phi and betaR and betaF
  void derivative_he () {
    

    he_s_u_mat.setZero();

    // he_s_mat = [he_alpha_f_mat, he_alpha_f_phi_mat, he_alpha_f_log_theta_vec, he_alpha_f_log_smoothing_f_vec, 0;
    //             he_alpha_f_phi_mat.transpose(), he_phi_mat, he_phi_log_theta_vec, 0, he_phi_log_smoothing_w_vec;
    //             he_alpha_f_log_theta_vec.transpose(), he_phi_log_theta_vec.transpose(), he_log_theta_scalar, 0, 0;
    //             he_alpha_f_log_smoothing_f_vec.transpose(), 0, 0, he_log_smoothing_f_scalar, 0;
    //             0, he_phi_log_smoothing_w_vec.transpose(), 0, 0, he_log_smoothing_w_scalar]
    he_s_u_mat.block(0, 0, kE, kE)  = he_alpha_f_mat;
    he_s_u_mat.block(0, kE, kE, kwopt) = he_alpha_f_phi_mat;
    he_s_u_mat.block(kE, 0, kwopt, kE) = he_alpha_f_phi_mat.transpose();
    he_s_u_mat.block(kE, kE, kwopt, kwopt) = he_phi_mat;


    he_s_u_mat.block(kE+kwopt, kE+kwopt, kbetaR, kbetaR) = he_betaR_mat;
    he_s_u_mat.block(kE+kwopt+kbetaR, kE+kwopt+kbetaR, kbetaF, kbetaF) = he_betaF_mat;

    he_s_u_mat.block(0,kE+kwopt,kE,kbetaR) = he_alpha_f_betaR_mat;
    he_s_u_mat.block(kE,kE+kwopt,kwopt,kbetaR) = he_phi_betaR_mat;
    he_s_u_mat.block(0,kE+kwopt+kbetaR,kE,kbetaF) = he_alpha_f_betaF_mat;
    he_s_u_mat.block(kE,kE+kwopt+kbetaR,kwopt,kbetaF) = he_phi_betaF_mat;
    he_s_u_mat.block(kE+kwopt, kE+kwopt+kbetaR, kbetaR, kbetaF) = he_betaR_betaF_mat;

    he_s_u_mat.block(kE+kwopt,0,kbetaR,kE) = he_alpha_f_betaR_mat.transpose();
    he_s_u_mat.block(kE+kwopt,kE,kbetaR,kwopt) = he_phi_betaR_mat.transpose();
    he_s_u_mat.block(kE+kwopt+kbetaR,0,kbetaF,kE) = he_alpha_f_betaF_mat.transpose();
    he_s_u_mat.block(kE+kwopt+kbetaR,kE,kbetaF,kwopt) = he_phi_betaF_mat.transpose();
    he_s_u_mat.block(kE+kwopt+kbetaR, kE+kwopt, kbetaF, kbetaR) = he_betaR_betaF_mat.transpose();

    // make it symmetric. Comment out ...
    // he_s_u_mat = (he_s_u_mat + he_s_u_mat.transpose())/2.0;
  }

  



  // ********* Derivatives *************

  // FUNCTIONS
  // 1. density function
  // d log(exponential family density) / d mu
  Vec dlogdensity_dmu () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = y(i) / mu(i) - (theta + y(i)) / (theta + mu(i));
    }
    return out;
  }
  // d^2 log(exponential family density) / d mu^2
  Vec d2logdensity_dmudmu () {
    Vec out(n);
    for (int i = 0; i < n; i++) {
      out(i) = - y(i) / pow(mu(i), 2) + (theta + y(i)) / pow(theta + mu(i), 2);
    }
    return out;
  }
 


  
  // 2. mean model
  // d log(mu) / d alpha_w
  Mat dlogmu_dw () {
    Mat out(n, kw);
    out.setZero();
    Vec Bf1st;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < M; j++) {
        Bf1st = BsplinevecCon1st(E(i,j), knots_f_list.at(j), 4, Zf_list.at(j));
        out.block(i, j*kwp, 1, kwp) = B_inner_list.at(j).row(i) * (Bf1st.dot(alpha_f.segment(j*kEp, kEp)));
      }
    }
    return out;
  }
  // d mu / d alpha_w
  Mat dmu_dw () {
    Mat out(n, kw);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dw_mat.row(i) * mu(i);
    }
    return out;
  }



  // 3. Re-parameterization
  // d alpha_w / d phi
  Mat dw_dphi () {
    // deriv_g <- diag(1/as.numeric(sqrt(t(phi_long) %*% Dw %*% phi_long)), kw) - phi_long %*% (t(phi_long) %*% Dw) * (as.numeric(t(phi_long) %*% Dw %*% phi_long)^(-3/2))
    // deriv_g <- deriv_g[1:(kw), 2:kw]

    // alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long));

    // diag(1/as.numeric(sqrt(t(phi_long) %*% Dw %*% phi_long)), kw)
    Mat out(kw, kwopt);
    out.setZero();
    Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D(kwp);
    Mat Ddense;
    Mat deriv_g2;
    Mat deriv_g;

    for (int i = 0; i < M; i++) {
      alpha_w_C_denominator = sqrt(phi_long.segment(i*kwp, kwp).dot(Dw * phi_long.segment(i*kwp, kwp)));

      D.diagonal().setConstant(1.0/alpha_w_C_denominator);
      Ddense = D.toDenseMatrix();

      // phi_long %*% (t(phi_long) %*% Dw) * (as.numeric(t(phi_long) %*% Dw %*% phi_long)^(-3/2))
      deriv_g2 = (phi_long.segment(i*kwp, kwp) * (Dw * phi_long.segment(i*kwp, kwp)).transpose()) * pow(alpha_w_C_denominator, -3);

      deriv_g = Ddense - deriv_g2;
      // Remove the first column
      out.block(i*kwp, i*kwpopt, kwp, kwpopt) = deriv_g.block(0, 1, kwp, kwp - 1) * K.block(i*(kwp-1), i*kwpopt, kwp-1, kwpopt); 
    }
    return out;
  }


  std::vector<Mat> d2w_dphidphi () {
    std::vector<Mat> out;
    // alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long)) = tmp^(1/2);
    Mat outtmp(kwopt, kwopt);
    
    Vec Dwphi;
    Mat outlarge(kwp, kwp);
    Scalar tmp1;
    Scalar tmp2;
    Mat m1(kwp, kwp);
    Mat m2(kwp, kwp);
    for (int i = 0; i < M; i++) {
      alpha_w_C_denominator = sqrt(phi_long.segment(i*kwp, kwp).dot(Dw * phi_long.segment(i*kwp, kwp)));
      tmp1 = pow(alpha_w_C_denominator, -3); // pow(tmp, -1.5)
      tmp2 = pow(alpha_w_C_denominator, -5); // pow(tmp, -2.5)
      Dwphi = Dw * phi_long.segment(i*kwp, kwp);
      
      for (int s = 0; s < kwp; s++) {
        if (s == 0) {
          outlarge = -1.0 * tmp1 * Dw + 3.0 * tmp2 * Dwphi * Dwphi.transpose();
        } else {
          m1.setZero();
          m1.row(s) = Dwphi.transpose()*tmp1;
          m1.col(s) = m1.col(s) + Dwphi*tmp1;
          m2 = -1.0 * tmp1 * Dw + 3.0 * tmp2 * Dwphi * Dwphi.transpose();
          outlarge = -1.0 * m1 + m2; 
        }
        outtmp.setZero();
        outtmp.block(i*kwpopt, i*kwpopt, kwpopt, kwpopt) = K.block(i*(kwp-1), i*kwpopt, kwp-1, kwpopt).transpose() * outlarge.block(1, 1, kwp-1, kwp-1) * K.block(i*(kwp-1), i*kwpopt, kwp-1, kwpopt);
        out.push_back(outtmp); // dim = (kwopt * kwopt)
      }
    }
    
    return out;
  }


  // *** Hessian ***
  Mat he_alpha_f () {
    Mat out1(kE, kE);
    Mat out2(kE, kE);
    Mat outpen(kE, kE);
    out1.setZero();
    out2.setZero();
    outpen.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Bf_matrix.row(i)).transpose() * (Bf_matrix.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Bf_matrix.row(i).transpose() * Bf_matrix.row(i));
    }

    for (int i = 0; i < M; i++) outpen.block(i*kEp, i*kEp, kEp, kEp) = smoothing_f(i)*Sf_list.at(i);
    
    return - out1 - out2 + outpen;
  }


  Mat he_betaR () {
    Mat out1(kbetaR, kbetaR);
    Mat out2(kbetaR, kbetaR);
    Mat Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    out1.setZero();
    out2.setZero();
    Ones.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Xrand.row(i)).transpose() * (Xrand.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Xrand.row(i).transpose() * Xrand.row(i));
    }
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      double ri = CppAD::Value(r(i));
      int ki = static_cast<int>(ri);
      for (int j = 0; j < ki; j++) Ones(begin + j, begin + j) = smoothing(i);
      begin += ki;
    }
    return - out1 - out2 + Ones;
  }


  Mat he_betaF () {
    Mat out1(kbetaF, kbetaF);
    Mat out2(kbetaF, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Xfix.row(i)).transpose() * (Xfix.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Xfix.row(i).transpose() * Xfix.row(i));
    }
    return - out1 - out2;
  }
  
  Mat he_alpha_w () {
    Mat out1(kw, kw);
    Mat out2(kw, kw);
    Mat outpen(kw, kw);
    out1.setZero();
    out2.setZero();
    outpen.setZero();
    Mat tmp(kw, kw);
    
    for (int i = 0; i < n; i++) {
      tmp.setZero();
      for (int j = 0; j < M; j++) tmp.block(j*kwp, j*kwp, kwp, kwp) = (BsplinevecCon2nd(E(i, j), knots_f_list.at(j), 4, Zf_list.at(j)).dot(alpha_f.segment(j*kEp, kEp))) * B_inner_list.at(j).row(i).transpose() * B_inner_list.at(j).row(i);
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dw_mat.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i) * tmp);
    }
    for (int j = 0; j < M; j++) outpen.block(j*kwp, j*kwp, kwp, kwp) = smoothing_w(j)*Sw_large;
    return - out1 - out2 + outpen;
  }



  Mat he_phi () {
    Mat out1 = dw_dphi_mat.transpose() * he_alpha_w_mat * dw_dphi_mat;
    Mat out2(kwopt, kwopt);
    out2.setZero();

    Vec gr_pen_w(kw);
    for (int i = 0; i < M; i++) gr_pen_w.segment(i*kwp, kwp) = smoothing_w(i) * Sw_large * alpha_w_C.segment(i*kwp, kwp);

    Vec gr_alpha_w_vec = - dmu_dw_mat.transpose() * dlogdensity_dmu_vec + gr_pen_w;

    for (int s = 0; s < kw; s++) {
      out2 = out2 + gr_alpha_w_vec(s) * d2w_dphidphi_list.at(s);
    }
    return out1 + out2;
  }


  Mat he_alpha_f_alpha_w () {
    Mat out1(kE, kw);
    Mat out2(kE, kw);
    out1.setZero();
    out2.setZero();
    Mat tmp(kE, kw);
    
    for (int i = 0; i < n; i++) {
      tmp.setZero();
      for (int j = 0; j < M; j++) tmp.block(j*kEp, j*kwp, kEp, kwp) = BsplinevecCon1st(E(i, j), knots_f_list.at(j), 4, Zf_list.at(j)) * B_inner_list.at(j).row(i);
      out1 += d2logdensity_dmudmu_vec(i) * (Bf_matrix.row(i) * mu(i)).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * Bf_matrix.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i)*tmp);
    }
    return - out1 - out2;
  }

  Mat he_alpha_f_phi () {
    Mat out = he_alpha_f_alpha_w_mat * dw_dphi_mat;
    return out;
  }
  Mat he_alpha_f_betaF () {
    Mat out1(kE, kbetaF);
    Mat out2(kE, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Bf_matrix.row(i)).transpose() * (Xfix.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * Bf_matrix.row(i).transpose() * (Xfix.row(i) * mu(i));
    }
    return - out1 - out2;
  }
  Mat he_alpha_f_betaR () {
    Mat out1(kE, kbetaR);
    Mat out2(kE, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Bf_matrix.row(i)).transpose() * (Xrand.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * Bf_matrix.row(i).transpose() * (Xrand.row(i) * mu(i));
    }
    return - out1 - out2;
  }
  Mat he_phi_betaR () {
    Mat out1(kwopt, kbetaR);
    Mat out2(kwopt, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dw_dphi_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * (Xrand.row(i) * mu(i)));
      out2 += dw_dphi_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dw_mat.row(i).transpose() * (Xrand.row(i) * mu(i)));
    }
    return - out1 - out2;
  }
  Mat he_phi_betaF () {
    Mat out1(kwopt, kbetaF);
    Mat out2(kwopt, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dw_dphi_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * (Xfix.row(i) * mu(i)));
      out2 += dw_dphi_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dw_mat.row(i).transpose() * (Xfix.row(i) * mu(i)));
    }
    return - out1 - out2;
  }
  Mat he_betaR_betaF () {
    Mat out1(kbetaR, kbetaF);
    Mat out2(kbetaR, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * (Xrand.row(i)).transpose() * (Xfix.row(i)) * mu(i) * mu(i);
      out2 += dlogdensity_dmu_vec(i) * Xrand.row(i).transpose() * (Xfix.row(i) * mu(i));
    }
    return - out1 - out2;
  }





  // *********** LAML ***********
  Scalar logdetH05() {
    Scalar logdetH05 = 0.0;
    Eigen::FullPivLU<Mat> lu(he_s_u_mat);
    Mat LU = lu.matrixLU();
    // Scalar c = lu.permutationP().determinant(); // -1 or 1
    Scalar lii;
    for (int i = 0; i < (kwopt+kE+kbetaR+kbetaF); i++) {
      lii = LU(i,i);
      // std::cout << "lii : " << CppAD::Value(lii) << std::endl;
      // std::cout << "c : " << CppAD::Value(c) << std::endl;
      // if (lii < 0.0) c *= -1;

      logdetH05 += log(CppAD::abs(lii));
    }
    // logdetH05 += log(c);



    return logdetH05/2.0;

  }
};







#endif