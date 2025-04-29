#ifndef ACEDLNMEIGENCLASS_HPP
#define ACEDLNMEIGENCLASS_HPP




#include <lambda_lanczos.hpp>
using lambda_lanczos::LambdaLanczos;


#include <cmath>

#include <iostream>
using namespace std;



// ************** PART 1: Define some functions **********************
void choleskyAD(Eigen::MatrixXd& L) {
  // Eigen::MatrixXd will be overwritten; lower triangle will be its Cholesky. Only the lower triangle is computed/stored
  int s = L.cols();
  for (int k = 0; k < s; k++) {
    // (a) define pivot
    L(k,k) = sqrt(L(k,k));

    // (b) adjust lead column
    for (int j = k+1; j < s; j++) L(j,k) /= L(k,k);
    for (int j = k+1; j < s; j++)
      for (int i = j; i < s; i++) L(i,j) -= L(i,k) * L(j,k);
  }
}

Eigen::MatrixXd invertL(Eigen::MatrixXd &L) {
  // inverse of a lower triangular matrix
  int n = L.cols();
  Eigen::MatrixXd M(n, n);
  M.setZero();
  for (int i = 0; i < n; i++)
  {
    M(i,i) = 1.0 / L(i,i);
    for (int j = 0; j < i; j++)
    {
      for (int k = j; k < i; k++) M(i,j) += L(i,k) * M(k,j);
      M(i,j) = -M(i,j) / L(i,i);
    }
  }
  return M;
}
// check whether there is nan in the input vector
bool hasNaN(Eigen::VectorXd vec) {
    for (int i = 0; i < vec.size(); i++) {
      if( std::isnan(vec(i))) return true; // has nan
    }
    return false; // No nan
}

// TODO: the bspline function evaluated at the points outside the boundaries are incorret!
// Bspline(l=0) = 0,0,0,0... It should be a linear function of l, not always equal to 0.
int knotindexEigen(double x,const Eigen::VectorXd t) {
  int q = t.size();
  int k=0;
  if (x < t(0)) return -1;
  while(x>=t(k)){
    k++;
    if (k >= q) break;
  }

  return k-1;
}

double weight(double x, const Eigen::VectorXd& t,int i,int k) {
  if (t(i+k-1) != t(i-1))
    return((x - t(i-1))/(t(i+k-1)-t(i-1)));
  return 0.;
}


double Bspline(double x, int j, const Eigen::VectorXd& t,int p) {
  // Evaluate the jth B-spline
  // B_p(x) of order p (degree p-1) at x
  if (p==1)
    return(x>=t(j-1) && x<t(j+1-1));

  double w1 = weight(x,t,j,p-1);
  double w2 = weight(x,t,j+1,p-1);
  double b1 = Bspline(x,j,t,p-1);
  double b2 = Bspline(x,j+1,t,p-1);

  return w1*b1 + (1.-w2)*b2;
}

double Bspline1st(double x, int j, const Eigen::VectorXd& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  double bb1 = Bspline(x,j+1,t,p-1);
  double bb2 = Bspline(x,j,t,p-1);

  double ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  double ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

double Bspline2nd(double x, int j, const Eigen::VectorXd& t,int p) {
  // 1st derivative of the jth B-spline
  // https://stats.stackexchange.com/questions/258786/what-is-the-second-derivative-of-a-b-spline
  double bb1 = Bspline1st(x,j+1,t,p-1);
  double bb2 = Bspline1st(x,j,t,p-1);

  double ww1 = 0.0;
  if (t(j+p-1) != t(j))
    ww1 = -1.0/(t(j+p-1)-t(j));
  double ww2 = 0.0;
  if (t(j+p-2) != t(j-1))
    ww2 = 1.0/(t(j+p-2)-t(j-1));

  return (p-1.0) * (ww1 * bb1 + ww2 * bb2);
}

Eigen::VectorXd Bsplinevec(double x, const Eigen::VectorXd& t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  // int k = knotindexEigen(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return b;
}

Eigen::VectorXd BsplinevecCon(double x, const Eigen::VectorXd& t, int p, Eigen::MatrixXd Z) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Eigen::VectorXd Bsplinevec1st(double x, const Eigen::VectorXd& t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  // int k = knotindexEigen(x,t);
  // for (int i=(k-(p-1));i<k+1;i++)
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return b;
}

// [[Rcpp::export]]
Eigen::VectorXd BsplinevecCon1st(double x, const Eigen::VectorXd& t, int p, Eigen::MatrixXd Z) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline1st(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}

Eigen::VectorXd Bsplinevec2nd(double x, const Eigen::VectorXd& t,int p) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return b;
}

// [[Rcpp::export]]
Eigen::VectorXd BsplinevecCon2nd(double x, const Eigen::VectorXd& t,int p, Eigen::MatrixXd Z) {
  int m = t.size() - p;
  Eigen::VectorXd b(m);
  b.setZero();
  for (int i=0;i<m;i++)
    b(i) = Bspline2nd(x,i+1,t,p);
  return Z.transpose()*b;
  // return b.transpose()*Z;
}





// https://github.com/tminka/lightspeed/blob/master/digamma.m
double lgamma1st (double x) {
  const double pi = 3.141592653589793238462643383279;
  const double large = 9.5;
  const double d1 = -0.5772156649015328606065121;
  const double d2 = pi*pi/6.0;
  const double small = 1e-6;
  const double s3 = 1.0/12.0;
  const double s4 = 1.0/120.0;
  const double s5 = 1.0/252.0;
  const double s6 = 1.0/240.0;
  const double s7 = 1.0/132.0;
  const double s8 = 691.0/32760.0;
  const double s9 = 1.0/12.0;
  const double s10 = 3617.0/8160.0;

  // Use de Moivre's expansion if x >= large = 9.5
  // calculate lgamma1st(x+10)
  double xplus10 = x + 10.0;
  double y = 0.0;
  double r = 1.0 / xplus10;
  y += log(xplus10) - 0.5 * r;
  r = r * r;
  y = y - r * ( s3 - r * ( s4 - r * (s5 - r * (s6 - r * s7))));

  // lgamma1st(x+10) = (1/x + 1/(x+1) + ... + 1/(x+9)) + lgamma1st(x)
  y = y - 1.0/x - 1.0/(x+1.0) - 1.0/(x+2.0) - 1.0/(x+3.0) - 1.0/(x+4.0) - 1.0/(x+5.0) - 1.0/(x+6.0) - 1.0/(x+7.0) - 1.0/(x+8.0) - 1/(x+9);

  return y;
}



// https://github.com/tminka/lightspeed/blob/master/trigamma.m

// double lgamma2nd (double x) {
//   const double pi = 3.141592653589793238462643383279;
//   const double c = pi*pi/6;
//   const double c1 = -2.404113806319188570799476;
//   const double b2 =  1.0/6.0;
//   const double b4 = -1.0/30.0;
//   const double b6 =  1.0/42.0;
//   const double b8 = -1.0/30.0;
//   const double b10 = 5.0/66.0;

//   // TO DO: % Reduce to trigamma(x+n) where ( X + N ) >= large.

//   double z = 1./(x*x);
//   double y = 0.5*z + (1.0 + z*(b2 + z*(b4 + z*(b6 + z*(b8 + z*b10))))) / x;

//   // std::cout << "x" << CppAD::Value(x) << std::endl;
//   // std::cout << "trigamma" << CppAD::Value(y) << std::endl;
//   return y;
// }

// **************** PART 2: g(mu) = DL term + linear term + smooth term *************************
class Model {
  // The DLNM model

private:
  // DATA
  const Eigen::VectorXd y; // Response
  const std::vector<Eigen::MatrixXd> Sf_list; // penalty matrix for f(E)
  const Eigen::MatrixXd Sw_large;
  const std::vector<Eigen::MatrixXd> B_inner_list;
  const std::vector<Eigen::VectorXd> knots_f_list; // knots for f(E) B-spline
  

  const Eigen::MatrixXd Xfix; // fixed effects
  const Eigen::MatrixXd Xrand; // random effects
  const Eigen::VectorXd r; // rank of each smooth
  const std::vector<Eigen::MatrixXd> Zf_list;

  const Eigen::VectorXd Xoffset; // offset
  
public:
  // DATA
  const Eigen::MatrixXd Dw;  // \int w(l)^2 dl = 1
  const Eigen::MatrixXd K;
  const Eigen::VectorXd a; 

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
  int M; // number of exposure. Additive distributed lag terms in this model 

  // PARAMETERS
  Eigen::VectorXd alpha_f;
  Eigen::VectorXd phi;
  Eigen::VectorXd phiKa; // phitKa = K * phi + a
  double log_theta;
  Eigen::VectorXd log_smoothing_f;
  Eigen::VectorXd log_smoothing_w;

  Eigen::VectorXd betaF; // parameters for fixed effects
  Eigen::VectorXd betaR; // parameters for random effects
  Eigen::VectorXd logsmoothing; // log smoothing parameters for random effects

  // Components generated
  double theta;
  Eigen::VectorXd smoothing_f;
  Eigen::VectorXd smoothing_w;
  Eigen::VectorXd smoothing;
  Eigen::VectorXd phi_long;
  double alpha_w_C_denominator;
  Eigen::VectorXd alpha_w_C;
  // Eigen::VectorXd alpha_w_C_pen;
  Eigen::MatrixXd Bf_matrix;
  // Eigen::VectorXd Bf;
  Eigen::MatrixXd E;
  Eigen::VectorXd eta;
  Eigen::VectorXd eta_remaining; // remaining terms = Xfix * betaF + Xrand * betaR
  Eigen::VectorXd mu; // log(mu) = eta + eta_remaining + Xoffset
  double NegLogL; // NegativeLogLikelihood value


  // Components for derivatives

  Eigen::MatrixXd dlogmu_df_mat;
  Eigen::MatrixXd dlogmu_dbetaR_mat;
  Eigen::MatrixXd dlogmu_dbetaF_mat;
  Eigen::MatrixXd dlogmu_dw_mat;
  Eigen::VectorXd dlogdensity_dmu_vec;
  Eigen::MatrixXd dmu_df_mat;
  Eigen::MatrixXd dmu_dbetaR_mat;
  Eigen::MatrixXd dmu_dbetaF_mat;
  Eigen::MatrixXd dmu_dw_mat;
  Eigen::MatrixXd dw_dphi_mat;
  Eigen::VectorXd gr_alpha_w_vec;
  Eigen::VectorXd d2logdensity_dmudmu_vec;
  // std::vector<Eigen::MatrixXd> d2mu_dfdf_list;
  // std::vector<Eigen::MatrixXd> d2mu_dbetaRdbetaR_list;
  // std::vector<Eigen::MatrixXd> d2mu_dbetaFdbetaF_list;
  // std::vector<Eigen::MatrixXd> d2logmu_dwdw_list;
  // std::vector<Eigen::MatrixXd> d2mu_dwdw_list;
  std::vector<Eigen::MatrixXd> d2w_dphidphi_list;
  // std::vector<Eigen::MatrixXd> d2logmu_dfdw_list;
  // std::vector<Eigen::MatrixXd> d2mu_dfdw_list;
  Eigen::MatrixXd he_alpha_w_mat;
  Eigen::MatrixXd he_alpha_f_alpha_w_mat;
  double dlogdensity_dtheta_scalar;
  // double d2logdensity_dthetadtheta_scalar;
  Eigen::VectorXd d2logdensity_dmudtheta_vec;


  // gradient and hessian for updating alpha_f, betaR and betaF
  Eigen::VectorXd gr_inner_vec;
  Eigen::MatrixXd he_inner_mat;

  // full gradient
  Eigen::VectorXd gr_alpha_f_vec;
  Eigen::VectorXd gr_betaR_vec;
  Eigen::VectorXd gr_betaF_vec;
  Eigen::VectorXd gr_phi_vec;
  Eigen::VectorXd gr_log_smoothing_f_vec;
  Eigen::VectorXd gr_log_smoothing_w_vec;
  double gr_log_theta_scalar;
  Eigen::VectorXd gr_logsmoothing_vec;

  Eigen::VectorXd gr_s_u_vec;
  Eigen::VectorXd gr_s_par_vec;

  // full hessian
  Eigen::MatrixXd he_alpha_f_mat;
  Eigen::MatrixXd he_betaR_mat;
  Eigen::MatrixXd he_betaF_mat;
  Eigen::MatrixXd he_phi_mat;
  Eigen::MatrixXd he_alpha_f_phi_mat;
  Eigen::MatrixXd he_alpha_f_betaF_mat;
  Eigen::MatrixXd he_alpha_f_betaR_mat;
  Eigen::MatrixXd he_phi_betaF_mat;
  Eigen::MatrixXd he_phi_betaR_mat;
  Eigen::MatrixXd he_betaR_betaF_mat;
  // double he_log_smoothing_f_scalar;
  // double he_log_smoothing_w_scalar;
  // double he_log_theta_scalar;
  Eigen::MatrixXd he_alpha_f_log_smoothing_f_mat;
  Eigen::MatrixXd he_betaR_logsmoothing_mat;
  Eigen::MatrixXd he_phi_log_smoothing_w_mat;
  Eigen::VectorXd he_alpha_f_log_theta_vec;
  Eigen::VectorXd he_phi_log_theta_vec;
  Eigen::VectorXd he_betaR_log_theta_vec;
  Eigen::VectorXd he_betaF_log_theta_vec;

  Eigen::MatrixXd he_s_u_mat;
  Eigen::MatrixXd he_s_par_u_mat;

  // To compute AIC
  double NegLogL_l; // NegativeLogLikelihood without penalty
  // matrix for I (hessian of log likelihood without penalty)
  Eigen::MatrixXd I_alpha_f_mat;
  Eigen::MatrixXd I_betaR_mat;
  Eigen::MatrixXd I_phi_mat;
  Eigen::MatrixXd I_alpha_w_mat;
  Eigen::MatrixXd I_mat; 

  // results for profile likelihood
  Eigen::VectorXd PL_gradient;
  Eigen::MatrixXd PL_hessian;
  int converge; // 0: converge. 99: not converge

  
  // Constructor
  Model(const Eigen::VectorXd& y_,
        const std::vector<Eigen::MatrixXd>& B_inner_list_,
        const std::vector<Eigen::VectorXd>& knots_f_list_,
        const Eigen::MatrixXd& Sw_large_,
        const std::vector<Eigen::MatrixXd>& Sf_list_,
        const Eigen::MatrixXd& Dw_,
        const Eigen::MatrixXd& Xrand_,
        const Eigen::MatrixXd& Xfix_,
        const std::vector<Eigen::MatrixXd>& Zf_list_,
        const Eigen::VectorXd& Xoffset_,
        const Eigen::VectorXd& r_,
        const Eigen::MatrixXd& K_,
        const Eigen::VectorXd& a_,
        Eigen::VectorXd& alpha_f_,
        Eigen::VectorXd& phi_,
        double log_theta_,
        const Eigen::VectorXd& log_smoothing_f_,
        const Eigen::VectorXd& log_smoothing_w_,
        Eigen::VectorXd& betaR_,
        Eigen::VectorXd& betaF_,
        Eigen::VectorXd& logsmoothing_) :
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
          Bf_matrix.block(i, j*kEp, 1, kEp) = BsplinevecCon(E(i,j), knots_f_list.at(j), 4, Zf_list.at(j)).transpose();
        }
        eta(i) = Bf_matrix.row(i).dot(alpha_f);
        eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
        mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
      }

      // Initialize the derivative components and NegativeLogLikelihood
      dw_dphi_mat = dw_dphi(); // d alpha_w / d phi
      d2w_dphidphi_list = d2w_dphidphi(); // d^2 alpha_w / d phi d phi
      gr_s_u_vec.resize(kwopt+kE+kbetaR+kbetaF);
      he_s_u_mat.resize(kwopt+kE+kbetaR+kbetaF, kwopt+kE+kbetaR+kbetaF);
      gr_s_par_vec.resize(1+2*M+p);
      he_s_par_u_mat.resize(1+2*M+p, kwopt+kE+kbetaR+kbetaF);
      gr_inner_vec.resize(kE+kbetaR+kbetaF);
      he_inner_mat.resize(kE+kbetaR+kbetaF, kE+kbetaR+kbetaF);

      derivative_coef(); // comment out for debug
      derivative_he(); // comment out for debug
      derivative_full(); // comment out for debug
      NegativeLogLikelihood(); // comment out for debug

      // Initialize PL
      PL_gradient.resize(kwopt);
      PL_hessian.resize(kwopt, kwopt);
    }

  // Functions to set parameters
  void setAlphaF(const Eigen::VectorXd alpha_f_) {
    alpha_f = alpha_f_;

    for (int i = 0; i < n; i++) {
      eta(i) = Bf_matrix.row(i).dot(alpha_f);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }

  void setPhi(const Eigen::VectorXd phi_) {
    phi = phi_;

    phiKa = K * phi + a;
    
    for (int i = 0; i < M; i++) {
      phi_long(i*kwp) = 1.0;
      phi_long.segment(i*kwp + 1, kwp-1) = phiKa.segment(i*(kwp-1), (kwp-1));
    }

    for (int i = 0; i < M; i++) {
      alpha_w_C_denominator = sqrt(phi_long.segment(i*kwp, kwp).dot(Dw * phi_long.segment(i*kwp, kwp)));
      alpha_w_C.segment(i*kwp, kwp) = phi_long.segment(i*kwp, kwp) / alpha_w_C_denominator;
      // alpha_w_C_pen.segment(i*(kwp-1), (kwp-1)) = phi.segment(i*(kwp-1), (kwp-1)) / alpha_w_C_denominator;
    }
    
    for (int i = 0; i < M; i++) {
      E.col(i) = B_inner_list.at(i) * alpha_w_C.segment(i*kwp, kwp);
    }

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < M; j++) {
        Bf_matrix.block(i, j*kEp, 1, kEp) = BsplinevecCon(E(i,j), knots_f_list.at(j), 4, Zf_list.at(j)).transpose();
      }
      eta(i) = Bf_matrix.row(i).dot(alpha_f);
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }

    dw_dphi_mat = dw_dphi(); // d alpha_w / d phi
    d2w_dphidphi_list = d2w_dphidphi(); // d^2 alpha_w / d phi d phi

  }
  void setBetaF(const Eigen::VectorXd betaF_) {
    betaF = betaF_;
    for (int i = 0; i < n; i++) {
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }
  void setBetaR(const Eigen::VectorXd betaR_) {
    betaR = betaR_;
    for (int i = 0; i < n; i++) {
      eta_remaining(i) = Xfix.row(i).dot(betaF) + Xrand.row(i).dot(betaR);
      mu(i) = exp(eta(i) + eta_remaining(i) + Xoffset(i));
    }
  }
  void setLogTheta(const double log_theta_) {
    log_theta = log_theta_;
    theta = exp(log_theta);
  }

  void setLogSmoothingF(const Eigen::VectorXd log_smoothing_f_) {
    log_smoothing_f = log_smoothing_f_;
    for (int i = 0; i < M; i++) smoothing_f(i) = exp(log_smoothing_f(i));
  }

  void setLogSmoothingW(const Eigen::VectorXd log_smoothing_w_) {
    log_smoothing_w = log_smoothing_w_;
    for (int i = 0; i < M; i++) smoothing_w(i) = exp(log_smoothing_w(i));
  }
  void setLogsmoothing(const Eigen::VectorXd logsmoothing_) { // log smoothing parameters for remaining terms
    logsmoothing = logsmoothing_;
    for (int i = 0; i < p; i++) smoothing(i) = exp(logsmoothing(i));
  }

  // get private members
  std::vector<Eigen::MatrixXd> getB_inner_list () {
    return B_inner_list;
  }
  Eigen::MatrixXd getDw () {
    return Dw;
  }

  std::vector<Eigen::VectorXd> getknots_f_list () {
    return knots_f_list;
  }
  std::vector<Eigen::MatrixXd> getZf_list () {
    return Zf_list;
  }
  Eigen::MatrixXd getXfix () {
    return Xfix;
  }
  Eigen::MatrixXd getXrand () {
    return Xrand;
  }
  Eigen::VectorXd getXoffset() {
    return Xoffset;
  }

  // Function to update derivatives.
  // RUN the function derivative_coef(), derivative_he() and derivative_full() after update parameters.
  // update derivatives related to spline coefficients alpha_f and phi, and betaR and betaF
  void derivative_coef() {
    dlogmu_dw_mat = dlogmu_dw();
    dlogmu_df_mat = dlogmu_df();
    dlogmu_dbetaR_mat = dlogmu_dbetaR();
    dlogmu_dbetaF_mat = dlogmu_dbetaF();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_df_mat = dmu_df();
    dmu_dbetaR_mat = dmu_dbetaR();
    dmu_dbetaF_mat = dmu_dbetaF();
    dmu_dw_mat = dmu_dw();
    gr_alpha_w_vec = gr_alpha_w();
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();
    // d2mu_dfdf_list = d2mu_dfdf();
    // d2mu_dbetaRdbetaR_list = d2mu_dbetaRdbetaR();
    // d2mu_dbetaFdbetaF_list = d2mu_dbetaFdbetaF();
    // d2logmu_dwdw_list = d2logmu_dwdw();
    // d2mu_dwdw_list = d2mu_dwdw();
    he_alpha_w_mat = he_alpha_w();
    // d2logmu_dfdw_list = d2logmu_dfdw();
    // d2mu_dfdw_list = d2mu_dfdw();
    he_alpha_f_alpha_w_mat = he_alpha_f_alpha_w();
    dlogdensity_dtheta_scalar = dlogdensity_dtheta();
    // d2logdensity_dthetadtheta_scalar = d2logdensity_dthetadtheta();
    d2logdensity_dmudtheta_vec = d2logdensity_dmudtheta();

    // obtain gradient
    gr_alpha_f_vec = gr_alpha_f();
    gr_betaR_vec = gr_betaR();
    gr_betaF_vec = gr_betaF();
    gr_phi_vec = gr_phi();
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
    gr_s_u_vec << gr_alpha_f_vec, gr_phi_vec, gr_betaR_vec, gr_betaF_vec;


    he_s_u_mat.setZero();


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

  // update derivatives related to overdispersion and smoothing parameters
  // Full derivative for LAML
  void derivative_full () {
    // obtain full gradient

    gr_log_smoothing_f_vec = gr_log_smoothing_f();
    gr_log_smoothing_w_vec = gr_log_smoothing_w();
    gr_log_theta_scalar = gr_log_theta();
    gr_logsmoothing_vec = gr_logsmoothing();



    // u represents spline coefficient alpha_f and phi, and betaR and betaF
    // par represents overdispersion and smoothing parameters

    gr_s_par_vec << gr_log_theta_scalar, gr_log_smoothing_f_vec, gr_log_smoothing_w_vec, gr_logsmoothing_vec;


    // obtain full hessian
    // he_log_smoothing_f_scalar = he_log_smoothing_f();
    // he_log_smoothing_w_scalar = he_log_smoothing_w();
    // he_log_theta_scalar = he_log_theta();
    he_alpha_f_log_smoothing_f_mat = he_alpha_f_log_smoothing_f();
    he_phi_log_smoothing_w_mat = he_phi_log_smoothing_w();
    he_betaR_logsmoothing_mat = he_betaR_logsmoothing();
    he_alpha_f_log_theta_vec = he_alpha_f_log_theta();
    he_phi_log_theta_vec = he_phi_log_theta();
    he_betaR_log_theta_vec = he_betaR_log_theta();
    he_betaF_log_theta_vec = he_betaF_log_theta();


    he_s_par_u_mat.setZero();



    he_s_par_u_mat.row(0) << he_alpha_f_log_theta_vec.transpose(), he_phi_log_theta_vec.transpose(), he_betaR_log_theta_vec.transpose(), he_betaF_log_theta_vec.transpose();
    he_s_par_u_mat.block(1, 0, M, kE) = he_alpha_f_log_smoothing_f_mat.transpose();
    he_s_par_u_mat.block(1+M, kE, M, kwopt) = he_phi_log_smoothing_w_mat.transpose();
    he_s_par_u_mat.block(1+2*M, kE+kwopt, p, kbetaR) = he_betaR_logsmoothing_mat.transpose();
  }

  // update variables related to alpha_f, betaR and betaF.
  // Used only in updating alpha_f, betaR and betaF. .
  void derivative_f () {
    dlogmu_df_mat = dlogmu_df();
    dlogmu_dbetaR_mat = dlogmu_dbetaR();
    dlogmu_dbetaF_mat = dlogmu_dbetaF();
    dlogdensity_dmu_vec = dlogdensity_dmu();
    dmu_df_mat = dmu_df();
    dmu_dbetaR_mat = dmu_dbetaR();
    dmu_dbetaF_mat = dmu_dbetaF();
    d2logdensity_dmudmu_vec = d2logdensity_dmudmu();
    // d2mu_dfdf_list = d2mu_dfdf();
    // d2mu_dbetaRdbetaR_list = d2mu_dbetaRdbetaR();
    // d2mu_dbetaFdbetaF_list = d2mu_dbetaFdbetaF();

    gr_alpha_f_vec = gr_alpha_f();
    gr_betaR_vec = gr_betaR();
    gr_betaF_vec = gr_betaF();

    he_alpha_f_mat = he_alpha_f();
    he_betaR_mat = he_betaR();
    he_betaF_mat = he_betaF();
    he_alpha_f_betaF_mat = he_alpha_f_betaF();
    he_alpha_f_betaR_mat = he_alpha_f_betaR();
    he_betaR_betaF_mat = he_betaR_betaF();

    gr_inner_vec << gr_alpha_f_vec, gr_betaR_vec, gr_betaF_vec;

    he_inner_mat.setZero();
    he_inner_mat.block(0,0,kE,kE) = he_alpha_f_mat;
    he_inner_mat.block(kE,kE,kbetaR,kbetaR) = he_betaR_mat;
    he_inner_mat.block(kE+kbetaR,kE+kbetaR,kbetaF,kbetaF) = he_betaF_mat;

    he_inner_mat.block(0,kE,kE,kbetaR) = he_alpha_f_betaR_mat;
    he_inner_mat.block(0,kE+kbetaR,kE,kbetaF) = he_alpha_f_betaF_mat;
    he_inner_mat.block(kE,kE+kbetaR,kbetaR,kbetaF) = he_betaR_betaF_mat;

    he_inner_mat.block(kE,0,kbetaR,kE) = he_alpha_f_betaR_mat.transpose();
    he_inner_mat.block(kE+kbetaR,0,kbetaF,kE) = he_alpha_f_betaF_mat.transpose();
    he_inner_mat.block(kE+kbetaR,kE,kbetaF,kbetaR) = he_betaR_betaF_mat.transpose();

  }


  // functions for NegativeLogLikelihood
  void NegativeLogLikelihood() {

  

    double loglik = 0;
    for (int i = 0; i < n; i++) {
      loglik += lgamma(y(i) + theta) - lgamma(theta) - lgamma(y(i) + 1) -
                                    theta * log(1 + mu(i)/theta) +
                                    y(i)*( eta(i) + eta_remaining(i) + Xoffset(i) - log_theta - log(1 + mu(i)/theta) );
    }

    // part 1: DLNM
    for (int i = 0; i < M; i++) {
      // Smooth Penalty
      loglik += -0.5 * smoothing_w(i) * alpha_w_C.segment(i*kwp, kwp).dot(Sw_large * alpha_w_C.segment(i*kwp, kwp)) - 0.5 * smoothing_f(i) * alpha_f.segment(i*kEp, kEp).dot(Sf_list.at(i) * alpha_f.segment(i*kEp, kEp));
      // Scale
      loglik += (kwp-1-1) / 2.0 * log_smoothing_w(i) + (kEp-1) / 2.0 * log_smoothing_f(i);
    }
    
    
    // part 2: Remaining smooth terms
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      Eigen::VectorXd betaRi(ki);
      for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      loglik += -0.5 * smoothing(i) * betaRi.dot(betaRi); // smooth penalty
      loglik += ki/2.0 * logsmoothing(i); // scale

      begin += ki;
    }


    NegLogL = -1.0 * loglik; // NEGATIVE log-likelihood
  }

  // functions for NegativeLogLikelihood WITHOUT penalty for AIC
  void NegativeLogLikelihood_l() {
    double loglik = 0;
    for (int i = 0; i < n; i++) {
      loglik += lgamma(y(i) + theta) - lgamma(theta) - lgamma(y(i) + 1) -
                                    theta * log(1 + mu(i)/theta) +
                                    y(i)*( eta(i) + eta_remaining(i) + Xoffset(i) - log_theta - log(1 + mu(i)/theta) );
    }
    NegLogL_l = -1.0 * loglik; // NEGATIVE log-likelihood
  }

  void prepare_AIC () {
    NegativeLogLikelihood_l();
    // hessian of log likelihood without penalty
    I_alpha_f_mat = I_alpha_f();
    I_betaR_mat = I_betaR();
    I_alpha_w_mat = I_alpha_w();
    I_phi_mat = I_phi();
    I_mat = he_s_u_mat;
    I_mat.block(0, 0, kE, kE)  = I_alpha_f_mat;
    I_mat.block(kE, kE, kwopt, kwopt) = I_phi_mat;
    I_mat.block(kE+kwopt, kE+kwopt, kbetaR, kbetaR) = I_betaR_mat;
  }


  // ********* Derivatives *************

  // FUNCTIONS
  // 1. density function
  // d log(exponential family density) / d mu
  Eigen::VectorXd dlogdensity_dmu () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = y(i) / mu(i) - (theta + y(i)) / (theta + mu(i));
    }
    return out;
  }
  // d^2 log(exponential family density) / d mu^2
  Eigen::VectorXd d2logdensity_dmudmu () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = - y(i) / pow(mu(i), 2) + (theta + y(i)) / pow(theta + mu(i), 2);
    }
    return out;
  }
  // d log(exponential family density) / d theta
  double dlogdensity_dtheta () {
    double out = 0.0;
    // std::cout << "x" << 3.5 << std::endl;

    // TO DO: optimize it. Use property of gamma function...
    for (int i = 0; i < n; i++) {
      out += log_theta - log(theta + mu(i)) + (mu(i) - y(i))/(theta+mu(i)) + lgamma1st(theta+y(i)) - lgamma1st(theta);
    }
    return out;
  }
  // d^2 log(exponential family density) / d theta^2
  // double d2logdensity_dthetadtheta () {
  //   double out = 0.0;

  //   for (int i = 0; i < n; i++) {
  //     out += 1/theta - 1/(theta + mu(i)) - (mu(i) - y(i)) / ((theta + mu(i))*(theta + mu(i))) + lgamma2nd(y(i) + theta) - lgamma2nd(theta);
  //   }
  //   return out;
  // }
  Eigen::VectorXd d2logdensity_dmudtheta () {
    Eigen::VectorXd out(n);
    for (int i = 0; i < n; i++) {
      out(i) = (y(i) - mu(i)) / pow(theta+mu(i), 2);
    }
    return out;
  }



  // 2. mean model
  // d log(mu) / d alpha_f
  Eigen::MatrixXd dlogmu_df () {
    return Bf_matrix;
  }
  // d mu / d alpha_f
  Eigen::MatrixXd dmu_df () {
    Eigen::MatrixXd out(n, kE);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_df_mat.row(i) * mu(i);
    }
    return out;
  }
  // d log(mu) / d alpha_w
  Eigen::MatrixXd dlogmu_dw () {
    Eigen::MatrixXd out(n, kw);
    out.setZero();
    Eigen::VectorXd Bf1st;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < M; j++) {
        Bf1st = BsplinevecCon1st(E(i,j), knots_f_list.at(j), 4, Zf_list.at(j));
        out.block(i, j*kwp, 1, kwp) = B_inner_list.at(j).row(i) * (Bf1st.dot(alpha_f.segment(j*kEp, kEp)));
      }
    }
    return out;
  }
  // d mu / d alpha_w
  Eigen::MatrixXd dmu_dw () {
    Eigen::MatrixXd out(n, kw);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dw_mat.row(i) * mu(i);
    }
    return out;
  }
  // d log(mu) / d betaR
  Eigen::MatrixXd dlogmu_dbetaR () {
    return Xrand;
  }
  // d mu / d betaR
  Eigen::MatrixXd dmu_dbetaR () {
    Eigen::MatrixXd out(n, kbetaR);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dbetaR_mat.row(i) * mu(i);
    }
    return out;
  }
  // d log(mu) / d betaF
  Eigen::MatrixXd dlogmu_dbetaF () {
    return Xfix;
  }
  // d mu / d betaR
  Eigen::MatrixXd dmu_dbetaF () {
    Eigen::MatrixXd out(n, kbetaF);
    for (int i = 0; i < n; i++) {
      out.row(i) = dlogmu_dbetaF_mat.row(i) * mu(i);
    }
    return out;
  }



  // 3. Re-parameterization
  // d alpha_w / d phi
  Eigen::MatrixXd dw_dphi () {
    // deriv_g <- diag(1/as.numeric(sqrt(t(phi_long) %*% Dw %*% phi_long)), kw) - phi_long %*% (t(phi_long) %*% Dw) * (as.numeric(t(phi_long) %*% Dw %*% phi_long)^(-3/2))
    // deriv_g <- deriv_g[1:(kw), 2:kw]

    // alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long));

    // diag(1/as.numeric(sqrt(t(phi_long) %*% Dw %*% phi_long)), kw)
    Eigen::MatrixXd out(kw, kwopt);
    out.setZero();
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> D(kwp);
    Eigen::MatrixXd Ddense;
    Eigen::MatrixXd deriv_g2;
    Eigen::MatrixXd deriv_g;

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


  std::vector<Eigen::MatrixXd> d2w_dphidphi () {
    std::vector<Eigen::MatrixXd> out;
    // alpha_w_C_denominator = sqrt(phi_long.dot(Dw * phi_long)) = tmp^(1/2);
    Eigen::MatrixXd outtmp(kwopt, kwopt);
    
    Eigen::VectorXd Dwphi;
    Eigen::MatrixXd outlarge(kwp, kwp);
    double tmp1;
    double tmp2;
    Eigen::MatrixXd m1(kwp, kwp);
    Eigen::MatrixXd m2(kwp, kwp);
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

  // *** GRADIENT ***
  Eigen::VectorXd gr_alpha_f () {
    Eigen::VectorXd tmp(kE);
    for (int i = 0; i < M; i++) tmp.segment(i*kEp, kEp) = smoothing_f(i) * Sf_list.at(i) * alpha_f.segment(i*kEp, kEp);
    Eigen::VectorXd out = - dmu_df_mat.transpose() * dlogdensity_dmu_vec + tmp;
    return out;
  }
  Eigen::VectorXd gr_betaR () {
    Eigen::VectorXd out = - dmu_dbetaR_mat.transpose() * dlogdensity_dmu_vec; // + smoothing * betaR;
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      // for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      // out += smoothing(i) * betaRi;
      for (int j = 0; j < ki; j++) out(begin + j) += smoothing(i) * betaR(begin + j);
      begin += ki;
    }
    return out;
  }
  Eigen::VectorXd gr_betaF () {
    Eigen::VectorXd out = - dmu_dbetaF_mat.transpose() * dlogdensity_dmu_vec;
    return out;
  }

  Eigen::VectorXd gr_alpha_w () {
    Eigen::VectorXd gr_pen_w(kw);
    for (int i = 0; i < M; i++) gr_pen_w.segment(i*kwp, kwp) = smoothing_w(i) * Sw_large * alpha_w_C.segment(i*kwp, kwp);
    Eigen::VectorXd out = - dmu_dw_mat.transpose() * dlogdensity_dmu_vec + gr_pen_w;
    return out;
  }
  Eigen::VectorXd gr_phi () {
    Eigen::VectorXd out = dw_dphi_mat.transpose() * gr_alpha_w_vec;
    return out;
  }
  Eigen::VectorXd gr_log_smoothing_f () {
    Eigen::VectorXd out(M);
    for (int i = 0; i < M; i++) out(i) = 0.5 * smoothing_f(i) * alpha_f.segment(i*kEp, kEp).dot(Sf_list.at(i) * alpha_f.segment(i*kEp, kEp)) - 0.5 * (kEp-1);
    return out;
  }
  Eigen::VectorXd gr_log_smoothing_w () {
    Eigen::VectorXd out(M);
    for (int i = 0; i < M; i++) out(i) = 0.5 * smoothing_w(i) * alpha_w_C.segment(i*kwp, kwp).dot(Sw_large * alpha_w_C.segment(i*kwp, kwp)) - 0.5 * (kwp-1-1);
    return out;
  }
  double gr_log_theta () {
    return -1.0 * theta * dlogdensity_dtheta_scalar;
  }
  Eigen::VectorXd gr_logsmoothing () {
    Eigen::VectorXd out(p);
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      Eigen::VectorXd betaRi(ki);
      for (int j = 0; j < ki; j++) betaRi(j) = betaR(begin + j);
      out(i) = 0.5 * smoothing(i) * betaRi.dot(betaRi) - 0.5*ki;
      begin += ki;
    }
    return out;
  }


  // *** Hessian ***
  Eigen::MatrixXd he_alpha_f () {
    Eigen::MatrixXd out1(kE, kE);
    Eigen::MatrixXd out2(kE, kE);
    Eigen::MatrixXd outpen(kE, kE);
    out1.setZero();
    out2.setZero();
    outpen.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));
    }

    for (int i = 0; i < M; i++) outpen.block(i*kEp, i*kEp, kEp, kEp) = smoothing_f(i)*Sf_list.at(i);
    
    return - out1 - out2 + outpen;
  }

  Eigen::MatrixXd I_alpha_f () { // hessian of negative likelihood without penalty 
    Eigen::MatrixXd out1(kE, kE);
    Eigen::MatrixXd out2(kE, kE);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_df_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_df_mat.row(i));
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd he_betaR () {
    Eigen::MatrixXd out1(kbetaR, kbetaR);
    Eigen::MatrixXd out2(kbetaR, kbetaR);
    Eigen::MatrixXd Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    out1.setZero();
    out2.setZero();
    Ones.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i));
    }
    int begin = 0;
    for (int i = 0; i < p; i++) {
      // Smooth Penalty
      int ki = static_cast<int>(r(i));
      for (int j = 0; j < ki; j++) Ones(begin + j, begin + j) = smoothing(i);
      begin += ki;
    }
    return - out1 - out2 + Ones;
  }

  Eigen::MatrixXd I_betaR () {  // hessian of negative likelihood without penalty 
    Eigen::MatrixXd out1(kbetaR, kbetaR);
    Eigen::MatrixXd out2(kbetaR, kbetaR);
    Eigen::MatrixXd Ones(kbetaR, kbetaR); // identity matrix with diagonal smoothing
    out1.setZero();
    out2.setZero();
    Ones.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaR_mat.row(i).transpose() * dlogmu_dbetaR_mat.row(i));
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd he_betaF () {
    Eigen::MatrixXd out1(kbetaF, kbetaF);
    Eigen::MatrixXd out2(kbetaF, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaF_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dbetaF_mat.row(i).transpose() * dlogmu_dbetaF_mat.row(i));
    }
    return - out1 - out2;
  }
  Eigen::MatrixXd he_alpha_w () {
    Eigen::MatrixXd out1(kw, kw);
    Eigen::MatrixXd out2(kw, kw);
    Eigen::MatrixXd outpen(kw, kw);
    out1.setZero();
    out2.setZero();
    outpen.setZero();
    Eigen::MatrixXd tmp(kw, kw);
    
    for (int i = 0; i < n; i++) {
      tmp.setZero();
      for (int j = 0; j < M; j++) tmp.block(j*kwp, j*kwp, kwp, kwp) = (BsplinevecCon2nd(E(i, j), knots_f_list.at(j), 4, Zf_list.at(j)).dot(alpha_f.segment(j*kEp, kEp))) * B_inner_list.at(j).row(i).transpose() * B_inner_list.at(j).row(i);
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dw_mat.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i) * tmp);
    }
    for (int j = 0; j < M; j++) outpen.block(j*kwp, j*kwp, kwp, kwp) = smoothing_w(j)*Sw_large;
    return - out1 - out2 + outpen;
  }



  Eigen::MatrixXd he_phi () {
    Eigen::MatrixXd out1 = dw_dphi_mat.transpose() * he_alpha_w_mat * dw_dphi_mat;
    Eigen::MatrixXd out2(kwopt, kwopt);
    out2.setZero();
    for (int s = 0; s < kw; s++) {
      out2 = out2 + gr_alpha_w_vec(s) * d2w_dphidphi_list.at(s);
    }
    return out1 + out2;
  }

  Eigen::MatrixXd I_alpha_w () { // hessian of negative likelihood without penalty  
    Eigen::MatrixXd out1(kw, kw);
    Eigen::MatrixXd out2(kw, kw);
    out1.setZero();
    out2.setZero();
    Eigen::MatrixXd tmp(kw, kw);
    
    for (int i = 0; i < n; i++) {
      tmp.setZero();
      for (int j = 0; j < M; j++) tmp.block(j*kwp, j*kwp, kwp, kwp) = (BsplinevecCon2nd(E(i, j), knots_f_list.at(j), 4, Zf_list.at(j)).dot(alpha_f.segment(j*kEp, kEp))) * B_inner_list.at(j).row(i).transpose() * B_inner_list.at(j).row(i);
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_dw_mat.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i) * tmp);
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd I_phi () { // hessian of negative likelihood without penalty  
    Eigen::MatrixXd out1 = dw_dphi_mat.transpose() * I_alpha_w_mat * dw_dphi_mat;
    Eigen::MatrixXd out2(kwopt, kwopt);
    out2.setZero();
    for (int s = 0; s < kw; s++) {
      out2 = out2 + gr_alpha_w_vec(s) * d2w_dphidphi_list.at(s);
    }
    return out1 + out2;
  }

  Eigen::MatrixXd he_alpha_f_alpha_w () {
    Eigen::MatrixXd out1(kE, kw);
    Eigen::MatrixXd out2(kE, kw);
    out1.setZero();
    out2.setZero();
    Eigen::MatrixXd tmp(kE, kw);
    
    for (int i = 0; i < n; i++) {
      tmp.setZero();
      for (int j = 0; j < M; j++) tmp.block(j*kEp, j*kwp, kEp, kwp) = BsplinevecCon1st(E(i, j), knots_f_list.at(j), 4, Zf_list.at(j)) * B_inner_list.at(j).row(i);
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dw_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * (mu(i) * dlogmu_df_mat.row(i).transpose() * dlogmu_dw_mat.row(i) + mu(i)*tmp);
    }
    return - out1 - out2;
  }

  Eigen::MatrixXd he_alpha_f_phi () {
    Eigen::MatrixXd out = he_alpha_f_alpha_w_mat * dw_dphi_mat;
    return out;
  }
  Eigen::MatrixXd he_alpha_f_betaF () {
    Eigen::MatrixXd out1(kE, kbetaF);
    Eigen::MatrixXd out2(kE, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    }
    return - out1 - out2;
  }
  Eigen::MatrixXd he_alpha_f_betaR () {
    Eigen::MatrixXd out1(kE, kbetaR);
    Eigen::MatrixXd out2(kE, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_df_mat.row(i).transpose() * dmu_dbetaR_mat.row(i);
    }
    return - out1 - out2;
  }
  Eigen::MatrixXd he_phi_betaR () {
    Eigen::MatrixXd out1(kwopt, kbetaR);
    Eigen::MatrixXd out2(kwopt, kbetaR);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dw_dphi_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dbetaR_mat.row(i));
      out2 += dw_dphi_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dw_mat.row(i).transpose() * dmu_dbetaR_mat.row(i));
    }
    return - out1 - out2;
  }
  Eigen::MatrixXd he_phi_betaF () {
    Eigen::MatrixXd out1(kwopt, kbetaF);
    Eigen::MatrixXd out2(kwopt, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += dw_dphi_mat.transpose() * (d2logdensity_dmudmu_vec(i) * dmu_dw_mat.row(i).transpose() * dmu_dbetaF_mat.row(i));
      out2 += dw_dphi_mat.transpose() * (dlogdensity_dmu_vec(i) * dlogmu_dw_mat.row(i).transpose() * dmu_dbetaF_mat.row(i));
    }
    return - out1 - out2;
  }
  Eigen::MatrixXd he_betaR_betaF () {
    Eigen::MatrixXd out1(kbetaR, kbetaF);
    Eigen::MatrixXd out2(kbetaR, kbetaF);
    out1.setZero();
    out2.setZero();
    for (int i = 0; i < n; i++) {
      out1 += d2logdensity_dmudmu_vec(i) * dmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
      out2 += dlogdensity_dmu_vec(i) * dlogmu_dbetaR_mat.row(i).transpose() * dmu_dbetaF_mat.row(i);
    }
    return - out1 - out2;
  }


  Eigen::MatrixXd he_alpha_f_log_smoothing_f () {
    Eigen::MatrixXd out(kE, M);
    out.setZero();
    Eigen::VectorXd tmp(kE);
    for (int i = 0; i < M; i++) {
      tmp.setZero();
      tmp.segment(i*kEp, kEp) = smoothing_f(i) * Sf_list.at(i) * alpha_f.segment(i*kEp, kEp);
      out.col(i) = tmp;
    }
    
    return out;
  }

  Eigen::MatrixXd he_phi_log_smoothing_w () {
    Eigen::MatrixXd out(kwopt, M);

    Eigen::VectorXd tmp(kw);
    
    for (int i = 0; i < M; i++) {
      tmp.setZero();
      tmp.segment(i*kwp, kwp) = smoothing_w(i) * Sw_large * alpha_w_C.segment(i*kwp, kwp);
      out.col(i) =  dw_dphi_mat.transpose() * tmp;
    }
    
    return out;
  }

  Eigen::MatrixXd he_betaR_logsmoothing () {
    Eigen::MatrixXd out(kbetaR, p);
    out.setZero();
    int begin = 0;
    for (int i = 0; i < p; i++) {
      int ki = static_cast<int>(r(i));
      for (int j = 0; j < ki; j++) out(begin + j, i) = smoothing(i) * betaR(begin + j);
      begin += ki;
    }
    return out;
  }
  Eigen::VectorXd he_alpha_f_log_theta () {
    // he_alpha_f_theta = dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec;
    return -1.0*dmu_df_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }
  Eigen::VectorXd he_phi_log_theta () {
    // he_alpha_w_theta = dmu_dw_mat.transpose() * d2logdensity_dmudtheta_vec;
    return -1.0*theta * dw_dphi_mat.transpose() * ( dmu_dw_mat.transpose() * d2logdensity_dmudtheta_vec );
  }
  Eigen::VectorXd he_betaR_log_theta () {
    return -1.0*dmu_dbetaR_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }
  Eigen::VectorXd he_betaF_log_theta () {
    return -1.0*dmu_dbetaF_mat.transpose() * d2logdensity_dmudtheta_vec * theta;
  }





  // *********** LAML ***********
  double logdetH05() {
    // double out = 0.5 * log(he_s_u_mat.determinant());
    // return out;

    double logdetH05 = 0.0;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(he_s_u_mat);
    Eigen::MatrixXd LU = lu.matrixLU();
    // double c = lu.permutationP().determinant(); // -1 or 1
    double lii;
    for (int i = 0; i < LU.rows(); i++) {
      lii = LU(i,i);

      logdetH05 += log(abs(lii));
    }
    // logdetH05 += log(c);
    return logdetH05/2.0;

  }
};





// optimize alpha_f and betaF for a given phi, log_smoothing_f, log_smoothing_w, and log_theta
// Use newton method with eigenvalue modification
void PL(Model& modelobj, bool verbose){
    int maxitr = 50, itr = 0;
    const double eps = 1e-05;
    double mineig = 1e-03; // minimum eigenvalue of Hessian, to ensure it is PD
    int maxstephalve = 50, stephalve = 0;
    int resetitr = 0, maxreset = 1; // if step is nan, reset coefficients as 0. only once.
    int additr = 50; // allow further iterations after resetting coefficients as 0

    Eigen::VectorXd alpha_f = modelobj.alpha_f;
    Eigen::VectorXd betaR = modelobj.betaR;
    Eigen::VectorXd betaF = modelobj.betaF;
    
    int kE = modelobj.kE;
    int kw = modelobj.kw;
    int kwopt = modelobj.kwopt;
    int kbetaR = modelobj.kbetaR;
    int kbetaF = modelobj.kbetaF;
    int converge = 0;

    int paraSize = kE+kbetaR+kbetaF;
    // Optimize ALPHA_F
    double u;
    double u_tmp;

  

    // update steps
    Eigen::VectorXd step(paraSize);
    step.setZero();

    Eigen::MatrixXd H(paraSize, paraSize);
    Eigen::VectorXd g(paraSize);

    g.setZero();
    H.setZero();

    // START DEFINE lanczos algorithm for smallest eigenvalue
    // Code from https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample4_use_Eigen_library.cpp
    // the matrix-vector multiplication routine
    auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
      auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
      auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

      // eigen_out = H * eigen_in; // Easy version
      eigen_out.noalias() += H * eigen_in; // Efficient version
    };

    LambdaLanczos<double> engine(mv_mul, paraSize, false, 1); // Find 1 minimum eigenvalue
    std::vector<double> smallest_eigenvalues;
    std::vector<std::vector<double>> smallest_eigenvectors;
    double smallest_eigval; // smallest eigenvalue
    // END DEFINE lanczos for smallest eigenvalue
    
    // eigen decomposition
    Eigen::VectorXd eigvals(paraSize);
    eigvals.setZero();
    Eigen::VectorXd invabseigvals(paraSize);
    invabseigvals.setZero();
    // double eigval;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(H,false); // Only values, not vectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(H,true); // Both values and vectors

    double delta = 1.0; // for New Q-Newton method
    double g_norm;

    for (int i = 0; i < paraSize; i++) g(i) = 1. + eps;

    while (itr < maxitr) {
      itr++;
      modelobj.derivative_f();

      g = modelobj.gr_inner_vec;
      g_norm = g.norm();
      if (g_norm < eps) break;
      modelobj.NegativeLogLikelihood();
      u = modelobj.NegLogL;

      H = modelobj.he_inner_mat;
      

      engine.run(smallest_eigenvalues, smallest_eigenvectors);
      smallest_eigval = smallest_eigenvalues[0]; // the smallest eigenvalue
      if ((smallest_eigval < 1e-2) || std::isnan(smallest_eigval)) {
        // Do Q-Newton's Step
        eigvec.compute(H); // Compute eigenvalues and vectors
        eigvals = eigvec.eigenvalues().array();

        if (abs(eigvals.prod()) < 1e-3) {
          for (int iii = 0; iii < paraSize; iii++) eigvals(iii) += delta*g_norm;
        }

        // for (int i = 0; i < paraSize; i++) invabseigvals(i) = 1. / max(abs(eigvals(i)), mineig); // flip signs
        for (int i = 0; i < paraSize; i++) invabseigvals(i) = 1. / abs(eigvals(i)); // flip signs
        step = eigvec.eigenvectors() * (invabseigvals.asDiagonal()) * (eigvec.eigenvectors().transpose()) * g;
      } else {
        // smallest eigenvalue > 1e-3
        // regular Newton's step
        // step = H.llt().solve(g);
        step = H.ldlt().solve(g);
      }

      // check nan in step
      // Really needed
      if(hasNaN(step)){
        if (resetitr < maxreset){
          resetitr++;
          alpha_f.setZero(); // reset alpha_f
          betaR.setZero();
          betaF.setZero();
          modelobj.setAlphaF(alpha_f);
          modelobj.setBetaR(betaR);
          modelobj.setBetaF(betaF);
          if(verbose) std::cout << "reset alpha_f and betaF as 0" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          continue; // do next iteration
        } else {
          converge = 99;
          break;
        }
      }

      alpha_f -= step.segment(0, kE);
      betaR -= step.segment(kE, kbetaR);
      betaF -= step.segment(kE+kbetaR, kbetaF);

      modelobj.setAlphaF(alpha_f);
      modelobj.setBetaR(betaR);
      modelobj.setBetaF(betaF);

      modelobj.NegativeLogLikelihood();

      u_tmp = modelobj.NegLogL;

      // halving if objective function increase.
      stephalve = 0;
      while ((u_tmp > u + 1e-8) & (stephalve < maxstephalve)){
        stephalve++;
        step /= 2.;

        alpha_f += step.segment(0, kE);
        betaR += step.segment(kE, kbetaR);
        betaF += step.segment(kE+kbetaR, kbetaF);
        modelobj.setAlphaF(alpha_f);
        modelobj.setBetaR(betaR);
        modelobj.setBetaF(betaF);

        modelobj.NegativeLogLikelihood();
        u_tmp = modelobj.NegLogL;
      }


      stephalve = 0;
      // Check feasibility of step. If u is nan then we went too far;
      // halve the step and try again
      while (std::isnan(u_tmp) & (stephalve < maxstephalve)) {
        stephalve++;
        step /= 2.; // This is still the step from the previous iteration

        alpha_f += step.segment(0, kE);
        betaR += step.segment(kE, kbetaR);
        betaF += step.segment(kE+kbetaR, kbetaF);
        modelobj.setAlphaF(alpha_f);
        modelobj.setBetaR(betaR);
        modelobj.setBetaF(betaF);
        modelobj.NegativeLogLikelihood();
        u_tmp = modelobj.NegLogL;
      }

      stephalve = 0;
      // if (stephalve > 0) std::cout << "Performed " << stephalve << " iterations of step-halving." << std::endl;
      if (std::isnan(u_tmp)) {
        // Step-halving didn't work
        // std::cout << "AlphaF: Step-halving failed with nan function value. Returning failure." << std::endl;
        converge = 99;
        break;
      }
    }
    if(itr == maxitr){
      // std::cout << "AlphaF: Newton method for updating alpha fails" << std::endl;
      converge = 99;
    }

    // if(verbose) std::cout << "-- AlphaF Gradient Max: " << g.maxCoeff() << std::endl;

    modelobj.derivative_coef();
    modelobj.NegativeLogLikelihood();

    Eigen::VectorXd gr_PL = modelobj.gr_phi_vec;
    Eigen::MatrixXd he_PL(kwopt, kwopt);
    
    // OLD: 
    // Eigen::MatrixXd mat1 = modelobj.he_alpha_f_mat.ldlt().solve(modelobj.he_alpha_f_phi_mat);
    // he_PL = modelobj.he_phi_mat - mat1.transpose() * modelobj.he_alpha_f_phi_mat;

    // NEW: 
    modelobj.derivative_f();
    Eigen::MatrixXd mat_tmp_PL(kE+kbetaF+kbetaR,kwopt);
    mat_tmp_PL.block(0, 0, kE, kwopt) = modelobj.he_alpha_f_phi_mat;
    mat_tmp_PL.block(kE, 0, kbetaR, kwopt) = modelobj.he_phi_betaR_mat.transpose();
    mat_tmp_PL.block(kE+kbetaR, 0, kbetaF, kwopt) = modelobj.he_phi_betaF_mat.transpose();
    Eigen::MatrixXd mat1 = modelobj.he_inner_mat.ldlt().solve(mat_tmp_PL);
    he_PL = modelobj.he_phi_mat - mat1.transpose() * mat_tmp_PL;

    modelobj.PL_gradient = gr_PL;
    modelobj.PL_hessian = he_PL;
    modelobj.converge = converge; // 0: converge. 99: not converge
}



void Inner(Model& modelobj, bool verbose) {
    // newton method
    int maxitr = 50, itr = 0;
    const double eps = 1e-03;
    const double largereps = 1e-02;
    double mineig = 1e-03; // minimum eigenvalue of Hessian, to ensure it is PD
    // double mineig = 1e-02; // minimum eigenvalue of Hessian, to ensure it is PD
    int maxstephalve = 50, stephalve = 0;
    int maxErangehalve = 20;
    int stephalve_inner = 0;
    int resetitr = 0, maxreset = 1; // if step is nan or always diverge, reset coefficients as 0. Only reset once
    int additr = 50; // allow further 20 iteraions after resetting coeffiicents as 0.

    // check non-moving step following https://github.com/awstringer1/varcomptest/blob/main/src/reml-ad.cpp
    int maxnonmovingsteps = 5; // Maximum number of iterations for which we will tolerate no movement. 5 in awstringer1/varcomptest
    double stepeps = 1e-12; // If max(abs(step)) < stepeps then we say the iteration resulted in no movement.
    int stepcounter = 0; // Count the number of non-moving steps

    Eigen::VectorXd phi = modelobj.phi;

    int kE = modelobj.kE;
    int kw = modelobj.kw;
    int kwopt = modelobj.kwopt;
    int M = modelobj.M;
    int kwp = modelobj.kwp;
    int kwpopt = modelobj.kwpopt;
    int n = modelobj.n;

    // catch double PL.fn
    double s;
    double s_tmp;

    // update steps
    Eigen::VectorXd step(kwopt);
    step.setZero();

    
    Eigen::MatrixXd H(kwopt, kwopt);
    Eigen::VectorXd g(kwopt);

    g.setZero();
    H.setZero();


    // START DEFINE lanczos algorithm for smallest eigenvalue
    // Code from https://github.com/mrcdr/lambda-lanczos/blob/master/src/samples/sample4_use_Eigen_library.cpp
    // the matrix-vector multiplication routine
    auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
      auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
      auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

      // eigen_out = H * eigen_in; // Easy version
      eigen_out.noalias() += H * eigen_in; // Efficient version
    };

    LambdaLanczos<double> engine(mv_mul, kwopt, false, 1); // Find 1 minimum eigenvalue
    std::vector<double> smallest_eigenvalues;
    std::vector<std::vector<double>> smallest_eigenvectors;
    double smallest_eigval; // smallest eigenvalue
    // END DEFINE lanczos for smalles eigenvalue
    

    // eigen decomposition
    Eigen::VectorXd eigvals(kwopt);
    eigvals.setZero();
    Eigen::VectorXd invabseigvals(kwopt);
    invabseigvals.setZero();
    // double eigval;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(H,false); // Only values, not vectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigvec(H,true); // Both values and vectors

    // check range of E
    Eigen::VectorXd phi_long(kw);
    Eigen::VectorXd phiKa(kw-M);
    Eigen::VectorXd alpha_w_C(kw);
    Eigen::MatrixXd E(n ,M);
    int kmin;
    int kmax;
    Eigen::VectorXd krange(M);
    int krangemin = 4;
    // int krangemin = 3;
    bool krangewarning = false;
    std::vector<Eigen::MatrixXd> B_inner_list = modelobj.getB_inner_list();
    Eigen::MatrixXd Dw = modelobj.getDw();
    std::vector<Eigen::VectorXd> knots_f_list = modelobj.getknots_f_list();

    Eigen::MatrixXd K = modelobj.K;
    Eigen::VectorXd a = modelobj.a;
    double alpha_w_C_denominator;
    
    double delta = 1.0; // for New Q-Newton method
    double g_norm;

    // initialize gradient
    for (int i = 0; i < kwopt; i++) g(i) = 1. + eps;

    if(verbose) std::cout << "* Start optimize profile likelihood" << std::endl;

    // start newton's method
    while (itr < maxitr) {
      itr++;
      // modelobj.NegativeLogLikelihood();

      resetcon_label: // reset phi as all zero. If the current phi always leads to the divergence in updating alpha_f
      // update alpha_f
      PL(modelobj, verbose);
      g = modelobj.PL_gradient;
      g_norm = g.norm();

      if (g_norm < eps) break;
      // std::cout << "iter" << itr << std::endl;
      // std::cout << "g_norm" << g_norm << std::endl;

      s = modelobj.NegLogL;
      // std::cout << "s: " << s << std::endl;
      H = modelobj.PL_hessian;

      engine.run(smallest_eigenvalues, smallest_eigenvectors);
      smallest_eigval = smallest_eigenvalues[0]; // the smallest eigenvalue
      if ((smallest_eigval < 1e-2) || std::isnan(smallest_eigval)) {
        // Do Q-Newton's Step
        eigvec.compute(H); // Compute eigenvalues and vectors
        eigvals = eigvec.eigenvalues().array();
        if (abs(eigvals.prod()) < 1e-3) {
          for (int iii = 0; iii < kwopt; iii++) eigvals(iii) += delta*g_norm;
        }
        // for (int i = 0; i < kwopt; i++) invabseigvals(i) = 1. / max(abs(eigvals(i)), mineig); // flip signs
        for (int i = 0; i < kwopt; i++) invabseigvals(i) = 1. / abs(eigvals(i)); // flip signs
        // std::cout << "invabseigvals max" << invabseigvals.maxCoeff() << std::endl;
        step = eigvec.eigenvectors() * (invabseigvals.asDiagonal()) * (eigvec.eigenvectors().transpose()) * g;
      } else {
        // smallest eigenvalue > 1e-3
        // regular Newton's step
        // step = H.llt().solve(g);
        step = H.ldlt().solve(g);
      }
      // check nan in step
      // NOT really needed. checking here to align with alpha_f.
      if(hasNaN(step)){
        if (resetitr < maxreset){
          resetitr++;
          phi.setZero(); // reset alpha_f
          modelobj.setPhi(phi);
          if(verbose) std::cout << "reset phi as 0 because of nan step" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          continue; // do next iteration
        } else {
          break;
        }
      }

      phi -= step;

      // ****** start checking range of E
      phiKa = K*phi + a;
      for (int i = 0; i < M; i++) {
        phi_long(i*kwp) = 1.0;
        phi_long.segment(i*kwp + 1, kwp-1) = phiKa.segment(i*(kwp-1), (kwp-1));
      }
      for (int i = 0; i < M; i++) {
        alpha_w_C_denominator = sqrt(phi_long.segment(i*kwp, kwp).dot(Dw * phi_long.segment(i*kwp, kwp)));
        alpha_w_C.segment(i*kwp, kwp) = phi_long.segment(i*kwp, kwp) / alpha_w_C_denominator;
      }
      for (int i = 0; i < M; i++) {
        E.col(i) = B_inner_list.at(i) * alpha_w_C.segment(i*kwp, kwp);
      }

      for (int i = 0; i < M; i++) {
        kmin = knotindexEigen(E.col(i).minCoeff(), knots_f_list.at(i));
        kmax = knotindexEigen(E.col(i).maxCoeff(), knots_f_list.at(i));  
        krange(i) = kmax - kmin;
      }
      
      stephalve = 0;

      while ((krange.minCoeff() < krangemin) & (stephalve < maxErangehalve)){
        stephalve++;
        // std::cout << "krange halving" << std::endl;
        for (int mi = 0; mi < M; mi++) {
          if(krange(mi) < krangemin) {
            step.segment(mi*kwpopt, kwpopt) = step.segment(mi*kwpopt, kwpopt)/2.;
            phi.segment(mi*kwpopt, kwpopt) = phi.segment(mi*kwpopt, kwpopt) + step.segment(mi*kwpopt, kwpopt); 
            // std::cout << "krange halving: " << mi << std::endl;
          }
        }
        
        // step /= 2.;
        // phi += step;

        phiKa = K*phi + a;
        for (int i = 0; i < M; i++) {
          phi_long(i*kwp) = 1.0;
          phi_long.segment(i*kwp + 1, kwp-1) = phiKa.segment(i*(kwp-1), (kwp-1));
        }
        for (int i = 0; i < M; i++) {
          alpha_w_C_denominator = sqrt(phi_long.segment(i*kwp, kwp).dot(Dw * phi_long.segment(i*kwp, kwp)));
          alpha_w_C.segment(i*kwp, kwp) = phi_long.segment(i*kwp, kwp) / alpha_w_C_denominator;
        }
        
        for (int i = 0; i < M; i++) {
          E.col(i) = B_inner_list.at(i) * alpha_w_C.segment(i*kwp, kwp);
        }

        for (int i = 0; i < M; i++) {
          kmin = knotindexEigen(E.col(i).minCoeff(), knots_f_list.at(i));
          kmax = knotindexEigen(E.col(i).maxCoeff(), knots_f_list.at(i));  
          krange(i) = kmax - kmin;
        }


        // if(verbose) std::cout << "E range krange" << krange.minCoeff() << std::endl;
        // std::cout << "halving" << std::endl;
        // std::cout << "alpha_f" << convertToDouble(modelobj.alpha_f) << std::endl;
        // Example: getLAML(c(3, 3, 3)) ## fails when kE = kw = 20. Nt = 1000. wl <- function(l) dnorm(l, mean = 10, sd = 10)/wl_de
      }
      if (krange.minCoeff() < krangemin) {
        krangewarning = true;
        if(verbose) std::cout << "E range krange: " << krange.transpose() << " < " << krangemin << std::endl;
        // std::cout << "Range of weighted exposure is small. Consider increasing kE and resetting starting values." << std::endl;
      }
      // finish checking ******
      modelobj.setPhi(phi);

      PL(modelobj, verbose);

      // halving if the optimization for alpha_f fails
      stephalve = 0;
      while ((modelobj.converge != 0) & (stephalve < maxstephalve)){
        stephalve++;
        // std::cout << "divergence halving" << std::endl;
        step /= 2.;
        phi += step;
        modelobj.setPhi(phi);
        PL(modelobj, verbose);
      }
      if (modelobj.converge != 0) {
        if (resetitr < maxreset){
          resetitr++;
          phi.setZero(); // reset alpha_f
          modelobj.setPhi(phi);
          if(verbose) std::cout << "reset phi as 0 because of divergence" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          goto resetcon_label; // do next iteration
        } else {
          std::cout << "Optimization for alpha_f fails" << std::endl;
          break;
        }
      }

      s_tmp = modelobj.NegLogL;

      // halving if objective function increase.
      stephalve = 0;
      while ((s_tmp > s + 1e-8) & (stephalve < maxstephalve)){
        stephalve++;
        // std::cout << "increase halving" << std::endl;
        step /= 2.;
        phi += step;
        modelobj.setPhi(phi);
        PL(modelobj, verbose);
        // when dealing with increase: halving if the optimization for alpha_f fails
        stephalve_inner = 0;
        while ((modelobj.converge != 0) & (stephalve_inner < maxstephalve)){
          stephalve_inner++;
          step /= 2.;
          phi += step;
          modelobj.setPhi(phi);
          PL(modelobj, verbose);
        }
        if (modelobj.converge != 0) {
          if (resetitr < maxreset){
            resetitr++;
            phi.setZero(); // reset alpha_f
            modelobj.setPhi(phi);
            if(verbose) std::cout << "reset phi as 0 because of divergence" << std::endl;
            itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
            goto resetcon_label; // do next iteration
          } else {
            std::cout << "Optimization for alpha_f fails" << std::endl;
            break;
          }
        }
        s_tmp = modelobj.NegLogL;
      }

      stephalve = 0;
      // Check feasibility of step. If u is nan then we went too far;
      // halve the step and try again
      while (std::isnan(s_tmp) & (stephalve < maxstephalve)) {
        stephalve++;
        step /= 2.; // This is still the step from the previous iteration
        phi += step;
        modelobj.setPhi(phi);
        PL(modelobj, verbose);
        // when dealing with NaN: halving if the optimization for alpha_f fails
        stephalve_inner = 0;
        while ((modelobj.converge != 0) & (stephalve_inner < maxstephalve)){
          stephalve_inner++;
          step /= 2.;
          phi += step;
          modelobj.setPhi(phi);
          PL(modelobj, verbose);
        }
        if (modelobj.converge != 0) {
          if (resetitr < maxreset){
            resetitr++;
            phi.setZero(); // reset alpha_f
            modelobj.setPhi(phi);
            if(verbose) std::cout << "reset phi as 0 because of divergence" << std::endl;
            itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
            goto resetcon_label; // do next iteration
          } else {
            break; // break the current innner while. The next code to run: if (std::isnan(s_tmp)) {...}
          }
        }
        s_tmp = modelobj.NegLogL;
      }

      if (std::isnan(s_tmp)) {
        if (resetitr < maxreset){
          resetitr++;
          phi.setZero(); // reset alpha_f
          modelobj.setPhi(phi);
          if(verbose) std::cout << "reset phi as 0 because of nan function" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          goto resetcon_label; // do next iteration
        } else {
          // Step-halving didn't work
          std::cout << "Phi: Step-halving failed with nan function value. Returning failure." << std::endl;
          break;
        }
      }
      // Count the number of iterations where we didn't move; if too many, we got stuck.
      // a part of the code follows https://github.com/awstringer1/varcomptest/blob/main/src/reml-ad.cpp
      if (step.lpNorm<Eigen::Infinity>() < stepeps) {
        stepcounter++;
        if (stepcounter > maxnonmovingsteps) {
          if (resetitr < maxreset){
            resetitr++;
            phi.setZero(); // reset alpha_f
            modelobj.setPhi(phi);
            if(verbose) std::cout << "reset phi as 0 because of non-moving" << std::endl;
            itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
            goto resetcon_label; // do next iteration
          } else {
            std::cout << "The algorithm hasn't moved for " << stepcounter << " steps; terminating. Please check the answer." << std::endl;
            break;
          }
        }
      } else {
        stepcounter = 0; // if move
      }

      // The last reset
      if((itr == maxitr) & (g.norm() >= largereps)) {
        if (resetitr < maxreset){
          resetitr++;
          phi.setZero(); // reset alpha_f
          modelobj.setPhi(phi);
          if(verbose) std::cout << "reset phi as 0. The last one" << std::endl;
          itr = std::max(0, itr - additr); // itr - additr if itr > additr. allow further additr iterations after resetting.
          goto resetcon_label; // do next iteration
        } else {
          if (verbose) {
            std::cout << "Newton method for updating weight function might fail. Profile Likelihood Gradient Max: " << g.maxCoeff() << std::endl;
          } else {
            // report less information if verbose is false
            if (g.maxCoeff() >= 1e-2) {
             std::cout << "Newton method for updating weight function might fail. Profile Likelihood Gradient Max: " << g.maxCoeff() << std::endl;
            }
          }
          break;
        }
      }
      // initilized count
      stephalve = 0;
    }
    if(krangewarning) std::cout << "Range of weighted exposure is small. Consider increasing kE or resetting starting values. Profile Likelihood Gradient Max: " << g.maxCoeff() << std::endl;
    if(verbose) std::cout << "* Finish middle opt. Profile Likelihood Gradient Max: " << g.maxCoeff() << std::endl;
}



#endif