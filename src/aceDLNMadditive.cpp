/** Include **/
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

#include <random> // for generating samples from standard normal distribution


#include "aceDLNMcppadClass.hpp"
#include "aceDLNMeigenClass.hpp"

// #include <LBFGSB.h>
// using namespace LBFGSpp;
// TODO: use https://github.com/yixuan/LBFGSpp


std::vector<Mat> convertListToVectorMat(List rList) {

  int n = rList.size();

  std::vector<Mat> matrixVector;


  for (int i = 0; i < n; ++i) {
    NumericMatrix rMatrix = as<NumericMatrix>(rList[i]);
    Eigen::Map<Eigen::MatrixXd> eigenMatrix(as<Eigen::Map<Eigen::MatrixXd>>(rMatrix));
    matrixVector.push_back(eigenMatrix.cast<Scalar>());
  }

  return matrixVector;
}

std::vector<Eigen::MatrixXd> convertListToVectorEigenMatrixXd(List rList) {

  int n = rList.size();

  std::vector<Eigen::MatrixXd> matrixVector;


  for (int i = 0; i < n; ++i) {
    NumericMatrix rMatrix = as<NumericMatrix>(rList[i]);
    Eigen::Map<Eigen::MatrixXd> eigenMatrix(as<Eigen::Map<Eigen::MatrixXd>>(rMatrix));
    matrixVector.push_back(eigenMatrix);
  }

  return matrixVector;
}

std::vector<Eigen::VectorXd> convertListToVectorEigenVectorXd(List rList) {

  int n = rList.size();
  std::vector<Eigen::VectorXd> vectorList;

  for (int i = 0; i < n; ++i) {
    NumericVector rVector = as<NumericVector>(rList[i]);
    Eigen::Map<Eigen::VectorXd> eigenVector(as<Eigen::Map<Eigen::VectorXd>>(rVector));
    vectorList.push_back(eigenVector);
  }

  return vectorList;
}


std::vector<Vec> convertListToVectorVec(List rList) {

  int n = rList.size();
  std::vector<Vec> vectorList;

  for (int i = 0; i < n; ++i) {
    NumericVector rVector = as<NumericVector>(rList[i]);
    Eigen::Map<Eigen::VectorXd> eigenVector(as<Eigen::Map<Eigen::VectorXd>>(rVector));
    vectorList.push_back(eigenVector.cast<Scalar>());
  }

  return vectorList;
}


struct LAMLResult {
    double fn;
    Eigen::VectorXd gradient;
};


CppAD::ADFun<double> LAMLTape(Model& modelobj, Modelcppad& modelcppadobj){
    // if(!modelcppadobj.ifhastape) {
      int kE = modelobj.kE;
      int kw = modelobj.kw;
      int kwopt = modelobj.kwopt;
      int kbetaR = modelobj.kbetaR;
      int kbetaF = modelobj.kbetaF;
      int p = modelobj.p;
      int M = modelobj.M;

      Eigen::VectorXd R_alpha_f = modelobj.alpha_f;
      Eigen::VectorXd R_betaR = modelobj.betaR;
      Eigen::VectorXd R_betaF = modelobj.betaF;
      Eigen::VectorXd R_phi = modelobj.phi;

      double R_log_theta = modelobj.log_theta;
      Eigen::VectorXd R_log_smoothing_f = modelobj.log_smoothing_f;
      Eigen::VectorXd R_log_smoothing_w = modelobj.log_smoothing_w;
      Eigen::VectorXd R_logsmoothing = modelobj.logsmoothing;

      Vec alpha_f = R_alpha_f.cast<Scalar>();
      Vec betaR = R_betaR.cast<Scalar>();
      Vec betaF = R_betaF.cast<Scalar>();
      Vec phi = R_phi.cast<Scalar>();

      Scalar log_theta = R_log_theta;
      Vec log_smoothing_f = R_log_smoothing_f.cast<Scalar>();
      Vec log_smoothing_w = R_log_smoothing_w.cast<Scalar>();
      Vec logsmoothing = R_logsmoothing.cast<Scalar>();

      // alpha_f, phi, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing


      Vec at(kE+kbetaR+kbetaF + kwopt + 1+2*M + p);
      at << alpha_f, phi, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing;


      // start reverse mode
      Vec result(1);
      CppAD::ADFun<double> gr;

      CppAD::Independent(at);
      modelcppadobj.setAlphaF(at.segment(0, kE));
      modelcppadobj.setPhi(at.segment(kE, kwopt));
      modelcppadobj.setBetaR(at.segment(kE+kwopt, kbetaR));
      modelcppadobj.setBetaF(at.segment(kE+kwopt+kbetaR, kbetaF));
      modelcppadobj.setLogTheta(at(kE+kbetaR+kbetaF+kwopt));
      modelcppadobj.setLogSmoothingF(at.segment(kE+kbetaR+kbetaF+kwopt+1, M));
      modelcppadobj.setLogSmoothingW(at.segment(kE+kbetaR+kbetaF+kwopt+1+M, M));
      modelcppadobj.setLogsmoothing(at.segment(kE+kbetaR+kbetaF+kwopt+1+2*M, p));


      modelcppadobj.derivative_coef();
      modelcppadobj.derivative_he();
      result(0) = modelcppadobj.logdetH05();

      gr.Dependent(at, result);
      // gr.optimize();
      // modelcppadobj.gr = gr;
      return gr;
      // END reverse



      // Forward mode 1
      // std::vector<Scalar> at_std(at.data(), at.data() + (kE+kbetaR+kbetaF + kw-1 + 3 + p));

      // CppAD::Independent(at_std);

      // std::vector<Scalar> result(1);
      // Vec at_eigen = Eigen::Map<Vec>(at_std.data(), (kE+kbetaR+kbetaF + kw-1 + 3 + p));
      // modelcppadobj.setAlphaF(at_eigen.segment(0, kE));
      // modelcppadobj.setPhi(at_eigen.segment(kE, kw-1));
      // modelcppadobj.setBetaR(at_eigen.segment(kE+kw-1, kbetaR));
      // modelcppadobj.setBetaF(at_eigen.segment(kE+kw-1+kbetaR, kbetaF));
      // modelcppadobj.setLogTheta(at_eigen(kE+kbetaR+kbetaF+kw-1));
      // modelcppadobj.setLogSmoothingF(at_eigen(kE+kbetaR+kbetaF+kw-1+1));
      // modelcppadobj.setLogSmoothingW(at_eigen(kE+kbetaR+kbetaF+kw-1+2));
      // modelcppadobj.setLogsmoothing(at_eigen.segment(kE+kbetaR+kbetaF+kw-1+3, p));

      // modelcppadobj.derivative_coef();
      // modelcppadobj.derivative_he();
      // result[0] = modelcppadobj.logdetH05();

      // CppAD::ADFun<double> gr(at_std, result);
      // END forward 1


      // Forward mode 2
      // CppAD::Independent(at);

      // Vec result(1);
      // modelcppadobj.setAlphaF(at.segment(0, kE));
      // modelcppadobj.setPhi(at.segment(kE, kw-1));
      // modelcppadobj.setBetaR(at.segment(kE+kw-1, kbetaR));
      // modelcppadobj.setBetaF(at.segment(kE+kw-1+kbetaR, kbetaF));
      // modelcppadobj.setLogTheta(at(kE+kbetaR+kbetaF+kw-1));
      // modelcppadobj.setLogSmoothingF(at(kE+kbetaR+kbetaF+kw-1+1));
      // modelcppadobj.setLogSmoothingW(at(kE+kbetaR+kbetaF+kw-1+2));
      // modelcppadobj.setLogsmoothing(at.segment(kE+kbetaR+kbetaF+kw-1+3, p));


      // modelcppadobj.derivative_coef();
      // modelcppadobj.derivative_he();
      // result(0) = modelcppadobj.logdetH05();

      // CppAD::ADFun<double> gr(at, result);
      // END forward 2


      // modelcppadobj.ifhastape = true;





    // }
}

LAMLResult LAML(Model& modelobj, Modelcppad& modelcppadobj) {
    LAMLResult result;
    double u_LAML;
    int kE = modelobj.kE;
    int kw = modelobj.kw;
    int kwopt = modelobj.kwopt;
    int kbetaR = modelobj.kbetaR;
    int kbetaF = modelobj.kbetaF;
    int p = modelobj.p;
    int M = modelobj.M;

    modelobj.derivative_coef();
    modelobj.derivative_he();
    modelobj.derivative_full();
    modelobj.NegativeLogLikelihood();

    Eigen::VectorXd alpha_f = modelobj.alpha_f;
    Eigen::VectorXd phi = modelobj.phi;
    Eigen::VectorXd betaR = modelobj.betaR;
    Eigen::VectorXd betaF = modelobj.betaF;

    Eigen::VectorXd gr_s_u_vec = modelobj.gr_s_u_vec;
    Eigen::VectorXd gr_s_par_vec = modelobj.gr_s_par_vec;
    Eigen::MatrixXd he_s_u_mat = modelobj.he_s_u_mat;
    Eigen::MatrixXd he_s_par_u_mat = modelobj.he_s_par_u_mat;


    double log_theta = modelobj.log_theta;
    Eigen::VectorXd log_smoothing_f = modelobj.log_smoothing_f;
    Eigen::VectorXd log_smoothing_w = modelobj.log_smoothing_w;
    Eigen::VectorXd logsmoothing = modelobj.logsmoothing;

    // First derivative of LAML

    CppAD::ADFun<double> cppadgr = LAMLTape(modelobj, modelcppadobj);



    Eigen::VectorXd at0(kE+kbetaR+kbetaF + kwopt + 1 + 2*M + p);
    at0 << alpha_f, phi, betaR, betaF, log_theta, log_smoothing_f, log_smoothing_w, logsmoothing;

    // reverse mode
    Eigen::VectorXd g_LAML(kE + kwopt + 1 + 2*M + kbetaR+kbetaF + p);
    g_LAML.setZero();
    g_LAML = cppadgr.Jacobian(at0);
    // END reverse mode


    // forward mode 1, 2
    // std::vector<Scalar> at0_std(at0.data(), at0.data() + (kE+kbetaR+kbetaF + kw-1 + 3 + p));
    // std::vector<double> dx((kE+kbetaR+kbetaF + kw-1 + 3 + p), 0.0);
    // std::vector<double> g_LAML_std((kE+kbetaR+kbetaF + kw-1 + 3 + p));
    //  std::vector<double> dy(1);
    // for (size_t i = 0; i < (kE+kbetaR+kbetaF + kw-1 + 3 + p); ++i) {
    //     dx[i] = 1.0;
    //     dy = modelcppadobj.gr.Forward(1, dx);
    //     g_LAML_std[i] = dy[0];
    //     dx[i] = 0.0;
    // }

    // Eigen::VectorXd g_LAML = Eigen::Map<Eigen::VectorXd>(g_LAML_std.data(), g_LAML_std.size());
    // END forward mode 1, 2




    u_LAML = modelobj.logdetH05() + modelobj.NegLogL - modelobj.n/2.0 * log(2*3.141592653589793238462643383279);


    // In R: grad[-(1:(kE+kw-1))] - H.full[-(1:(kE+kw-1)),(1:(kE+kw-1))] %*% as.vector(solve(H.alpha, grad[(1:(kE+kw-1))]))
    Eigen::VectorXd g1 = g_LAML.segment(0, kE+kwopt+kbetaR+kbetaF) + gr_s_u_vec;
    Eigen::VectorXd g2 = g_LAML.segment(kE+kwopt+kbetaR+kbetaF, 1+2*M+p) + gr_s_par_vec;
    Eigen::VectorXd gr = g2 - he_s_par_u_mat * he_s_u_mat.ldlt().solve(g1);

    result.fn = u_LAML;
    result.gradient = gr;
    return result;
}




// [[Rcpp::export]]
List aceDLNMadditivebuild(const Eigen::VectorXd R_y,
                   const List R_B_inner_list,
                   const List R_knots_f_list, //  const Eigen::VectorXd R_knots_f,
                   const Eigen::MatrixXd R_Sw_large,
                   const List R_Sf_list, //const Eigen::MatrixXd R_Sf,
                   const Eigen::MatrixXd R_Dw,
                   const Eigen::MatrixXd R_Xrand,
                   const Eigen::MatrixXd R_Xfix,
                   const List R_Zf_list, // const Eigen::MatrixXd R_Zf,
                   const Eigen::VectorXd R_Xoffset,
                   const Eigen::VectorXd R_r,
                   const Eigen::MatrixXd R_K,
                   const Eigen::VectorXd R_a,
                   Eigen::VectorXd R_alpha_f,
                   Eigen::VectorXd R_phi,
                   double R_log_theta,
                   Eigen::VectorXd R_log_smoothing_f,
                   Eigen::VectorXd R_log_smoothing_w,
                   Eigen::VectorXd R_betaR,
                   Eigen::VectorXd R_betaF,
                   Eigen::VectorXd R_logsmoothing) {
    // convert
    Vec y = R_y.cast<Scalar>();
    std::vector<Mat> B_inner_list = convertListToVectorMat(R_B_inner_list);
    std::vector<Eigen::MatrixXd> R_B_inner_Eigen_list = convertListToVectorEigenMatrixXd(R_B_inner_list);

    // Vec knots_f = R_knots_f.cast<Scalar>();
    std::vector<Vec> knots_f_list = convertListToVectorVec(R_knots_f_list);
    std::vector<Eigen::VectorXd> R_knots_f_Eigen_list = convertListToVectorEigenVectorXd(R_knots_f_list);
    Mat Sw_large = R_Sw_large.cast<Scalar>();
    // Mat Sf = R_Sf.cast<Scalar>();
    std::vector<Mat> Sf_list = convertListToVectorMat(R_Sf_list);
    std::vector<Eigen::MatrixXd> R_Sf_Eigen_list = convertListToVectorEigenMatrixXd(R_Sf_list);
    Mat Dw = R_Dw.cast<Scalar>();
    Mat Xrand = R_Xrand.cast<Scalar>();
    Mat Xfix = R_Xfix.cast<Scalar>();
    // Mat Zf = R_Zf.cast<Scalar>();
    std::vector<Mat> Zf_list = convertListToVectorMat(R_Zf_list);
    std::vector<Eigen::MatrixXd> R_Zf_Eigen_list = convertListToVectorEigenMatrixXd(R_Zf_list);
    Vec Xoffset = R_Xoffset.cast<Scalar>();
    Vec r = R_r.cast<Scalar>();
    Mat K = R_K.cast<Scalar>();
    Vec a = R_a.cast<Scalar>();
    Vec alpha_f = R_alpha_f.cast<Scalar>();
    Vec phi = R_phi.cast<Scalar>();
    Scalar log_theta = R_log_theta;
    Vec log_smoothing_f = R_log_smoothing_f.cast<Scalar>();
    Vec log_smoothing_w = R_log_smoothing_w.cast<Scalar>();
    Vec betaR = R_betaR.cast<Scalar>();
    Vec betaF = R_betaF.cast<Scalar>();
    Vec logsmoothing = R_logsmoothing.cast<Scalar>();


    Modelcppad* modelcppadobj_ptr = new Modelcppad(y, B_inner_list, knots_f_list, Sw_large, Sf_list, Dw,
                                                   Xrand, Xfix, Zf_list, Xoffset, r, K, a,
                                                   alpha_f, phi, log_theta, log_smoothing_f, log_smoothing_w,
                                                   betaR, betaF, logsmoothing);

    Rcpp::XPtr<Modelcppad> ptrcppad(modelcppadobj_ptr);

    Model* modelobj_ptr = new Model(R_y, R_B_inner_Eigen_list, R_knots_f_Eigen_list, R_Sw_large, R_Sf_Eigen_list, R_Dw,
                                    R_Xrand, R_Xfix, R_Zf_Eigen_list, R_Xoffset, R_r, R_K, R_a,
                                    R_alpha_f, R_phi, R_log_theta, R_log_smoothing_f, R_log_smoothing_w,
                                    R_betaR, R_betaF, R_logsmoothing);
    Rcpp::XPtr<Model> ptr(modelobj_ptr);


    return List::create(Named("address.eigen") = ptr,
                        Named("address.cppad") = ptrcppad);
}






// [[Rcpp::export]]
List aceDLNMadditiveopt(SEXP ptr,
                SEXP ptrcppad,
                Eigen::VectorXd R_alpha_f,
                Eigen::VectorXd R_phi,
                double R_log_theta,
                Eigen::VectorXd R_log_smoothing_f,
                Eigen::VectorXd R_log_smoothing_w,
                Eigen::VectorXd R_betaR,
                Eigen::VectorXd R_betaF,
                Eigen::VectorXd R_logsmoothing,
                bool verbose) {


    Rcpp::XPtr<Model> modelobj_ptr(ptr);
    Rcpp::XPtr<Modelcppad> modelcppadobj_ptr(ptrcppad);

    Model& modelobj = *modelobj_ptr;
    Modelcppad& modelcppadobj = *modelcppadobj_ptr;

    modelobj.setAlphaF(R_alpha_f);
    modelobj.setPhi(R_phi);
    modelobj.setBetaR(R_betaR);
    modelobj.setBetaF(R_betaF);
    modelobj.setLogTheta(R_log_theta);
    modelobj.setLogSmoothingF(R_log_smoothing_f);
    modelobj.setLogSmoothingW(R_log_smoothing_w);
    modelobj.setLogsmoothing(R_logsmoothing);



    // Vec alpha_f = R_alpha_f.cast<Scalar>();
    // Vec phi = R_phi.cast<Scalar>();
    // Scalar log_theta = R_log_theta;
    // Scalar log_smoothing_f = R_log_smoothing_f;
    // Scalar log_smoothing_w = R_log_smoothing_w;
    // Vec betaR = R_betaR.cast<Scalar>();
    // Vec betaF = R_betaF.cast<Scalar>();
    // Vec logsmoothing = R_logsmoothing.cast<Scalar>();

    // modelobj.setAlphaF(alpha_f);
    // modelobj.setPhi(phi);
    // modelobj.setBetaR(betaR);
    // modelobj.setBetaF(betaF);
    // modelobj.setLogTheta(log_theta);
    // modelobj.setLogSmoothingF(log_smoothing_f);
    // modelobj.setLogSmoothingW(log_smoothing_w);
    // modelobj.setLogsmoothing(logsmoothing);

    // Inner opt
    Inner(modelobj, verbose);
    // PL(modelobj, verbose);
    // get gr of LAML
    LAMLResult LAMLresult;
    LAMLresult = LAML(modelobj, modelcppadobj); // comment out for debug


    // if(hasNaN(LAMLresult.gradient)) {
    //   // rebuild the tape
    //   modelcppadobj.ifhastape = false;
    //   LAMLresult = LAML(modelobj, modelcppadobj);
    // }
    return List::create(Named("LAML.fn") = LAMLresult.fn,
                        Named("LAML.gradient") = LAMLresult.gradient,
                        Named("alpha_f.mod") = modelobj.alpha_f,
                        Named("phi.mod") = modelobj.phi,
                        Named("betaR.mod") = modelobj.betaR,
                        Named("betaF.mod") = modelobj.betaF,
                        Named("address") = modelobj_ptr
                        );
}




// [[Rcpp::export]]
List aceDLNMadditiveCI(SEXP ptr,
               const int Rci,
               const int rseed,
               bool ifeta,
               bool delta,
               bool verbose) {

  Rcpp::XPtr<Model> modelobj_ptr(ptr);
  Model& modelobj = *modelobj_ptr;


  Eigen::VectorXd R_alpha_f = modelobj.alpha_f;
  Eigen::VectorXd R_phi = modelobj.phi;
  Eigen::VectorXd R_betaR = modelobj.betaR;
  Eigen::VectorXd R_betaF = modelobj.betaF;

  Eigen::MatrixXd R_Dw = modelobj.Dw;
  Eigen::MatrixXd R_K = modelobj.K;
  Eigen::VectorXd R_a = modelobj.a;

  int kw = modelobj.kw;
  int kwopt = modelobj.kwopt;
  int kE = modelobj.kE;
  int kbetaR = modelobj.kbetaR;
  int kbetaF = modelobj.kbetaF;
  int M = modelobj.M;
  int kwp = modelobj.kwp;
  int kEp = modelobj.kEp;

  int paraSize = kE+kwopt+kbetaR+kbetaF;
  int paraSizefull;

  // hessian
  Eigen::MatrixXd R_he;
  Eigen::VectorXd R_alpha_w(kw);

  // Vectors for sampling
  Eigen::VectorXd R_phi_sample(kwopt);
  Eigen::VectorXd R_phiKa_sample(kw-1);
  Eigen::VectorXd R_alpha_w_sample(kw);
  Eigen::VectorXd R_alpha_f_sample(kE);
  Eigen::VectorXd R_betaR_sample(kbetaR);
  Eigen::VectorXd R_betaF_sample(kbetaF);

  // Matrices to save results
  Eigen::MatrixXd phi_sample_mat(Rci, kwopt);
  Eigen::MatrixXd alpha_w_sample_mat(Rci, kw);
  Eigen::MatrixXd alpha_f_sample_mat(Rci, kE);
  Eigen::MatrixXd betaR_sample_mat(Rci, kbetaR);
  Eigen::MatrixXd betaF_sample_mat(Rci, kbetaF);

  int n = modelobj.n;
  
  // components for eta
  Eigen::MatrixXd R_E;
  Eigen::VectorXd eta_sample;
  Eigen::MatrixXd eta_sample_mat;
  Eigen::VectorXd Esurface_sample;
  Eigen::MatrixXd Esurface_sample_mat;
  std::vector<Eigen::MatrixXd> R_B_inner_list = modelobj.getB_inner_list();
  std::vector<Eigen::VectorXd>  R_knots_f_list = modelobj.getknots_f_list();
  std::vector<Eigen::MatrixXd> R_Zf_list = modelobj.getZf_list();
  Eigen::MatrixXd R_Xfix = modelobj.getXfix();
  Eigen::MatrixXd R_Xrand = modelobj.getXrand();
  Eigen::VectorXd R_Xoffset = modelobj.getXoffset();
  Eigen::VectorXd R_Bf(kE);

  if(ifeta) {
    R_E.resize(n, M);
    eta_sample.resize(n);
    eta_sample_mat.resize(Rci, n);
    Esurface_sample.resize(n);
    Esurface_sample_mat.resize(Rci, n);
  }

  // Mode of phi
  Eigen::VectorXd R_phi_mod = R_phi;
  // Generate phi
  double R_alpha_w_C_denominator;
  Eigen::VectorXd R_phi_long(kw);
  Eigen::VectorXd R_phiKa = R_K * R_phi + R_a;


  paraSizefull = paraSize;

  // Joint
  R_he = modelobj.he_s_u_mat;
  Eigen::VectorXd R_u_mod(paraSizefull);
  // Hessian
  // cholesky of inverse Hessian
  Eigen::MatrixXd R_he_u_L(paraSize, paraSize);
  Eigen::MatrixXd R_he_u_L_inv(paraSizefull, paraSize);
  Eigen::VectorXd zjoint(paraSize);
  Eigen::VectorXd samplejoint(paraSizefull);


  R_u_mod << R_alpha_f, R_phi, R_betaR, R_betaF;

  // cholesky of inverse Hessian
  R_he_u_L = R_he.llt().matrixL();

  R_he_u_L_inv = (invertL(R_he_u_L)).transpose();

  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(rseed);
  std::normal_distribution<> dist(0, 1);


  for(int i = 0; i < Rci; i++)
  {
    // Jointly sample
    for (int j = 0; j < paraSize; j++) {
      zjoint(j) = dist(gen);
    }
    samplejoint = R_u_mod + R_he_u_L_inv * zjoint;
    // get alpha_f
    R_alpha_f_sample = samplejoint.segment(0, kE);


    // get phi
    R_phi_sample = samplejoint.segment(kE, kwopt);
    R_phiKa_sample = R_K * R_phi_sample + R_a;

    // get alpha_w
    // for (int j = 0; j < (kw - 1); j++) {
    //   R_phi_long(j + 1) = R_phiKa_sample(j);
    // }
    // R_alpha_w_C_denominator = sqrt(R_phi_long.dot(R_Dw * R_phi_long));
    // R_alpha_w_sample = R_phi_long / R_alpha_w_C_denominator;

    for (int j = 0; j < M; j++) {
      R_phi_long(j*kwp) = 1.0;
      R_phi_long.segment(j*kwp + 1, kwp-1) = R_phiKa_sample.segment(j*(kwp-1), (kwp-1));
    }
    for (int j = 0; j < M; j++) {
      R_alpha_w_C_denominator = sqrt(R_phi_long.segment(j*kwp, kwp).dot(R_Dw * R_phi_long.segment(j*kwp, kwp)));
      R_alpha_w_sample.segment(j*kwp, kwp) = R_phi_long.segment(j*kwp, kwp) / R_alpha_w_C_denominator;
    }

    // get betaR
    R_betaR_sample = samplejoint.segment(kE+kwopt, kbetaR);
    // get betaF
    R_betaF_sample = samplejoint.segment(kE+kwopt+kbetaR, kbetaF);

    // save
    phi_sample_mat.row(i) = R_phi_sample.transpose();


    if(ifeta) {
      for (int j = 0; j < M; j++) {
        R_E.col(j) = R_B_inner_list.at(j) * R_alpha_w_sample.segment(j*kwp, kwp);
      }
      for (int ii = 0; ii < n; ii++) {
        for (int j = 0; j < M; j++) {
          R_Bf.segment(j*kEp, kEp) = BsplinevecCon(R_E(ii,j), R_knots_f_list.at(j), 4, R_Zf_list.at(j));
        }
        Esurface_sample(ii) = R_Bf.dot(R_alpha_f_sample);
        eta_sample(ii) = Esurface_sample(ii) + R_Xfix.row(ii).dot(R_betaF_sample) + R_Xrand.row(ii).dot(R_betaR_sample) + R_Xoffset(ii);
      }
      Esurface_sample_mat.row(i) = Esurface_sample.transpose() - Esurface_sample.mean()*Eigen::VectorXd::Ones(n);
      eta_sample_mat.row(i) = eta_sample.transpose();
    }

    // save
    alpha_w_sample_mat.row(i) = R_alpha_w_sample.transpose();
    alpha_f_sample_mat.row(i) = R_alpha_f_sample.transpose();
    betaR_sample_mat.row(i) = R_betaR_sample.transpose();
    betaF_sample_mat.row(i) = R_betaF_sample.transpose();
  }

  return List::create(Named("phi_sample") = phi_sample_mat,
                    Named("alpha_w_sample") = alpha_w_sample_mat,
                    Named("alpha_f_sample") = alpha_f_sample_mat,
                    Named("betaR_sample") = betaR_sample_mat,
                    Named("betaF_sample") = betaF_sample_mat,
                    Named("eta_sample_mat") = eta_sample_mat,
                    Named("Esurface_sample_mat") = Esurface_sample_mat,
                    Named("Hessian_inner") = R_he);

}



// [[Rcpp::export]]
List ConditionalAICaceDLNMadditive(SEXP ptr) {
  
  Rcpp::XPtr<Model> modelobj_ptr(ptr);
  Model& modelobj = *modelobj_ptr;
  
  modelobj.prepare_AIC();
  // hessian
  Eigen::MatrixXd R_he;
  R_he = modelobj.he_s_u_mat;
  // I 
  Eigen::MatrixXd R_I;
  R_I = modelobj.I_mat;
  // 
  Eigen::MatrixXd mat_AIC = R_he.ldlt().solve(R_I);

  double l = modelobj.NegLogL_l;
  // double edf1 = (2 * mat_AIC - mat_AIC * mat_AIC).trace();
  double edf = mat_AIC.trace(); 
  // the widely used version of conditional AIC proposed by Hastie and Tibshirani (1990). 
  // See Wood et al. 2016 JASA
  double AIC = 2.0*l + 2.0*edf;

  return List::create(Named("AIC") = AIC,
                      Named("l") = -1.0*l,
                      Named("edf") = edf);
}

