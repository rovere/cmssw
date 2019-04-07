#define _USE_MATH_DEFINES

#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <memory>  // unique_ptr
#include<chrono>

#include <TFile.h>
#include <TH1F.h>
#include <TMath.h> // Chi2 probabilities

//#define USE_BL

#ifdef USE_BL
#include "RecoPixelVertexing/PixelTrackFitting/interface/BrokenLine.h"
#else
#include "RecoPixelVertexing/PixelTrackFitting/interface/RiemannFit.h"
#endif

using namespace std;
using namespace Eigen;
using namespace Rfit;
using std::unique_ptr;

namespace Rfit {
using Vector3i = Eigen::Matrix<int, 3, 1>;
using Vector4i = Eigen::Matrix<int, 4, 1>;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Vector8d = Eigen::Matrix<double, 8, 1>;
};  // namespace Rfit

// quadruplets...
struct hits_gen {
  Matrix3xNd<4> hits;
  Eigen::Matrix<float,6,4> hits_ge;
  Vector5d true_par;
};

struct geometry {
  Vector8d barrel;
  Vector4i barrel_2;
  Vector8d R_err;
  Vector8d Rp_err;
  Vector8d z_err;
  Vector6d hand;
  Vector3i hand_2;
  Vector6d xy_err;
  Vector6d zh_err;
  double z_max;
  double r_max;
};

void test_helix_fit();

constexpr int c_speed = 299792458;
constexpr double pi = M_PI;
default_random_engine generator(1);

/**
 * Returns the points smeared according to the specified errors in the radial
 * direction (dist_R), in the direction perpendicular to R (dist_Rp) and in Z.
 * The smearing is customized for the barrel and endcap case differently. The
 * transition between the 2 region is set by the caller function and passed
 * here as a boolean paramter.
*/


void smearing(const Vector5d& err, const bool& isbarrel, double& x, double& y, double& z) {
  normal_distribution<double> dist_R(0., err[0]);
  normal_distribution<double> dist_Rp(0., err[1]);
  normal_distribution<double> dist_z(0., err[2]);
  normal_distribution<double> dist_xyh(0., err[3]);
  normal_distribution<double> dist_zh(0., err[4]);
  if (isbarrel) {
    double dev_Rp = dist_Rp(generator);
    double dev_R = dist_R(generator);
    double R = sqrt(Rfit::sqr(x) + Rfit::sqr(y));
//    x += dev_Rp * +y / R + dev_R * -x / R;
//    y += dev_Rp * -x / R + dev_R * -y / R;
    // Original signs seem wrong to me, unless I swapped the direction of the
    // R correction.
    x += dev_Rp * +y / R + dev_R * x / R;
    y += dev_Rp * -x / R + dev_R * y / R;
    z += dist_z(generator);
  } else {
    x += dist_xyh(generator);
    y += dist_xyh(generator);
    z += dist_zh(generator);
  }
}


template<int N>
void Hits_cov(Eigen::Matrix<float,6,4> & V, const unsigned int& i,
    const unsigned int& n, const Matrix3xNd<N>& hits,
    const Vector5d& err, bool isbarrel) {
  if (isbarrel) {
    /**
     * Index of the covariance matrix:
     *
     * Errors computed starting from (in the barrel case):
     * x ==> x + delta_r*x/R + delta_rp*y/R = x + delta_x
     * y ==> y + delta_r*y/R - delta_rp*x/R = y + delta_y
     * z ==> z + delta_z
     *
     * 0 = Cov(x,x) = (delta_r^2*x^2 + delta_rp^2*y^2) / R^2
     * 1 = Cov(x,y) = Cov(y,x) = (delta_r^2 - delta_rp^2) * (x*y)/R^2
     * 2 = Cov(y,y) = (delta_r^2*y^2 + delta_rp^2*x^2) / R^2
     * 5 = Cov(z,z)
     * All other terms are set to 0
     *
     * Forward case:
     * 0 = Cov(x,x) = delta_rp**2
     * 2 = Cov(y,y) = delta_rp**2
     * 5 = Cov(z,z) = delta_z**2
     *
     * All other terms are set to 0
     *
     */
    double R2 = Rfit::sqr(hits(0, i)) + Rfit::sqr(hits(1, i));
    V.col(i)[0] =
        (Rfit::sqr(err[1]) * Rfit::sqr(hits(1, i)) + Rfit::sqr(err[0]) * Rfit::sqr(hits(0, i))) /
        R2;
    V.col(i)[1] =
        (Rfit::sqr(err[0]) - Rfit::sqr(err[1])) * hits(1, i) * hits(0, i) / R2;
    V.col(i)[2] =
        (Rfit::sqr(err[1]) * Rfit::sqr(hits(0, i)) + Rfit::sqr(err[0]) * Rfit::sqr(hits(1, i))) /
        R2;
    V.col(i)[5] = Rfit::sqr(err[2]);
  } else {
    V.col(i)[0] = Rfit::sqr(err[3]);
    V.col(i)[2] = Rfit::sqr(err[3]);
    V.col(i)[5] = Rfit::sqr(err[4]);
  }
}

hits_gen Hits_gen(const unsigned int& n, const Matrix<double, 6, 1>& gen_par) {
  hits_gen gen;
  gen.hits = MatrixXd::Zero(3, n);
  gen.hits_ge = Eigen::Matrix<float,6,4>::Zero();
  // err /= 10000.;
  constexpr double rad[8] = {2.95, 6.8, 10.9, 16., 3.1, 7., 11., 16.2};
  // constexpr double R_err[8] = {5./10000, 5./10000, 5./10000, 5./10000, 5./10000,
  // 5./10000, 5./10000, 5./10000};  constexpr double Rp_err[8] = {35./10000, 18./10000,
  // 15./10000, 34./10000, 35./10000, 18./10000, 15./10000, 34./10000};  constexpr double z_err[8] =
  // {72./10000, 38./10000, 25./10000, 56./10000, 72./10000, 38./10000, 25./10000, 56./10000};
  constexpr double R_err[8] = {10. / 10000, 10. / 10000, 10. / 10000, 10. / 10000,
                               10. / 10000, 10. / 10000, 10. / 10000, 10. / 10000};
  constexpr double Rp_err[8] = {35. / 10000, 35. / 10000, 35. / 10000, 35. / 10000,
                                35. / 10000, 35. / 10000, 35. / 10000, 35. / 10000};
  constexpr double z_err[8] = {70. / 10000, 70. / 10000, 70. / 10000, 70. / 10000,
                               70. / 10000, 70. / 10000, 70. / 10000, 70. / 10000};
  // (x2, y2) is the center of the circle in the transverse plane
  const double x2 = gen_par(0) + gen_par(4) * cos(gen_par(3) * pi / 180);
  const double y2 = gen_par(1) + gen_par(4) * sin(gen_par(3) * pi / 180);
  // alpha is the angle between the X-Axis and the center of the circle in the
  // transverse plane. It is identical to gen_par(3), in radians.
  const double alpha = atan2(y2, x2);

  for (unsigned int i = 0; i < n; ++i) {
    const double a = gen_par(4);
    const double b = rad[i];
    const double c = sqrt(Rfit::sqr(x2) + Rfit::sqr(y2));
    const double beta = acos((Rfit::sqr(a) - Rfit::sqr(b) - Rfit::sqr(c)) / (-2. * b * c));
    const double gamma = alpha + beta;
    gen.hits(0, i) = rad[i] * cos(gamma);
    gen.hits(1, i) = rad[i] * sin(gamma);
    gen.hits(2, i) = gen_par(2) + 1 / tan(gen_par(5) * pi / 180) * 2. *
                                      asin(sqrt(Rfit::sqr((gen_par(0) - gen.hits(0, i))) +
                                                Rfit::sqr((gen_par(1) - gen.hits(1, i)))) /
                                           (2. * gen_par(4))) *
                                      gen_par(4);
    // We ideally stop being in the barrel if our abs(z) position is beyond
    // 27cm. Quite arbitrary but not far from reality for CMS.
    bool isbarrel = (std::abs(gen.hits(2, i)) < 27);
    Vector5d err;
    // The first 3 errors are used in the barrel case and represent the error
    // in the radial direction (usually very small), the error in the direction
    // orthogonal to the radial one, usually on the plane of the detector (tens
    // of microns) and along the Z axis (several tens of microns, since that's
    // the least precise coordinate for the barrel case). The last 2 errors are
    // used in the forward case: the fourth in the plane of the forward
    // detector (tens of microns) and the fifth is along Z (few microns, since
    // this is the most precise coordinates in the forward region).
    err << R_err[i], Rp_err[i], z_err[i], Rp_err[i], R_err[i];
    smearing(err, isbarrel, gen.hits(0, i), gen.hits(1, i), gen.hits(2, i));
    Hits_cov(gen.hits_ge, i, n, gen.hits, err, isbarrel);
  }

  return gen;
}

/**
 * Take the input 6 parameters: x,y,z,phi,R,theta and generate the parameters
 * that are actually fit using the helix and line fit. In particular:
 *
 * X0 and Y0 are the center of the circle in the tansverse plane
 *
 * The radius of the circle in the transverse plane, in cm, is left unchanged.
 *
 * Finally, the three parameters (X0, Y0, R) are translated back into (phi,
 * Tip, p_t). Phi, in this case, is the angle between p_t and the X-Axis, as on
 * page 5 of the CMS Note. Tip is the signed distance of the circle with
 * respect to the origin. P_t is in GeV/c.
 *
 * The fourth parameter becomes the cot(Theta), in radians.
 *
 * The fifth parameters should become Zip (to be checked, the aritmetic I saw
 * should add up to 0 correction, but with lots of useless computations???)
 *
 * This function therefore returns:
 *
 * (Phi, Tip(signed), p_t, cot(theta)[radians], Zip).
 */

Vector5d True_par(const Matrix<double, 6, 1>& gen_par, const int& charge, const double& B_field) {
  Vector5d true_par;
  const double x0 = gen_par(0) + gen_par(4) * cos(gen_par(3) * pi / 180);
  const double y0 = gen_par(1) + gen_par(4) * sin(gen_par(3) * pi / 180);
  circle_fit circle;
  circle.par << x0, y0, gen_par(4);
  circle.q = 1;
  Rfit::par_uvrtopak(circle, B_field, false);
  true_par.block(0, 0, 3, 1) = circle.par;
  true_par(3) = 1 / tan(gen_par(5) * pi / 180);
  const int dir = ((gen_par(0) - cos(true_par(0) - pi / 2) * true_par(1)) * (gen_par(1) - y0) -
                       (gen_par(1) - sin(true_par(0) - pi / 2) * true_par(1)) * (gen_par(0) - x0) >
                   0)
                      ? -1
                      : 1;
  true_par(4) = gen_par(2) +
                1 / tan(gen_par(5) * pi / 180) * dir * 2.f *
                    asin(sqrt(Rfit::sqr((gen_par(0) - cos(true_par(0) - pi / 2) * true_par(1))) +
                              Rfit::sqr((gen_par(1) - sin(true_par(0) - pi / 2) * true_par(1)))) /
                         (2.f * gen_par(4))) *
                    gen_par(4);
  return true_par;
}

/**
 * Take the input parameters: x,y,z,phi,p_t and eta and transform them in the
 * format that is most useful to generate the points along the helix.
 *
 * The x, y and z floats are left unchanges and represent the d0 and Zip of the
 * helix at the reference point, which is the point of closest approach to the
 * beamline.
 *
 * Phi is, instead, translated by 90 degrees according to the charge of the
 * particle: if negative we add 90 degrees, otherwise we subtract 90.  This is
 * very much evident from the usual right-hand rule, the Lorentz force and the
 * drawing of the circle in the transverse plane, considering B along the
 * positive Z axis. Phi is basically the angle between the origin and the
 * center of the circle in the transverse plane.
 *
 * The fifth parameter, p_t, is translated into the radius of the circle in the
 * transverse plane using the usual convention:
 *
 * p_t = 0.299 * R[m] * B[T] =
 *
 * so that:
 *
 * R[cm] = p_t/b_field, where:
 *
 * b_field = c[m/s] * 10^-9 * 10^-2(m->cm)
 *
 * The sixth parameter, eta, is translated into the corresponding theta angle,
 * in degrees.
 *
 */

Matrix<double, 6, 1> New_par(const Matrix<double, 6, 1>& gen_par, const int& charge,
                             const double& B_field) {
  Matrix<double, 6, 1> new_par;
  new_par.block(0, 0, 3, 1) = gen_par.block(0, 0, 3, 1);
  new_par(3) = gen_par(3) - charge * 90;
  new_par(4) = gen_par(4) / B_field;
//  new_par(5) = atan(sinh(gen_par(5))) * 180 / pi;
  new_par(5) = 2.*atan(exp(-gen_par(5))) * 180 / pi;
  return new_par;
}

template<typename Fit, size_t N>
void computePull(std::array<Fit, N> & fit, const char * label,
    int n_, int iteration, const Vector5d & true_par) {
  Eigen::Matrix<double, 41, Eigen::Dynamic, 1> score(41, iteration);

  std::string histo_name("Phi Pull");
  histo_name += label;
  TH1F phi_pull(histo_name.data(), histo_name.data(), 100, -10., 10.);
  histo_name = "dxy Pull ";
  histo_name += label;
  TH1F dxy_pull(histo_name.data(), histo_name.data(), 100, -10., 10.);
  histo_name = "dz Pull ";
  histo_name += label;
  TH1F dz_pull(histo_name.data(), histo_name.data(), 100, -10., 10.);
  histo_name = "Theta Pull ";
  histo_name += label;
  TH1F theta_pull(histo_name.data(), histo_name.data(), 100, -10., 10.);
  histo_name = "Pt Pull ";
  histo_name += label;
  TH1F pt_pull(histo_name.data(), histo_name.data(), 100, -10., 10.);
  histo_name = "Phi Error ";
  histo_name += label;
  TH1F phi_error(histo_name.data(), histo_name.data(), 100, 0., 0.1);
  histo_name = "dxy error ";
  histo_name += label;
  TH1F dxy_error(histo_name.data(), histo_name.data(), 100, 0., 0.1);
  histo_name = "dz error ";
  histo_name += label;
  TH1F dz_error(histo_name.data(), histo_name.data(), 100, 0., 0.1);
  histo_name = "Theta error ";
  histo_name += label;
  TH1F theta_error(histo_name.data(), histo_name.data(), 100, 0., 0.1);
  histo_name = "Pt error ";
  histo_name += label;
  TH1F pt_error(histo_name.data(), histo_name.data(), 100, 0., 0.1);
  histo_name = "Chi2 Line";
  histo_name += label;
  TH1F chi2_line(histo_name.data(), histo_name.data(), 100, 0., 10.);
  histo_name = "Chi2 Circle";
  histo_name += label;
  TH1F chi2_circle(histo_name.data(), histo_name.data(), 100, 0., 10.);
  histo_name = "Chi2 Line Prob";
  histo_name += label;
  TH1F chi2_lineProb(histo_name.data(), histo_name.data(), 100, 0., 1.);
  histo_name = "Chi2 Circle Prob";
  histo_name += label;
  TH1F chi2_circleProb(histo_name.data(), histo_name.data(), 100, 0., 1.);
  for (int x = 0; x < iteration; x++) {
    // Compute PULLS information
    score(0, x) = (fit[x].par(0) - true_par(0)) / sqrt(fit[x].cov(0, 0));
    score(1, x) = (fit[x].par(1) - true_par(1)) / sqrt(fit[x].cov(1, 1));
    score(2, x) = (fit[x].par(2) - true_par(2)) / sqrt(fit[x].cov(2, 2));
    score(3, x) = (fit[x].par(3) - true_par(3)) / sqrt(fit[x].cov(3, 3));
    score(4, x) = (fit[x].par(4) - true_par(4)) / sqrt(fit[x].cov(4, 4));
    phi_pull.Fill(score(0, x));
    dxy_pull.Fill(score(1, x));
    pt_pull.Fill(score(2, x));
    theta_pull.Fill(score(3, x));
    dz_pull.Fill(score(4, x));
    phi_error.Fill(sqrt(fit[x].cov(0, 0)));
    dxy_error.Fill(sqrt(fit[x].cov(1, 1)));
    pt_error.Fill(sqrt(fit[x].cov(2, 2)));
    theta_error.Fill(sqrt(fit[x].cov(3, 3)));
    dz_error.Fill(sqrt(fit[x].cov(4, 4)));
    chi2_circle.Fill(fit[x].chi2_circle);
    chi2_line.Fill(fit[x].chi2_line);
    chi2_circleProb.Fill(TMath::Prob(fit[x].chi2_circle, n_ - 3));
    chi2_lineProb.Fill(TMath::Prob(fit[x].chi2_line, n_ - 2));
    score(5, x) =
      (fit[x].par(0) - true_par(0)) * (fit[x].par(1) - true_par(1)) / (fit[x].cov(0, 1));
    score(6, x) =
      (fit[x].par(0) - true_par(0)) * (fit[x].par(2) - true_par(2)) / (fit[x].cov(0, 2));
    score(7, x) =
      (fit[x].par(1) - true_par(1)) * (fit[x].par(2) - true_par(2)) / (fit[x].cov(1, 2));
    score(8, x) =
      (fit[x].par(3) - true_par(3)) * (fit[x].par(4) - true_par(4)) / (fit[x].cov(3, 4));
    score(9, x) = fit[x].chi2_circle;
    score(25, x) = fit[x].chi2_line;
    score(10, x) = sqrt(fit[x].cov(0, 0)) / fit[x].par(0) * 100;
    score(13, x) = sqrt(fit[x].cov(3, 3)) / fit[x].par(3) * 100;
    score(14, x) = sqrt(fit[x].cov(4, 4)) / fit[x].par(4) * 100;
    score(15, x) = (fit[x].par(0) - true_par(0)) * (fit[x].par(3) - true_par(3)) /
      sqrt(fit[x].cov(0, 0)) / sqrt(fit[x].cov(3, 3));
    score(16, x) = (fit[x].par(1) - true_par(1)) * (fit[x].par(3) - true_par(3)) /
      sqrt(fit[x].cov(1, 1)) / sqrt(fit[x].cov(3, 3));
    score(17, x) = (fit[x].par(2) - true_par(2)) * (fit[x].par(3) - true_par(3)) /
      sqrt(fit[x].cov(2, 2)) / sqrt(fit[x].cov(3, 3));
    score(18, x) = (fit[x].par(0) - true_par(0)) * (fit[x].par(4) - true_par(4)) /
      sqrt(fit[x].cov(0, 0)) / sqrt(fit[x].cov(4, 4));
    score(19, x) = (fit[x].par(1) - true_par(1)) * (fit[x].par(4) - true_par(4)) /
      sqrt(fit[x].cov(1, 1)) / sqrt(fit[x].cov(4, 4));
    score(20, x) = (fit[x].par(2) - true_par(2)) * (fit[x].par(4) - true_par(4)) /
      sqrt(fit[x].cov(2, 2)) / sqrt(fit[x].cov(4, 4));
    score(21, x) = (fit[x].par(0) - true_par(0)) * (fit[x].par(1) - true_par(1)) /
      sqrt(fit[x].cov(0, 0)) / sqrt(fit[x].cov(1, 1));
    score(22, x) = (fit[x].par(0) - true_par(0)) * (fit[x].par(2) - true_par(2)) /
      sqrt(fit[x].cov(0, 0)) / sqrt(fit[x].cov(2, 2));
    score(23, x) = (fit[x].par(1) - true_par(1)) * (fit[x].par(2) - true_par(2)) /
      sqrt(fit[x].cov(1, 1)) / sqrt(fit[x].cov(2, 2));
    score(24, x) = (fit[x].par(3) - true_par(3)) * (fit[x].par(4) - true_par(4)) /
      sqrt(fit[x].cov(3, 3)) / sqrt(fit[x].cov(4, 4));
    score(30, x) = fit[x].par(0);
    score(31, x) = fit[x].par(1);
    score(32, x) = fit[x].par(2);
    score(33, x) = fit[x].par(3);
    score(34, x) = fit[x].par(4);
    score(35, x) = sqrt(fit[x].cov(0,0));
    score(36, x) = sqrt(fit[x].cov(1,1));
    score(37, x) = sqrt(fit[x].cov(2,2));
    score(38, x) = sqrt(fit[x].cov(3,3));
    score(39, x) = sqrt(fit[x].cov(4,4));
  }

  double phi_ = score.row(0).mean();
  double a_ = score.row(1).mean();
  double pt_ = score.row(2).mean();
  double coT_ = score.row(3).mean();
  double Zip_ = score.row(4).mean();
  std::cout << std::setprecision(5) << std::scientific << label << " AVERAGE FITTED VALUES: \n"
    << "phi: " << score.row(30).mean() << " +/- " << score.row(35).mean() << " [+/-] " << sqrt(score.row(35).array().abs2().mean() - score.row(35).mean()*score.row(35).mean()) << std::endl
    << "d0:  " << score.row(31).mean() << " +/- " << score.row(36).mean() << " [+/-] " << sqrt(score.row(36).array().abs2().mean() - score.row(36).mean()*score.row(36).mean()) << std::endl
    << "pt:  " << score.row(32).mean() << " +/- " << score.row(37).mean() << " [+/-] " << sqrt(score.row(37).array().abs2().mean() - score.row(37).mean()*score.row(37).mean()) << std::endl
    << "coT: " << score.row(33).mean() << " +/- " << score.row(38).mean() << " [+/-] " << sqrt(score.row(38).array().abs2().mean() - score.row(38).mean()*score.row(38).mean()) << std::endl
    << "Zip: " << score.row(34).mean() << " +/- " << score.row(39).mean() << " [+/-] " << sqrt(score.row(39).array().abs2().mean() - score.row(39).mean()*score.row(39).mean()) << std::endl;

  Matrix5d correlation;
  correlation << 1., score.row(21).mean(), score.row(22).mean(), score.row(15).mean(),
              score.row(20).mean(), score.row(21).mean(), 1., score.row(23).mean(), score.row(16).mean(),
              score.row(19).mean(), score.row(22).mean(), score.row(23).mean(), 1., score.row(17).mean(),
              score.row(20).mean(), score.row(15).mean(), score.row(16).mean(), score.row(17).mean(), 1.,
              score.row(24).mean(), score.row(18).mean(), score.row(19).mean(), score.row(20).mean(),
              score.row(24).mean(), 1.;

  cout << "\n" << label << " PULLS (mean, sigma, relative_error):\n"
    << "phi:  " << phi_ << "     "
    << sqrt((score.row(0).array() - phi_).square().sum() / (iteration - 1)) << "   "
    << abs(score.row(10).mean()) << "%\n"
    << "a0 :  " << a_ << "     "
    << sqrt((score.row(1).array() - a_).square().sum() / (iteration - 1)) << "   "
    << abs(score.row(11).mean()) << "%\n"
    << "pt :  " << pt_ << "     "
    << sqrt((score.row(2).array() - pt_).square().sum() / (iteration - 1)) << "   "
    << abs(score.row(12).mean()) << "%\n"
    << "coT:  " << coT_ << "     "
    << sqrt((score.row(3).array() - coT_).square().sum() / (iteration - 1)) << "   "
    << abs(score.row(13).mean()) << "%\n"
    << "Zip:  " << Zip_ << "     "
    << sqrt((score.row(4).array() - Zip_).square().sum() / (iteration - 1)) << "   "
    << abs(score.row(14).mean()) << "%\n\n"
    << "cov(phi,a0)_:  " << score.row(5).mean() << "\n"
    << "cov(phi,pt)_:  " << score.row(6).mean() << "\n"
    << "cov(a0,pt)_:   " << score.row(7).mean() << "\n"
    << "cov(coT,Zip)_: " << score.row(8).mean() << "\n\n"
    << "chi2_circle:  " << score.row(9).mean() << " vs " << n_ - 3 << "\n"
    << "chi2_line:    " << score.row(25).mean() << " vs " << n_ - 2 << "\n\n"
    << "correlation matrix:\n"
    << correlation << "\n\n"
    << endl;

  phi_pull.Fit("gaus", "Q");
  dxy_pull.Fit("gaus", "Q");
  dz_pull.Fit("gaus", "Q");
  theta_pull.Fit("gaus", "Q");
  pt_pull.Fit("gaus", "Q");
  phi_pull.Write();
  dxy_pull.Write();
  dz_pull.Write();
  theta_pull.Write();
  pt_pull.Write();
  phi_error.Write();
  dxy_error.Write();
  dz_error.Write();
  theta_error.Write();
  pt_error.Write();
  chi2_circle.Write();
  chi2_line.Write();
  chi2_circleProb.Write();
  chi2_lineProb.Write();
}


void test_helix_fit(bool getcin) {
  int n_;
  const double B_field = 3.8 * c_speed / pow(10, 9) / 100;
  Matrix<double, 6, 1> gen_par;
  Vector5d true_par;
  Vector5d err;
  generator.seed(1);
  std::cout << std::setprecision(6);
  cout << "_________________________________________________________________________\n";
  cout << "n x(cm) y(cm) z(cm) phi(grad) R(Gev/c) eta iteration debug" << endl;
  if (getcin) {
    cout << "hits: ";
    cin  >> n_;
    cout << "x: ";
    cin  >> gen_par(0);
    cout << "y: ";
    cin  >> gen_par(1);
    cout << "z: ";
    cin  >> gen_par(2);
    cout << "phi: ";
    cin  >> gen_par(3);
    cout << "p_t: ";
    cin  >> gen_par(4);
    cout << "eta: ";
    cin  >> gen_par(5);
  } else {
     n_ = 4;
     gen_par(0) = -0.1;  // x
     gen_par(1) = 0.1;   // y
     gen_par(2) = -1.;   // z
     gen_par(3) = 45.;   // phi
     gen_par(4) = 10.;   // R (p_t)
     gen_par(5) = 1.;    // eta
  }

  std::cout << "\nInput parameters:\n" << std::setprecision(4)
    << "x:        " << gen_par(0) << "\n"
    << "y:        " << gen_par(1) << "\n"
    << "z:        " << gen_par(2) << "\n"
    << "phi:      " << gen_par(3) << "\n"
    << "phi(rad): " << gen_par(3)*M_PI/180. << "\n"
    << "p_t:      " << gen_par(4) << "\n"
    << "eta:      " << gen_par(5) << std::endl;

  const int iteration = 5000;
  gen_par = New_par(gen_par, 1, B_field);
  true_par = True_par(gen_par, 1, B_field);
  std::array<helix_fit, iteration> helixRiemann_fit;

  std::cout << "\nTransformed Input parameters:\n" << std::setprecision(4)
    << "x:     " << gen_par(0) << "\n"
    << "y:     " << gen_par(1) << "\n"
    << "z:     " << gen_par(2) << "\n"
    << "phi_c: " << gen_par(3) << "\n"
    << "R:     " << gen_par(4) << "\n"
    << "theta: " << gen_par(5) << "\n"
    << "dxy:   " << std::sqrt(gen_par(0)*gen_par(0) + gen_par(1)*gen_par(1)) << "\n"
    << std::endl;
  std::cout << "\nTrue parameters:\n" << std::setprecision(4)
    << "phi:  " << true_par(0) << "\n"
    << "dxy:  " << true_par(1) << "\n"
    << "pt:   " << true_par(2) << "\n"
    << "CotT: " << true_par(3) << "\n"
    << "Zip:  " << true_par(4) << "\n"
    << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  auto delta = start-start;
  for (int i = 0; i < 100*iteration; i++) {
    hits_gen gen;
    gen = Hits_gen(n_, gen_par);
    //      gen.hits = MatrixXd::Zero(3, 4);
    //      gen.hits_cov = MatrixXd::Zero(3 * 4, 3 * 4);
    //      gen.hits.col(0) << 1.82917642593, 2.0411875248, 7.18495464325;
    //      gen.hits.col(1) << 4.47041416168, 4.82704305649, 18.6394691467;
    //      gen.hits.col(2) << 7.25991010666, 7.74653434753, 30.6931324005;
    //      gen.hits.col(3) << 8.99161434174, 9.54262828827, 38.1338043213;
    delta -= std::chrono::high_resolution_clock::now()-start;
    helixRiemann_fit[i%iteration] =
#ifdef USE_BL
      BrokenLine::BL_Helix_fit(gen.hits, gen.hits_ge, B_field);
#else
      Rfit::Helix_fit(gen.hits, gen.hits_ge, B_field, true);
#endif
    delta += std::chrono::high_resolution_clock::now()-start;

    if (helixRiemann_fit[i%iteration].par(0)>10.) std::cout << "error" << std::endl;
    if (0==i)
      cout << std::setprecision(6)
        << "phi:  " << helixRiemann_fit[i].par(0) << " +/- " << sqrt(helixRiemann_fit[i].cov(0, 0)) << " vs "
        << true_par(0) << endl
        << "Tip:  " << helixRiemann_fit[i].par(1) << " +/- " << sqrt(helixRiemann_fit[i].cov(1, 1)) << " vs "
        << true_par(1) << endl
        << "p_t:  " << helixRiemann_fit[i].par(2) << " +/- " << sqrt(helixRiemann_fit[i].cov(2, 2)) << " vs "
        << true_par(2) << endl
        << "theta:" << helixRiemann_fit[i].par(3) << " +/- " << sqrt(helixRiemann_fit[i].cov(3, 3)) << " vs "
        << true_par(3) << endl
        << "Zip:  " << helixRiemann_fit[i].par(4) << " +/- " << sqrt(helixRiemann_fit[i].cov(4, 4)) << " vs "
        << true_par(4) << endl
        << "charge:" << helixRiemann_fit[i].q << " vs 1" << endl
        << "covariance matrix:" << endl
        << helixRiemann_fit[i].cov << endl
        << "Initial hits:\n" << gen.hits << endl
        << "Initial Covariance:\n" << gen.hits_ge << endl;

  }
  std::cout << "elapsted time " << double(std::chrono::duration_cast<std::chrono::nanoseconds>(delta).count())/1.e6 << std::endl;
  computePull(helixRiemann_fit, "Riemann", n_, iteration, true_par);
}

int main(int nargs, char**) {
  TFile f("TestFitResults.root", "RECREATE");
  test_helix_fit(nargs>1);
  f.Close();
  return 0;
}

