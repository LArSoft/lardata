#include "lardata/RecoObjects/TrackStatePropagator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

using namespace recob::tracking;

namespace trkf {

  TrackStatePropagator::TrackStatePropagator(double minStep,
                                             double maxElossFrac,
                                             int maxNit,
                                             double tcut,
                                             double wrongDirDistTolerance,
                                             bool propPinvErr)
    : fMinStep(minStep)
    , fMaxElossFrac(maxElossFrac)
    , fMaxNit(maxNit)
    , fTcut(tcut)
    , fWrongDirDistTolerance(wrongDirDistTolerance)
    , fPropPinvErr(propPinvErr)
  {
    larprop = lar::providerFrom<detinfo::LArPropertiesService>();
  }

  using PropDirection = TrackStatePropagator::PropDirection;

  TrackState TrackStatePropagator::propagateToPlane(bool& success,
                                                    const detinfo::DetectorPropertiesData& detProp,
                                                    const TrackState& origin,
                                                    const Plane& target,
                                                    bool dodedx,
                                                    bool domcs,
                                                    PropDirection dir) const
  {
    //
    // 1- find distance to target plane
    auto [distance, sperp] = distancePairToPlane(success, origin, target);
    //
    if ((distance < -fWrongDirDistTolerance && dir == FORWARD) ||
        (distance > fWrongDirDistTolerance && dir == BACKWARD)) {
      success = false;
      return origin;
    }
    //
    // 2- propagate 3d position by distance, form propagated state on plane parallel to origin plane
    Point_t p = propagatedPosByDistance(
      origin.position(), origin.momentum() * origin.parameters()[4], distance);
    TrackState tmpState(
      SVector5(0., 0., origin.parameters()[2], origin.parameters()[3], origin.parameters()[4]),
      origin.covariance(),
      Plane(p, origin.plane().direction()),
      origin.isTrackAlongPlaneDir(),
      origin.pID());
    //
    // 3- rotate state to target plane
    double dw2dw1 = 0;
    tmpState = rotateToPlane(success, tmpState, target, dw2dw1);
    SVector5 par5d = tmpState.parameters();
    SMatrixSym55 cov5d = tmpState.covariance();
    //
    // 4- compute jacobian to propagate uncertainties
    SMatrix55 pm = ROOT::Math::SMatrixIdentity(); //diagonal elements are 1
    pm(0, 2) = sperp;                             // du2/d(dudw1);
    pm(1, 3) = sperp;                             // dv2/d(dvdw1);
    //
    // 5- apply material effects, performing more iterations if the distance is long
    bool arrived = false;
    int nit = 0; // Iteration count.
    double deriv = 1.;
    SMatrixSym55 noise_matrix;
    while (!arrived) {
      ++nit;
      if (nit > fMaxNit) {
        success = false;
        return origin;
      }
      // Estimate maximum step distance, such that fMaxElossFrac of initial energy is lost by dedx
      const double mass = origin.mass();
      const double p = 1. / par5d[4];
      const double e = std::hypot(p, mass);
      const double t = e - mass;
      const double dedx = 0.001 * detProp.Eloss(std::abs(p), mass, fTcut);
      const double range = t / dedx;
      const double smax = std::max(fMinStep, fMaxElossFrac * range);
      double s = distance;
      if (domcs && smax > 0 && std::abs(s) > smax) {
        if (fMaxNit == 1) {
          success = false;
          return origin;
        }
        s = (s > 0 ? smax : -smax);
        distance -= s;
      }
      else
        arrived = true;
      // now apply material effects
      if (domcs) {
        bool flip = false;
        if (origin.isTrackAlongPlaneDir() == true && dw2dw1 < 0.) flip = true;
        if (origin.isTrackAlongPlaneDir() == false && dw2dw1 > 0.) flip = true;
        bool ok = apply_mcs(detProp,
                            par5d[2],
                            par5d[3],
                            par5d[4],
                            origin.mass(),
                            s,
                            range,
                            p,
                            e * e,
                            flip,
                            noise_matrix);
        if (!ok) {
          success = false;
          return origin;
        }
      }
      if (dodedx) { apply_dedx(par5d(4), detProp, dedx, e, origin.mass(), s, deriv); }
    }
    if (fPropPinvErr) pm(4, 4) *= deriv;
    //
    // 6- create final track state
    cov5d = ROOT::Math::Similarity(pm, cov5d); //*rj
    cov5d = cov5d + noise_matrix;
    TrackState trackState(
      par5d, cov5d, target, origin.momentum().Dot(target.direction()) > 0, origin.pID());
    return trackState;
  }

  TrackState TrackStatePropagator::rotateToPlane(bool& success,
                                                 const TrackState& origin,
                                                 const Plane& target,
                                                 double& dw2dw1) const
  {
    const bool isTrackAlongPlaneDir = origin.momentum().Dot(target.direction()) > 0;
    //
    SVector5 par5 = origin.parameters();
    const double sinA1 = origin.plane().sinAlpha();
    const double cosA1 = origin.plane().cosAlpha();
    const double sinA2 = target.sinAlpha();
    const double cosA2 = target.cosAlpha();
    const double sinB1 = origin.plane().sinBeta();
    const double cosB1 = origin.plane().cosBeta();
    const double sinB2 = target.sinBeta();
    const double cosB2 = target.cosBeta();
    const double sindB = -sinB1 * cosB2 + cosB1 * sinB2;
    const double cosdB = cosB1 * cosB2 + sinB1 * sinB2;
    const double ruu = cosA1 * cosA2 + sinA1 * sinA2 * cosdB;
    const double ruv = sinA2 * sindB;
    const double ruw = sinA1 * cosA2 - cosA1 * sinA2 * cosdB;
    const double rvu = -sinA1 * sindB;
    const double rvv = cosdB;
    const double rvw = cosA1 * sindB;
    const double rwu = cosA1 * sinA2 - sinA1 * cosA2 * cosdB;
    const double rwv = -cosA2 * sindB;
    const double rww = sinA1 * sinA2 + cosA1 * cosA2 * cosdB;
    dw2dw1 = par5[2] * rwu + par5[3] * rwv + rww;
    if (dw2dw1 == 0.) {
      success = false;
      return origin;
    }
    const double dudw2 = (par5[2] * ruu + par5[3] * ruv + ruw) / dw2dw1;
    const double dvdw2 = (par5[2] * rvu + par5[3] * rvv + rvw) / dw2dw1;
    SMatrix55 pm;
    //
    pm(0, 0) = ruu - dudw2 * rwu; // du2/du1
    pm(1, 0) = rvu - dvdw2 * rwu; // dv2/du1
    pm(2, 0) = 0.;                // d(dudw2)/du1
    pm(3, 0) = 0.;                // d(dvdw2)/du1
    pm(4, 0) = 0.;                // d(pinv2)/du1
    //
    pm(0, 1) = ruv - dudw2 * rwv; // du2/dv1
    pm(1, 1) = rvv - dvdw2 * rwv; // dv2/dv1
    pm(2, 1) = 0.;                // d(dudw2)/dv1
    pm(3, 1) = 0.;                // d(dvdw2)/dv1
    pm(4, 1) = 0.;                // d(pinv2)/dv1
    //
    pm(0, 2) = 0.;                           // du2/d(dudw1);
    pm(1, 2) = 0.;                           // dv2/d(dudw1);
    pm(2, 2) = (ruu - dudw2 * rwu) / dw2dw1; // d(dudw2)/d(dudw1);
    pm(3, 2) = (rvu - dvdw2 * rwu) / dw2dw1; // d(dvdw2)/d(dudw1);
    pm(4, 2) = 0.;                           // d(pinv2)/d(dudw1);
    //
    pm(0, 3) = 0.;                           // du2/d(dvdw1);
    pm(1, 3) = 0.;                           // dv2/d(dvdw1);
    pm(2, 3) = (ruv - dudw2 * rwv) / dw2dw1; // d(dudw2)/d(dvdw1);
    pm(3, 3) = (rvv - dvdw2 * rwv) / dw2dw1; // d(dvdw2)/d(dvdw1);
    pm(4, 3) = 0.;                           // d(pinv2)/d(dvdw1);
    //
    pm(0, 4) = 0.; // du2/d(pinv1);
    pm(1, 4) = 0.; // dv2/d(pinv1);
    pm(2, 4) = 0.; // d(dudw2)/d(pinv1);
    pm(3, 4) = 0.; // d(dvdw2)/d(pinv1);
    pm(4, 4) = 1.; // d(pinv2)/d(pinv1);
    //
    par5[0] = (origin.position().X() - target.position().X()) * cosA2 +
              (origin.position().Y() - target.position().Y()) * sinA2 * sinB2 -
              (origin.position().Z() - target.position().Z()) * sinA2 * cosB2;
    par5[1] = (origin.position().Y() - target.position().Y()) * cosB2 +
              (origin.position().Z() - target.position().Z()) * sinB2;
    par5[2] = dudw2;
    par5[3] = dvdw2;
    //
    success = true;
    return TrackState(par5,
                      ROOT::Math::Similarity(pm, origin.covariance()),
                      Plane(origin.position(), target.direction()),
                      isTrackAlongPlaneDir,
                      origin.pID());
  }

  double TrackStatePropagator::distanceToPlane(bool& success,
                                               const Point_t& origpos,
                                               const Vector_t& origmom,
                                               const Plane& target) const
  {
    const Point_t& targpos = target.position();
    const Vector_t& targdir = target.direction();
    //check that origmom is not along the plane, i.e. targdir.Dot(origmom.Unit())=0
    if (targdir.Dot(origmom.Unit()) == 0) {
      success = false;
      return DBL_MAX;
    }
    //distance along track direction
    double s = targdir.Dot(targpos - origpos) / targdir.Dot(origmom.Unit());
    success = true;
    return s;
  }

  double TrackStatePropagator::perpDistanceToPlane(bool& success,
                                                   const Point_t& origpos,
                                                   const Plane& target) const
  {
    const Point_t& targpos = target.position();
    const Vector_t& targdir = target.direction();
    //point-plane distance projected along direction orthogonal to the plane
    double sperp = targdir.Dot(targpos - origpos);
    success = true;
    return sperp;
  }

  std::pair<double, double> TrackStatePropagator::distancePairToPlane(bool& success,
                                                                      const Point_t& origpos,
                                                                      const Vector_t& origmom,
                                                                      const Plane& target) const
  {
    const Point_t& targpos = target.position();
    const Vector_t& targdir = target.direction();
    //check that origmom is not along the plane, i.e. targdir.Dot(origmom.Unit())=0
    if (targdir.Dot(origmom.Unit()) == 0) {
      success = false;
      return std::pair<double, double>(DBL_MAX, DBL_MAX);
    }
    //point-plane distance projected along direction orthogonal to the plane
    double sperp = targdir.Dot(targpos - origpos);
    //distance along track direction
    double s = sperp / targdir.Dot(origmom.Unit());
    success = true;
    return std::pair<double, double>(s, sperp);
  }

  void TrackStatePropagator::apply_dedx(double& pinv,
                                        detinfo::DetectorPropertiesData const& detProp,
                                        double dedx,
                                        double e1,
                                        double mass,
                                        double s,
                                        double& deriv) const
  {
    // For infinite initial momentum, return with infinite momentum.
    if (pinv == 0.) return;
    //
    const double emid = e1 - 0.5 * s * dedx;
    if (emid > mass) {
      const double pmid = std::sqrt(emid * emid - mass * mass);
      const double e2 = e1 - 0.001 * s * detProp.Eloss(pmid, mass, fTcut);
      if (e2 > mass) {
        const double p2 = std::sqrt(e2 * e2 - mass * mass);
        double pinv2 = 1. / p2;
        if (pinv < 0.) pinv2 = -pinv2;
        // derivative
        deriv = pinv2 * pinv2 * pinv2 * e2 / (pinv * pinv * pinv * e1);
        // update result.
        pinv = pinv2;
      }
    }
    return;
  }

  bool TrackStatePropagator::apply_mcs(detinfo::DetectorPropertiesData const& detProp,
                                       double dudw,
                                       double dvdw,
                                       double pinv,
                                       double mass,
                                       double s,
                                       double range,
                                       double p,
                                       double e2,
                                       bool flipSign,
                                       SMatrixSym55& noise_matrix) const
  {
    // If distance is zero, or momentum is infinite, return zero noise.

    if (pinv == 0. || s == 0.) return true;

    // Use crude estimate of the range of the track.
    if (range > 100.) range = 100.;
    const double p2 = p * p;

    // Calculate the radiation length in cm.
    const double x0 = larprop->RadiationLength() / detProp.Density();

    // Calculate projected rms scattering angle.
    // Use the estimted range in the logarithm factor.
    // Use the incremental propagation distance in the square root factor.
    const double betainv = std::sqrt(1. + pinv * pinv * mass * mass);
    const double theta_fact = (0.0136 * pinv * betainv) * (1. + 0.038 * std::log(range / x0));
    const double theta02 = theta_fact * theta_fact * std::abs(s / x0);

    // Calculate some common factors needed for multiple scattering.
    const double ufact2 = 1. + dudw * dudw;
    const double vfact2 = 1. + dvdw * dvdw;
    const double uvfact2 = 1. + dudw * dudw + dvdw * dvdw;
    const double uvfact = std::sqrt(uvfact2);
    const double uv = dudw * dvdw;
    const double dist2_3 = s * s / 3.;
    double dist_2 = std::abs(s) / 2.;
    if (flipSign) dist_2 = -dist_2;

    // Calculate energy loss fluctuations.

    const double evar = 1.e-6 * detProp.ElossVar(p, mass) * std::abs(s); // E variance (GeV^2).
    const double pinvvar = evar * e2 / (p2 * p2 * p2); // Inv. p variance (1/GeV^2)

    // Update elements of noise matrix.

    // Position submatrix.
    noise_matrix(0, 0) += dist2_3 * theta02 * ufact2; // sigma^2(u,u)
    noise_matrix(1, 0) += dist2_3 * theta02 * uv;     // sigma^2(u,v)
    noise_matrix(1, 1) += dist2_3 * theta02 * vfact2; // sigma^2(v,v)

    // Slope submatrix.
    noise_matrix(2, 2) += theta02 * uvfact2 * ufact2; // sigma^2(u', u')
    noise_matrix(3, 2) += theta02 * uvfact2 * uv;     // sigma^2(v', u')
    noise_matrix(3, 3) += theta02 * uvfact2 * vfact2; // sigma^2(v', v')

    // Same-view position-slope correlations.
    noise_matrix(2, 0) += dist_2 * theta02 * uvfact * ufact2; // sigma^2(u', u)
    noise_matrix(3, 1) += dist_2 * theta02 * uvfact * vfact2; // sigma^2(v', v)

    // Opposite-view position-slope correlations.
    noise_matrix(2, 1) += dist_2 * theta02 * uvfact * uv; // sigma^2(u', v)
    noise_matrix(3, 0) += dist_2 * theta02 * uvfact * uv; // sigma^2(v', u)

    // Momentum correlations (zero).
    // noise_matrix(4,0) += 0.;                                   // sigma^2(pinv, u)
    // noise_matrix(4,1) += 0.;                                   // sigma^2(pinv, v)
    // noise_matrix(4,2) += 0.;                                   // sigma^2(pinv, u')
    // noise_matrix(4,3) += 0.;                                   // sigma^2(pinv, v')

    // Energy loss fluctuations.
    if (fPropPinvErr) noise_matrix(4, 4) += pinvvar; // sigma^2(pinv, pinv)

    // Done (success).
    return true;
  }

} // end namespace trkf
