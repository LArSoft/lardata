///////////////////////////////////////////////////////////////////////
///
/// \file   PropYZLine.cxx
///
/// \brief  Propagate to SurfYZLine surface.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/PropYZLine.h"
#include "cetlib_except/exception.h"
#include "lardata/RecoObjects/InteractGeneral.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "lardata/RecoObjects/SurfYZLine.h"
#include "lardata/RecoObjects/SurfYZPlane.h"
#include <cmath>

namespace trkf {

  /// Constructor.
  ///
  /// Arguments.
  ///
  /// tcut   - Delta ray energy cutoff for calculating dE/dx.
  /// doDedx - dE/dx enable flag.
  ///
  PropYZLine::PropYZLine(detinfo::DetectorPropertiesData const& detProp, double tcut, bool doDedx)
    : Propagator(detProp,
                 tcut,
                 doDedx,
                 (tcut >= 0. ? std::make_shared<const InteractGeneral>(detProp, tcut) :
                               std::shared_ptr<const Interactor>{}))
  {}

  /// Propagate without error.
  /// Optionally return propagation matrix and noise matrix.
  ///
  /// Arguments:
  ///
  /// trk   - Track to propagate.
  /// psurf - Destination surface.
  /// dir   - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  /// doDedx - dE/dx enable/disable flag.
  /// prop_matrix - Pointer to optional propagation matrix.
  /// noise_matrix - Pointer to optional noise matrix.
  ///
  /// Returned value: propagation distance + success flag.
  ///
  std::optional<double> PropYZLine::short_vec_prop(KTrack& trk,
                                                   const std::shared_ptr<const Surface>& psurf,
                                                   Propagator::PropDirection dir,
                                                   bool doDedx,
                                                   TrackMatrix* prop_matrix,
                                                   TrackError* noise_matrix) const
  {
    // Set the default return value to be unitialized with value 0.

    std::optional<double> result{std::nullopt};

    // Get destination surface and surface parameters.
    // Return failure if wrong surface type.

    const SurfYZLine* to = dynamic_cast<const SurfYZLine*>(&*psurf);
    if (to == 0) return result;
    double x02 = to->x0();
    double y02 = to->y0();
    double z02 = to->z0();
    double phi2 = to->phi();

    // Remember starting track.

    KTrack trk0(trk);

    // Get track position.

    double xyz[3];
    trk.getPosition(xyz);
    double x01 = xyz[0];
    double y01 = xyz[1];
    double z01 = xyz[2];

    // Propagate to origin surface.

    TrackMatrix local_prop_matrix;
    TrackMatrix* plocal_prop_matrix = (prop_matrix == 0 ? 0 : &local_prop_matrix);
    std::optional<double> result1 = origin_vec_prop(trk, psurf, plocal_prop_matrix);
    if (!result1) return result1;

    // Get the intermediate track state vector and track parameters.

    const TrackVector& vec = trk.getVector();
    if (vec.size() != 5)
      throw cet::exception("PropYZLine")
        << "Track state vector has wrong size" << vec.size() << "\n";
    double r1 = vec(0);
    double v1 = vec(1);
    double phid1 = vec(2);
    double eta1 = vec(3);
    double pinv = vec(4);

    // Calculate transcendental functions.

    double sinphid1 = std::sin(phid1);
    double cosphid1 = std::cos(phid1);
    double sh1 = std::sinh(eta1);
    double ch1 = std::cosh(eta1);
    double sinphi2 = std::sin(phi2);
    double cosphi2 = std::cos(phi2);

    // Calculate the initial position in the intermediate coordinate system.

    double u1 = -r1 * sinphid1;
    double w1 = r1 * cosphid1;

    // Calculate initial position in the destination coordinate system.

    double u2 = x01 - x02 + u1;
    double v2 = (y01 - y02) * cosphi2 + (z01 - z02) * sinphi2 + v1;
    double w2 = -(y01 - y02) * sinphi2 + (z01 - z02) * cosphi2 + w1;

    // Calculate the impact parameter in the destination coordinate system.

    double r2 = w2 * cosphid1 - u2 * sinphid1;

    // Calculate the perpendicular propagation distance.

    double d2 = -(w2 * sinphid1 + u2 * cosphid1);

    // Calculate the final position in the destination coordinate system.

    //double u2p = -r2 * sinphid1;
    double v2p = v2 + d2 * sh1;
    //double w2p = r2 * cosphid1;

    // Calculate the signed propagation distance.

    double s = d2 * ch1;

    // Check if propagation was in the right direction.
    // (Compare sign of s with requested direction).

    bool sok = (dir == Propagator::UNKNOWN || (dir == Propagator::FORWARD && s >= 0.) ||
                (dir == Propagator::BACKWARD && s <= 0.));

    // If wrong direction, return failure without updating the track
    // or propagation matrix.

    if (!sok) {
      trk = trk0;
      return result;
    }

    // Find final momentum.

    double deriv = 1.;
    auto pinv2 = std::make_optional(pinv);
    if (getDoDedx() && doDedx && s != 0.) {
      double* pderiv = (prop_matrix != 0 ? &deriv : 0);
      pinv2 = dedx_prop(pinv, trk.Mass(), s, pderiv);
    }

    // Return failure in case of range out.

    if (!pinv2) {
      trk = trk0;
      return result;
    }

    // Update default result to success and store propagation distance.

    result = std::make_optional(s);

    // Update propagation matrix (if requested).

    if (prop_matrix != 0) {
      TrackMatrix pm;
      pm.resize(vec.size(), vec.size(), false);

      // Calculate partial derivatives.

      pm(0, 0) = 1.; // dr2/dr1
      pm(1, 0) = 0.; // dv2/dr1
      pm(2, 0) = 0.; // d(phi2)/dr1
      pm(3, 0) = 0.; // d(eta2)/dr1
      pm(4, 0) = 0.; // d(pinv2)/dr1

      pm(0, 1) = 0.; // dr2/dv1
      pm(1, 1) = 1.; // dv2/dv1
      pm(2, 1) = 0.; // d(phi2)/dv1
      pm(3, 1) = 0.; // d(eta2)/dv1
      pm(4, 1) = 0.; // d(pinv2)/dv1

      pm(0, 2) = d2;        // dr2/d(phi1);
      pm(1, 2) = -r2 * sh1; // dv2/d(phi1);
      pm(2, 2) = 1.;        // d(phi2)/d(phi1);
      pm(3, 2) = 0.;        // d(eta2)/d(phi1);
      pm(4, 2) = 0.;        // d(pinv2)/d(phi1);

      pm(0, 3) = 0.;       // dr2/d(eta1);
      pm(1, 3) = d2 * ch1; // dv2/d(eta1);
      pm(2, 3) = 0.;       // d(phi2)/d(eta1);
      pm(3, 3) = 1.;       // d(eta2)/d(eta1);
      pm(4, 3) = 0.;       // d(pinv2)/d(eta1);

      pm(0, 4) = 0.;    // dr2/d(pinv1);
      pm(1, 4) = 0.;    // dv2/d(pinv1);
      pm(2, 4) = 0.;    // d(phi2)/d(pinv1);
      pm(3, 4) = 0.;    // d(eta2)/d(pinv1);
      pm(4, 4) = deriv; // d(pinv2)/d(pinv1);

      // Compose the final propagation matrix from zero-distance propagation and
      // parallel surface propagation.

      *prop_matrix = prod(pm, *plocal_prop_matrix);
    }

    // Update noise matrix (if requested).

    if (noise_matrix != 0) {
      noise_matrix->resize(vec.size(), vec.size(), false);
      if (getInteractor().get() != 0) {
        bool ok = getInteractor()->noise(trk, s, *noise_matrix);
        if (!ok) {
          trk = trk0;
          return std::nullopt;
        }
      }
      else
        noise_matrix->clear();
    }

    // Construct track vector at destination surface.

    TrackVector vec2(vec.size());
    vec2(0) = r2;
    vec2(1) = v2p;
    vec2(2) = phid1;
    vec2(3) = eta1;
    vec2(4) = *pinv2;

    // Update track.

    trk.setSurface(psurf);
    trk.setVector(vec2);

    // Done.

    return result;
  }

  /// Propagate without error to dynamically generated origin surface.
  /// Optionally return propagation matrix.
  ///
  /// Arguments:
  ///
  /// trk - Track to propagate.
  /// porient - Orientation surface.
  /// prop_matrix - Pointer to optional propagation matrix.
  ///
  /// Returned value: propagation distance + success flag.
  ///
  /// Propagation distance is always zero after successful propagation.
  ///
  std::optional<double> PropYZLine::origin_vec_prop(KTrack& trk,
                                                    const std::shared_ptr<const Surface>& porient,
                                                    TrackMatrix* prop_matrix) const
  {
    // Set the default return value to be unitialized with value 0.

    std::optional<double> result{std::nullopt};

    // Remember starting track.

    KTrack trk0(trk);

    // Get initial track parameters and direction.
    // Note the initial track can be on any type of surface.

    TrackVector vec = trk.getVector(); // Modifiable copy.
    if (vec.size() != 5)
      throw cet::exception("PropYZPlane")
        << "Track state vector has wrong size" << vec.size() << "\n";
    Surface::TrackDirection dir = trk.getDirection();

    // Get track position.

    double xyz[3];
    trk.getPosition(xyz);
    double x02 = xyz[0];
    double y02 = xyz[1];
    double z02 = xyz[2];

    // Generate the origin surface, which will be the destination surface.
    // Return failure if orientation surface is the wrong type.

    const SurfYZLine* orient = dynamic_cast<const SurfYZLine*>(&*porient);
    if (orient == 0) return result;
    double phi2 = orient->phi();
    std::shared_ptr<const Surface> porigin(new SurfYZLine(x02, y02, z02, phi2));

    // Test initial surface types.

    if (const SurfYZLine* from = dynamic_cast<const SurfYZLine*>(&*trk.getSurface())) {

      // Initial surface is SurfYZLine.
      // Get surface paramters.

      double phi1 = from->phi();

      // Transform track to origin surface.

      bool ok = transformYZLine(phi1, phi2, vec, dir, prop_matrix);
      result = std::make_optional(0.);
      if (!ok) return std::nullopt;
    }
    else if (const SurfYZPlane* from = dynamic_cast<const SurfYZPlane*>(&*trk.getSurface())) {

      // Initial surface is SurfYZPlane.
      // Get surface paramters.

      double phi1 = from->phi();

      // Transform track to origin surface.

      bool ok = transformYZPlane(phi1, phi2, vec, dir, prop_matrix);
      result = std::make_optional(0.);
      if (!ok) return std::nullopt;
    }
    else if (const SurfXYZPlane* from = dynamic_cast<const SurfXYZPlane*>(&*trk.getSurface())) {

      // Initial surface is SurfXYZPlane.
      // Get surface paramters.

      double theta1 = from->theta();
      double phi1 = from->phi();

      // Transform track to origin surface.

      bool ok = transformXYZPlane(theta1, phi1, phi2, vec, dir, prop_matrix);
      result = std::make_optional(0.);
      if (!ok) return std::nullopt;
    }

    // Update track.

    trk.setSurface(porigin);
    trk.setVector(vec);
    trk.setDirection(dir);

    // Final validity check.

    if (!trk.isValid()) {
      trk = trk0;
      result = std::nullopt;
    }

    // Done.

    return result;
  }

  // Transform track parameters from SurfYZLine to SurfYZLine.

  bool PropYZLine::transformYZLine(double phi1,
                                   double phi2,
                                   TrackVector& vec,
                                   Surface::TrackDirection& dir,
                                   TrackMatrix* prop_matrix) const
  {
    // Calculate surface transcendental functions.

    double sindphi = std::sin(phi2 - phi1);
    double cosdphi = std::cos(phi2 - phi1);

    // Get the initial track parameters.

    double r1 = vec(0);
    double phid1 = vec(2);
    double eta1 = vec(3);

    // Calculate elements of rotation matrix from initial coordinate
    // system to destination coordinte system.

    double rvv = cosdphi;
    double rvw = sindphi;

    double rwv = -sindphi;
    double rww = cosdphi;

    // Calculate track transcendental functions.

    double sinphid1 = std::sin(phid1);
    double cosphid1 = std::cos(phid1);
    double sh1 = 1. / std::cosh(eta1); // sech(eta1)
    double th1 = std::tanh(eta1);

    // Calculate initial position in Cartesian coordinates.

    double u1 = -r1 * sinphid1;
    double w1 = r1 * cosphid1;

    // Calculate direction in destination coordinate system.

    double du2 = sh1 * cosphid1;
    double dv2 = th1 * cosdphi + sh1 * sinphid1 * sindphi;
    double dw2 = -th1 * sindphi + sh1 * sinphid1 * cosdphi;
    double duw2 = std::hypot(du2, dw2);

    // Calculate final direction track parameters.

    double phid2 = atan2(dw2, du2);
    double eta2 = std::asinh(dv2 / duw2);

    // Update propagation matrix (if requested).

    if (prop_matrix != 0) {
      TrackMatrix& pm = *prop_matrix;
      pm.resize(vec.size(), vec.size(), false);

      // Calculate partial derivatives.

      // Partials of initial positions and directions wrt initial t.p.'s.

      double du1dr1 = -sinphid1;
      double du1dphi1 = -w1;

      double dw1dr1 = cosphid1;
      double dw1dphi1 = u1;

      double ddu1dphi1 = -sinphid1 * sh1;
      double ddu1deta1 = -cosphid1 * sh1 * th1;

      double ddv1deta1 = sh1 * sh1;

      double ddw1dphi1 = cosphid1 * sh1;
      double ddw1deta1 = -sinphid1 * sh1 * th1;

      // Rotate partials to destination coordinate system.

      double du2dr1 = du1dr1;
      double dv2dr1 = rvw * dw1dr1;
      double dw2dr1 = rww * dw1dr1;

      double dv2dv1 = rvv;
      double dw2dv1 = rwv;

      double du2dphi1 = du1dphi1;
      double dv2dphi1 = rvw * dw1dphi1;
      double dw2dphi1 = rww * dw1dphi1;

      double ddu2dphi1 = ddu1dphi1;
      double ddv2dphi1 = rvw * ddw1dphi1;
      double ddw2dphi1 = rww * ddw1dphi1;

      double ddu2deta1 = ddu1deta1;
      double ddv2deta1 = rvv * ddv1deta1 + rvw * ddw1deta1;
      double ddw2deta1 = rwv * ddv1deta1 + rww * ddw1deta1;

      // Partials of final t.p. wrt final position and direction.

      double dr2du2 = -dw2 / duw2;
      double dr2dw2 = du2 / duw2;

      double dphi2ddu2 = -dw2 / (duw2 * duw2);
      double dphi2ddw2 = du2 / (duw2 * duw2);

      double deta2ddv2 = 1. / (duw2 * duw2);

      // Partials of final t.p. wrt initial t.p.

      double dr2dr1 = dr2du2 * du2dr1 + dr2dw2 * dw2dr1;
      double dr2dv1 = dr2dw2 * dw2dv1;
      double dr2dphi1 = dr2du2 * du2dphi1 + dr2dw2 * dw2dphi1;

      double dphi2dphi1 = dphi2ddu2 * ddu2dphi1 + dphi2ddw2 * ddw2dphi1;
      double dphi2deta1 = dphi2ddu2 * ddu2deta1 + dphi2ddw2 * ddw2deta1;

      double deta2dphi1 = deta2ddv2 * ddv2dphi1;
      double deta2deta1 = deta2ddv2 * ddv2deta1;

      // We still need to calculate the corretion due to the dependence of the
      // propagation distance on the initial track parameters.  This correction is
      // needed even though the actual propagation distance is zero.

      // This correction only effects the v track parameter, since the v parameter
      // the only parameter that actually dependes on the propagation distance.

      // Partials of propagation distance wrt position and direction in the destination
      // coordinate system.

      double dsdu2 = -du2 / (duw2 * duw2);
      double dsdw2 = -dw2 / (duw2 * duw2);

      // Partials of propagation distance wrt initial t.p.

      double dsdr1 = dsdu2 * du2dr1 + dsdw2 * dw2dr1;
      double dsdv1 = dsdw2 * dw2dv1;
      double dsdphi1 = dsdu2 * du2dphi1 + dsdw2 * dw2dphi1;

      // Calculate correction to v parameter partials wrt initial t.p. due to path length.

      dv2dr1 += dv2 * dsdr1;
      dv2dv1 += dv2 * dsdv1;
      dv2dphi1 += dv2 * dsdphi1;

      // Fill derivative matrix.

      pm(0, 0) = dr2dr1; // dr2/dr1
      pm(1, 0) = dv2dr1; // dv2/dr1
      pm(2, 0) = 0.;     // d(phi2)/dr1
      pm(3, 0) = 0.;     // d(eta2)/dr1
      pm(4, 0) = 0.;     // d(pinv2)/dr1

      pm(0, 1) = dr2dv1; // dr2/dv1
      pm(1, 1) = dv2dv1; // dv2/dv1
      pm(2, 1) = 0.;     // d(phi2)/dv1
      pm(3, 1) = 0.;     // d(eta2)/dv1
      pm(4, 1) = 0.;     // d(pinv2)/dv1

      pm(0, 2) = dr2dphi1;   // dr2/d(phi1);
      pm(1, 2) = dv2dphi1;   // dv2/d(phi1);
      pm(2, 2) = dphi2dphi1; // d(phi2)/d(phi1);
      pm(3, 2) = deta2dphi1; // d(eta2)/d(phi1);
      pm(4, 2) = 0.;         // d(pinv2)/d(phi1);

      pm(0, 3) = 0.;         // dr2/d(eta1);
      pm(1, 3) = 0.;         // dv2/d(eta1);
      pm(2, 3) = dphi2deta1; // d(phi2)/d(eta1);
      pm(3, 3) = deta2deta1; // d(eta2)/d(eta1);
      pm(4, 3) = 0.;         // d(pinv2)/d(eta1);

      pm(0, 4) = 0.; // dr2/d(pinv1);
      pm(1, 4) = 0.; // dv2/d(pinv1);
      pm(2, 4) = 0.; // d(phi2)/d(pinv1);
      pm(3, 4) = 0.; // d(eta2)/d(pinv1);
      pm(4, 4) = 1.; // d(pinv2)/d(pinv1);
    }

    // Update track vector.

    vec(0) = 0.;
    vec(1) = 0.;
    vec(2) = phid2;
    vec(3) = eta2;

    // Done (success).

    return true;
  }

  // Transform track parameters from SurfYZPlane to SurfYZLine.

  bool PropYZLine::transformYZPlane(double phi1,
                                    double phi2,
                                    TrackVector& vec,
                                    Surface::TrackDirection& dir,
                                    TrackMatrix* prop_matrix) const
  {
    // Calculate surface transcendental functions.

    double sindphi = std::sin(phi2 - phi1);
    double cosdphi = std::cos(phi2 - phi1);

    // Get the initial track parameters.

    double dudw1 = vec(2);
    double dvdw1 = vec(3);

    // Make sure initial track has a valid direction.

    double dirf = 1.;
    if (dir == Surface::BACKWARD)
      dirf = -1.;
    else if (dir != Surface::FORWARD)
      return false;

    // Calculate elements of rotation matrix from initial coordinate
    // system to destination coordinte system.

    double rvv = cosdphi;
    double rvw = sindphi;

    double rwv = -sindphi;
    double rww = cosdphi;

    // Calculate direction in the starting coordinate system.

    double dw1 = dirf / std::sqrt(1. + dudw1 * dudw1 + dvdw1 * dvdw1);
    double du1 = dudw1 * dw1;
    double dv1 = dvdw1 * dw1;

    // Rotate direction vector into destination coordinate system.

    double du2 = du1;
    double dv2 = rvv * dv1 + rvw * dw1;
    double dw2 = rwv * dv1 + rww * dw1;
    double duw2 = std::hypot(du2, dw2);

    // Calculate final direction track parameters.

    double phid2 = atan2(dw2, du2);
    double eta2 = std::asinh(dv2 / duw2);

    // Update propagation matrix (if requested).

    if (prop_matrix != 0) {
      TrackMatrix& pm = *prop_matrix;
      pm.resize(vec.size(), vec.size(), false);

      // Calculate partial derivatives.

      // Partials of initial positions and directions wrt initial t.p.'s.

      double ddu1ddudw1 = (1. + dvdw1 * dvdw1) * dw1 * dw1 * dw1;
      double ddu1ddvdw1 = -dudw1 * dvdw1 * dw1 * dw1 * dw1;

      double ddv1ddudw1 = -dudw1 * dvdw1 * dw1 * dw1 * dw1;
      double ddv1ddvdw1 = (1. + dudw1 * dudw1) * dw1 * dw1 * dw1;

      double ddw1ddudw1 = -dudw1 * dw1 * dw1 * dw1;
      double ddw1ddvdw1 = -dvdw1 * dw1 * dw1 * dw1;

      // Rotate partials to destination coordinate system.

      double dv2dv1 = rvv;
      double dw2dv1 = rwv;

      double ddu2ddudw1 = ddu1ddudw1;
      double ddv2ddudw1 = rvv * ddv1ddudw1 + rvw * ddw1ddudw1;
      double ddw2ddudw1 = rwv * ddv1ddudw1 + rww * ddw1ddudw1;

      double ddu2ddvdw1 = ddu1ddvdw1;
      double ddv2ddvdw1 = rvv * ddv1ddvdw1 + rvw * ddw1ddvdw1;
      double ddw2ddvdw1 = rwv * ddv1ddvdw1 + rww * ddw1ddvdw1;

      // Partials of final t.p. wrt final position and direction.

      double dr2du2 = -dw2 / duw2;
      double dr2dw2 = du2 / duw2;

      double dphi2ddu2 = -dw2 / (duw2 * duw2);
      double dphi2ddw2 = du2 / (duw2 * duw2);

      double deta2ddv2 = 1. / (duw2 * duw2);

      // Partials of final t.p. wrt initial t.p.

      double dr2du1 = dr2du2;
      double dr2dv1 = dr2dw2 * dw2dv1;

      double dphi2ddudw1 = dphi2ddu2 * ddu2ddudw1 + dphi2ddw2 * ddw2ddudw1;
      double dphi2ddvdw1 = dphi2ddu2 * ddu2ddvdw1 + dphi2ddw2 * ddw2ddvdw1;

      double deta2ddudw1 = deta2ddv2 * ddv2ddudw1;
      double deta2ddvdw1 = deta2ddv2 * ddv2ddvdw1;

      // We still need to calculate the corretion due to the dependence of the
      // propagation distance on the initial track parameters.  This correction is
      // needed even though the actual propagation distance is zero.

      // This correction only effects the v track parameter, since the v parameter
      // the only parameter that actually dependes on the propagation distance.

      // Partials of propagation distance wrt position and direction in the destination
      // coordinate system.

      double dsdu2 = -du2 / (duw2 * duw2);
      double dsdw2 = -dw2 / (duw2 * duw2);

      // Partials of propagation distance wrt initial t.p.

      double dsdu1 = dsdu2;
      double dsdv1 = dsdw2 * dw2dv1;

      // Calculate correction to v parameter partials wrt initial t.p. due to path length.

      double dv2du1 = dv2 * dsdu1;
      dv2dv1 += dv2 * dsdv1;

      // Fill matrix.

      pm(0, 0) = dr2du1; // dr2/du1
      pm(1, 0) = dv2du1; // dv2/du1
      pm(2, 0) = 0.;     // d(phi2)/du1
      pm(3, 0) = 0.;     // d(eta2)/du1
      pm(4, 0) = 0.;     // d(pinv2)/du1

      pm(0, 1) = dr2dv1; // dr2/dv1
      pm(1, 1) = dv2dv1; // dv2/dv1
      pm(2, 1) = 0.;     // d(phi2)/dv1
      pm(3, 1) = 0.;     // d(eta2)/dv1
      pm(4, 1) = 0.;     // d(pinv2)/dv1

      pm(0, 2) = 0.;          // dr2/d(dudw1);
      pm(1, 2) = 0.;          // dv2/d(dudw1);
      pm(2, 2) = dphi2ddudw1; // d(dudw2)/d(dudw1);
      pm(3, 2) = deta2ddudw1; // d(eta2)/d(dudw1);
      pm(4, 2) = 0.;          // d(pinv2)/d(dudw1);

      pm(0, 3) = 0.;          // dr2/d(dvdw1);
      pm(1, 3) = 0.;          // dv2/d(dvdw1);
      pm(2, 3) = dphi2ddvdw1; // d(phi2)/d(dvdw1);
      pm(3, 3) = deta2ddvdw1; // d(eta2)/d(dvdw1);
      pm(4, 3) = 0.;          // d(pinv2)/d(dvdw1);

      pm(0, 4) = 0.; // dr2/d(pinv1);
      pm(1, 4) = 0.; // dv2/d(pinv1);
      pm(2, 4) = 0.; // d(phi2)/d(pinv1);
      pm(3, 4) = 0.; // d(eta2)/d(pinv1);
      pm(4, 4) = 1.; // d(pinv2)/d(pinv1);
    }

    // Update track vector.

    vec(0) = 0.;
    vec(1) = 0.;
    vec(2) = phid2;
    vec(3) = eta2;

    // Done (success).

    return true;
  }

  // Transform track parameters from SurfXYZPlane to SurfYZLine.

  bool PropYZLine::transformXYZPlane(double theta1,
                                     double phi1,
                                     double phi2,
                                     TrackVector& vec,
                                     Surface::TrackDirection& dir,
                                     TrackMatrix* prop_matrix) const
  {
    // Calculate surface transcendental functions.

    double sinth1 = std::sin(theta1);
    double costh1 = std::cos(theta1);

    double sindphi = std::sin(phi2 - phi1);
    double cosdphi = std::cos(phi2 - phi1);

    // Get the initial track parameters.

    double dudw1 = vec(2);
    double dvdw1 = vec(3);

    // Make sure initial track has a valid direction.

    double dirf = 1.;
    if (dir == Surface::BACKWARD)
      dirf = -1.;
    else if (dir != Surface::FORWARD)
      return false;

    // Calculate elements of rotation matrix from initial coordinate
    // system to destination coordinte system.

    double ruu = costh1;
    double ruw = sinth1;

    double rvu = -sinth1 * sindphi;
    double rvv = cosdphi;
    double rvw = costh1 * sindphi;

    double rwu = -sinth1 * cosdphi;
    double rwv = -sindphi;
    double rww = costh1 * cosdphi;

    // Calculate direction in the starting coordinate system.

    double dw1 = dirf / std::sqrt(1. + dudw1 * dudw1 + dvdw1 * dvdw1);
    double du1 = dudw1 * dw1;
    double dv1 = dvdw1 * dw1;

    // Rotate direction vector into destination coordinate system.

    double du2 = ruu * du1 + ruw * dw1;
    double dv2 = rvu * du1 + rvv * dv1 + rvw * dw1;
    double dw2 = rwu * du1 + rwv * dv1 + rww * dw1;
    double duw2 = std::hypot(du2, dw2);

    // Calculate final direction track parameters.

    double phid2 = atan2(dw2, du2);
    double eta2 = std::asinh(dv2 / duw2);

    // Update propagation matrix (if requested).

    if (prop_matrix != 0) {
      TrackMatrix& pm = *prop_matrix;
      pm.resize(vec.size(), vec.size(), false);

      // Calculate partial derivatives.

      // Partials of initial positions and directions wrt initial t.p.'s.

      double ddu1ddudw1 = (1. + dvdw1 * dvdw1) * dw1 * dw1 * dw1;
      double ddu1ddvdw1 = -dudw1 * dvdw1 * dw1 * dw1 * dw1;

      double ddv1ddudw1 = -dudw1 * dvdw1 * dw1 * dw1 * dw1;
      double ddv1ddvdw1 = (1. + dudw1 * dudw1) * dw1 * dw1 * dw1;

      double ddw1ddudw1 = -dudw1 * dw1 * dw1 * dw1;
      double ddw1ddvdw1 = -dvdw1 * dw1 * dw1 * dw1;

      // Rotate partials to destination coordinate system.

      double du2du1 = ruu;
      double dv2du1 = rvu;
      double dw2du1 = rwu;

      double dv2dv1 = rvv;
      double dw2dv1 = rwv;

      double ddu2ddudw1 = ruu * ddu1ddudw1 + ruw * ddw1ddudw1;
      double ddv2ddudw1 = rvu * ddu1ddudw1 + rvv * ddv1ddudw1 + rvw * ddw1ddudw1;
      double ddw2ddudw1 = rwu * ddu1ddudw1 + rwv * ddv1ddudw1 + rww * ddw1ddudw1;

      double ddu2ddvdw1 = ruu * ddu1ddvdw1 + ruw * ddw1ddvdw1;
      double ddv2ddvdw1 = rvu * ddu1ddvdw1 + rvv * ddv1ddvdw1 + rvw * ddw1ddvdw1;
      double ddw2ddvdw1 = rwu * ddu1ddvdw1 + rwv * ddv1ddvdw1 + rww * ddw1ddvdw1;

      // Partials of final t.p. wrt final position and direction.

      double dr2du2 = -dw2 / duw2;
      double dr2dw2 = du2 / duw2;

      double dphi2ddu2 = -dw2 / (duw2 * duw2);
      double dphi2ddw2 = du2 / (duw2 * duw2);

      double deta2ddv2 = 1. / (duw2 * duw2);

      // Partials of final t.p. wrt initial t.p.

      double dr2du1 = dr2du2 * du2du1 + dr2dw2 * dw2du1;
      double dr2dv1 = dr2dw2 * dw2dv1;

      double dphi2ddudw1 = dphi2ddu2 * ddu2ddudw1 + dphi2ddw2 * ddw2ddudw1;
      double dphi2ddvdw1 = dphi2ddu2 * ddu2ddvdw1 + dphi2ddw2 * ddw2ddvdw1;

      double deta2ddudw1 = deta2ddv2 * ddv2ddudw1;
      double deta2ddvdw1 = deta2ddv2 * ddv2ddvdw1;

      // We still need to calculate the corretion due to the dependence of the
      // propagation distance on the initial track parameters.  This correction is
      // needed even though the actual propagation distance is zero.

      // This correction only effects the v track parameter, since the v parameter
      // the only parameter that actually dependes on the propagation distance.

      // Partials of propagation distance wrt position and direction in the destination
      // coordinate system.

      double dsdu2 = -du2 / (duw2 * duw2);
      double dsdw2 = -dw2 / (duw2 * duw2);

      // Partials of propagation distance wrt initial t.p.

      double dsdu1 = dsdu2 * du2du1 + dsdw2 * dw2du1;
      double dsdv1 = dsdw2 * dw2dv1;

      // Calculate correction to v parameter partials wrt initial t.p. due to path length.

      dv2du1 += dv2 * dsdu1;
      dv2dv1 += dv2 * dsdv1;

      // Fill matrix.

      pm(0, 0) = dr2du1; // dr2/du1
      pm(1, 0) = dv2du1; // dv2/du1
      pm(2, 0) = 0.;     // d(phi2)/du1
      pm(3, 0) = 0.;     // d(eta2)/du1
      pm(4, 0) = 0.;     // d(pinv2)/du1

      pm(0, 1) = dr2dv1; // dr2/dv1
      pm(1, 1) = dv2dv1; // dv2/dv1
      pm(2, 1) = 0.;     // d(phi2)/dv1
      pm(3, 1) = 0.;     // d(eta2)/dv1
      pm(4, 1) = 0.;     // d(pinv2)/dv1

      pm(0, 2) = 0.;          // dr2/d(dudw1);
      pm(1, 2) = 0.;          // dv2/d(dudw1);
      pm(2, 2) = dphi2ddudw1; // d(dudw2)/d(dudw1);
      pm(3, 2) = deta2ddudw1; // d(eta2)/d(dudw1);
      pm(4, 2) = 0.;          // d(pinv2)/d(dudw1);

      pm(0, 3) = 0.;          // dr2/d(dvdw1);
      pm(1, 3) = 0.;          // dv2/d(dvdw1);
      pm(2, 3) = dphi2ddvdw1; // d(phi2)/d(dvdw1);
      pm(3, 3) = deta2ddvdw1; // d(eta2)/d(dvdw1);
      pm(4, 3) = 0.;          // d(pinv2)/d(dvdw1);

      pm(0, 4) = 0.; // dr2/d(pinv1);
      pm(1, 4) = 0.; // dv2/d(pinv1);
      pm(2, 4) = 0.; // d(phi2)/d(pinv1);
      pm(3, 4) = 0.; // d(eta2)/d(pinv1);
      pm(4, 4) = 1.; // d(pinv2)/d(pinv1);
    }

    // Update track vector.

    vec(0) = 0.;
    vec(1) = 0.;
    vec(2) = phid2;
    vec(3) = eta2;

    // Done (success).

    return true;
  }
} // end namespace trkf
