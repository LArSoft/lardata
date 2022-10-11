///////////////////////////////////////////////////////////////////////
///
/// \file   InteractPlane.cxx
///
/// \brief  Interactor for planar surfaces.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "cetlib_except/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/RecoObjects/InteractPlane.h"
#include "lardata/RecoObjects/SurfPlane.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

namespace trkf {

  /// Constructor.
  ///
  /// Arguments:
  ///
  /// tcut - Maximum delta ray energy.
  ///
  InteractPlane::InteractPlane(detinfo::DetectorPropertiesData const& detProp, double tcut)
    : Interactor(tcut), fDetProp{detProp}
  {}

  /// Calculate noise matrix.
  ///
  /// Arguments:
  ///
  /// trk          - Original track.
  /// s            - Path distance.
  /// noise_matrix - Resultant noise matrix.
  ///
  /// Returns: True if success.
  ///
  /// Currently calculate noise from multiple scattering only.
  ///
  /// Note about multiple scattering calculation:
  ///
  /// In the case of normal incident track (u' = v' = 0), the multiple
  /// scattering calculations used in this class reduce to the
  /// familiar small-angle formulas found in the particle data book
  /// for a thick scatterer (meaning multiple scattering modifies both
  /// position and slope errors).  However, the distance used in the
  /// logarithm factor of the rms scattering angle is an estimate of
  /// the total track length, rather than the incremental track
  /// length.
  ///
  /// For non-normal incident track, the error ellipse is elongated in
  /// the radial direciton of the uv plane by factor sqrt(1 + u'^2 +
  /// v'^2).  This is equivalent to expansion in the u direction by
  /// factor sqrt(1 + u'^2), and expansion in the v direction by
  /// sqrt(1 + v'^2), with uv correlation u' v' / sqrt((1 + u'^2)(1 +
  /// v'^2)).
  ///
  /// Correlation between position and slope in the same view is
  /// sqrt(3)/2 regardless of normal incidence.
  ///
  /// Correlation between position and slope in the opposite view is
  /// (sqrt(3)/2) u' v' / sqrt((1 + u'^2)(1 + v'^2))
  ///
  bool InteractPlane::noise(const KTrack& trk, double s, TrackError& noise_matrix) const
  {
    // Get LAr service.

    auto const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();

    // Make sure we are on a plane surface (throw exception if not).

    const SurfPlane* psurf = dynamic_cast<const SurfPlane*>(&*trk.getSurface());
    if (psurf == nullptr)
      throw cet::exception("InteractPlane") << "InteractPlane called for non-planar surface.\n";

    // Clear noise matrix.

    noise_matrix.clear();

    // Unpack track parameters.

    const TrackVector& vec = trk.getVector();
    double dudw = vec[2];
    double dvdw = vec[3];
    double pinv = vec[4];
    double mass = trk.Mass();

    // If distance is zero, or momentum is infinite, return zero noise.

    if (pinv == 0. || s == 0.) return true;

    // Make a crude estimate of the range of the track.

    double p = 1. / std::abs(pinv);
    double p2 = p * p;
    double e2 = p2 + mass * mass;
    double e = std::sqrt(e2);
    double t = e - mass;
    double dedx = 0.001 * fDetProp.Eloss(p, mass, getTcut());
    double range = t / dedx;
    if (range > 100.) range = 100.;

    // Calculate the radiation length in cm.

    double x0 = larprop->RadiationLength() / fDetProp.Density();

    // Calculate projected rms scattering angle.
    // Use the estimted range in the logarithm factor.
    // Use the incremental propagation distance in the square root factor.

    double betainv = std::sqrt(1. + pinv * pinv * mass * mass);
    double theta_fact = (0.0136 * pinv * betainv) * (1. + 0.038 * std::log(range / x0));
    double theta02 = theta_fact * theta_fact * std::abs(s / x0);

    // Calculate some sommon factors needed for multiple scattering.

    double ufact2 = 1. + dudw * dudw;
    double vfact2 = 1. + dvdw * dvdw;
    double uvfact2 = 1. + dudw * dudw + dvdw * dvdw;
    double uvfact = std::sqrt(uvfact2);
    double uv = dudw * dvdw;
    double dist2_3 = s * s / 3.;
    double dist_2 = std::abs(s) / 2.;
    if (trk.getDirection() == Surface::BACKWARD) dist_2 = -dist_2;

    // Calculate energy loss fluctuations.

    double evar = 1.e-6 * fDetProp.ElossVar(p, mass) * std::abs(s); // E variance (GeV^2).
    double pinvvar = evar * e2 / (p2 * p2 * p2);                    // Inv. p variance (1/GeV^2)

    // Fill elements of noise matrix.

    // Position submatrix.

    noise_matrix(0, 0) = dist2_3 * theta02 * ufact2; // sigma^2(u,u)
    noise_matrix(1, 0) = dist2_3 * theta02 * uv;     // sigma^2(u,v)
    noise_matrix(1, 1) = dist2_3 * theta02 * vfact2; // sigma^2(v,v)

    // Slope submatrix.

    noise_matrix(2, 2) = theta02 * uvfact2 * ufact2; // sigma^2(u', u')
    noise_matrix(3, 2) = theta02 * uvfact2 * uv;     // sigma^2(v', u')
    noise_matrix(3, 3) = theta02 * uvfact2 * vfact2; // sigma^2(v', v')

    // Same-view position-slope correlations.

    noise_matrix(2, 0) = dist_2 * theta02 * uvfact * ufact2; // sigma^2(u', u)
    noise_matrix(3, 1) = dist_2 * theta02 * uvfact * vfact2; // sigma^2(v', v)

    // Opposite-view position-slope correlations.

    noise_matrix(2, 1) = dist_2 * theta02 * uvfact * uv; // sigma^2(u', v)
    noise_matrix(3, 0) = dist_2 * theta02 * uvfact * uv; // sigma^2(v', u)

    // Momentum correlations (zero).

    noise_matrix(4, 0) = 0.; // sigma^2(pinv, u)
    noise_matrix(4, 1) = 0.; // sigma^2(pinv, v)
    noise_matrix(4, 2) = 0.; // sigma^2(pinv, u')
    noise_matrix(4, 3) = 0.; // sigma^2(pinv, v')

    // Energy loss fluctuations.

    noise_matrix(4, 4) = pinvvar; // sigma^2(pinv, pinv)

    // Done (success).

    return true;
  }

} // end namespace trkf
