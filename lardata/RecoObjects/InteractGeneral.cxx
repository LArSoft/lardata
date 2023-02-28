///////////////////////////////////////////////////////////////////////
///
/// \file   InteractGeneral.cxx
///
/// \brief  Interactor for general surfaces.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/InteractGeneral.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include <cmath>

namespace trkf {

  /// Constructor.
  ///
  /// Arguments:
  ///
  /// tcut - Maximum delta ray energy.
  ///
  InteractGeneral::InteractGeneral(detinfo::DetectorPropertiesData const& detProp, double tcut)
    : Interactor(tcut), fInteract(detProp, tcut), fProp(detProp, -1., false)
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
  /// We make a zero distance propagation to a plane surface
  /// (SurfXYZPlane) that is normal to the track direction.
  /// Then calculate the noise matrix on that surface and
  /// transform back to the original surface.
  ///
  bool InteractGeneral::noise(const KTrack& trk, double s, TrackError& noise_matrix) const
  {
    // Get track position and direction.

    double xyz[3];
    double mom[3];
    trk.getPosition(xyz);
    trk.getMomentum(mom);

    // Generate a SurfXYZPlane with origin at current track position, and
    // normal to current track direction.

    std::shared_ptr<Surface> psurf(
      new SurfXYZPlane(xyz[0], xyz[1], xyz[2], mom[0], mom[1], mom[2]));

    // Propagate track to newly created surface (zero-distance propagation).

    TrackMatrix prop_matrix;
    KTrack temp_trk = trk;
    std::optional<double> result =
      fProp.short_vec_prop(temp_trk, psurf, Propagator::UNKNOWN, false, &prop_matrix);

    // Return failure if propagation did not succeed.

    if (!result) return false;

    // Calculate noise on plane surface.

    TrackError plane_noise(5);
    fInteract.noise(temp_trk, s, plane_noise);

    // Transform noise matrix to original surface using inverse of propagation matrix.

    invert(prop_matrix);
    TrackMatrix temp = prod(plane_noise, trans(prop_matrix));
    TrackMatrix temp2 = prod(prop_matrix, temp);
    noise_matrix = ublas::symmetric_adaptor<TrackMatrix>(temp2);

    // Done (success).

    return true;
  }
} // end namespace trkf
