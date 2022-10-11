///////////////////////////////////////////////////////////////////////
///
/// \file   PropAny.cxx
///
/// \brief  Propagate between any two surfaces.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/PropAny.h"
#include "cetlib_except/exception.h"
#include "lardata/RecoObjects/InteractPlane.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "lardata/RecoObjects/SurfYZLine.h"
#include "lardata/RecoObjects/SurfYZPlane.h"

namespace trkf {

  /// Constructor.
  ///
  /// Arguments.
  ///
  /// tcut   - Delta ray energy cutoff for calculating dE/dx.
  /// doDedx - dE/dx enable flag.
  ///
  PropAny::PropAny(detinfo::DetectorPropertiesData const& detProp, double tcut, bool doDedx)
    : Propagator(detProp,
                 tcut,
                 doDedx,
                 (tcut >= 0 ? std::make_shared<InteractPlane const>(detProp, tcut) :
                              std::shared_ptr<Interactor const>{}))
    , fPropYZLine(detProp, tcut, doDedx)
    , fPropYZPlane(detProp, tcut, doDedx)
    , fPropXYZPlane(detProp, tcut, doDedx)
  {}

  /// Propagate without error.
  /// Optionally return propagation matrix and noise matrix.
  /// This method tests the type of the destination surface, and calls
  /// the corresponding typed propagator.
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
  std::optional<double> PropAny::short_vec_prop(KTrack& trk,
                                                const std::shared_ptr<const Surface>& psurf,
                                                Propagator::PropDirection dir,
                                                bool doDedx,
                                                TrackMatrix* prop_matrix,
                                                TrackError* noise_matrix) const
  {
    // Default result.

    std::optional<double> result{std::nullopt};

    // Test the type of the destination surface.

    if (dynamic_cast<const SurfYZLine*>(&*psurf))
      result = fPropYZLine.short_vec_prop(trk, psurf, dir, doDedx, prop_matrix, noise_matrix);
    else if (dynamic_cast<const SurfYZPlane*>(&*psurf))
      result = fPropYZPlane.short_vec_prop(trk, psurf, dir, doDedx, prop_matrix, noise_matrix);
    else if (dynamic_cast<const SurfXYZPlane*>(&*psurf))
      result = fPropXYZPlane.short_vec_prop(trk, psurf, dir, doDedx, prop_matrix, noise_matrix);
    else
      throw cet::exception("PropAny") << "Destination surface has unknown type.\n";

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
  std::optional<double> PropAny::origin_vec_prop(KTrack& trk,
                                                 const std::shared_ptr<const Surface>& porient,
                                                 TrackMatrix* prop_matrix) const
  {
    // Default result.

    std::optional<double> result{std::nullopt};

    // Test the type of the destination surface.

    if (dynamic_cast<const SurfYZLine*>(&*porient))
      result = fPropYZLine.origin_vec_prop(trk, porient, prop_matrix);
    else if (dynamic_cast<const SurfYZPlane*>(&*porient))
      result = fPropYZPlane.origin_vec_prop(trk, porient, prop_matrix);
    else if (dynamic_cast<const SurfXYZPlane*>(&*porient))
      result = fPropXYZPlane.origin_vec_prop(trk, porient, prop_matrix);
    else
      throw cet::exception("PropAny") << "Destination surface has unknown type.\n";

    // Done.

    return result;
  }

} // end namespace trkf
