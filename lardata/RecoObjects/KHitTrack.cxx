///////////////////////////////////////////////////////////////////////
///
/// \file   KHitTrack.cxx
///
/// \brief  Basic Kalman filter track class, with measurements.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KHitTrack.h"

namespace trkf {

  /// Default constructor.
  KHitTrack::KHitTrack() {}

  /// Initializing constructor - KFitTrack + measurement.
  ///
  /// Arguments:
  ///
  /// trf - KFitTrack.
  /// hit - Measurement.
  ///
  KHitTrack::KHitTrack(const KFitTrack& trf, const std::shared_ptr<const KHitBase>& hit)
    : KFitTrack(trf), fHit(hit)
  {}

  /// Initializing constructor - KETrack.
  ///
  /// Arguments:
  ///
  /// tre - KETrack.
  ///
  KHitTrack::KHitTrack(const KETrack& tre) : KFitTrack(tre) {}

  /// Destructor.
  KHitTrack::~KHitTrack() {}

  /// Printout
  std::ostream& KHitTrack::Print(std::ostream& out, bool doTitle) const
  {
    if (doTitle) out << "KHitTrack:\n";

    // Print base class.

    KFitTrack::Print(out, false);

    // Print information specific to this class.

    if (fHit.get() != 0) out << "  " << *fHit;
    return out;
  }

} // end namespace trkf
