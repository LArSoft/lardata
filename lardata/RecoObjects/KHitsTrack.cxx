///////////////////////////////////////////////////////////////////////
///
/// \file   KHitsTrack.cxx
///
/// \brief  Basic Kalman filter track class, with measurements.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KHitsTrack.h"

namespace trkf {

  /// Default constructor.
  KHitsTrack::KHitsTrack() {}

  /// Initializing constructor - KFitTrack.
  ///
  /// Arguments:
  ///
  /// trf - KFitTrack.
  ///
  KHitsTrack::KHitsTrack(const KFitTrack& trf) : KFitTrack(trf) {}

  /// Initializing constructor - KETrack.
  ///
  /// Arguments:
  ///
  /// tre - KETrack.
  ///
  KHitsTrack::KHitsTrack(const KETrack& tre) : KFitTrack(tre) {}

  /// Destructor.
  KHitsTrack::~KHitsTrack() {}

  /// Printout
  std::ostream& KHitsTrack::Print(std::ostream& out, bool doTitle) const
  {
    if (doTitle) out << "KHitsTrack:\n";

    // Print information specific to this class.

    out << "  " << fHits.size() << " hits.\n";

    // Print base class.

    KFitTrack::Print(out, false);
    return out;
  }

} // end namespace trkf
