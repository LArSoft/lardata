///////////////////////////////////////////////////////////////////////
///
/// \file   KHitBase.cxx
///
/// \brief  Base class for Kalman filter measurement.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KHitBase.h"

#include <ostream>

namespace trkf {

  /// Default Constructor.
  KHitBase::KHitBase() : fPredDist(0.), fID(0), fMeasPlane(-1) {}

  /// Initializing Constructor.
  ///
  /// Arguments:
  ///
  /// psurf - Measurement surface pointer.
  ///
  KHitBase::KHitBase(const std::shared_ptr<const Surface>& psurf, int plane)
    : fPredDist(0.), fID(0), fMeasSurf(psurf), fMeasPlane(plane)
  {}

  /// Printout
  std::ostream& KHitBase::Print(std::ostream& out, bool doTitle) const
  {
    if (doTitle) out << "KHitBase:\n";
    out << "  Measurement Surface: " << *fMeasSurf << "\n";
    out << "  Measurement Plane: " << fMeasPlane << "\n";
    if (fPredSurf.get() != 0) {
      out << "  Prediction Surface: " << *fPredSurf << "\n";
      out << "  Prediction Distance: " << fPredDist << "\n";
    }
    return out;
  }

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const KHitBase& trk) { return trk.Print(out); }

} // end namespace trkf
