///////////////////////////////////////////////////////////////////////
///
/// \file   Surface.cxx
///
/// \brief  Base class for Kalman filter surface.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/Surface.h"

namespace trkf {

  /// Default constructor.
  Surface::Surface() {}

  /// Destructor.
  Surface::~Surface() {}

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const Surface& surf)
  {
    return surf.Print(out);
  }

  /// Calculate difference of two track parameter vectors.
  /// This method has a default implementation which is just the numeric difference.
  /// Surfaces that require a more sophisticated difference (e.g. phi-wrap difference)
  /// should override this method.
  TrackVector Surface::getDiff(const TrackVector& vec1, const TrackVector& vec2) const
  {
    return vec1 - vec2;
  }

} // end namespace trkf
