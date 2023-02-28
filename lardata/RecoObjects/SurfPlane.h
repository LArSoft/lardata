////////////////////////////////////////////////////////////////////////
///
/// \file   SurfPlane.h
///
/// \brief  Base class for Kalman filter planar surfaces.
///
/// \author H. Greenlee
///
/// This class acts as an intermediate layer between abstract surfaces
/// (base class Surface), and concrete planar surfaces (like
/// SurfYZPlane and SurfXYZPlane).  It does not include any data
/// members of its own.  However, it guarantees the existence of a
/// local Cartesian coordinate system (u,v,w), in which this plane is
/// located at w=0, and the track parameters are (u, v, du/dw, dv/dw,
/// q/p).  This class does not specify how the local Cartesian
/// coordinate system is related to the global coordinate system.
///
////////////////////////////////////////////////////////////////////////

#ifndef SURFPLANE_H
#define SURFPLANE_H

#include "lardata/RecoObjects/KalmanLinearAlgebra.h"
#include "lardata/RecoObjects/Surface.h"

namespace trkf {

  class SurfPlane : public Surface {
  public:
    /// Default constructor.
    SurfPlane();

    /// Destructor.
    virtual ~SurfPlane();

    // Overrides.

    /// Get pointing error of track.
    double PointingError(const TrackVector& vec, const TrackError& err) const;

    /// Get starting error matrix for Kalman filter.
    void getStartingError(TrackError& err) const;
  };
}

#endif
