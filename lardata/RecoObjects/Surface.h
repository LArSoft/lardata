////////////////////////////////////////////////////////////////////////
///
/// \file   Surface.h
///
/// \brief  Base class for Kalman filter surface.
///
/// \author H. Greenlee
///
/// Surfaces may have the following distinct uses in the context
/// of the Kalman filter.
///
/// 1.  Destination for track propagation.
/// 2.  Define a local 3D coordinate systems.
/// 3.  Have a standard set of track parameters, which are meaningful
///     in the local coordinate system.
///
/// Larsoft global coordinates are (x,y,z).
///
/// Surface local coordinates are called (u,v,w), where w=0 is the
/// surface, and (u,v) are the local coordinates within the surface.
///
/// Notes about track and surface directions:
///
/// 1.  Surfaces are in general orientable, which means that a track
///     that is located at a surface can be propagating in the forward
///     direction with respect to the surface (dw/ds > 0) or in the
///     backward direction (dw/ds < 0).
/// 2.  For some kinds surfaces and track parameters, the track
///     direction is implied by the track parameters themselves.
///     For others it isn't, and must be supplied externally.
/// 3.  A surface can be queried to find the track direction implied
///     by track parameters, (via method Surface::getDirection), with
///     result returned via enum TrackDirection (possible values FORWARD,
///     BACKWARD, UNKNOWN).
/// 4.  In all situations, a direction implied by track parameters
///     has precedence over an externally supplied one.
///
/// This class doesn't have any attributes of its own, but it provides
/// several virtual methods that derived classes can or must override.
///
////////////////////////////////////////////////////////////////////////

#ifndef SURFACE_H
#define SURFACE_H

#include "lardata/RecoObjects/KalmanLinearAlgebra.h"
#include <iosfwd>

namespace trkf {

  class Surface {
  public:
    /// Track direction enum.
    enum TrackDirection { FORWARD, BACKWARD, UNKNOWN };

    /// Default constructor.
    Surface();

    /// Destructor.
    virtual ~Surface();

    // Virtual methods.

    /// Clone method.
    virtual Surface* clone() const = 0;

    /// Surface-specific tests of validity of track parameters.
    virtual bool isTrackValid(const TrackVector& vec) const = 0;

    /// Transform global to local coordinates.
    virtual void toLocal(const double xyz[3], double uvw[3]) const = 0;

    /// Transform local to global coordinates.
    virtual void toGlobal(const double uvw[3], double xyz[3]) const = 0;

    /// Calculate difference of two track parameter vectors.
    virtual TrackVector getDiff(const TrackVector& vec1, const TrackVector& vec2) const;

    /// Get position of track.
    virtual void getPosition(const TrackVector& vec, double xyz[3]) const = 0;

    /// Get direction of track (default UNKNOWN).
    virtual TrackDirection getDirection(const TrackVector& /* vec */,
                                        TrackDirection dir = UNKNOWN) const
    {
      return dir;
    }

    /// Get momentum vector of track.
    virtual void getMomentum(const TrackVector& vec,
                             double mom[3],
                             TrackDirection dir = UNKNOWN) const = 0;

    /// Get pointing error of track.
    virtual double PointingError(const TrackVector& vec, const TrackError& err) const = 0;

    /// Get starting error matrix for Kalman filter.
    virtual void getStartingError(TrackError& err) const = 0;

    /// Test whether two surfaces are parallel, within tolerance.
    virtual bool isParallel(const Surface& surf) const = 0;

    /// Find perpendicular forward distance to a parallel surface
    virtual double distanceTo(const Surface& surf) const = 0;

    /// Test two surfaces for equality, within tolerance.
    virtual bool isEqual(const Surface& surf) const = 0;

    /// Printout
    virtual std::ostream& Print(std::ostream& out) const = 0;
  };

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const Surface& surf);
}

#endif
