///////////////////////////////////////////////////////////////////////
///
/// \file   SurfPlane.cxx
///
/// \brief  Base class for Kalman filter planar surfaces.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/SurfPlane.h"
#include <cmath>

namespace trkf {

  /// Default constructor.
  SurfPlane::SurfPlane() {}

  /// Destructor.
  SurfPlane::~SurfPlane() {}

  /// Get pointing error of track.
  ///
  /// Arguments:
  ///
  /// vec - Track parameters.
  /// err - Track error matrix.
  ///
  /// Returns: Pointing error.
  ///
  /// This method calculates the track pointing error based on the
  /// slope track paramers and errors (parameters 2 and 3).
  ///
  double SurfPlane::PointingError(const TrackVector& vec, const TrackError& err) const
  {
    // Get slope parameters and error matrix.

    double xp = vec(2);
    double yp = vec(3);
    double exx = err(2, 2);
    double eyy = err(3, 3);
    double exy = err(3, 2);

    // Calculate error matrix of pointing unit vector in some coordinate system.

    double den = 1. + xp * xp + yp * yp;
    double den3 = den * den * den;

    double vxx = ((1. + yp * yp) * (1. + yp * yp) * exx + xp * xp * yp * yp * eyy -
                  2. * xp * yp * (1. + yp * yp) * exy) /
                 den3;
    double vyy = (xp * xp * yp * yp * exx + (1. + xp * xp) * (1. + xp * xp) * eyy -
                  2. * xp * yp * (1. + xp * xp) * exy) /
                 den3;
    double vzz = (xp * xp * exx + yp * yp * eyy + 2. * xp * yp * exy) / den3;

    double vxy = (-xp * yp * (1. + yp * yp) * exx - xp * yp * (1. + xp * xp) * eyy +
                  (1. + xp * xp + yp * yp + 2. * xp * xp * yp * yp) * exy) /
                 den3;
    double vyz =
      (xp * xp * yp * exx - yp * (1. + xp * xp) * eyy - xp * (1. + xp * xp - yp * yp) * exy) / den3;
    double vxz =
      (-xp * (1. + yp * yp) * exx + xp * yp * yp * eyy - yp * (1. - xp * xp + yp * yp) * exy) /
      den3;

    // Calculate square root of the largest eigenvalue of error matrix.

    double ddd2 = vxx * vxx + vyy * vyy + vzz * vzz - 2. * vxx * vyy - 2. * vxx * vzz -
                  2. * vyy * vzz + 4. * vxy * vxy + 4. * vyz * vyz + 4. * vxz * vxz;
    double ddd = sqrt(ddd2 > 0. ? ddd2 : 0.);
    double lambda2 = 0.5 * (vxx + vyy + vzz + ddd);
    double lambda = sqrt(lambda2 > 0. ? lambda2 : 0.);

    return lambda;
  }

  /// Get starting error matrix for Kalman filter.
  ///
  /// Arguments:
  ///
  /// err - Error matrix.
  ///
  void SurfPlane::getStartingError(TrackError& err) const
  {
    err.resize(5, false);
    err.clear();
    err(0, 0) = 1000.;
    err(1, 1) = 1000.;
    err(2, 2) = 0.25;
    err(3, 3) = 0.25;
    err(4, 4) = 10.;
  }

} // end namespace trkf
