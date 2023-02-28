///////////////////////////////////////////////////////////////////////
///
/// \file   SurfLine.cxx
///
/// \brief  Base class for Kalman filter planar surfaces.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/SurfLine.h"
#include <cmath>

namespace trkf {

  /// Default constructor.
  SurfLine::SurfLine() {}

  /// Destructor.
  SurfLine::~SurfLine() {}

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
  double SurfLine::PointingError(const TrackVector& vec, const TrackError& err) const
  {
    // Get slope parameters and error matrix.

    double phi = vec(2);
    double eta = vec(3);
    double epp = err(2, 2); // sigma^2(phi,phi)
    double ehh = err(3, 3); // sigma^2(eta,eta)
    double ehp = err(3, 2); // sigma^2(eta,phi)

    // Calculate error matrix of pointing unit vector in some coordinate system.

    double sh = 1. / std::cosh(eta); // sech(eta)
    double sh2 = sh * sh;
    double sh3 = sh * sh2;
    double sh4 = sh * sh3;

    double th = std::tanh(eta);
    double th2 = th * th;

    double cphi = std::cos(phi);
    double cphi2 = cphi * cphi;

    double sphi = std::sin(phi);
    double sphi2 = sphi * sphi;

    double vxx = sh2 * th2 * cphi2 * ehh + sh2 * sphi2 * epp + 2. * sh2 * th * sphi * cphi * ehp;
    double vyy = sh2 * th2 * sphi2 * ehh + sh2 * cphi2 * epp - 2. * sh2 * th * sphi * cphi * ehp;
    double vzz = sh4 * epp;

    double vxy =
      sh2 * th2 * sphi * cphi * ehh - sh2 * sphi * cphi * epp - sh2 * th * (cphi2 - sphi2) * ehp;
    double vyz = -sh3 * th * sphi * ehh + sh3 * cphi * ehp;
    double vxz = -sh3 * th * cphi * ehh - sh3 * sphi * ehp;

    // For debugging.  The determinant of the error matrix should be zero.

    // double det = vxx*vyy*vzz + 2.*vxy*vyz*vxz - vxx*vyz*vyz - vyy*vxz*vxz - vzz*vxy*vxy;

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
  void SurfLine::getStartingError(TrackError& err) const
  {
    err.resize(5, false);
    err.clear();
    err(0, 0) = 1000.;
    err(1, 1) = 1000.;
    err(2, 2) = 10.;
    err(3, 3) = 10.;
    err(4, 4) = 10.;
  }

} // end namespace trkf
