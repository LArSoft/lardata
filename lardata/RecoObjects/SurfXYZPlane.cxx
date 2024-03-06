///////////////////////////////////////////////////////////////////////
///
/// \file   SurfXYZPlane.cxx
///
/// \brief  General planar surface.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "TVector2.h"
#include "cetlib_except/exception.h"
#include <cmath>

namespace trkf {

  // Static attributes.

  double SurfXYZPlane::fPhiTolerance = 1.e-10;
  double SurfXYZPlane::fThetaTolerance = 1.e-10;
  double SurfXYZPlane::fSepTolerance = 1.e-6;

  /// Default constructor.
  SurfXYZPlane::SurfXYZPlane() : fX0(0.), fY0(0.), fZ0(0.), fPhi(0.), fTheta(0.) {}

  /// Initializing constructor.
  ///
  /// Arguments:
  ///
  /// x0, y0, z0 - Global coordinates of local origin.
  /// phi - Rotation angle about x-axis (wire angle).
  /// theta - Rotation angle about y'-axis (projected Lorentz angle).
  ///
  SurfXYZPlane::SurfXYZPlane(double x0, double y0, double z0, double phi, double theta)
    : fX0(x0), fY0(y0), fZ0(z0), fPhi(phi), fTheta(theta)
  {}

  /// Initializing constructor (normal vector).
  ///
  /// Arguments:
  ///
  /// x0, y0, z0 - Global coordinates of local origin.
  /// nx, ny, nz - Normal vector in global coordinate system.
  ///
  SurfXYZPlane::SurfXYZPlane(double x0, double y0, double z0, double nx, double ny, double nz)
    : fX0(x0), fY0(y0), fZ0(z0), fPhi(0.), fTheta(0.)
  {
    // Calculate angles.

    double nyz = std::hypot(ny, nz);
    fTheta = atan2(nx, nyz);
    fPhi = 0.;
    if (nyz != 0.) fPhi = atan2(-ny, nz);
  }

  /// Destructor.
  SurfXYZPlane::~SurfXYZPlane() {}

  /// Clone method.
  Surface* SurfXYZPlane::clone() const
  {
    return new SurfXYZPlane(*this);
  }

  /// Surface-specific tests of validity of track parameters.
  bool SurfXYZPlane::isTrackValid(const TrackVector& vec) const
  {
    return true;
  }

  /// Transform global to local coordinates.
  ///
  /// Arguments:
  ///
  /// xyz - Cartesian coordinates in global coordinate system.
  /// uvw - Cartesian coordinates in local coordinate system.
  ///
  void SurfXYZPlane::toLocal(const double xyz[3], double uvw[3]) const
  {
    double sinth = std::sin(fTheta);
    double costh = std::cos(fTheta);
    double sinphi = std::sin(fPhi);
    double cosphi = std::cos(fPhi);

    // u = (x-x0)*cos(theta) + (y-y0)*sin(theta)*sin(phi) - (z-z0)*sin(theta)*cos(phi)
    uvw[0] =
      (xyz[0] - fX0) * costh + (xyz[1] - fY0) * sinth * sinphi - (xyz[2] - fZ0) * sinth * cosphi;

    // v =                     (y-y0)*cos(phi)            + (z-z0)*sin(phi)
    uvw[1] = (xyz[1] - fY0) * cosphi + (xyz[2] - fZ0) * sinphi;

    // w = (x-x0)*sin(theta) - (y-y0)*cos(theta)*sin(phi) + (z-z0)*cos(theta)*cos(phi)
    uvw[2] =
      (xyz[0] - fX0) * sinth - (xyz[1] - fY0) * costh * sinphi + (xyz[2] - fZ0) * costh * cosphi;
  }

  /// Transform local to global coordinates.
  ///
  /// Arguments:
  ///
  /// uvw - Cartesian coordinates in local coordinate system.
  /// xyz - Cartesian coordinates in global coordinate system.
  ///
  void SurfXYZPlane::toGlobal(const double uvw[3], double xyz[3]) const
  {
    double sinth = std::sin(fTheta);
    double costh = std::cos(fTheta);
    double sinphi = std::sin(fPhi);
    double cosphi = std::cos(fPhi);

    // x = x0 + u*cos(theta)                       + w*sin(theta)
    xyz[0] = fX0 + uvw[0] * costh + uvw[2] * sinth;

    // y = y0 + u*sin(theta)*sin(phi) + v*cos(phi) - w*cos(theta)*sin(phi)
    xyz[1] = fY0 + uvw[0] * sinth * sinphi + uvw[1] * cosphi - uvw[2] * costh * sinphi;

    // z = z0 - u*sin(theta)*cos(phi) + v*sin(phi) + w*cos(theta)*cos(phi)
    xyz[2] = fZ0 - uvw[0] * sinth * cosphi + uvw[1] * sinphi + uvw[2] * costh * cosphi;
  }

  /// Get position of track.
  ///
  /// Arguments:
  ///
  /// vec - Track state vector.
  /// xyz - Position in global coordinate system.
  ///
  void SurfXYZPlane::getPosition(const TrackVector& vec, double xyz[3]) const
  {
    // Get position in local coordinate system.

    double uvw[3];
    uvw[0] = vec(0);
    uvw[1] = vec(1);
    uvw[2] = 0.;

    // Transform to global coordinate system.

    toGlobal(uvw, xyz);
    return;
  }

  /// Get momentum vector of track.
  ///
  /// Arguments:
  ///
  /// vec - Track state vector.
  /// mom - Momentum vector in global coordinate system.
  /// dir - Track direction.
  ///
  void SurfXYZPlane::getMomentum(const TrackVector& vec, double mom[3], TrackDirection dir) const
  {

    // Get momentum.

    double invp = std::abs(vec(4));
    double p = 1. / std::max(invp, 1.e-3); // Capped at 1000. GeV/c.

    // Get track slope parameters.

    double dudw = vec(2);
    double dvdw = vec(3);

    // Calculate dw/ds.

    double dwds = 1. / std::sqrt(1. + dudw * dudw + dvdw * dvdw);
    TrackDirection realdir = getDirection(vec, dir); // Should be same as original direction.
    if (realdir == BACKWARD)
      dwds = -dwds;
    else if (realdir != FORWARD)
      throw cet::exception("SurfXYZPlane") << "Track direction not specified.\n";

    // Calculate momentum vector in local coordinate system.

    double pu = p * dudw * dwds;
    double pv = p * dvdw * dwds;
    double pw = p * dwds;

    // Rotate momentum to global coordinte system.

    double sinth = std::sin(fTheta);
    double costh = std::cos(fTheta);
    double sinphi = std::sin(fPhi);
    double cosphi = std::cos(fPhi);

    mom[0] = pu * costh + pw * sinth;
    mom[1] = pu * sinth * sinphi + pv * cosphi - pw * costh * sinphi;
    mom[2] = -pu * sinth * cosphi + pv * sinphi + pw * costh * cosphi;

    return;
  }

  /// Test whether two surfaces are parallel, within tolerance.
  /// This method will only return true if the other surface
  /// is a SurfXYZPlane.
  ///
  /// Arguments:
  ///
  /// surf - Other surface.
  ///
  /// Returned value: true if parallel.
  ///
  bool SurfXYZPlane::isParallel(const Surface& surf) const
  {
    bool result = false;

    // Test if the other surface is a SurfXYZPlane.

    const SurfXYZPlane* psurf = dynamic_cast<const SurfXYZPlane*>(&surf);
    if (psurf != 0) {

      // Test whether surface angle parameters are the same
      // with tolerance.

      double delta_phi = TVector2::Phi_mpi_pi(fPhi - psurf->phi());
      double delta_theta = fTheta - psurf->theta();
      if (std::abs(delta_phi) <= fPhiTolerance && std::abs(delta_theta) <= fThetaTolerance)
        result = true;
    }
    return result;
  }

  /// Find perpendicular forward distance to a parallel surface.
  ///
  /// Throw an exception if the other surface is not parallel.
  ///
  /// Assuming the other surface is parallel, the distance is
  /// simply the w-coordinate of the other surface, and is signed.
  ///
  /// Arguments:
  ///
  /// surf - Other surface.
  ///
  /// Returned value: Distance.
  ///
  double SurfXYZPlane::distanceTo(const Surface& surf) const
  {
    // Check if the other surface is parallel to this one.

    bool parallel = isParallel(surf);
    if (!parallel)
      throw cet::exception("SurfXYZPlane") << "Attempt to find distance to non-parallel surface.\n";

    // Find the origin of the other surface in global coordinates,
    // then convert to our local coordinates.

    double otheruvw[3] = {0., 0., 0.};
    double xyz[3];
    double myuvw[3];
    surf.toGlobal(otheruvw, xyz);
    toLocal(xyz, myuvw);

    // Distance is local w-coordinate of other surface origin.

    return myuvw[2];
  }

  /// Test two surfaces for equality, within tolerance.
  /// Here equal is defined as parallel and having zero separation,
  /// within tolerance.  Note that this definition of equality allows
  /// the two surfaces to have different origins.
  ///
  /// Arguments:
  ///
  /// surf - Other surface.
  ///
  /// Returned values: true if equal.
  ///
  bool SurfXYZPlane::isEqual(const Surface& surf) const
  {
    bool result = false;

    // Test if the other surface is parallel.

    bool parallel = isParallel(surf);
    if (parallel) {
      double dist = distanceTo(surf);
      if (std::abs(dist) <= fSepTolerance) result = true;
    }

    return result;
  }

  /// Printout
  std::ostream& SurfXYZPlane::Print(std::ostream& out) const
  {
    out << "SurfXYZPlane{ x0=" << fX0 << ", y0=" << fY0 << ", z0=" << fZ0 << ", phi=" << fPhi
        << ", theta=" << fTheta << "}";
    return out;
  }

} // end namespace trkf
