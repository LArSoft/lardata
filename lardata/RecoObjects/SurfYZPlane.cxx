///////////////////////////////////////////////////////////////////////
///
/// \file   SurfYZPlane.cxx
///
/// \brief  Planar surface parallel to x-axis.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/SurfYZPlane.h"
#include "TVector2.h"
#include "cetlib_except/exception.h"
#include <cmath>

namespace trkf {

  // Static attributes.

  double SurfYZPlane::fPhiTolerance = 1.e-10;
  double SurfYZPlane::fSepTolerance = 1.e-6;

  /// Default constructor.
  SurfYZPlane::SurfYZPlane() : fX0(0.), fY0(0.), fZ0(0.), fPhi(0.) {}

  /// Initializing constructor.
  ///
  /// Arguments:
  ///
  /// x0, y0, z0 - Global coordinates of local origin.
  /// phi - Rotation angle about x-axis.
  ///
  SurfYZPlane::SurfYZPlane(double x0, double y0, double z0, double phi)
    : fX0(x0), fY0(y0), fZ0(z0), fPhi(phi)
  {}

  /// Destructor.
  SurfYZPlane::~SurfYZPlane() {}

  /// Clone method.
  Surface* SurfYZPlane::clone() const { return new SurfYZPlane(*this); }

  /// Surface-specific tests of validity of track parameters.
  bool SurfYZPlane::isTrackValid(const TrackVector& vec) const { return true; }

  /// Transform global to local coordinates.
  ///
  /// Arguments:
  ///
  /// xyz - Cartesian coordinates in global coordinate system.
  /// uvw - Cartesian coordinates in local coordinate system.
  ///
  void SurfYZPlane::toLocal(const double xyz[3], double uvw[3]) const
  {
    double sinphi = std::sin(fPhi);
    double cosphi = std::cos(fPhi);

    // u = x-x0
    uvw[0] = xyz[0] - fX0;

    // v =  (y-y0)*cos(phi) + (z-z0)*sin(phi)
    uvw[1] = (xyz[1] - fY0) * cosphi + (xyz[2] - fZ0) * sinphi;

    // w = -(y-y0)*sin(phi) + (z-z0)*cos(phi)
    uvw[2] = -(xyz[1] - fY0) * sinphi + (xyz[2] - fZ0) * cosphi;
  }

  /// Transform local to global coordinates.
  ///
  /// Arguments:
  ///
  /// uvw - Cartesian coordinates in local coordinate system.
  /// xyz - Cartesian coordinates in global coordinate system.
  ///
  void SurfYZPlane::toGlobal(const double uvw[3], double xyz[3]) const
  {
    double sinphi = std::sin(fPhi);
    double cosphi = std::cos(fPhi);

    // x = x0 + u
    xyz[0] = fX0 + uvw[0];

    // y = y0 + v*cos(phi) - w*sin(phi)
    xyz[1] = fY0 + uvw[1] * cosphi - uvw[2] * sinphi;

    // z = z0 + v*sin(phi) + w*cos(phi)
    xyz[2] = fZ0 + uvw[1] * sinphi + uvw[2] * cosphi;
  }

  /// Get position of track.
  ///
  /// Arguments:
  ///
  /// vec - Track state vector.
  /// xyz - Position in global coordinate system.
  ///
  void SurfYZPlane::getPosition(const TrackVector& vec, double xyz[3]) const
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
  void SurfYZPlane::getMomentum(const TrackVector& vec, double mom[3], TrackDirection dir) const
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
      throw cet::exception("SurfYZPlane") << "Track direction not specified.\n";

    // Calculate momentum vector in local coordinate system.

    double pu = p * dudw * dwds;
    double pv = p * dvdw * dwds;
    double pw = p * dwds;

    // Rotate momentum to global coordinte system.

    double sinphi = std::sin(fPhi);
    double cosphi = std::cos(fPhi);

    mom[0] = pu;
    mom[1] = pv * cosphi - pw * sinphi;
    mom[2] = pv * sinphi + pw * cosphi;

    return;
  }

  /// Test whether two surfaces are parallel, within tolerance.
  /// This method will only return true if the other surface
  /// is a SurfYZPlane.
  ///
  /// Arguments:
  ///
  /// surf - Other surface.
  ///
  /// Returned value: true if parallel.
  ///
  bool SurfYZPlane::isParallel(const Surface& surf) const
  {
    bool result = false;

    // Test if the other surface is a SurfYZPlane.

    const SurfYZPlane* psurf = dynamic_cast<const SurfYZPlane*>(&surf);
    if (psurf != 0) {

      // Test whether surface angle parameters are the same
      // within tolerance.

      double delta_phi = TVector2::Phi_mpi_pi(fPhi - psurf->phi());
      if (std::abs(delta_phi) <= fPhiTolerance) result = true;
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
  double SurfYZPlane::distanceTo(const Surface& surf) const
  {
    // Check if the other surface is parallel to this one.

    bool parallel = isParallel(surf);
    if (!parallel)
      throw cet::exception("SurfYZPlane") << "Attempt to find distance to non-parallel surface.\n";

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
  /// Here equal is defined as having all surface parameters the same,
  /// not just having the surfaces coincide spatially, so that the
  /// local coordinate systems are the same between the two surfaces.
  ///
  /// Arguments:
  ///
  /// surf - Other surface.
  ///
  /// Returned values: true if equal.
  ///
  bool SurfYZPlane::isEqual(const Surface& surf) const
  {
    bool result = false;

    // Test if the other surface is a SurfYZPlane.

    const SurfYZPlane* psurf = dynamic_cast<const SurfYZPlane*>(&surf);
    if (psurf != 0) {

      // Test whether surface parameters are the same within tolerance.

      double delta_phi = TVector2::Phi_mpi_pi(fPhi - psurf->phi());
      double dx = fX0 - psurf->x0();
      double dy = fY0 - psurf->y0();
      double dz = fZ0 - psurf->z0();
      if (std::abs(delta_phi) <= fPhiTolerance && std::abs(dx) <= fSepTolerance &&
          std::abs(dy) <= fSepTolerance && std::abs(dz) <= fSepTolerance)
        result = true;
    }
    return result;
  }

  /// Printout
  std::ostream& SurfYZPlane::Print(std::ostream& out) const
  {
    out << "SurfYZPlane{ x0=" << fX0 << ", y0=" << fY0 << ", z0=" << fZ0 << ", phi=" << fPhi << "}";
    return out;
  }

} // end namespace trkf
