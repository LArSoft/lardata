///////////////////////////////////////////////////////////////////////
///
/// \file   SurfYZLine.cxx
///
/// \brief  Line surface perpendicular to x-axis.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/SurfYZLine.h"
#include "TMath.h"
#include "TVector2.h"
#include "cetlib_except/exception.h"
#include <cmath>

namespace trkf {

  // Static attributes.

  double SurfYZLine::fPhiTolerance = 1.e-10;
  double SurfYZLine::fSepTolerance = 1.e-6;

  /// Default constructor.
  SurfYZLine::SurfYZLine() : fX0(0.), fY0(0.), fZ0(0.), fPhi(0.) {}

  /// Initializing constructor.
  ///
  /// Arguments:
  ///
  /// x0, y0, z0 - Global coordinates of local origin.
  /// phi - Rotation angle about x-axis.
  ///
  SurfYZLine::SurfYZLine(double x0, double y0, double z0, double phi)
    : fX0(x0), fY0(y0), fZ0(z0), fPhi(phi)
  {}

  /// Destructor.
  SurfYZLine::~SurfYZLine() {}

  /// Clone method.
  Surface* SurfYZLine::clone() const
  {
    return new SurfYZLine(*this);
  }

  /// Surface-specific tests of validity of track parameters.
  bool SurfYZLine::isTrackValid(const TrackVector& vec) const
  {
    // Limit allowed range of eta parameter.

    return std::abs(vec(3)) < 10.;
  }

  /// Transform global to local coordinates.
  ///
  /// Arguments:
  ///
  /// xyz - Cartesian coordinates in global coordinate system.
  /// uvw - Cartesian coordinates in local coordinate system.
  ///
  void SurfYZLine::toLocal(const double xyz[3], double uvw[3]) const
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
  void SurfYZLine::toGlobal(const double uvw[3], double xyz[3]) const
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

  /// Calculate difference of two track parameter vectors, taking into account phi wrap.
  ///
  /// Arguments:
  ///
  /// vec1 - First vector.
  /// vec2 - Second vector.
  ///
  /// Returns: vec1 - vec2
  ///
  TrackVector SurfYZLine::getDiff(const TrackVector& vec1, const TrackVector& vec2) const
  {
    TrackVector result = vec1 - vec2;
    while (result(2) <= -TMath::Pi())
      result(2) += TMath::TwoPi();
    while (result(2) > TMath::Pi())
      result(2) -= TMath::TwoPi();
    return result;
  }

  /// Get position of track.
  ///
  /// Arguments:
  ///
  /// vec - Track state vector.
  /// xyz - Position in global coordinate system.
  ///
  void SurfYZLine::getPosition(const TrackVector& vec, double xyz[3]) const
  {
    // Get position in local coordinate system.

    double phi = vec(2);
    double sinphi = std::sin(phi);
    double cosphi = std::cos(phi);

    double uvw[3];
    uvw[0] = -vec(0) * sinphi;
    uvw[1] = vec(1);
    uvw[2] = vec(0) * cosphi;

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
  /// dir - Track direction (ignored).
  ///
  void SurfYZLine::getMomentum(const TrackVector& vec, double mom[3], TrackDirection dir) const
  {

    // Get momentum.

    double invp = std::abs(vec(4));
    double p = 1. / std::max(invp, 1.e-3); // Capped at 1000. GeV/c.

    // Get track direction parameters.

    double phi = vec(2);
    double eta = vec(3);

    double sinphi = std::sin(phi);
    double cosphi = std::cos(phi);
    double sh = 1. / std::cosh(eta); // sech(eta)
    double th = std::tanh(eta);

    // Calculate momentum vector in local coordinate system.

    double pu = p * cosphi * sh;
    double pv = p * th;
    double pw = p * sinphi * sh;

    // Rotate momentum to global coordinte system.

    double sinfphi = std::sin(fPhi);
    double cosfphi = std::cos(fPhi);

    mom[0] = pu;
    mom[1] = pv * cosfphi - pw * sinfphi;
    mom[2] = pv * sinfphi + pw * cosfphi;

    return;
  }

  /// Test whether two surfaces are parallel, within tolerance.
  /// This method will only return true if the other surface
  /// is a SurfYZLine.
  ///
  /// Arguments:
  ///
  /// surf - Other surface.
  ///
  /// Returned value: true if parallel.
  ///
  bool SurfYZLine::isParallel(const Surface& surf) const
  {
    bool result = false;

    // Test if the other surface is a SurfYZLine.

    const SurfYZLine* psurf = dynamic_cast<const SurfYZLine*>(&surf);
    if (psurf != 0) {

      // Test whether surface angle parameters are the same
      // within tolerance.

      double delta_phi = TVector2::Phi_mpi_pi(fPhi - psurf->phi());
      if (std::abs(delta_phi) <= fPhiTolerance) result = true;
    }
    return result;
  }

  /// Find perpendicular distance to a parallel surface.
  ///
  /// Throw an exception if the other surface is not parallel.
  ///
  /// Arguments:
  ///
  /// surf - Other surface.
  ///
  /// Returned value: Distance.
  ///
  double SurfYZLine::distanceTo(const Surface& surf) const
  {
    // Check if the other surface is parallel to this one.

    bool parallel = isParallel(surf);
    if (!parallel)
      throw cet::exception("SurfYZLine") << "Attempt to find distance to non-parallel surface.\n";

    // Find the origin of the other surface in global coordinates,
    // then convert to our local coordinates.

    double otheruvw[3] = {0., 0., 0.};
    double xyz[3];
    double myuvw[3];
    surf.toGlobal(otheruvw, xyz);
    toLocal(xyz, myuvw);

    // Distance of v-axis to other surface origin.

    return std::hypot(myuvw[0], myuvw[2]);
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
  bool SurfYZLine::isEqual(const Surface& surf) const
  {
    bool result = false;

    // Test if the other surface is a SurfYZLine.

    const SurfYZLine* psurf = dynamic_cast<const SurfYZLine*>(&surf);
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
  std::ostream& SurfYZLine::Print(std::ostream& out) const
  {
    out << "SurfYZLine{ x0=" << fX0 << ", y0=" << fY0 << ", z0=" << fZ0 << ", phi=" << fPhi << "}";
    return out;
  }

} // end namespace trkf
