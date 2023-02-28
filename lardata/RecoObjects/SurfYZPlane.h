////////////////////////////////////////////////////////////////////////
///
/// \file   SurfYZPlane.h
///
/// \brief  Planar surface parallel to x-axis.
///
/// \author H. Greenlee
///
/// This class represents a planar surface parallel to the global
/// x-axis, or equivalently, the normal vector is in the yz-plane.
///
/// This surface is defined by four parameters, which are,
/// (x0, y0, z0) - Local origin in yz-plane.
/// phi - Rotation angle around x-axis.
///
/// The local uvw coordinate system is related to the global xyz
/// coordinate system as follows.
///
/// u = x-x0
/// v =  (y-y0)*cos(phi) + (z-z0)*sin(phi)
/// w = -(y-y0)*sin(phi) + (z-z0)*cos(phi)
///
/// or inversely,
///
/// x = x0 + u
/// y = y0 + v*cos(phi) - w*sin(phi)
/// z = z0 + v*sin(phi) + w*cos(phi)
///
/// Track parameters on this type of surface are as follows.
///
/// 1. u
/// 2. v
/// 3. du/dw
/// 4. dv/dw
/// 5. 1/p (nonmagnetic) or q/p (magnetic)
///
////////////////////////////////////////////////////////////////////////

#ifndef SURFYZPLANE_H
#define SURFYZPLANE_H

#include "lardata/RecoObjects/SurfPlane.h"

namespace trkf {

  class SurfYZPlane : public SurfPlane {
  public:
    /// Default constructor.
    SurfYZPlane();

    /// Initializing constructor.
    SurfYZPlane(double x0, double y0, double z0, double phi);

    /// Destructor.
    virtual ~SurfYZPlane();

    // Accessors.
    double x0() const { return fX0; }   ///< X origin.
    double y0() const { return fY0; }   ///< Y origin.
    double z0() const { return fZ0; }   ///< Z origin.
    double phi() const { return fPhi; } ///< Rotation angle about x-axis.

    /// Clone method.
    virtual Surface* clone() const;

    /// Surface-specific tests of validity of track parameters.
    virtual bool isTrackValid(const TrackVector& vec) const;

    /// Transform global to local coordinates.
    virtual void toLocal(const double xyz[3], double uvw[3]) const;

    /// Transform local to global coordinates.
    virtual void toGlobal(const double uvw[3], double xyz[3]) const;

    /// Get position of track.
    virtual void getPosition(const TrackVector& vec, double xyz[3]) const;

    /// Get momentum vector of track.
    virtual void getMomentum(const TrackVector& vec,
                             double mom[3],
                             TrackDirection dir = UNKNOWN) const;

    /// Test whether two surfaces are parallel, within tolerance.
    virtual bool isParallel(const Surface& surf) const;

    /// Find perpendicular forward distance to a parallel surface
    virtual double distanceTo(const Surface& surf) const;

    /// Test two surfaces for equality, within tolerance.
    virtual bool isEqual(const Surface& surf) const;

    /// Printout
    virtual std::ostream& Print(std::ostream& out) const;

  private:
    // Static attributes.

    static double fPhiTolerance; ///< Phi tolerance for parallel.
    static double fSepTolerance; ///< Separation tolerance for equal.

    // Attributes.

    double fX0;  ///< X origin.
    double fY0;  ///< Y origin.
    double fZ0;  ///< Z origin.
    double fPhi; ///< Rotation angle about x-axis.
  };
}

#endif
