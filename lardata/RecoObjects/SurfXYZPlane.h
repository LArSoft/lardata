////////////////////////////////////////////////////////////////////////
///
/// \file   SurfXYZPlane.h
///
/// \brief  General planar surface.
///
/// \author H. Greenlee
///
/// This class represents an (almost) general planar surface.  It is
/// intended to represent the measurement surface defined by a wire
/// and the drift velocity in the case of a magnetic LAr TPC with a
/// nonzero Lorentz angle.
///
/// In contrast to a completely general planar surface, the rotation
/// part of the global to local coordinate transformation for this
/// surface is defined by two rotation angles (not three Euler
/// angles), which are the wire angle and the projected Lorentz angle.
/// What we calling the projected Lorentz angle is the dihedral angle
/// between the plane defined by the wire and the x-axis, and the
/// plane defined by the wire and the drift velocity, which can be
/// different for different views.
///
/// In the case of a nonzero Lorentz angle, the drift velocity is not
/// perpendicular to a readout wire, so the drift velocity vector will
/// not coincide with a coordinate axis in the local coordinate
/// system.  Nevertheless, the projected Lorentz angle provides
/// sufficient information to implement the prediction function for
/// this type of surface, provided the drift velocity perpendicular to
/// the readout planes (which is the same for every view) is available
/// externally.
///
/// This surface is defined by five parameters, which are,
/// (x0, y0, z0) - Local origin.
/// phi - Rotation angle around x-axis (wire angle).
/// theta - Rotation angle around y'-axis (projected Lorentz angle).
///
/// The local uvw coordinate system is related to the global xyz
/// coordinate system via an intermediate x'y'z' system as follows.
///
/// x' = x-x0
/// y' =  (y-y0)*cos(phi) + (z-z0)*sin(phi)
/// z' = -(y-y0)*sin(phi) + (z-z0)*cos(phi)
///
/// u = x'*cos(theta) - z'*sin(theta)
/// v = y'
/// w = x'*sin(theta) + z'*cos(theta)
///
/// u = (x-x0)*cos(theta) + (y-y0)*sin(theta)*sin(phi) - (z-z0)*sin(theta)*cos(phi)
/// v =                     (y-y0)*cos(phi)            + (z-z0)*sin(phi)
/// w = (x-x0)*sin(theta) - (y-y0)*cos(theta)*sin(phi) + (z-z0)*cos(theta)*cos(phi)
///
/// or inversely,
///
/// x' = u*cos(theta) + w*sin(theta)
/// y' = v
/// z' = -u*sin(theta) + w*cos(theta)
///
/// x = x0 + x'
/// y = y0 + y'*cos(phi) - z'*sin(phi)
/// z = z0 + y'*sin(phi) + z'*cos(phi)
///
/// x = x0 + u*cos(theta)                       + w*sin(theta)
/// y = y0 + u*sin(theta)*sin(phi) + v*cos(phi) - w*cos(theta)*sin(phi)
/// z = z0 - u*sin(theta)*cos(phi) + v*sin(phi) + w*cos(theta)*cos(phi)
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

#ifndef SURFXYZPLANE_H
#define SURFXYZPLANE_H

#include "lardata/RecoObjects/SurfPlane.h"

namespace trkf {

  class SurfXYZPlane : public SurfPlane {
  public:
    /// Default constructor.
    SurfXYZPlane();

    /// Initializing constructor (angles).
    SurfXYZPlane(double x0, double y0, double z0, double phi, double theta);

    /// Initializing constructor (normal vector).
    SurfXYZPlane(double x0, double y0, double z0, double nx, double ny, double nz);

    /// Destructor.
    virtual ~SurfXYZPlane();

    // Accessors.
    double x0() const { return fX0; }       ///< X origin.
    double y0() const { return fY0; }       ///< Y origin.
    double z0() const { return fZ0; }       ///< Z origin.
    double phi() const { return fPhi; }     ///< Rot. angle about x-axis (wire angle).
    double theta() const { return fTheta; } ///< Rot. angle about y'-axis (projected Lorentz angle).

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

    static double fPhiTolerance;   ///< Phi tolerance for parallel.
    static double fThetaTolerance; ///< Theta tolerance for parallel.
    static double fSepTolerance;   ///< Separation tolerance for equal.

    // Attributes.

    double fX0;    ///< X origin.
    double fY0;    ///< Y origin.
    double fZ0;    ///< Z origin.
    double fPhi;   ///< Rotation angle about x-axis (wire angle).
    double fTheta; ///< Rotation angle about y'-axis (projected Lorentz angle).
  };
}

#endif
