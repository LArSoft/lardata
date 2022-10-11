////////////////////////////////////////////////////////////////////////
///
/// \file   SurfYZLine.h
///
/// \brief  Line surface perpendicular to x-axis.
///
/// \author H. Greenlee
///
/// This class represents a line surface perpendicular to the global
/// x-axis, or equivalently, parallel to the yz-plane.
///
/// The surface parameters and local coordinate system of this
/// surface are the same as SurfYZPlane.
///
/// This surface is defined by four parameters, which are,
/// (x0, y0, z0) - Local origin in global coordinates.
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
/// Track parameters on this type of surface are:
///
/// 1. r
/// 2. v
/// 3. phi
/// 4. eta
/// 5. 1/p (nonmagnetic) or q/p (magnetic)
///
/// r = Signed impoact parameter.  Absolute value of r is the perpendicular
///     distance of the track to the v-axis at the point of closest
///     approach to v-axis.  Sign of r matches sign of L_v (v projection
///     of angular momentum).
/// v = V-coordinate of track at point of closest approach to v-axis.
/// phi = Direction of track in u-w plane (phi = arctan(w/u)).
/// eta = Pseudorapidity with respect to v-axis.
/// q/p or 1/p = Inverse momentum.
///
/// In terms of these parameters, the point of closest approach to the
/// v-axis is
///
/// u = -r sin(phi)
/// v = v
/// w = r cos(phi)
///
/// The unit direction vector is
///
/// du/ds = cos(phi) sech(eta)
/// dv/ds = tanh(eta)
/// dw/ds = sin(phi) sech(eta)
///
/// Inversely:
///
/// phi = atan(dw/du) = atan2(dw/ds, du/ds)
/// eta = atanh(dv/ds) = asinh(dv/duw)
/// r = w cos(phi) - u sin(phi)
///
///
////////////////////////////////////////////////////////////////////////

#ifndef SURFYZLINE_H
#define SURFYZLINE_H

#include "lardata/RecoObjects/SurfLine.h"

namespace trkf {

  class SurfYZLine : public SurfLine {
  public:
    /// Default constructor.
    SurfYZLine();

    /// Initializing constructor.
    SurfYZLine(double x0, double y0, double z0, double phi);

    /// Destructor.
    virtual ~SurfYZLine();

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

    /// Calculate difference of two track parameter vectors.
    virtual TrackVector getDiff(const TrackVector& vec1, const TrackVector& vec2) const;

    /// Get position of track.
    virtual void getPosition(const TrackVector& vec, double xyz[3]) const;

    /// Get momentum vector of track.
    virtual void getMomentum(const TrackVector& vec,
                             double mom[3],
                             TrackDirection dir = UNKNOWN) const;

    /// Test whether two surfaces are parallel, within tolerance.
    virtual bool isParallel(const Surface& surf) const;

    /// Find perpendicular distance to a parallel surface
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
