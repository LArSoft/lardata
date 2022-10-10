#ifndef TRACKINGPLANEHELPER_H
#define TRACKINGPLANEHELPER_H

#include "TMath.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardata/RecoObjects/SurfWireX.h"
#include "lardataobj/RecoBase/TrackingPlane.h"
#include "lardataobj/RecoBase/Trajectory.h"

namespace recob {
  namespace tracking {

    /// \file  lardata/RecoObjects/TrackingPlaneHelper.h
    ///
    /// \author  G. Cerati (FNAL, MicroBooNE)
    /// \date    2017
    /// \version 1.0

    /// helper function to construct a recob::tracking::Plane from a Point_t and a Vector_t; the point is on the plane, the vector is orthogonal to the plane.
    inline Plane makePlane(recob::tracking::Point_t const& pos,
                           recob::tracking::Vector_t const& dir)
    {
      return Plane(pos, dir);
    }

    /// helper function to construct a recob::tracking::Plane from a recob::Trajectory::TrajectoryPoint_t.
    inline Plane makePlane(recob::Trajectory::TrajectoryPoint_t const& s)
    {
      return Plane(s.position, s.direction());
    }

    /// helper function to construct a recob::tracking::Plane from a trkf::SurfWireX object.
    inline Plane makePlane(trkf::SurfWireX const& s)
    {
      return Plane(Point_t(s.x0(), s.y0(), s.z0()),
                   Vector_t(0, -std::sin(s.phi()), std::cos(s.phi())));
    }

    /// helper function to construct a recob::tracking::Plane from a geo::WireGeo object. The plane will contain the wire and the x axis, assumed to be the drift direction (to be generalized).
    inline Plane makePlane(geo::WireGeo const& wgeom)
    {
      auto xyz = wgeom.GetCenter<Point_t>();
      xyz.SetX(0.);
      double const phi = TMath::PiOver2() - wgeom.ThetaZ();
      return Plane{xyz, Vector_t{0, -std::sin(phi), std::cos(phi)}};
    }

  }
}

#endif
