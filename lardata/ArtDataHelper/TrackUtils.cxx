/**
 * @file   lardata/ArtDataHelper/TrackUtils.cxx
 * @brief  Utility functions to extract information from `recob::Track` - implementation
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 8th, 2016
 * @see    lardata/ArtDataHelper/TrackUtils.h
 *
 */

// our header
#include "lardata/ArtDataHelper/TrackUtils.h"

// LArSoft libraries
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Vector_t
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()

// framework libraries
#include "cetlib_except/exception.h"

// ROOT libraries

// C/C++ standard libraries
#include <cmath>



//------------------------------------------------------------------------------
double lar::util::TrackProjectedLength(recob::Track const& track, geo::View_t view) {

   if(view == geo::kUnknown) {
      throw cet::exception("TrackProjectedLength") << "cannot provide projected length for "
        << "unknown view\n";
   }
   double length = 0.;

   auto const* geom = lar::providerFrom<geo::Geometry>();
   double angleToVert = 0.;
   for(unsigned int i = 0; i < geom->Nplanes(); ++i){
      if(geom->Plane(i).View() == view){
         angleToVert = geom->Plane(i).Wire(0).ThetaZ(false) - 0.5*::util::pi<>();
         break;
      }
   }

   // now loop over all points in the trajectory and add the contribution to the
   // to the desired view

   for(size_t p = 1; p < track.NumberTrajectoryPoints(); ++p){
      const auto& pos_cur = track.LocationAtPoint(p);
      const auto& pos_prev = track.LocationAtPoint(p - 1);
      double dist = std::sqrt( std::pow(pos_cur.x() - pos_prev.x(), 2) +
         std::pow(pos_cur.y() - pos_prev.y(), 2) +
         std::pow(pos_cur.z() - pos_prev.z(), 2) );

      // (sin(angleToVert),cos(angleToVert)) is the direction perpendicular to wire
      // fDir[p-1] is the direction between the two relevant points
      const auto& dir_prev = track.DirectionAtPoint(p - 1);
      double cosgamma = std::abs(std::sin(angleToVert)*dir_prev.Y() +
         std::cos(angleToVert)*dir_prev.Z() );

      /// @todo is this right, or should it be dist*cosgamma???
      length += dist/cosgamma;
   } // end loop over distances between trajectory points

   return length;
} // lar::util::TrackProjectedLength()



//------------------------------------------------------------------------------
double lar::util::TrackPitchInView
  (recob::Track const& track, geo::View_t view, size_t trajectory_point /* = 0 */)
{
   /*
    * The plan:
    * 1. find the wire plane we are talking about
    *    (in the right TPC and with the right view)
    * 2. ask the plane the answer
    *
    */


   if(trajectory_point >= track.NumberTrajectoryPoints()) {
      cet::exception("TrackPitchInView") << "ERROR: Asking for trajectory point #"
         << trajectory_point << " when trajectory vector size is of size "
         << track.NumberTrajectoryPoints() << ".\n";
   }
   recob::Track::TrajectoryPoint_t const& point
     = track.TrajectoryPoint(trajectory_point);

   // this throws if the position is not in any TPC,
   // or if there is no view with specified plane
   auto const& geom = *(lar::providerFrom<geo::Geometry>());
   geo::PlaneGeo const& plane = geom.PositionToTPC(point.position).Plane(view);

#if 0 // this can be enabled after `geo::PlaneGeo::InterWireProjectedDistance()` becomes available in larcorealg
   double const d = plane.InterWireProjectedDistance(point.direction());

   // do we prefer to just return the value and let the caller check it?
   if (d > 50.0 * plane.WirePitch()) { // after many pitches track would scatter
      throw cet::exception("Track")
        << "track at point #" << trajectory_point
        << " is almost parallel to the wires in view "
        << geo::PlaneGeo::ViewName(view) << " (wire direction is "
        << plane.GetWireDirection<geo::Vector_t>() << "; track direction is "
        << point.direction()
        << ").\n";
   }
   return d;

#else // !0
   //
   // 2. project the direction of the track on that plane
   //
   // this is the projection of the direction of the track on the wire plane;
   // it is 2D and its second component ("Y()") is on wire coordinate direction;
   // remember that the projection modulus is smaller than 1 because it is
   // the 3D direction versor, deprived of its drift direction component
   auto const& proj = plane.Projection(point.direction());

   if (lar::util::RealComparisons(1e-4).zero(proj.Y())) {
      throw cet::exception("Track")
        << "track at point #" << trajectory_point
        << " is almost parallel to the wires in view "
        << geo::PlaneGeo::ViewName(view) << " (wire direction is "
        << plane.GetWireDirection<geo::Vector_t>() << "; track direction is "
        << point.direction()
        << ", its projection on plane " << plane.ID() << " is " << proj
        << ").\n";
   }

   //
   // 3. scale that projection so that it covers a wire pitch worth in the wire
   //    coordinate direction;
   //    WirePitch() is what gives this vector a physical size [cm]
   //
   return proj.R() / std::abs(proj.Y()) * plane.WirePitch();
#endif // ?0

} // lar::util::TrackPitchInView()

//------------------------------------------------------------------------------
