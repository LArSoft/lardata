/**
 * @file   lardata/ArtDataHelper/TrackUtils.h
 * @brief  Utility functions to extract information from `recob::Track`
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 8th, 2016
 * @see    TrackUtils.cxx
 *
 *
 * `lar::util::TrackProjectedLength()` and `lar::util::TrackPitchInView()` have
 * been factored out from `recob::Track`, from `recob::Track::ProjectedLength()`
 * and `recob::Track::PitchInView()` respectively.
 */

#ifndef LARDATA_ARTDATAHELPER_TRACKUTILS_H
#define LARDATA_ARTDATAHELPER_TRACKUTILS_H

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::View_t

namespace recob {
  class Track;
}

namespace lar::util {

  /**
   * @brief Returns the length of the projection of a track on a view
   * @param track the track to be projected on a view
   * @param view the view to project the track on
   * @return length of the projection, in centimetres
   *
   *
   * @todo CAREFUL: using view to determine projected length does not work for DUNE.
   * Need to think more about this.
   *
   */
  double TrackProjectedLength(recob::Track const& track, geo::View_t view);

  /**
   * @brief Returns the projected length of track on a wire pitch step [cm]
   * @param track the track to be projected on a view
   * @param view the view for track projection
   * @param trajectory_point at which point of the track to look for the pitch (default: `0`)
   * @return wire pitch along the track direction at its specified point [cm]
   * @throw cet::exception (category `"TrackPitchInView"`) if the
   *                       `trajectory_point` index is not valid in `track`
   * @throw cet::exception (category `"Geometry"`) if the point is in no TPC
   * @throw cet::exception (category `"TPCGeo"`) if the `view` is
   *                       unknown, not available or otherwise invalid
   * @throw cet::exception (category `"Track"`) if the track projection on
   *                       the wire plane is parallel to the wires (< 0.01%)
   *
   * This function returns the distance covered by the track between two
   * wires, projected on the wire plane.
   * The direction of the track is the one at the specified trajectory point
   * (the first one by default). That direction is projected on the wire
   * plane with the specified `view` within the TPC that contains that
   * point.
   *
   * The returned value is the distance, in centimeters, between two
   * consecutive wires on that projected direction. This is always a
   * positive number, regardless the direction of the track, and never
   * smaller than the wire pitch on the projection wire plane.
   */
  double TrackPitchInView(recob::Track const& track,
                          geo::View_t view,
                          size_t trajectory_point = 0U);

} // namespace lar::util

#endif // LARDATA_ARTDATAHELPER_TRACKUTILS_H
