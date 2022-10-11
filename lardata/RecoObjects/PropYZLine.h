////////////////////////////////////////////////////////////////////////
///
/// \file   PropYZLine.h
///
/// \brief  Propagate to SurfYZLine surface.
///
/// \author H. Greenlee
///
/// Class for propagating to a destionation SurfYZLine surface.
///
////////////////////////////////////////////////////////////////////////

#ifndef PROPYZLINE_H
#define PROPYZLINE_H

#include "lardata/RecoObjects/Propagator.h"

namespace detinfo {
  class DetectorPropertiesData;
}

namespace trkf {

  class PropYZLine : public trkf::Propagator {
  public:
    PropYZLine(detinfo::DetectorPropertiesData const& detProp, double tcut, bool doDedx);

    Propagator* clone() const override { return new PropYZLine(*this); }

    /// Propagate without error.
    std::optional<double> short_vec_prop(KTrack& trk,
                                         const std::shared_ptr<const Surface>& surf,
                                         Propagator::PropDirection dir,
                                         bool doDedx,
                                         TrackMatrix* prop_matrix = 0,
                                         TrackError* noise_matrix = 0) const override;

    /// Propagate without error to surface whose origin parameters coincide with track position.
    std::optional<double> origin_vec_prop(KTrack& trk,
                                          const std::shared_ptr<const Surface>& porient,
                                          TrackMatrix* prop_matrix = 0) const override;

  private:
    /// The following methods transform the track parameters from
    /// initial surface to SurfYZLine origin surface, and generate a
    /// propagation matrix.  The first group of function parameters
    /// are the orientation surface parameters of the initial surface.
    /// The second group of function parameters are the orientation
    /// parameters of the of the destination surface.  The origin
    /// parameters of the destination surface are assumed to match the
    /// position of the track.

    /// Transform yz line -> yz line.

    bool transformYZLine(double phi1,
                         double phi2,
                         TrackVector& vec,
                         Surface::TrackDirection& dir,
                         TrackMatrix* prop_matrix) const;

    /// Transform yz plane -> yz line.

    bool transformYZPlane(double phi1,
                          double phi2,
                          TrackVector& vec,
                          Surface::TrackDirection& dir,
                          TrackMatrix* prop_matrix) const;
    /// Transform xyz plane -> yz line.

    bool transformXYZPlane(double theta1,
                           double phi1,
                           double phi2,
                           TrackVector& vec,
                           Surface::TrackDirection& dir,
                           TrackMatrix* prop_matrix) const;
  };
}

#endif
