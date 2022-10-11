////////////////////////////////////////////////////////////////////////
///
/// \file   PropAny.h
///
/// \brief  Propagate between any two surfaces.
///
/// \author H. Greenlee
///
/// Class for propagating between any two surfaces.  This propagator
/// tests the type of the destination surface, and calls the
/// propagator.
///
////////////////////////////////////////////////////////////////////////

#ifndef PROPANY_H
#define PROPANY_H

#include "lardata/RecoObjects/PropXYZPlane.h"
#include "lardata/RecoObjects/PropYZLine.h"
#include "lardata/RecoObjects/PropYZPlane.h"
#include "lardata/RecoObjects/Propagator.h"

namespace trkf {

  class PropAny : public trkf::Propagator {
  public:
    /// Constructor.
    PropAny(detinfo::DetectorPropertiesData const& detProp, double tcut, bool doDedx);

    Propagator* clone() const override { return new PropAny(*this); }

    /// Propagate without error.
    std::optional<double> short_vec_prop(KTrack& trk,
                                         const std::shared_ptr<const Surface>& surf,
                                         Propagator::PropDirection dir,
                                         bool doDedx,
                                         TrackMatrix* prop_matrix = 0,
                                         TrackError* noise_matrix = 0) const override;

    /// Propagate without error to surface whose origin parameters coincide with track position.
    virtual std::optional<double> origin_vec_prop(KTrack& trk,
                                                  const std::shared_ptr<const Surface>& porient,
                                                  TrackMatrix* prop_matrix = 0) const override;

  private:
    /// Underlying propagators.

    PropYZLine fPropYZLine;
    PropYZPlane fPropYZPlane;
    PropXYZPlane fPropXYZPlane;
  };
}

#endif
