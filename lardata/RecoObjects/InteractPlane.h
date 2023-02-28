////////////////////////////////////////////////////////////////////////
///
/// \file   InteractPlane.h
///
/// \brief  Interactor for planar surfaces.
///
/// \author H. Greenlee
///
/// This class calculates propagation noise for tracks on planar
/// surfaces.  In particular, this class works for any surface that
/// has a local Cartesian coordinate system in which the track
/// parameters are (u, v, u'=du/dw, v'=dv/dw, q/p).
///
////////////////////////////////////////////////////////////////////////

#ifndef INTERACTPLANE_H
#define INTERACTPLANE_H

#include "lardata/RecoObjects/Interactor.h"

namespace detinfo {
  class DetectorPropertiesData;
}

namespace trkf {

  class InteractPlane : public trkf::Interactor {
  public:
    InteractPlane(detinfo::DetectorPropertiesData const& detProp, double tcut);

    Interactor* clone() const override { return new InteractPlane(*this); }

    /// Calculate noise matrix.
    bool noise(const KTrack& trk, double s, TrackError& noise_matrix) const override;

  private:
    detinfo::DetectorPropertiesData const& fDetProp;
  };
}

#endif
