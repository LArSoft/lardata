////////////////////////////////////////////////////////////////////////
///
/// \file   InteractGeneral.h
///
/// \brief  Interactor for planar surfaces.
///
/// \author H. Greenlee
///
/// This class calculates propagation noise for tracks on any surface.
/// This class works by transforming tracks to a planar surface that is
/// normal to the track and calculating the noise matrix on that surface,
/// then transforming the noise matrix back to the original surface.
///
////////////////////////////////////////////////////////////////////////

#ifndef INTERACTGENERAL_H
#define INTERACTGENERAL_H

#include "lardata/RecoObjects/InteractPlane.h"
#include "lardata/RecoObjects/PropAny.h"

namespace trkf {

  class InteractGeneral : public trkf::Interactor {
  public:
    explicit InteractGeneral(detinfo::DetectorPropertiesData const& detProp, double tcut);

    Interactor* clone() const override { return new InteractGeneral(*this); }
    bool noise(const KTrack& trk, double s, TrackError& noise_matrix) const override;

  private:
    InteractPlane fInteract;
    PropAny fProp;
  };
}

#endif
