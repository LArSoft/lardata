////////////////////////////////////////////////////////////////////////
///
/// \file   KHitWireX.h
///
/// \brief  Kalman filter wire-time measurement on a SurfWireX surface.
///
/// \author H. Greenlee
///
/// This class is a type of one-dimensional Kalman filter measurement
/// reprsenting a single wire-time hit on a surface parallel to the
/// x-axis (approprite for nonmagnetic LAr tpc).
///
/// This class derives from base class KHit<1>, which is the general
/// one-dimensional measurement base class.  This class has a constructor
/// from an art::Ptr<recob::Hit>, which pointer is retained as a data
/// product.  This class overrides virtual base class method subpredict.
///
///
////////////////////////////////////////////////////////////////////////

#ifndef KHITWIREX_H
#define KHITWIREX_H

#include "canvas/Persistency/Common/Ptr.h"
#include "lardata/RecoObjects/KHit.h"
#include "lardataobj/RecoBase/Hit.h"

namespace detinfo {
  class DetectorPropertiesData;
}

namespace trkf {

  class KHitWireX : public KHit<1> {
  public:
    /// Constructor from Hit.
    KHitWireX(const detinfo::DetectorPropertiesData& detProp,
              const art::Ptr<recob::Hit>& hit,
              const std::shared_ptr<const Surface>& psurf);

    /// Constructor from wire id (mainly for testing).
    KHitWireX(const geo::WireID& wireid, double x, double xerr);

    // Accessors.

    /// Get original hit.
    const art::Ptr<recob::Hit>& getHit() const { return fHit; }

    // Prediction method.
    bool subpredict(const KETrack& tre,
                    KVector<1>::type& pvec,
                    KSymMatrix<1>::type& perr,
                    KHMatrix<1>::type& hmatrix) const override;

  private:
    // Attributes.

    art::Ptr<recob::Hit> fHit;
  };
}

#endif
