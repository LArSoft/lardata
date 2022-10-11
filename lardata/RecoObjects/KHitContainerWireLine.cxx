///////////////////////////////////////////////////////////////////////
///
/// \file   KHitContainerWireLine.cxx
///
/// \brief  A KHitContainer for KHitWireLine type measurements.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KHitContainerWireLine.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoObjects/KHitWireLine.h"

namespace trkf {

  /// Fill container.
  ///
  /// Arguments:
  ///
  /// hits       - RecoBase/Hit collection.
  /// only_plane - Choose hits from this plane if >= 0.
  ///
  /// This method converts the hits in the input collection into
  /// KHitWireLine objects and inserts them into the base class.
  ///
  void KHitContainerWireLine::fill(const detinfo::DetectorPropertiesData& detProp,
                                   const art::PtrVector<recob::Hit>& hits,
                                   int only_plane)
  {
    // Get services.

    art::ServiceHandle<geo::Geometry const> geom;

    // Loop over hits.

    for (art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end();
         ++ihit) {
      const recob::Hit& hit = **ihit;

      // Extract the wireid from the Hit.
      geo::WireID hitWireID = hit.WireID();

      uint32_t channel = hit.Channel();

      // Choose plane.
      if (only_plane >= 0 && hitWireID.Plane != (unsigned int)(only_plane)) continue;

      // Make a new KHitGroup for each hit.

      getUnsorted().push_back(KHitGroup());
      KHitGroup* pgr = &(getUnsorted().back());
      if (!pgr) {
        throw cet::exception("KHitContainerWireLine")
          << __func__ << ": no group map for channel " << channel << "\n";
      }

      pgr->addHit(std::make_shared<KHitWireLine>(detProp, *ihit, pgr->getSurface()));
    }
  }

} // end namespace trkf
