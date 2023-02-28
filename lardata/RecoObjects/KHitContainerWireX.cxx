///////////////////////////////////////////////////////////////////////
///
/// \file   KHitContainerWireX.cxx
///
/// \brief  A KHitContainer for KHitWireX type measurements.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <map>

#include "cetlib_except/exception.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoObjects/KHitContainerWireX.h"
#include "lardata/RecoObjects/KHitWireX.h"

namespace trkf {

  /// Fill container.
  ///
  /// Arguments:
  ///
  /// hits       - RecoBase/Hit collection.
  /// only_plane - Choose hits from this plane if >= 0.
  ///
  /// This method converts the hits in the input collection into
  /// KHitWireX objects and inserts them into the base class.  Hits
  /// corresponding to the same readout wire are grouped together as
  /// KHitGroup objects.
  ///
  void KHitContainerWireX::fill(const detinfo::DetectorPropertiesData& detProp,
                                const art::PtrVector<recob::Hit>& hits,
                                int only_plane)
  {
    // Get services.

    art::ServiceHandle<geo::Geometry const> geom;

    // Make a temporary map from channel number to KHitGroup objects.
    // The KHitGroup pointers are borrowed references to KHitGroup
    // objects stored by value in the base class.

    std::map<unsigned int, KHitGroup*> group_map;

    // Loop over hits.

    for (art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin(); ihit != hits.end();
         ++ihit) {
      const recob::Hit& hit = **ihit;

      // Extract the wire id from the Hit.
      geo::WireID hitWireID = hit.WireID();

      uint32_t channel = hit.Channel();

      // Choose plane.
      if (only_plane >= 0 && hitWireID.Plane != (unsigned int)(only_plane)) continue;

      // See if we need to make a new KHitGroup.

      KHitGroup* pgr = 0;
      if (group_map.count(channel) == 0) {
        getUnsorted().push_back(KHitGroup());
        pgr = &(getUnsorted().back());
        group_map[channel] = pgr;
      }
      else
        pgr = group_map[channel];
      if (!pgr) {
        throw cet::exception("KHitContainerWireX")
          << __func__ << ": no group map for channel " << channel << "\n";
      }

      pgr->addHit(std::make_shared<KHitWireX>(detProp, *ihit, pgr->getSurface()));
    }
  }

} // end namespace trkf
