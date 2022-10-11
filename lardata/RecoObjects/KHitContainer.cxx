///////////////////////////////////////////////////////////////////////
///
/// \file   KHitContainer.cxx
///
/// \brief  A collection of KHitGroups.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KHitContainer.h"

#include "cetlib_except/exception.h"

namespace trkf {

  void fill(const art::PtrVector<recob::Hit>& hits, int only_plane) {}

  /// Clear all lists.
  void KHitContainer::clear()
  {
    fSorted.clear();
    fUnsorted.clear();
    fUnused.clear();
  }

  /// Move all objects to unsorted list (from sorted and unused lists).
  void KHitContainer::reset()
  {
    fUnsorted.splice(fUnsorted.end(), fSorted);
    fUnsorted.splice(fUnsorted.end(), fUnused);
  }

  /// (Re)sort objects in unsorted and sorted lists.
  ///
  /// Arguments:
  ///
  /// trk         - Track to be propagated.
  /// addUnsorted - If true, include unsorted objects in sort.
  /// prop        - Propagator.
  /// dir         - Propagation direction.
  ///
  void KHitContainer::sort(const KTrack& trk,
                           bool addUnsorted,
                           const Propagator& prop,
                           Propagator::PropDirection dir)
  {
    // Maybe transfer all objects in unsorted list to the sorted list.

    if (addUnsorted) fSorted.splice(fSorted.end(), fUnsorted);

    // Loop over objects in sorted list.

    for (std::list<KHitGroup>::iterator igr = fSorted.begin(); igr != fSorted.end();) {

      KHitGroup& gr = *igr;

      // Get destination surface.

      const std::shared_ptr<const Surface>& psurf = gr.getSurface();

      // Make a fresh copy of the track and propagate it to
      // the destination surface.

      KTrack trkp = trk;
      std::optional<double> dist = prop.vec_prop(trkp, psurf, dir, false, 0, 0);
      if (!dist) {

        // If propagation failed, reset the path flag for this surface
        // and move the KHitGroup to the unsorted list.  Be careful to
        // keep the list iterator valid.

        gr.setPath(false, 0.);
        std::list<KHitGroup>::iterator it = igr;
        ++igr;
        fUnsorted.splice(fUnsorted.end(), fSorted, it);
      }
      else {

        // Otherwise (if propagation succeeded), set the path distance
        // and advance the iterator to the next KHitGroup.

        gr.setPath(true, *dist);
        ++igr;
      }
    }

    // Finally, sort the sorted list in order of path distance.

    fSorted.sort();
  }

  /// Return the plane with the most KHitGroups in the unsorted list.
  unsigned int KHitContainer::getPreferredPlane() const
  {
    // Count hits in each plane.

    std::vector<unsigned int> planehits(3, 0);

    // Loop over KHitGroups in the unsorted list.

    for (std::list<KHitGroup>::const_iterator igr = fUnsorted.begin(); igr != fUnsorted.end();
         ++igr) {

      const KHitGroup& gr = *igr;

      // Get plane of this KHitGroup.

      int plane = gr.getPlane();
      ++planehits.at(plane);
    }

    // Figure out which plane has the most hits.

    unsigned int prefplane = 0;
    for (unsigned int i = 0; i < planehits.size(); ++i) {
      if (planehits[i] >= planehits[prefplane]) prefplane = i;
    }
    return prefplane;
  }

} // end namespace trkf
