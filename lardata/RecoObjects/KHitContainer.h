////////////////////////////////////////////////////////////////////////
///
/// \file   KHitContainer.h
///
/// \brief  A collection of KHitGroups.
///
/// \author H. Greenlee
///
/// This class internally maintains three STL lists of KHitGroup objects.
///
/// 1.  Sorted KHitGroup objects (have path length).
/// 2.  Unsorted KHitGroup objects (don't currently have path length).
/// 3.  Unused KHitGroup objects.
///
/// The following methods are provided.
///
/// 1.  Sort
///
/// A KTrack object and propagation direction are passed as arguments.
/// The track is propagated without error to each object on the sorted
/// and (maybe) the unsorted list.  Reachable objects have their path
/// length updated, are moved to the sorted list, and are eventually
/// sorted.  Unreachable objects are moved to the unsorted list.
///
/// Here are the envisioned use cases of this class.
///
/// 1.  At the beginning of the event, a set of candidate measurements
///     are loaded into the unsorted list.
/// 2.  Candidate measurements are sorted using a seed track.
/// 3.  During the progress of the Kalman filter, candidate measurements
///     are visited in order from the sorted list.
/// 4.  If necessary, candidate measurements can be resorted during the
///     progress of the Kalman filter using the updated track.
/// 5.  After candidate measurements are disposed of (added to track or
///     rejected), they are moved to the unused list.
/// 6.  The Kalman filter can be repeated using a new seed track by
///     moving all objects to the unsorted list.
///
/// Most of these use cases involve transfering objects among the three
/// lists.  These kinds of operations can be accomplished using STL
/// list splice method without copying the objects.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHITCONTAINER_H
#define KHITCONTAINER_H

#include "canvas/Persistency/Common/PtrVector.h"
#include "lardata/RecoObjects/KHitGroup.h"
#include "lardata/RecoObjects/KTrack.h"
#include "lardata/RecoObjects/Propagator.h"
#include "lardataobj/RecoBase/Hit.h"
#include <list>

namespace detinfo {
  class DetectorPropertiesData;
}

namespace trkf {

  class KHitContainer {
  public:
    virtual ~KHitContainer() = default;

    virtual void fill(detinfo::DetectorPropertiesData const& clock_data,
                      const art::PtrVector<recob::Hit>& hits,
                      int only_plane) = 0;

    const std::list<KHitGroup>& getSorted() const { return fSorted; }
    const std::list<KHitGroup>& getUnsorted() const { return fUnsorted; }
    const std::list<KHitGroup>& getUnused() const { return fUnused; }

    // Non-const Accessors.

    std::list<KHitGroup>& getSorted() { return fSorted; }     ///< Sorted list.
    std::list<KHitGroup>& getUnsorted() { return fUnsorted; } ///< Unsorted list.
    std::list<KHitGroup>& getUnused() { return fUnused; }     ///< Unused list.

    /// Clear all lists.
    void clear();

    /// Move all objects to unsorted list (from sorted and unused lists).
    void reset();

    /// (Re)sort objects in unsorted and sorted lists.
    void sort(const KTrack& trk,
              bool addUnsorted,
              const Propagator& prop,
              Propagator::PropDirection dir = Propagator::UNKNOWN);

    /// Return the plane with the most KHitGroups in the unsorted list.
    unsigned int getPreferredPlane() const;

  private:
    // Attributes.

    std::list<KHitGroup> fSorted;   ///< Sorted KHitGroup objects.
    std::list<KHitGroup> fUnsorted; ///< Unsorted KHitGroup objects.
    std::list<KHitGroup> fUnused;   ///< Unused KHitGroup objects.
  };
}

#endif
