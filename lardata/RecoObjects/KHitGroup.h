////////////////////////////////////////////////////////////////////////
///
/// \file   KHitGroup.h
///
/// \brief  A collection of measurements on the same surface.
///
/// \author H. Greenlee
///
/// This class represents a collection of measurements on a common
/// surface.  The measurements are polymorphic, stored as pointers
/// to KHitBase base class.
///
/// The idea behind this class is that the contained measurement are
/// mutually exclusive for inclusion in a single track.
///
/// This class includes the following attributes.
///
/// 1.  Pointer to common surface.
/// 2.  Plane index.
/// 3.  Measurement collection.
/// 4.  Estimated path flag.
/// 5.  Estimated path distance.
///
/// The last two attributes is included as an aid in sorting
/// measurements for inclusion in tracks.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHITGROUP_H
#define KHITGROUP_H

#include "lardata/RecoObjects/KHitBase.h"
#include <vector>

namespace trkf {

  class KHitGroup {
  public:
    /// Default constructor.
    KHitGroup(bool has_path = false, double path = 0.);

    /// Destructor.
    virtual ~KHitGroup();

    // Accessors.

    /// Surface accessor.
    const std::shared_ptr<const Surface>& getSurface() const { return fSurf; }

    /// Plane index.
    int getPlane() const { return fPlane; }

    /// Measurement collection accessor.
    const std::vector<std::shared_ptr<const KHitBase>>& getHits() const { return fHits; }

    /// Path flag.
    bool getHasPath() const { return fHasPath; }

    /// Estimated path distance.
    double getPath() const { return fPath; }

    // Modifiers.

    /// Clear the collection.
    void clear() { fHits.clear(); }

    /// Add a mesaurement into the colleciton.
    void addHit(const std::shared_ptr<const KHitBase>& hit);

    /// Set path flag and estimated path distance.
    void setPath(bool has_path, double path)
    {
      fHasPath = has_path;
      fPath = path;
    }

    // Relational operators, sort by estimated path distance.

    bool operator==(const KHitGroup& obj) const; ///< Equivalance operator.
    bool operator<(const KHitGroup& obj) const;  ///< Less than operator.

  private:
    // Attributes.

    std::shared_ptr<const Surface> fSurf;               ///< Common surface.
    int fPlane;                                         ///< Plane index of measurements.
    std::vector<std::shared_ptr<const KHitBase>> fHits; ///< Measuement collection.
    bool fHasPath;                                      ///< Path flag.
    double fPath;                                       ///< Estimated path distance.
  };
}

#endif
