///////////////////////////////////////////////////////////////////////
///
/// \file   KHitGroup.cxx
///
/// \brief  A collection of measurements on the same surface.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KHitGroup.h"
#include "cetlib_except/exception.h"

namespace trkf {

  /// Default Constructor.
  ///
  /// Arguments:
  ///
  /// has_path - Path flag (optional).
  /// path     - Estimated path distance (optional).
  ///
  KHitGroup::KHitGroup(bool has_path, double path) : fPlane(-1), fHasPath(has_path), fPath(path) {}

  /// Destructor.
  KHitGroup::~KHitGroup() {}

  /// Add a mesaurement into the colleciton.
  ///
  /// Arguments:
  ///
  /// hit - Measurement to add.
  ///
  void KHitGroup::addHit(const std::shared_ptr<const KHitBase>& hit)
  {
    // Make sure that the measurement pointer is not null (throw exception if null).

    if (hit.get() == 0) throw cet::exception("KHitGroup") << "Attempt to add null measurement.\n";

    // If the stored common surface pointer has not yet been initialized,
    // initialize it now.
    // Otherwise, make sure that the new measurement surface matches
    // the common surface (pointer to same surface object).
    // Throw exception if the new surface doesn't match.

    if (hit->getMeasPlane() < 0)
      throw cet::exception("KHitGroup")
        << __func__ << ": invalid hit plane " << hit->getMeasPlane() << "\n";

    if (fSurf.get() == 0) {
      fSurf = hit->getMeasSurface();
      fPlane = hit->getMeasPlane();
    }
    else {
      if (fSurf.get() != hit->getMeasSurface().get())
        throw cet::exception("KHitGroup") << "Attempt to add non-matching measurement.\n";
      if (hit->getMeasPlane() != fPlane) {
        throw cet::exception("KHitGroup") << __func__ << ": hit plane mismatch, "
                                          << hit->getMeasPlane() << " vs. " << fPlane << "\n";
      }
      if (fPlane < 0)
        throw cet::exception("KHitGroup") << __func__ << ": invalid plane " << fPlane << "\n";
    }

    // Everything OK.  Add the measurement.

    fHits.push_back(hit);
  }

  /// Equivalance operator.
  ///
  /// Objects with path flag false compare equal.
  /// Objects with path flag true compare according to estimated path distance.
  ///
  /// Arguments:
  ///
  /// g - Comparison object.
  ///
  /// Returned value: True if equal.
  ///
  bool KHitGroup::operator==(const KHitGroup& obj) const
  {
    bool result = (!fHasPath && !obj.fHasPath) || (fHasPath && obj.fHasPath && fPath == obj.fPath);
    return result;
  }

  /// Less than operator.
  ///
  /// Objects with path flag false compare equal (return false).
  /// Objects with path flag true compare according to estimated path distance.
  /// Objects with path flag false are not comparable to objects
  /// with path flag true (throw exception).
  ///
  /// Arguments:
  ///
  /// g - Comparison object.
  ///
  /// Returned value: True if less than.
  ///
  bool KHitGroup::operator<(const KHitGroup& obj) const
  {
    bool result = false;
    if (fHasPath != obj.fHasPath)
      throw cet::exception("KHitGroup") << "Attempt to compare incomparable objects.\n";
    if (fHasPath) result = fPath < obj.fPath;
    return result;
  }

} // end namespace trkf
