////////////////////////////////////////////////////////////////////////
///
/// \file   KHitsTrack.h
///
/// \brief  Basic Kalman filter track class, with measurements.
///
/// \author H. Greenlee
///
/// This class inherits the following attributes from KFitTrack.
///
/// 1. Surface.
/// 2. Track state vector.
/// 3. Track direction parameter.
/// 4. Track error matrix.
/// 5. Propagation distance.
/// 6. Fit chisquare.
/// 7. Fit status.
///
/// This class adds the following attributes of its own.
///
/// 8. A collection of measurements, in the form of a vector
///    of shared pointers to KHitBase.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHITSTRACK_H
#define KHITSTRACK_H

#include "lardata/RecoObjects/KFitTrack.h"
#include "lardata/RecoObjects/KHitBase.h"
#include <vector>

namespace trkf {

  class KHitsTrack : public KFitTrack {
  public:
    /// Default constructor.
    KHitsTrack();

    /// Initializing constructor - KFitTrack.
    KHitsTrack(const KFitTrack& trf);

    /// Initializing constructor - KETrack.
    KHitsTrack(const KETrack& tre);

    /// Destructor.
    virtual ~KHitsTrack();

    // Accessor.

    /// Measurement collection.
    const std::vector<std::shared_ptr<const KHitBase>>& getHits() { return fHits; }

    // Modifiers.

    /// Add a measurement.
    void addHit(const std::shared_ptr<const KHitBase>& hit) { fHits.push_back(hit); }

    /// Printout
    virtual std::ostream& Print(std::ostream& out, bool doTitle = true) const;

  private:
    // Attributes.

    /// Measurement collection.
    std::vector<std::shared_ptr<const KHitBase>> fHits;
  };
}

#endif
