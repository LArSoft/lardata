////////////////////////////////////////////////////////////////////////
///
/// \file   KETrack.h
///
/// \brief  Basic Kalman filter track class, with error.
///
/// \author H. Greenlee
///
/// This class inherits the following attributes from KTrack.
///
/// 1. Surface.
/// 2. Track state vector.
/// 3. Track direction parameter.
///
/// This class adds the following attributes of its own.
///
/// 4. Track error matrix.
///
////////////////////////////////////////////////////////////////////////

#ifndef KETRACK_H
#define KETRACK_H

#include "lardata/RecoObjects/KTrack.h"
#include <optional>

namespace trkf {

  class KETrack : public KTrack {
  public:
    /// Default constructor.
    KETrack();

    /// Constructor - specify surface only.
    KETrack(const std::shared_ptr<const Surface>& psurf);

    /// Constructor - surface + track parameters + error matrix.
    KETrack(const std::shared_ptr<const Surface>& psurf,
            const TrackVector& vec,
            const TrackError& err,
            Surface::TrackDirection dir = Surface::UNKNOWN,
            int pdg = 0);

    /// Constructor - KTrack + error matrix.
    KETrack(const KTrack& trk, const TrackError& err);

    /// Destructor.
    virtual ~KETrack();

    // Accessors.

    const TrackError& getError() const { return fErr; } ///< Track error matrix.
    double PointingError() const;                       ///< Pointing error (radians).

    // Modifiers.

    TrackError& getError() { return fErr; }              ///< Modifiable error matrix.
    void setError(const TrackError& err) { fErr = err; } ///< Set error matrix.

    /// Combine two tracks.
    std::optional<double> combineTrack(const KETrack& tre);

    /// Printout
    virtual std::ostream& Print(std::ostream& out, bool doTitle = true) const;

  private:
    // Attributes.

    TrackError fErr; ///< Track error matrix.
  };
}

#endif
