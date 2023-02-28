////////////////////////////////////////////////////////////////////////
///
/// \file   KTrack.h
///
/// \brief  Basic Kalman filter track class, without error.
///
/// \author H. Greenlee
///
/// This class include the following attributes.
///
/// 1. Surface.
/// 2. Track state vector.
/// 3. Track direction parameter.
/// 4. Particle id hypothesis.
///
/// The surface attribute is polymorphic, and is held via
/// std::shared_ptr type of smart pointer, which handles memory
/// management using reference-counted shared ownership.
///
////////////////////////////////////////////////////////////////////////

#ifndef KTRACK_H
#define KTRACK_H

#include "lardata/RecoObjects/KalmanLinearAlgebra.h"
#include "lardata/RecoObjects/Surface.h"
#include <iosfwd>
#include <memory>

namespace trkf {

  class KTrack {
  public:
    /// Enum

    /// Default constructor.
    KTrack();

    /// Constructor - specify surface only.
    KTrack(const std::shared_ptr<const Surface>& psurf);

    /// Constructor - surface + track parameters.
    KTrack(std::shared_ptr<const Surface> psurf,
           const TrackVector& vec,
           Surface::TrackDirection dir = Surface::UNKNOWN,
           int pdg = 0);

    /// Destructor.
    virtual ~KTrack();

    // Accessors.

    const std::shared_ptr<const Surface>& getSurface() const { return fSurf; } ///< Surface.
    const TrackVector& getVector() const { return fVec; } ///< Track state vector.
    Surface::TrackDirection getDirection() const;         ///< Track direction.
    int PdgCode() const { return fPdgCode; }              ///< Pdg code.
    double Mass() const;                                  ///< Based on pdg code.

    // Modifiers.

    TrackVector& getVector() { return fVec; } ///< Modifiable state vector.

    /// Set surface.
    void setSurface(const std::shared_ptr<const Surface>& psurf) { fSurf = psurf; }
    void setVector(const TrackVector& vec) { fVec = vec; }         ///< Set state vector.
    void setDirection(Surface::TrackDirection dir) { fDir = dir; } ///< Set direction.
    void setPdgCode(int pdg) { fPdgCode = pdg; }                   ///< Set pdg code.

    /// Test if track is valid.
    bool isValid() const;

    /// Get position of track.
    void getPosition(double xyz[3]) const;

    /// Get momentum vector of track.
    void getMomentum(double mom[3]) const;

    /// Get x-latitude.
    double XLatitude() const;

    /// Get x-longitude.
    double XLongitude() const;

    /// Printout
    virtual std::ostream& Print(std::ostream& out, bool doTitle = true) const;

  private:
    // Attributes.

    std::shared_ptr<const Surface> fSurf; ///< Track surface.
    TrackVector fVec;                     ///< Track state vector.
    Surface::TrackDirection fDir;         ///< Track direction.
    int fPdgCode;                         ///< Pdg id. hypothesis.
  };

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const KTrack& trk);
}

#endif
