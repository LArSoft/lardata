#ifndef KFTRACKSTATE_H
#define KFTRACKSTATE_H

#include "lardata/RecoObjects/TrackState.h"

namespace trkf {

  /// \file  lardata/RecoObjects/KFTrackState.h
  /// \class KFTrackState
  ///
  /// \brief Extension of a TrackState to perform KalmanFilter calculations.
  ///
  /// \author  G. Cerati (FNAL, MicroBooNE)
  /// \date    2017
  /// \version 1.0
  ///
  /// This class extends the concept of TrackState, providing functionalities for Kalman Filter-specific calculations such as 'update' and 'combine'.
  /// It holds a TrackState by value, which is modified in place when updating or combining.
  ///

  class KFTrackState {
  public:
    //
    KFTrackState(const SVector5& trackStatePar,
                 const SMatrixSym55& trackStateCov,
                 const Plane& plane,
                 bool trackAlongPlaneDir,
                 int pid)
      : fTrackState(trackStatePar, trackStateCov, plane, trackAlongPlaneDir, pid)
    {}
    KFTrackState(TrackState&& trackState) : fTrackState(std::move(trackState)) {}

    /// Update the TrackState given a HitState (they need to be on the same plane)
    bool updateWithHitState(const HitState& hitstate);

    /// Combine the TrackState given another TrackState (they need to be on the same plane)
    bool combineWithTrackState(const TrackState& trackstate);

    /// Get the (const reference to the) TrackState
    const TrackState& trackState() const { return fTrackState; }
    /// Set the TrackState
    void setTrackState(TrackState&& s) { fTrackState = std::move(s); }
    //
    //@{
    /// This function calls the homonymous function of the stored TrackState
    const SVector5& parameters() const { return fTrackState.parameters(); }
    const SMatrixSym55& covariance() const { return fTrackState.covariance(); }
    const Plane& plane() const { return fTrackState.plane(); }
    const Point_t& position() const { return fTrackState.position(); }
    const Vector_t& momentum() const { return fTrackState.momentum(); }
    int pID() const { return fTrackState.pID(); }
    double mass() const { return fTrackState.mass(); }
    const SVector6 parameters6D() const { return fTrackState.parameters6D(); }
    bool isTrackAlongPlaneDir() const { return fTrackState.isTrackAlongPlaneDir(); }
    //
    double residual(const HitState& hitstate) const { return fTrackState.residual(hitstate); }
    double combinedError2(const HitState& hitstate) const
    {
      return fTrackState.combinedError2(hitstate);
    }
    double combinedError(const HitState& hitstate) const
    {
      return fTrackState.combinedError(hitstate);
    }
    double chi2(const HitState& hitstate) const { return fTrackState.chi2(hitstate); }
    //
    void setCovariance(const SMatrixSym55& trackStateCov)
    {
      fTrackState.setCovariance(trackStateCov);
    }
    void setParameters(const SVector5& trackStatePar) { fTrackState.setParameters(trackStatePar); }
    //@}
    //
    /// Printout information
    std::ostream& dump(std::ostream& out = std::cout) const
    {
      out << "KFTrackState with pID=" << pID() << " mass=" << mass() << "\npars=" << parameters()
          << " position=" << position() << " momentum=" << momentum() << "\ncov=\n"
          << covariance() << "\non plane with pos=" << plane().position()
          << " and dir=" << plane().direction() << " along=" << isTrackAlongPlaneDir() << "\n";
      return out;
    }
    //
  private:
    TrackState fTrackState;
  };

}
#endif
