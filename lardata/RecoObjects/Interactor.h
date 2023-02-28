////////////////////////////////////////////////////////////////////////
///
/// \file   Interactor.h
///
/// \brief  Base class for Kalman filter track interactor.
///
/// \author H. Greenlee
///
/// This class defined the general interface for calculating
/// propagation noise.
///
/// This class defined a single pure virtual method called "noise"
/// whose purpose is to calculate the propagation noise matrix
/// associated with the propagation of a track over a specified
/// distance.
///
////////////////////////////////////////////////////////////////////////

#ifndef INTERACTOR_H
#define INTERACTOR_H

#include "lardata/RecoObjects/KTrack.h"
#include "lardata/RecoObjects/KalmanLinearAlgebra.h"

namespace trkf {

  class Interactor {
  public:
    explicit Interactor(double tcut);
    virtual ~Interactor();

    double getTcut() const { return fTcut; }

    /// Clone method.
    virtual Interactor* clone() const = 0;

    /// Calculate noise matrix.
    virtual bool noise(const KTrack& trk, double s, TrackError& noise_matrix) const = 0;

  private:
    double fTcut; ///< Maximum delta ray energy for dE/dx.
  };
}

#endif
