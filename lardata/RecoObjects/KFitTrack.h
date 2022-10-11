////////////////////////////////////////////////////////////////////////
///
/// \file   KFitTrack.h
///
/// \brief  Basic Kalman filter track class, with fit information.
///
/// \author H. Greenlee
///
/// This class inherits the following attributes from KETrack.
///
/// 1. Surface.
/// 2. Track state vector.
/// 3. Track direction parameter.
/// 4. Track error matrix.
///
/// This class adds the following attributes of its own.
///
/// 5. Propagation distance.
/// 6. Fit chisquare.
/// 7. Fit status.
///
/// The fit status attribute is specified using enum FitStatus.
/// The meanings of the FitStatus enum are given below.
/// In general, any track that is intended to be used for physics
/// analysis should have status OPTIMAL.
///
////////////////////////////////////////////////////////////////////////

#ifndef KFITTRACK_H
#define KFITTRACK_H

#include "lardata/RecoObjects/KETrack.h"

namespace trkf {

  class KFitTrack : public KETrack {
  public:
    /// Fit status enum.
    enum FitStatus {
      INVALID,            // No valid fit information.
      UNKNOWN,            // Unknown.
      FORWARD,            // Fit based on past measurements, including current surface.
      FORWARD_PREDICTED,  // Fit based on past measurements, not including current surface.
      BACKWARD,           // Fit based on future measurements, including current surface.
      BACKWARD_PREDICTED, // Fit based on future measurements, not including current surface.
      OPTIMAL,            // Fit based on all measurements, including current surface.
      OPTIMAL_PREDICTED   // Fit based on all measurements, except current surface.
    };

    /// Default constructor.
    KFitTrack();

    /// Initializing constructor.
    KFitTrack(const KETrack& tre, double s = 0., double chisq = 0., FitStatus stat = INVALID);

    /// Destructor.
    virtual ~KFitTrack();

    // Accessors.

    double getPath() const { return fPath; }    ///< Propagation distance.
    double getChisq() const { return fChisq; }  ///< Fit chisquare.
    FitStatus getStat() const { return fStat; } ///< Fit status.

    // Modifiers.

    void setPath(double path) { fPath = path; }     ///< Set propagation distance.
    void setChisq(double chisq) { fChisq = chisq; } ///< Set chisquare.
    void setStat(FitStatus stat) { fStat = stat; }  ///< Set fit status.

    /// Combine two tracks.
    bool combineFit(const KFitTrack& trf);

    /// Printout
    virtual std::ostream& Print(std::ostream& out, bool doTitle = true) const;

  private:
    // Attributes.

    double fPath;    ///< Propagation distance.
    double fChisq;   ///< Fit chisquare.
    FitStatus fStat; ///< Fit status.
  };
}

#endif
