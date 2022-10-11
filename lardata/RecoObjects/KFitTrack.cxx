///////////////////////////////////////////////////////////////////////
///
/// \file   KFitTrack.cxx
///
/// \brief  Basic Kalman filter track class, with fit information.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KFitTrack.h"
#include "cetlib_except/exception.h"

namespace trkf {

  /// Default constructor.
  KFitTrack::KFitTrack() : fPath(0.), fChisq(0.), fStat(INVALID) {}

  /// Initializing constructor.
  ///
  /// Arguments:
  ///
  /// tre   - KETrack.
  /// s     - Path distance.
  /// chisq - Fit chisquare.
  /// stat  - Fit status.
  ///
  KFitTrack::KFitTrack(const KETrack& tre, double s, double chisq, FitStatus stat)
    : KETrack(tre), fPath(s), fChisq(chisq), fStat(stat)
  {}

  /// Destructor.
  KFitTrack::~KFitTrack() {}

  /// Combine two tracks.
  ///
  /// Arguments:
  ///
  /// trf - Another track.
  ///
  /// Returns: True if success.
  ///
  /// This method updates the current track to be the weighted average
  /// of itself and another track.  Track parameters and error matrix
  /// are updated by method KETrack::combineTrack.  The updated
  /// chissquare is the sum of the original chisquare, the chisquare
  /// from the other track, and the chisquare of the combination.  The
  /// path distance is not updated.
  ///
  /// The updated status is derived from the input
  /// status according to the following table.
  ///
  /// FORWARD + BACKWARD_PREDICTED = OPTIMAL
  /// FORWARD_PREDICTED + BACKWARD = OPTIMAL
  /// FORWARD_PREDICTED + BACKWARD_PREDICTED = OPTIMAL_PREDICTED
  ///
  /// Any other combination of input statuses will cause this method
  /// to return failure.  This method can also return failure because
  /// of singular error matrices.  In case of failure, the original
  /// track is not updated.
  ///
  /// The two tracks being combined should be on the same surface
  /// (throw exception if not).
  ///
  bool KFitTrack::combineFit(const KFitTrack& trf)
  {
    // Make sure that the two track surfaces are the same.
    // Throw an exception if they are not.

    if (!getSurface()->isEqual(*trf.getSurface()))
      throw cet::exception("KFitTrack") << "Track combination surfaces are not the same.\n";

    // Default result failure.

    bool result = false;

    // Make sure input statuses are compatible and figure out the
    // updated status.

    FitStatus stat1 = getStat();
    FitStatus stat2 = trf.getStat();
    FitStatus statu = UNKNOWN;
    if ((stat1 == FORWARD && stat2 == BACKWARD_PREDICTED) ||
        (stat1 == FORWARD_PREDICTED && stat2 == BACKWARD) ||
        (stat1 == BACKWARD && stat2 == FORWARD_PREDICTED) ||
        (stat1 == BACKWARD_PREDICTED && stat2 == FORWARD))
      statu = OPTIMAL;
    else if (stat1 == FORWARD_PREDICTED && stat2 == BACKWARD_PREDICTED)
      statu = OPTIMAL_PREDICTED;
    else
      statu = UNKNOWN;
    if (statu != UNKNOWN) {

      // Update track parameters and error matrix.

      std::optional<double> chisq = combineTrack(trf);

      // Update status and chisquare.

      if (!!chisq) {
        result = true;
        fStat = statu;
        fChisq = fChisq + trf.getChisq() + *chisq;
      }
    }

    // Done.

    return result;
  }

  /// Printout
  std::ostream& KFitTrack::Print(std::ostream& out, bool doTitle) const
  {
    if (doTitle) out << "KFitTrack:\n";

    // Print information specific to this class.

    out << "  Distance = " << fPath << "\n";
    out << "  Chisquare = " << fChisq << "\n";
    out << "  Status = "
        << (fStat == INVALID ?
              "INVALID" :
              (fStat == UNKNOWN ?
                 "UNKNOWN" :
                 (fStat == FORWARD ?
                    "FORWARD" :
                    (fStat == FORWARD_PREDICTED ?
                       "FORWARD_PREDICTED" :
                       (fStat == BACKWARD ?
                          "BACKWARD" :
                          (fStat == BACKWARD_PREDICTED ?
                             "BACKWARD_PREDICTED" :
                             (fStat == OPTIMAL ? "OPTIMAL" : "OPTIMAL_PREDICTED")))))))
        << "\n";

    // Print base class.

    KETrack::Print(out, false);
    return out;
  }

} // end namespace trkf
