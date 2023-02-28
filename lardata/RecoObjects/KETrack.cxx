///////////////////////////////////////////////////////////////////////
///
/// \file   KETrack.cxx
///
/// \brief  Basic Kalman filter track class, with error.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KETrack.h"
#include "cetlib_except/exception.h"
#include <cmath>

namespace trkf {

  /// Default constructor.
  KETrack::KETrack() {}

  /// Constructor - specify surface only.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  ///
  KETrack::KETrack(const std::shared_ptr<const Surface>& psurf) : KTrack(psurf) {}

  /// Constructor - surface + track parameters + error matrix.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  /// vec   - Track state vector.
  /// err   - Track error matrix.
  /// dir   - Track direction.
  /// pdg   - Pdg code.
  ///
  KETrack::KETrack(const std::shared_ptr<const Surface>& psurf,
                   const TrackVector& vec,
                   const TrackError& err,
                   Surface::TrackDirection dir,
                   int pdg)
    : KTrack(psurf, vec, dir, pdg), fErr(err)
  {}

  /// Constructor - KTrack + error matrix.
  ///
  /// Arguments:
  ///
  /// trk - KTrack.
  /// err - Track error matrix.
  ///
  KETrack::KETrack(const KTrack& trk, const TrackError& err) : KTrack(trk), fErr(err) {}

  /// Destructor.
  KETrack::~KETrack() {}

  /// Calculate track pointing error (sigma, in radians).
  ///
  /// This method calculates a single pointing error (sigma, in
  /// radians) based on the track parameters and error matrix.  The
  /// actual calculation is done by the similarly names method of the
  /// surface class, since this class doesn't know what the track
  /// parameters mean.
  ///
  double KETrack::PointingError() const
  {
    if (!isValid())
      throw cet::exception("KETrack") << "Pointing error requested for invalid track.\n";
    return getSurface()->PointingError(getVector(), fErr);
  }

  /// Combine two tracks.
  ///
  /// Arguments:
  ///
  /// tre - Another track.
  ///
  /// Returns: Chisquare + success flag.
  ///
  /// This method updates the current track to be the weighted average
  /// of itself and another track.  The chisquare of the combination
  /// is returned as the result value.  The combination can fail
  /// because the sum of the two error matrices is singular, in which
  /// case the success flag embedded in the return value is false.
  ///
  std::optional<double> KETrack::combineTrack(const KETrack& tre)
  {
    // Make sure that the two track surfaces are the same.
    // Throw an exception if they are not.

    if (!getSurface()->isEqual(*tre.getSurface()))
      throw cet::exception("KETrack") << "Track combination surfaces are not the same.\n";

    // Default result is failure.

    std::optional<double> result{std::nullopt};

    // We will use asymmetric versions of the updating formulas, such
    // that the result is calculated as a perturbation on the
    // better-measured track.  We define the better measured track as
    // the one with the smaller error matrix trace.

    // Extract the two state vectors and error matrices as pointers.

    const TrackVector* vec1 = &getVector();
    const TrackError* err1 = &getError();
    const TrackVector* vec2 = &tre.getVector();
    const TrackError* err2 = &tre.getError();

    // Calculate the traces of the error matrices.

    double tr1 = 0;
    for (unsigned int i = 0; i < err1->size1(); ++i)
      tr1 += (*err1)(i, i);

    double tr2 = 0;
    for (unsigned int i = 0; i < err2->size1(); ++i)
      tr2 += (*err2)(i, i);

    // Define vec1, err1 as belong to the better measured track.
    // Swap if necessary.

    if (tr1 > tr2) {
      const TrackVector* tvec = vec1;
      vec1 = vec2;
      vec2 = tvec;
      const TrackError* terr = err1;
      err1 = err2;
      err2 = terr;
    }

    // Calculate the difference vector and difference error matrix.

    TrackVector dvec = *vec1 - *vec2;
    TrackError derr = *err1 + *err2;

    // Invert the difference error matrix.
    // This is the only place where a detectable failure can occur.

    bool ok = syminvert(derr);
    if (ok) {

      // Calculate updated state vector.
      // vec1 = vec1 - err1 * derr * dvec

      TrackVector tvec1 = prod(derr, dvec);
      TrackVector tvec2 = prod(*err1, tvec1);
      TrackVector tvec3 = *vec1 - tvec2;
      setVector(tvec3);

      // Calculate updated error matrix.
      // err1 = err1 - err1 * derr * err1

      TrackMatrix terr1 = prod(derr, *err1);
      TrackMatrix terr2 = prod(*err1, terr1);
      TrackError terr2s = ublas::symmetric_adaptor<TrackMatrix>(terr2);
      TrackError terr3 = *err1 - terr2s;
      setError(terr3);

      // Calculate chisquare.
      // chisq = dvec^T * derr * dvec

      TrackVector dvec1 = prod(derr, dvec);
      double chisq = inner_prod(dvec, dvec1);
      result = std::make_optional(chisq);
    }

    // Final validity check.

    if (!isValid()) result = std::nullopt;

    // Done.

    return result;
  }

  /// Printout
  std::ostream& KETrack::Print(std::ostream& out, bool doTitle) const
  {
    if (doTitle) out << "KETrack:\n";

    // Print base class.

    KTrack::Print(out, false);

    // Print diagonal errors.

    out << "  Diagonal errors:\n"
        << "  [";
    for (unsigned int i = 0; i < fErr.size1(); ++i) {
      if (i != 0) out << ", ";
      double err = fErr(i, i);
      err = (err >= 0. ? std::sqrt(err) : -std::sqrt(-err));
      out << err;
    }
    out << "]\n";

    // Print correlations.

    out << "  Correlation matrix:";
    for (unsigned int i = 0; i < fErr.size1(); ++i) {
      if (i == 0)
        out << "\n  [";
      else
        out << "\n   ";
      for (unsigned int j = 0; j <= i; ++j) {
        if (j != 0) out << ", ";
        if (i == j)
          out << 1.;
        else {
          double eiijj = fErr(i, i) * fErr(j, j);
          double eij = fErr(i, j);
          if (eiijj != 0.)
            eij /= std::sqrt(std::abs(eiijj));
          else
            eij = 0.;
          out << eij;
        }
      }
    }
    out << "]\n";
    out << "  Pointing error = " << PointingError() << "\n";
    return out;
  }

} // end namespace trkf
