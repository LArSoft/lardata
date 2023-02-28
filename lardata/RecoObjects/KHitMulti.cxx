///////////////////////////////////////////////////////////////////////
///
/// \file   KHitMulti.cxx
///
/// \brief  Compound Kalman Filter measurement.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KHitMulti.h"
#include "cetlib_except/exception.h"

namespace trkf {

  /// Default Constructor.
  KHitMulti::KHitMulti() : fMeasDim(0), fChisq(0.) {}

  /// Initializing Constructor.
  ///
  /// Arguments:
  ///
  /// psurf - Measurement surface pointer.
  ///
  KHitMulti::KHitMulti(const std::shared_ptr<const Surface>& psurf)
    : KHitBase(psurf), fMeasDim(0), fChisq(0.)
  {}

  /// Destructor.
  KHitMulti::~KHitMulti() {}

  /// Add a measurement.
  ///
  /// Arguments:
  ///
  /// pmeas - Measurement.
  ///
  /// This method tries to dynamic cast the measurement to a supported
  /// type.  If the dynamic cast fails, throw an exception.
  ///
  void KHitMulti::addMeas(const std::shared_ptr<const KHitBase>& pmeas)
  {
    // It is an error to pass in a null pointer.

    if (pmeas.get() == 0)
      throw cet::exception("KHitMulti") << "Attempt to add null measurement pointer.\n";

    // Do the dynamic cast.

    std::shared_ptr<const KHit<1>> pmeas1 =
      std::dynamic_pointer_cast<const KHit<1>, const KHitBase>(pmeas);

    // Throw an exception if dynamic cast failed.

    if (pmeas1.get() == 0)
      throw cet::exception("KHitMulti") << "Dynamic cast for KHitBase pointer failed.\n";
    addMeas(pmeas1);
  }

  /// Add a measurement.
  ///
  /// Arguments:
  ///
  /// pmeas - Measurement.
  ///
  void KHitMulti::addMeas(const std::shared_ptr<const KHit<1>>& pmeas)
  {
    // It is an error to pass in a null pointer.

    if (pmeas.get() == 0)
      throw cet::exception("KHitMulti") << "Attempt to add null measurement pointer.\n";

    // Add the measurement.

    ++fMeasDim;
    fMeasVec.push_back(pmeas);
  }

  /// Prediction method (return false if fail).
  ///
  /// Arguments:
  ///
  /// tre  - Track hypothesis.
  /// prop - Propagator.
  /// ref  - Reference track.
  ///
  /// Returns: True if success.
  ///
  /// This class calls the predict method of each underlying
  /// measurement and updates the combined prediction attributes.
  ///
  bool KHitMulti::predict(const KETrack& tre, const Propagator& prop, const KTrack* ref) const
  {
    // Resize and clear all linear algebra objects.

    fMvec.resize(fMeasDim, false);
    fMvec.clear();

    fMerr.resize(fMeasDim, false);
    fMerr.clear();

    fPvec.resize(fMeasDim, false);
    fPvec.clear();

    fPerr.resize(fMeasDim, false);
    fPerr.clear();

    fRvec.resize(fMeasDim, false);
    fRvec.clear();

    fRerr.resize(fMeasDim, false);
    fRerr.clear();

    fRinv.resize(fMeasDim, false);
    fRinv.clear();

    fH.resize(fMeasDim, tre.getVector().size());
    fH.clear();

    // Update the prediction surface to be the track surface.

    fPredSurf = tre.getSurface();
    fPredDist = 0.;

    // Result.

    bool ok = true;

    // Loop over one-dimensional measurements.

    for (unsigned int im = 0; ok && im < fMeasVec.size(); ++im) {
      const KHit<1>& meas = *(fMeasVec[im]);

      // Update prediction for this measurement.

      ok = meas.predict(tre, prop, ref);
      if (!ok) break;

      //

      // Update objects that are concatenations of underlying measurements.

      fMvec(im) = meas.getMeasVector()(0);       // Measurement vector.
      fMerr(im, im) = meas.getMeasError()(0, 0); // Measurement error matrix.
      fPvec(im) = meas.getPredVector()(0);       // Prediction vector.

      // H-matrix.

      for (unsigned int j = 0; j < meas.getH().size2(); ++j)
        fH(im, j) = meas.getH()(0, j);
    }
    if (ok) {

      // Calculate prediction error matrix.
      // T = H C H^T.

      ublas::matrix<double> temp(fH.size2(), fMeasDim);
      ublas::matrix<double> temp2(fMeasDim, fMeasDim);
      temp = prod(tre.getError(), trans(fH));
      temp2 = prod(fH, temp);
      fPerr = ublas::symmetric_adaptor<ublas::matrix<double>>(temp2);

      // Update residual

      fRvec = fMvec - fPvec;
      fRerr = fMerr + fPerr;
      fRinv = fRerr;
      ok = syminvert(fRinv);
      if (ok) {

        // Calculate incremental chisquare.

        ublas::vector<double> rtemp = prod(fRinv, fRvec);
        fChisq = inner_prod(fRvec, rtemp);
      }
    }

    // If a problem occured at any step, clear the prediction surface pointer.

    if (!ok) {
      fPredSurf.reset();
      fPredDist = 0.;
    }

    // Done.

    return ok;
  }

  /// Update track method.
  ///
  /// Arguments:
  ///
  /// tre - Track to be updated.
  ///
  /// This method is almost an exact copy of the update method in KHit<N>.
  ///
  void KHitMulti::update(KETrack& tre) const
  {
    // Make sure that the track surface and the prediction surface are the same.
    // Throw an exception if they are not.

    if (!getPredSurface()->isEqual(*tre.getSurface()))
      throw cet::exception("KHitMulti") << "Track surface not the same as prediction surface.\n";

    const TrackVector& tvec = tre.getVector();
    const TrackError& terr = tre.getError();
    TrackVector::size_type size = tvec.size();

    // Calculate gain matrix.

    ublas::matrix<double> temp(size, fMeasDim);
    ublas::matrix<double> gain(size, fMeasDim);
    temp = prod(trans(fH), fRinv);
    gain = prod(terr, temp);

    // Calculate updated track state.

    TrackVector newvec = tre.getVector() + prod(gain, fRvec);

    // Calculate updated error matrix.

    TrackMatrix fact = ublas::identity_matrix<TrackVector::value_type>(size);
    fact -= prod(gain, fH);
    TrackMatrix errtemp1 = prod(terr, trans(fact));
    TrackMatrix errtemp2 = prod(fact, errtemp1);
    TrackError errtemp2s = ublas::symmetric_adaptor<TrackMatrix>(errtemp2);
    ublas::matrix<double> errtemp3 = prod(fMerr, trans(gain));
    TrackMatrix errtemp4 = prod(gain, errtemp3);
    TrackError errtemp4s = ublas::symmetric_adaptor<TrackMatrix>(errtemp4);
    TrackError newerr = errtemp2s + errtemp4s;

    // Update track.

    tre.setVector(newvec);
    tre.setError(newerr);
  }

  /// Printout
  std::ostream& KHitMulti::Print(std::ostream& out, bool doTitle) const
  {
    if (doTitle) out << "KHitMulti:\n";
    return out;
  }

} // end namespace trkf
