////////////////////////////////////////////////////////////////////////
///
/// \file   KHit.h
///
/// \brief  Kalman filter measurement class template.
///
/// \author H. Greenlee
///
/// Class KHit represents a general measurement on a surface.  It is
/// specialized compared to base class KHitBase by specifying the
/// dimension of the measurement vector by an integer template
/// parameter N.
///
/// KHit class inherits the following attribute from KHitBase.
///
/// 1.  Measurement surface.
/// 2.  Prediction surface.
///
/// KHit adds the following attributes.
///
/// 3.  Measurement vector.
/// 4.  Measurement error matrix.
/// 5.  Prediction vector.
/// 6.  Prediction error matrix.
/// 7.  Residual vector.
/// 8.  Residual error matrix.
/// 9.  Inverse of residual error matrix.
/// 10.  Kalman H-matrix.
/// 11.  Incremental chisquare.
///
/// The first two attributes (measurement vector + error matrix) are
/// filled during construction, and the remaining attributes are left
/// in a default state.
///
/// The remaining attributes (also prediction surface attribute of
/// base class) are filled by the prediction method.  The attributes
/// that are filled in by the prediction method are mutable, and
/// prediction method is const.  The actual calculation of the
/// prediction vector, prediction error matrix, and H-matrix are left
/// to be implemented by a derived class, which must override the pure
/// virtual method subpredict.  The remaining unfilled attributes are
/// calculated locally in this class in the prediction method
/// inherited from KHitBase.
///
/// The measurement and prediction surfaces are not required to be the
/// same.  If they are different, the method predict is required to
/// make an internal propagation from the prediction surface to the
/// measurement surface, which propagation influeces the calculated
/// H-matrix, as well as the prediction vector and error.
///
/// The intended use case is as follows.
///
/// 1.  Track (KETrack) is propagated to measurement surface.
///     Methods predict and update will throw an exception of track
///     surface does not match measurement surface.
/// 2.  Prediction is updated by calling method KHit::predict.
/// 3.  At this point the calling program can make a cut on the
///     incremental chisquare, returned by method KHit::chisq.
/// 4.  If the chisquare cut passes, update track by calling method
///     KHit::update.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHIT_H
#define KHIT_H

#include "cetlib_except/exception.h"
#include "lardata/RecoObjects/KHitBase.h"
#include "lardata/RecoObjects/Propagator.h"

namespace trkf {

  template <int N>
  class KHit : public KHitBase {
  public:
    /// Default constructor.
    KHit();

    /// Initializing Constructor -- surface only.
    KHit(const std::shared_ptr<const Surface>& psurf);

    /// Fully Initializing Constructor.
    KHit(const std::shared_ptr<const Surface>& psurf,
         const typename KVector<N>::type& mvec,
         const typename KSymMatrix<N>::type& merr);

    /// Destructor.
    virtual ~KHit();

    // Modifiers.

    /// Set measurement vector.
    void setMeasVector(const typename KVector<N>::type& mvec) { fMvec = mvec; }

    /// Set measurement error.
    void setMeasError(const typename KSymMatrix<N>::type& merr) { fMerr = merr; }

    // Accessors.

    /// Measurement vector.
    const typename KVector<N>::type& getMeasVector() const { return fMvec; }

    /// Measurement error matrix.
    const typename KSymMatrix<N>::type& getMeasError() const { return fMerr; }

    /// Prediction vector.
    const typename KVector<N>::type& getPredVector() const { return fPvec; }

    /// Prediction matrix.
    const typename KSymMatrix<N>::type& getPredError() const { return fPerr; }

    /// Residual vector.
    const typename KVector<N>::type& getResVector() const { return fRvec; }

    /// Residual error matrix.
    const typename KSymMatrix<N>::type& getResError() const { return fRerr; }

    /// Residual inv. error matrix.
    const typename KSymMatrix<N>::type& getResInvError() const { return fRinv; }

    /// Kalman H-matrix.
    const typename KHMatrix<N>::type& getH() const { return fH; }

    /// Incremental chisquare.
    double getChisq() const { return fChisq; }

    // Overrides.
    // Implementation of overrides is found at the bottom of this header.

    /// Prediction method (return false if fail).
    bool predict(const KETrack& tre, const Propagator& prop, const KTrack* ref = 0) const;

    /// Update track method.
    void update(KETrack& tre) const;

    // Pure virtual methods.

    /// Calculate prediction function (return via arguments).
    virtual bool subpredict(const KETrack& tre,
                            typename KVector<N>::type& pvec,
                            typename KSymMatrix<N>::type& perr,
                            typename KHMatrix<N>::type& hmatrix) const = 0;

    /// Printout
    virtual std::ostream& Print(std::ostream& out, bool doTitle = true) const;

  private:
    // Attributes.

    typename KVector<N>::type fMvec;            ///< Measurement vector.
    typename KSymMatrix<N>::type fMerr;         ///< Measurement error matrix.
    mutable typename KVector<N>::type fPvec;    ///< Prediction vector.
    mutable typename KSymMatrix<N>::type fPerr; ///< Prediction  error matrix.
    mutable typename KVector<N>::type fRvec;    ///< Residual vector.
    mutable typename KSymMatrix<N>::type fRerr; ///< Residual error matrix.
    mutable typename KSymMatrix<N>::type fRinv; ///< Residual inverse error matrix.
    mutable typename KHMatrix<N>::type fH;      ///< Kalman H-matrix.
    mutable double fChisq;                      ///< Incremental chisquare.
  };

  // Method implementations.

  /// Default constructor.
  template <int N>
  KHit<N>::KHit() : fChisq(0.)
  {}

  /// Initializing Constructor -- surface only.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  ///
  template <int N>
  KHit<N>::KHit(const std::shared_ptr<const Surface>& psurf) : KHitBase(psurf), fChisq(0.)
  {}

  /// Fully Initializing Constructor.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  /// mvec  - Measurement vector.
  /// merr  - Measurement error.
  ///
  template <int N>
  KHit<N>::KHit(const std::shared_ptr<const Surface>& psurf,
                const typename KVector<N>::type& mvec,
                const typename KSymMatrix<N>::type& merr)
    : KHitBase(psurf), fMvec(mvec), fMerr(merr), fChisq(0.)
  {}

  /// Destructor.
  template <int N>
  KHit<N>::~KHit()
  {}

  /// Prediction method.
  ///
  /// Arguments;
  ///
  /// tre  - Track prediction.
  /// prop - Propagator.
  /// ref  - Reference track.
  ///
  template <int N>
  bool KHit<N>::predict(const KETrack& tre, const Propagator& prop, const KTrack* ref) const
  {
    // Update the prediction surface to be the track surface.

    fPredSurf = tre.getSurface();
    fPredDist = 0.;

    // Default result.

    bool ok = false;

    // Update prediction vector, error matrox, and H-matrix.

    // First test whether the prediction surface matches the
    // measurement surface.

    if (getMeasSurface()->isEqual(*fPredSurf)) {

      // Prediction and measurement surfaces agree.
      //Just call subpredict method (don't do propagation).

      ok = subpredict(tre, fPvec, fPerr, fH);
    }
    else {

      // Track surface does not match the prediction surface, so an
      // internal propagation will be required.  Throw an exception if
      // no propagator was provided.

      // First make a copy of the original KETrack.

      KETrack treprop(tre);

      // Also make a copy of the reference track (if specified).

      KTrack refprop;
      KTrack* prefprop = 0;
      if (ref != 0) {
        refprop = *ref;
        prefprop = &refprop;
      }

      // Make a no-noise, no-dE/dx propagation to the measurement
      // surface.  But do calculate the propagation matrix, which we
      // will use to update the H-matrix calculated in the derived
      // class.

      TrackMatrix prop_matrix;
      std::optional<double> dist = prop.err_prop(
        treprop, getMeasSurface(), Propagator::UNKNOWN, false, prefprop, &prop_matrix);
      if (dist) {

        // Update prediction distance.

        fPredDist = *dist;

        // Now we are ready to calculate the prediction on the
        // measurement surface.

        typename KHMatrix<N>::type hmatrix;
        ok = subpredict(treprop, fPvec, fPerr, hmatrix);
        if (ok) {

          // Use the propagation matrix to transform the H-matrix back
          // to the prediction surface.

          fH = prod(hmatrix, prop_matrix);
        }
      }
    }
    if (ok) {

      // Update residual

      fRvec = fMvec - fPvec;
      fRerr = fMerr + fPerr;
      fRinv = fRerr;
      ok = syminvert(fRinv);
      if (ok) {

        // Calculate incremental chisquare.

        // (workaround: if we use the copy constructor, gcc emits a spurious warning)
        //	typename KVector<N>::type rtemp = prod(fRinv, fRvec);
        fChisq = inner_prod(fRvec, prod(fRinv, fRvec));
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
  template <int N>
  void KHit<N>::update(KETrack& tre) const
  {
    // Make sure that the track surface and the prediction surface are the same.
    // Throw an exception if they are not.

    if (!getPredSurface()->isEqual(*tre.getSurface()))
      throw cet::exception("KHit") << "Track surface not the same as prediction surface.\n";

    const TrackVector& tvec = tre.getVector();
    const TrackError& terr = tre.getError();
    TrackVector::size_type size = tvec.size();

    // Calculate gain matrix.

    typename KGMatrix<N>::type temp(size, N);
    typename KGMatrix<N>::type gain(size, N);
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
    typename KHMatrix<N>::type errtemp3 = prod(fMerr, trans(gain));
    TrackMatrix errtemp4 = prod(gain, errtemp3);
    TrackError errtemp4s = ublas::symmetric_adaptor<TrackMatrix>(errtemp4);
    TrackError newerr = errtemp2s + errtemp4s;

    // Update track.

    tre.setVector(newvec);
    tre.setError(newerr);
  }

  /// Printout
  template <int N>
  std::ostream& KHit<N>::Print(std::ostream& out, bool doTitle) const
  {
    if (doTitle) out << "KHit<" << N << ">:\n";

    // Print base class.

    KHitBase::Print(out, false);

    // Print measurement vector.

    out << "  Measurement vector:\n"
        << "  [";
    for (unsigned int i = 0; i < fMvec.size(); ++i) {
      if (i != 0) out << ", ";
      out << fMvec(i);
    }
    out << "]\n";

    // Print diagonal measurement  errors.

    out << "  Diagonal measurement errors:\n"
        << "  [";
    for (unsigned int i = 0; i < fMerr.size1(); ++i) {
      if (i != 0) out << ", ";
      double err = fMerr(i, i);
      err = (err >= 0. ? std::sqrt(err) : -std::sqrt(-err));
      out << err;
    }
    out << "]\n";

    // Print measurement correlations.

    if (fMerr.size1() > 1) {
      out << "  Measurement correlation matrix:";
      for (unsigned int i = 0; i < fMerr.size1(); ++i) {
        if (i == 0)
          out << "\n  [";
        else
          out << "\n   ";
        for (unsigned int j = 0; j <= i; ++j) {
          if (j != 0) out << ", ";
          if (i == j)
            out << 1.;
          else {
            double eiijj = fMerr(i, i) * fMerr(j, j);
            double eij = fMerr(i, j);
            if (eiijj != 0.)
              eij /= std::sqrt(std::abs(eiijj));
            else
              eij = 0.;
            out << eij;
          }
        }
      }
      out << "]\n";
    }

    // Print prediction vector.

    out << "  Prediction vector:\n"
        << "  [";
    for (unsigned int i = 0; i < fPvec.size(); ++i) {
      if (i != 0) out << ", ";
      out << fPvec(i);
    }
    out << "]\n";

    // Print diagonal prediction  errors.

    out << "  Diagonal prediction errors:\n"
        << "  [";
    for (unsigned int i = 0; i < fPerr.size1(); ++i) {
      if (i != 0) out << ", ";
      double err = fPerr(i, i);
      err = (err >= 0. ? std::sqrt(err) : -std::sqrt(-err));
      out << err;
    }
    out << "]\n";

    // Print prediction correlations.

    if (fPerr.size1() > 1) {
      out << "  Prediction correlation matrix:";
      for (unsigned int i = 0; i < fPerr.size1(); ++i) {
        if (i == 0)
          out << "\n  [";
        else
          out << "\n   ";
        for (unsigned int j = 0; j <= i; ++j) {
          if (j != 0) out << ", ";
          if (i == j)
            out << 1.;
          else {
            double eiijj = fPerr(i, i) * fPerr(j, j);
            double eij = fPerr(i, j);
            if (eiijj != 0.)
              eij /= std::sqrt(std::abs(eiijj));
            else
              eij = 0.;
            out << eij;
          }
        }
      }
      out << "]\n";
    }

    // Print residual vector.

    out << "  Residual vector:\n"
        << "  [";
    for (unsigned int i = 0; i < fRvec.size(); ++i) {
      if (i != 0) out << ", ";
      out << fRvec(i);
    }
    out << "]\n";

    // Print diagonal residual  errors.

    out << "  Diagonal residual errors:\n"
        << "  [";
    for (unsigned int i = 0; i < fRerr.size1(); ++i) {
      if (i != 0) out << ", ";
      double err = fRerr(i, i);
      err = (err >= 0. ? std::sqrt(err) : -std::sqrt(-err));
      out << err;
    }
    out << "]\n";

    // Print residual correlations.

    if (fRerr.size1() > 1) {
      out << "  Residual correlation matrix:";
      for (unsigned int i = 0; i < fRerr.size1(); ++i) {
        if (i == 0)
          out << "\n  [";
        else
          out << "\n   ";
        for (unsigned int j = 0; j <= i; ++j) {
          if (j != 0) out << ", ";
          if (i == j)
            out << 1.;
          else {
            double eiijj = fRerr(i, i) * fRerr(j, j);
            double eij = fRerr(i, j);
            if (eiijj != 0.)
              eij /= std::sqrt(std::abs(eiijj));
            else
              eij = 0.;
            out << eij;
          }
        }
      }
      out << "]\n";
    }

    // Print incremental chisquare.

    out << "  Incremental chisquare = " << fChisq << "\n";

    // Done.

    return out;
  }
}

#endif
