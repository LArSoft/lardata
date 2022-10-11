////////////////////////////////////////////////////////////////////////
///
/// \file   KHitMulti.h
///
/// \brief  Compound Kalman Filter measurement.
///
/// \author H. Greenlee
///
/// This class allows a collection of Kalman Filter measurements to
/// function as a single measurement.  This class inherits from
/// KHitBase and fulfills the KHitBase interface.  The following
/// attributes are owned by base class KHitBase.
///
/// 1.  Measurement surface.
/// 2.  Prediction surface.
///
/// This class adds the following attributes of its own.
///
/// 1.  A collection of underlying measurements.
/// 2.  Measurement vector.
/// 3.  Measurement error matrix.
/// 4.  Prediction vector.
/// 5.  Prediction error matrix.
/// 6.  Residual vector.
/// 7.  Residual error matrix.
/// 8.  Inverse of residual error matrix.
/// 9.  Kalman H-matrix.
/// 10.  Incremental chisquare.
///
/// The only supported kind of underlying measurement is KHit<1> (more
/// types could be added if needed).  The dimension of the measurement
/// space is dynamic (potentially large) and specified at run time.
///
/// The measurement vector, measurement error matrix, prediction
/// vector, and H-matrix are simply the concatenation of the
/// corressponding quantities from the underlying measurements.  The
/// measurement error matrix is block-diagonal with no correlation
/// between the underlying measurements.  The prediction error matrix
/// is calculated from the error matrix of the track hypothesis and
/// the full H-matrix (there will be correlations between
/// measurements).  Residuals and chisquare are calculated in the
/// usual way.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHITMULTI_H
#define KHITMULTI_H

#include "lardata/RecoObjects/KHit.h"

namespace trkf {

  class KHitMulti : public KHitBase {
  public:
    /// Default constructor.
    KHitMulti();

    /// Initializing Constructor -- measurement surface only.
    KHitMulti(const std::shared_ptr<const Surface>& psurf);

    /// Destructor.
    virtual ~KHitMulti();

    // Accessors.

    /// Measurement space dimension.
    int getMeasDim() const { return fMeasDim; }

    /// Measurement collection.
    const std::vector<std::shared_ptr<const KHit<1>>>& getMeasVec() const { return fMeasVec; }

    /// Measurement vector.
    const ublas::vector<double>& getMeasVector() const { return fMvec; }

    /// Measurement error matrix.
    const ublas::symmetric_matrix<double>& getMeasError() const { return fMerr; }

    /// Prediction vector.
    const ublas::vector<double>& getPredVector() const { return fPvec; }

    /// Prediction matrix.
    const ublas::symmetric_matrix<double>& getPredError() const { return fPerr; }

    /// Residual vector.
    const ublas::vector<double>& getResVector() const { return fRvec; }

    /// Residual error matrix.
    const ublas::symmetric_matrix<double>& getResError() const { return fRerr; }

    /// Residual inv. error matrix.
    const ublas::symmetric_matrix<double>& getResInvError() const { return fRinv; }

    /// Kalman H-matrix.
    const ublas::matrix<double>& getH() const { return fH; }

    /// Incremental chisquare.
    double getChisq() const { return fChisq; }

    // Modifiers.

    /// Add a measurement of unknown type.
    void addMeas(const std::shared_ptr<const KHitBase>& pmeas);

    /// Add a one-dimensional measurement.
    void addMeas(const std::shared_ptr<const KHit<1>>& pmeas);

    // Overrides.

    /// Prediction method (return false if fail).
    bool predict(const KETrack& tre, const Propagator& prop, const KTrack* ref = 0) const;

    /// Update track method.
    void update(KETrack& tre) const;

    /// Printout
    virtual std::ostream& Print(std::ostream& out, bool doTitle = true) const;

  private:
    // Attributes.

    int fMeasDim;                                         ///< Measurement space dimension.
    std::vector<std::shared_ptr<const KHit<1>>> fMeasVec; ///< Underlying measurements.
    mutable ublas::vector<double> fMvec;                  ///< Measurement vector.
    mutable ublas::symmetric_matrix<double> fMerr;        ///< Measurement error matrix.
    mutable ublas::vector<double> fPvec;                  ///< Prediction vector.
    mutable ublas::symmetric_matrix<double> fPerr;        ///< Prediction  error matrix.
    mutable ublas::vector<double> fRvec;                  ///< Residual vector.
    mutable ublas::symmetric_matrix<double> fRerr;        ///< Residual error matrix.
    mutable ublas::symmetric_matrix<double> fRinv;        ///< Residual inverse error matrix.
    mutable ublas::matrix<double> fH;                     ///< Kalman H-matrix.
    mutable double fChisq;                                ///< Incremental chisquare.
  };
}

#endif
