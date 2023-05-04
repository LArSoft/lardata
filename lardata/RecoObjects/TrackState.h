#ifndef TRACKSTATE_H
#define TRACKSTATE_H

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/RecoObjects/TrackingPlaneHelper.h"
#include "lardataobj/RecoBase/TrackingPlane.h"
#include "lardataobj/RecoBase/TrackingTypes.h"

namespace trkf {

  using SVector5 = recob::tracking::SVector5;
  using SVector6 = recob::tracking::SVector6;
  using SMatrix55 = recob::tracking::SMatrix55;
  using SMatrixSym55 = recob::tracking::SMatrixSym55;
  using SMatrixSym66 = recob::tracking::SMatrixSym66;
  using Plane = recob::tracking::Plane;
  using Point_t = recob::tracking::Point_t;
  using Vector_t = recob::tracking::Vector_t;

  namespace {
    const double elmass = 0.000510998; // Electron
    const double mumass = 0.105658367; // Muon
    const double pimass = 0.13957;     // Charged pion
    const double kmass = 0.493677;     // Charged kaon
    const double pmass = 0.938272;     // Proton
  }

  /// \file  lardata/RecoObjects/TrackState.h
  /// \class HitState
  ///
  /// \brief Class for a measurement on a recob::tracking::Plane (plane defined by a wire and the drift direction).
  ///
  /// \author  G. Cerati (FNAL, MicroBooNE)
  /// \date    2017
  /// \version 1.0
  ///
  /// This class collects the measurement information from a Hit on wire.
  /// The information are the measured (1D) position, its error, and the measurement plane (defined by the wire and the drift direction)
  ///

  class HitState {
  public:
    HitState(double hitMeas,
             double hitMeasErr2,
             geo::WireID const& wireId,
             const geo::WireGeo& wgeom)
      : fHitMeas(hitMeas)
      , fHitMeasErr2(hitMeasErr2)
      , fWireId(wireId)
      , fPlane(recob::tracking::makePlane(wgeom))
    {}
    HitState(double hitMeas, double hitMeasErr2, geo::WireID&& wireId, const geo::WireGeo& wgeom)
      : fHitMeas(hitMeas)
      , fHitMeasErr2(hitMeasErr2)
      , fWireId(std::move(wireId))
      , fPlane(recob::tracking::makePlane(wgeom))
    {}
    double hitMeas() const { return fHitMeas; }
    double hitMeasErr2() const { return fHitMeasErr2; }
    const Plane& plane() const { return fPlane; }
    const geo::WireID& wireId() const { return fWireId; }
    std::ostream& dump(std::ostream& out = std::cout) const
    {
      out << "HitState with meas=" << hitMeas() << " err2=" << hitMeasErr2()
          << " plane=" << wireId().Plane << " wire=" << wireId().Wire
          << " on plane with pos=" << plane().position() << " and dir=" << plane().direction()
          << "\n";
      return out;
    }

  private:
    double fHitMeas;
    double fHitMeasErr2;
    const geo::WireID fWireId;
    Plane fPlane;
  };

  /// \file  lardata/RecoObjects/TrackState.h
  /// \class TrackState
  ///
  /// \brief Class for track parameters (and errors) defined on a recob::tracking::Plane.
  ///
  /// \author  G. Cerati (FNAL, MicroBooNE)
  /// \date    2017
  /// \version 1.0
  ///
  /// This class collects the track parameters (and errors) defined on a recob::tracking::Plane.
  /// It stores the 5d parameters and covariance, plus the global position and momentum.
  /// Given a HitState on the same plane, it provides easy access to functionalities like chi2 and residual.
  ///

  class TrackState {
  public:
    TrackState(const SVector5& trackStatePar,
               const SMatrixSym55& trackStateCov,
               const Plane& plane,
               bool trackAlongPlaneDir,
               int pid)
      : fTrackStatePar(trackStatePar), fTrackStateCov(trackStateCov), fPlane(plane), fPid(pid)
    {
      SVector6 par6d = fPlane.Local5DToGlobal6DParameters(fTrackStatePar, trackAlongPlaneDir);
      fPos = Point_t(par6d[0], par6d[1], par6d[2]);
      fMom = Point_t(par6d[3], par6d[4], par6d[5]);
    }
    //
    /// track parameters defined on the plane
    const SVector5& parameters() const { return fTrackStatePar; }
    /// track parameter covariance matrix on the plane
    const SMatrixSym55& covariance() const { return fTrackStateCov; }
    /// plane where the parameters are defined
    const Plane& plane() const { return fPlane; }
    /// position of the track
    const Point_t& position() const { return fPos; }
    /// momentum of the track
    const Vector_t& momentum() const { return fMom; }
    /// particle id hypthesis of the track
    int pID() const { return fPid; }
    /// mass hypthesis of the track
    double mass() const
    {
      if (abs(fPid) == 11) { return elmass; }
      if (abs(fPid) == 13) { return mumass; }
      if (abs(fPid) == 211) { return pimass; }
      if (abs(fPid) == 321) { return kmass; }
      if (abs(fPid) == 2212) { return pmass; }
      return util::kBogusD;
    }
    /// track parameters in global cartesian coordinates
    SVector6 parameters6D() const
    {
      return SVector6(fPos.X(), fPos.Y(), fPos.Z(), fMom.X(), fMom.Y(), fMom.Z());
    }
    /// track parameter covariance matrix in global cartesian coordinates
    SMatrixSym66 covariance6D() const
    {
      return fPlane.Local5DToGlobal6DCovariance(fTrackStateCov, true, fMom);
    }
    //
    /// is the track momentum along the plane direction?
    bool isTrackAlongPlaneDir() const { return fMom.Dot(fPlane.direction()) > 0; }
    //
    /// Printout information
    std::ostream& dump(std::ostream& out = std::cout) const
    {
      out << "TrackState with pID=" << pID() << " mass=" << mass() << "\npars=" << parameters()
          << " position=" << position() << " momentum=" << momentum() << "\ncov=\n"
          << covariance() << "\non plane with pos=" << plane().position()
          << " and dir=" << plane().direction() << " along=" << isTrackAlongPlaneDir() << "\n";
      return out;
    }
    //
    /// Residual of the TrackState with respect to a HitState. The two states must be on the same plane; it is responsibility of the user to enforce this.
    inline double residual(const HitState& hitstate) const
    {
      return hitstate.hitMeas() - fTrackStatePar(0);
    }

    /// Combined squared error of the TrackState with respect to a HitState. The two states must be on the same plane; it is responsibility of the user to enforce this.
    inline double combinedError2(const HitState& hitstate) const
    {
      return hitstate.hitMeasErr2() + fTrackStateCov(0, 0);
    }

    /// Combined error of the TrackState with respect to a HitState. The two states must be on the same plane; it is responsibility of the user to enforce this.
    inline double combinedError(const HitState& hitstate) const
    {
      return sqrt(combinedError2(hitstate));
    }

    /// Chi2 of the TrackState with respect to a HitState. The two states must be on the same plane; it is responsibility of the user to enforce this.
    inline double chi2(const HitState& hitstate) const
    {
      return residual(hitstate) * residual(hitstate) / combinedError2(hitstate);
    }

    /// Set the covariance matrix of the TrackState.
    void setCovariance(const SMatrixSym55& trackStateCov) { fTrackStateCov = trackStateCov; }

    /// Set the parameters of the TrackState; also update the global position and momentum accordingly.
    void setParameters(const SVector5& trackStatePar)
    {
      fTrackStatePar = trackStatePar;
      SVector6 par6d = fPlane.Local5DToGlobal6DParameters(trackStatePar, isTrackAlongPlaneDir());
      fPos = Point_t(par6d[0], par6d[1], par6d[2]);
      fMom = Vector_t(par6d[3], par6d[4], par6d[5]);
    }
    //
  private:
    SVector5 fTrackStatePar;     ///< track parameters defined on the plane
    SMatrixSym55 fTrackStateCov; ///< track parameter covariance matrix on the plane
    Plane fPlane;                ///< plane where the parameters are defined
    int fPid;                    ///< particle id hypthesis of the track
    Point_t fPos;                ///< position of the track (cached)
    Vector_t fMom;               ///< momentum of the track (cached)
  };

}

#endif
