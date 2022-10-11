#ifndef TRACKSTATEPROPAGATOR_H
#define TRACKSTATEPROPAGATOR_H

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "lardata/RecoObjects/TrackState.h"
#include "lardataobj/RecoBase/TrackingPlane.h"
#include "lardataobj/RecoBase/TrackingTypes.h"

#include <utility>

namespace detinfo {
  class DetectorPropertiesData;
  class LArProperties;
}

namespace trkf {

  /// \class TrackStatePropagator
  ///
  /// \brief Class for propagation of a trkf::TrackState to a recob::tracking::Plane
  ///
  /// \author  G. Cerati (FNAL, MicroBooNE)
  /// \date    2017
  /// \version 1.0
  ///
  /// This class holds the functionalities needed to propagate a trkf::TrackState to a recob::tracking::Plane.
  /// While the core physics is mainly duplicated from trkf::Propagator and its derived classes (kudos to H. Greenlee),
  /// the code and the interface are optimized for usage with classes based on SMatrix (e.g. TrackState) and
  /// for the needs of TrackKalmanFitter.
  ///
  /// While the propagated position can be directly computed, accounting for the material effects
  /// in the covariance matrix requires an iterative procedure in case of long propagations distances.
  ///
  /// For configuration options see TrackStatePropagator#Config
  ///

  class TrackStatePropagator {
  public:
    using Plane = recob::tracking::Plane;
    using Point_t = recob::tracking::Point_t;
    using Vector_t = recob::tracking::Vector_t;
    using SVector6 = recob::tracking::SVector6;

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<double> minStep{Name("minStep"),
                                  Comment("Minimum propagation step length guaranteed."),
                                  1.0};
      fhicl::Atom<double> maxElossFrac{
        Name("maxElossFrac"),
        Comment("Maximum propagation step length based on fraction of energy loss."),
        0.1};
      fhicl::Atom<int> maxNit{
        Name("maxNit"),
        Comment("Maximum number of iterations when applying material effects."),
        10};
      fhicl::Atom<double> tcut{Name("tcut"), Comment("Maximum delta ray energy for dE/dx."), 10.};
      fhicl::Atom<double> wrongDirDistTolerance{
        Name("wrongDirDistTolerance"),
        Comment("Allowed propagation distance in the wrong direction."),
        0.01};
      fhicl::Atom<bool> propPinvErr{
        Name("propPinvErr"),
        Comment("Propagate error on 1/p or not (in order to avoid infs, it should be set to false "
                "when 1/p not updated)."),
        false};
    };
    using Parameters = fhicl::Table<Config>;

    /// Propagation direction enum.
    enum PropDirection { FORWARD = 0, BACKWARD = 1, UNKNOWN = 2 };

    /// Constructor from parameter values.
    TrackStatePropagator(double minStep,
                         double maxElossFrac,
                         int maxNit,
                         double tcut,
                         double wrongDirDistTolerance,
                         bool propPinvErr);

    /// Constructor from Parameters (fhicl::Table<Config>).
    explicit TrackStatePropagator(Parameters const& p)
      : TrackStatePropagator(p().minStep(),
                             p().maxElossFrac(),
                             p().maxNit(),
                             p().tcut(),
                             p().wrongDirDistTolerance(),
                             p().propPinvErr())
    {}

    /// Main function for propagation of a TrackState to a Plane
    TrackState propagateToPlane(bool& success,
                                const detinfo::DetectorPropertiesData& detProp,
                                const TrackState& origin,
                                const Plane& target,
                                bool dodedx,
                                bool domcs,
                                PropDirection dir = FORWARD) const;

    /// Rotation of a TrackState to a Plane (zero distance propagation)
    TrackState rotateToPlane(bool& success, const TrackState& origin, const Plane& target) const
    {
      double dw2dw1 = 0.;
      return rotateToPlane(success, origin, target, dw2dw1);
    }

    /// Quick accesss to the propagated position given a distance
    Point_t propagatedPosByDistance(const Point_t& origpos,
                                    const Vector_t& origdir,
                                    double distance) const
    {
      return origpos + distance * origdir;
    }

    //@{
    /// Distance of a TrackState (Point and Vector) to a Plane, along the TrackState direction
    double distanceToPlane(bool& success,
                           const Point_t& origpos,
                           const Vector_t& origdir,
                           const Plane& target) const;
    double distanceToPlane(bool& success, const TrackState& origin, const Plane& target) const
    {
      return distanceToPlane(success, origin.position(), origin.momentum().Unit(), target);
    }
    //@}

    //@{
    /// Distance of a TrackState (Point) to a Plane along the direction orthogonal to the Plane
    double perpDistanceToPlane(bool& success, const Point_t& origpos, const Plane& target) const;
    double perpDistanceToPlane(bool& success, const TrackState& origin, const Plane& target) const
    {
      return perpDistanceToPlane(success, origin.position(), target);
    }
    //@}

    //@{
    /// Return both direction types in one go
    std::pair<double, double> distancePairToPlane(bool& success,
                                                  const Point_t& origpos,
                                                  const Vector_t& origdir,
                                                  const Plane& target) const;
    std::pair<double, double> distancePairToPlane(bool& success,
                                                  const TrackState& origin,
                                                  const Plane& target) const
    {
      return distancePairToPlane(success, origin.position(), origin.momentum().Unit(), target);
    }
    //@}

    /// Apply energy loss.
    void apply_dedx(double& pinv,
                    detinfo::DetectorPropertiesData const& detProp,
                    double dedx,
                    double e1,
                    double mass,
                    double s,
                    double& deriv) const;

    /// Apply multiple coulomb scattering.
    bool apply_mcs(detinfo::DetectorPropertiesData const& detProp,
                   double dudw,
                   double dvdw,
                   double pinv,
                   double mass,
                   double s,
                   double range,
                   double p,
                   double e2,
                   bool flipSign,
                   SMatrixSym55& noise_matrix) const;

    /// get Tcut parameter used in DetectorPropertiesService Eloss method
    double getTcut() const { return fTcut; }

  private:
    /// Rotation of a TrackState to a Plane (zero distance propagation), keeping track of dw2dw1 (needed by mcs)
    TrackState rotateToPlane(bool& success,
                             const TrackState& origin,
                             const Plane& target,
                             double& dw2dw1) const;

    double fMinStep;      ///< Minimum propagation step length guaranteed.
    double fMaxElossFrac; ///< Maximum propagation step length based on fraction of energy loss.
    int fMaxNit;          ///< Maximum number of iterations.
    double fTcut;         ///< Maximum delta ray energy for dE/dx.
    double fWrongDirDistTolerance; ///< Allowed propagation distance in the wrong direction.
    bool
      fPropPinvErr; ///< Propagate error on 1/p or not (in order to avoid infs, it should be set to false when 1/p not updated)
    const detinfo::LArProperties* larprop;
  };
}

#endif
