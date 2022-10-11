///////////////////////////////////////////////////////////////////////
///
/// \file   Propagator.cxx
///
/// \brief  Base class for Kalman filter propagator.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/Propagator.h"
#include "cetlib_except/exception.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

namespace trkf {

  /// Constructor.
  ///
  /// Arguments:
  ///
  /// tcut   - Maximum delta ray energy.
  /// doDedx - dE/dx enable flag.
  ///
  Propagator::Propagator(detinfo::DetectorPropertiesData const& detProp,
                         double tcut,
                         bool doDedx,
                         const std::shared_ptr<const Interactor>& interactor)
    : fDetProp{detProp}, fTcut(tcut), fDoDedx(doDedx), fInteractor(interactor)
  {}

  /// Destructor.
  Propagator::~Propagator() = default;

  /// Propagate without error (long distance).
  ///
  /// Arguments:
  ///
  /// trk          - Track to propagate.
  /// psurf        - Destination surface.
  /// dir          - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  /// doDedx       - dE/dx enable/disable flag.
  /// prop_matrix  - Return propagation matrix if not null.
  /// noise_matrix - Return noise matrix if not null.
  ///
  /// Returned value: Propagation distance & success flag.
  ///
  /// This method calls virtual method short_vec_prop in steps of some
  /// maximum size.
  ///
  std::optional<double> Propagator::vec_prop(KTrack& trk,
                                             const std::shared_ptr<const Surface>& psurf,
                                             PropDirection dir,
                                             bool doDedx,
                                             TrackMatrix* prop_matrix,
                                             TrackError* noise_matrix) const
  {
    std::optional<double> result{std::nullopt};

    // Get the inverse momentum (assumed to be track parameter four).

    double pinv = trk.getVector()(4);

    // If dE/dx is not requested, or if inverse momentum is zero, then
    // it is safe to propagate in one step.  In this case, just pass
    // the call to short_vec_prop with unlimited distance.

    bool dedx = getDoDedx() && doDedx;
    if (!dedx || pinv == 0.)
      result = short_vec_prop(trk, psurf, dir, dedx, prop_matrix, noise_matrix);

    else {

      // dE/dx is requested.  In this case we limit the maximum
      // propagation distance such that the kinetic energy of the
      // particle should not change by more thatn 10%.

      // Initialize propagation matrix to unit matrix (if specified).

      int nvec = trk.getVector().size();
      if (prop_matrix) *prop_matrix = ublas::identity_matrix<TrackVector::value_type>(nvec);

      // Initialize noise matrix to zero matrix (if specified).

      if (noise_matrix) {
        noise_matrix->resize(nvec, nvec, false);
        noise_matrix->clear();
      }

      // Remember the starting track.

      KTrack trk0(trk);

      // Make pointer variables pointing to local versions of the
      // propagation and noise matrices, or null if not specified.

      TrackMatrix local_prop_matrix;
      TrackMatrix* plocal_prop_matrix = (prop_matrix == 0 ? 0 : &local_prop_matrix);
      TrackError local_noise_matrix;
      TrackError* plocal_noise_matrix = (noise_matrix == 0 ? 0 : &local_noise_matrix);

      // Cumulative propagation distance.

      double s = 0.;

      // Begin stepping loop.
      // We put a maximum iteration count to prevent infinite loops caused by
      // floating point pathologies.  The iteration count is large enough to reach
      // any point in the tpc using the minimum step size (for a reasonable tpc).

      bool done = false;
      int nitmax = 10000; // Maximum number of iterations.
      int nit = 0;        // Iteration count.
      while (!done) {

        // If the iteration count exceeds the maximum, return failure.

        ++nit;
        if (nit > nitmax) {
          trk = trk0;
          result = std::nullopt;
          return result;
        }

        // Estimate maximum step distance according to the above
        // stated principle.

        pinv = trk.getVector()(4);
        double mass = trk.Mass();
        double p = 1. / std::abs(pinv);
        double e = std::hypot(p, mass);
        double t = p * p / (e + mass);
        double dedx = 0.001 * fDetProp.Eloss(p, mass, fTcut);
        double smax = 0.1 * t / dedx;
        if (smax <= 0.)
          throw cet::exception("Propagator") << __func__ << ": maximum step " << smax << "\n";

        // Always allow a step of at least 0.3 cm (about one wire spacing).

        if (smax < 0.3) smax = 0.3;

        // First do a test propagation (without dE/dx and errors) to
        // find the distance to the destination surface.

        KTrack trktest(trk);
        std::optional<double> dist = short_vec_prop(trktest, psurf, dir, false, 0, 0);

        // If the test propagation failed, return failure.

        if (!dist) {
          trk = trk0;
          return dist;
        }

        // Generate destionation surface for this step (either final
        // destination, or some intermediate surface).

        std::shared_ptr<const Surface> pstep;
        if (std::abs(*dist) <= smax) {
          done = true;
          pstep = psurf;
        }
        else {

          // Generate intermediate surface.
          // First get point where track will intersect intermediate surface.

          double xyz0[3]; // Starting point.
          trk.getPosition(xyz0);
          double xyz1[3]; // Destination point.
          trktest.getPosition(xyz1);
          double frac = smax / std::abs(*dist);
          double xyz[3]; // Intermediate point.
          xyz[0] = xyz0[0] + frac * (xyz1[0] - xyz0[0]);
          xyz[1] = xyz0[1] + frac * (xyz1[1] - xyz0[1]);
          xyz[2] = xyz0[2] + frac * (xyz1[2] - xyz0[2]);

          // Choose orientation of intermediate surface perpendicular
          // to track.

          double mom[3];
          trk.getMomentum(mom);

          // Make intermediate surface object.

          pstep = std::shared_ptr<const Surface>(
            new SurfXYZPlane(xyz[0], xyz[1], xyz[2], mom[0], mom[1], mom[2]));
        }

        // Do the actual step propagation.

        dist = short_vec_prop(trk, pstep, dir, doDedx, plocal_prop_matrix, plocal_noise_matrix);

        // If the step propagation failed, return failure.

        if (!dist) {
          trk = trk0;
          return dist;
        }

        // Update cumulative propagation distance.

        s += *dist;

        // Update cumulative propagation matrix (left-multiply).

        if (prop_matrix != 0) {
          TrackMatrix temp = prod(*plocal_prop_matrix, *prop_matrix);
          *prop_matrix = temp;
        }

        // Update cumulative noise matrix.

        if (noise_matrix != 0) {
          TrackMatrix temp = prod(*noise_matrix, trans(*plocal_prop_matrix));
          TrackMatrix temp2 = prod(*plocal_prop_matrix, temp);
          *noise_matrix = ublas::symmetric_adaptor<TrackMatrix>(temp2);
          *noise_matrix += *plocal_noise_matrix;
        }
      }

      // Set the final result (distance + success).

      result = std::make_optional(s);
    }

    // Done.

    return result;
  }

  /// Linearized propagate without error.
  ///
  /// Arguments:
  ///
  /// trk          - Track to propagate.
  /// psurf        - Destination surface.
  /// dir          - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  /// doDedx       - dE/dx enable/disable flag.
  /// ref          - Reference track (for linearized propagation).  Can be null.
  /// prop_matrix  - Return propagation matrix if not null.
  /// noise_matrix - Return noise matrix if not null.
  ///
  /// Returned value: Propagation distance & success flag.
  ///
  /// If the reference track is null, this method simply calls vec_prop.
  ///
  std::optional<double> Propagator::lin_prop(KTrack& trk,
                                             const std::shared_ptr<const Surface>& psurf,
                                             PropDirection dir,
                                             bool doDedx,
                                             KTrack* ref,
                                             TrackMatrix* prop_matrix,
                                             TrackError* noise_matrix) const
  {
    // Default result.

    std::optional<double> result;

    if (ref == 0)
      result = vec_prop(trk, psurf, dir, doDedx, prop_matrix, noise_matrix);
    else {

      // A reference track has been provided.

      // It is an error (throw exception) if the reference track and
      // the track to be propagted are not on the same surface.

      if (!trk.getSurface()->isEqual(*(ref->getSurface())))
        throw cet::exception("Propagator")
          << "Input track and reference track not on same surface.\n";

      // Remember the starting track and reference track.

      KTrack trk0(trk);
      KTrack ref0(*ref);

      // Propagate the reference track.  Make sure we calculate the
      // propagation matrix.

      TrackMatrix prop_temp;
      if (prop_matrix == 0) prop_matrix = &prop_temp;

      // Do the propgation.  The returned result will be the result of
      // this propagatrion.

      result = vec_prop(*ref, psurf, dir, doDedx, prop_matrix, noise_matrix);
      if (!!result) {

        // Propagation of reference track succeeded.  Update the track
        // state vector and surface of the track to be propagated.

        TrackVector diff = trk.getSurface()->getDiff(trk.getVector(), ref0.getVector());
        TrackVector newvec = ref->getVector() + prod(*prop_matrix, diff);

        // Store updated state vector and surface.

        trk.setVector(newvec);
        trk.setSurface(psurf);
        trk.setDirection(ref->getDirection());
        if (!trk.getSurface()->isEqual(*(ref->getSurface())))
          throw cet::exception("Propagator") << __func__ << ": surface mismatch";

        // Final validity check.  In case of failure, restore the track
        // and reference track to their starting values.

        if (!trk.isValid()) {
          result = std::nullopt;
          trk = trk0;
          *ref = ref0;
        }
      }
      else {

        // Propagation failed.
        // Restore the reference track to its starting value, so that we ensure
        // the reference track and the actual track remain on the same surface.

        trk = trk0;
        *ref = ref0;
      }
    }

    // Done.

    return result;
  }

  /// Propagate with error, but without noise (i.e. reversibly).
  ///
  /// Arguments:
  ///
  /// tre         - Track to propagate.
  /// psurf       - Destination surface.
  /// dir         - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  /// doDedx      - dE/dx enable/disable flag.
  /// ref         - Reference track (for linearized propagation).  Can be null.
  /// prop_matrix - Return propagation matrix if not null.
  ///
  /// Returned value: propagation distance + success flag.
  ///
  std::optional<double> Propagator::err_prop(KETrack& tre,
                                             const std::shared_ptr<const Surface>& psurf,
                                             PropDirection dir,
                                             bool doDedx,
                                             KTrack* ref,
                                             TrackMatrix* prop_matrix) const
  {
    // Propagate without error, get propagation matrix.

    TrackMatrix prop_temp;
    if (prop_matrix == 0) prop_matrix = &prop_temp;
    std::optional<double> result = lin_prop(tre, psurf, dir, doDedx, ref, prop_matrix, 0);

    // If propagation succeeded, update track error matrix.

    if (!!result) {
      TrackMatrix temp = prod(tre.getError(), trans(*prop_matrix));
      TrackMatrix temp2 = prod(*prop_matrix, temp);
      TrackError newerr = ublas::symmetric_adaptor<TrackMatrix>(temp2);
      tre.setError(newerr);
    }

    // Done.

    return result;
  }

  /// Propagate with error and noise.
  ///
  /// Arguments:
  ///
  /// tre    - Track to propagate.
  /// psurf  - Destination surface.
  /// dir    - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  /// doDedx - dE/dx enable/disable flag.
  /// ref    - Reference track (for linearized propagation).  Can be null.
  ///
  /// Returned value: propagation distance + success flag.
  ///
  std::optional<double> Propagator::noise_prop(KETrack& tre,
                                               const std::shared_ptr<const Surface>& psurf,
                                               PropDirection dir,
                                               bool doDedx,
                                               KTrack* ref) const
  {
    // Propagate without error, get propagation matrix and noise matrix.

    TrackMatrix prop_matrix;
    TrackError noise_matrix;
    std::optional<double> result =
      lin_prop(tre, psurf, dir, doDedx, ref, &prop_matrix, &noise_matrix);

    // If propagation succeeded, update track error matrix.

    if (!!result) {
      TrackMatrix temp = prod(tre.getError(), trans(prop_matrix));
      TrackMatrix temp2 = prod(prop_matrix, temp);
      TrackError newerr = ublas::symmetric_adaptor<TrackMatrix>(temp2);
      newerr += noise_matrix;
      tre.setError(newerr);
    }

    // Done.

    return result;
  }

  /// Method to calculate updated momentum due to dE/dx.
  ///
  /// Arguments:
  ///
  /// pinv  - Initial inverse momentum (units c/GeV).
  /// mass  - Particle mass (GeV/c^2).
  /// s     - Path distance.
  /// deriv - Pointer to store derivative d(pinv2)/d(pinv1) if nonzero.
  ///
  /// Returns: Final inverse momentum (pinv2) + success flag.
  ///
  /// Failure is returned in case of range out.
  ///
  /// Inverse momentum can be signed (q/p).  Returned inverse momentum
  /// has the same sign as the input.
  ///
  /// In this method, we are solving the differential equation in
  /// terms of energy.
  ///
  /// dE/dx = -f(E)
  ///
  /// where f(E) is the stopping power returned by method
  /// LArProperties::Eloss.
  ///
  /// We expect that this method will be called exclusively for short
  /// distance propagation.  The differential equation is solved using
  /// the midpoint method using a single step, which requires two
  /// evaluations of f(E).
  ///
  /// dE = -s*f(E1)
  /// E2 = E1 - s*f(E1 + 0.5*dE)
  ///
  /// The derivative is calculated assuming E2 = E1 + constant, giving
  ///
  /// d(pinv2)/d(pinv1) = pinv2^3 E2 / (pinv1^3 E1).
  ///
  ///
  std::optional<double> Propagator::dedx_prop(double pinv,
                                              double mass,
                                              double s,
                                              double* deriv) const
  {
    // For infinite initial momentum, return with success status,
    // still infinite momentum.

    if (pinv == 0.) return std::make_optional(0.);

    // Set the default return value to be uninitialized with value 0.

    std::optional<double> result{std::nullopt};

    // Calculate final energy.

    double p1 = 1. / std::abs(pinv);
    double e1 = std::hypot(p1, mass);
    double de = -0.001 * s * fDetProp.Eloss(p1, mass, fTcut);
    double emid = e1 + 0.5 * de;
    if (emid > mass) {
      double pmid = std::sqrt(emid * emid - mass * mass);
      double e2 = e1 - 0.001 * s * fDetProp.Eloss(pmid, mass, fTcut);
      if (e2 > mass) {
        double p2 = std::sqrt(e2 * e2 - mass * mass);
        double pinv2 = 1. / p2;
        if (pinv < 0.) pinv2 = -pinv2;

        // Calculation was successful, update result.

        result = std::make_optional(pinv2);

        // Also calculate derivative, if requested.

        if (deriv != 0) *deriv = pinv2 * pinv2 * pinv2 * e2 / (pinv * pinv * pinv * e1);
      }
    }

    // Done.

    return result;
  }
} // end namespace trkf
