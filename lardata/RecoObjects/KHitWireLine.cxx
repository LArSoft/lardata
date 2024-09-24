///////////////////////////////////////////////////////////////////////
///
/// \file   KHitWireLine.cxx
///
/// \brief  Kalman filter wire-time measurement on a SurfWireLine surface.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/KHitWireLine.h"
#include "cetlib_except/exception.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoObjects/SurfWireLine.h"

namespace trkf {

  /// Constructor.
  ///
  /// Arguments:
  ///
  /// hit   - Hit.
  /// psurf - Measurement surface (can be null).
  ///
  /// The measurement surface is only a suggestion.  It is allowed to
  /// be specified to allow measurements to whare surfaces to save
  /// memory.
  ///
  KHitWireLine::KHitWireLine(const detinfo::DetectorPropertiesData& detProp,
                             const art::Ptr<recob::Hit>& hit,
                             const std::shared_ptr<const Surface>& psurf)
    : KHit(psurf), fHit(hit)
  {
    // Extract wire id.
    geo::WireID wireid = hit->WireID();

    // Extract time information from hit.

    double t = hit->PeakTime();
    double terr = hit->SigmaPeakTime();

    // Don't let the time error be less than 1./sqrt(12.) ticks.
    // This should be removed when hit errors are fixed.

    if (terr < 1. / std::sqrt(12.)) terr = 1. / std::sqrt(12.);

    // Calculate position and error.

    double x = detProp.ConvertTicksToX(t, wireid.Plane, wireid.TPC, wireid.Cryostat);
    double xerr = terr * detProp.GetXTicksCoefficient();

    // Check the surface (determined by wire id + drift time).  If the
    // surface pointer is null, make a new SurfWireLine surface and
    // update the base class appropriately.  Otherwise, just check
    // that the specified surface agrees with the wire id + drift time.

    if (psurf.get() == 0) {
      std::shared_ptr<const Surface> new_psurf(new SurfWireLine(wireid, x));
      setMeasSurface(new_psurf);
    }
    else {
      SurfWireLine check_surf(wireid, x);
      if (!check_surf.isEqual(*psurf))
        throw cet::exception("KHitWireLine") << "Measurement surface doesn't match hit.\n";
    }

    setMeasPlane(wireid.Plane);

    // Update measurement vector and error matrix.

    trkf::KVector<1>::type mvec(1, 0.);
    setMeasVector(mvec);

    trkf::KSymMatrix<1>::type merr(1);
    merr(0, 0) = xerr * xerr;
    setMeasError(merr);

    // Set the unique id from a combination of the channel number and the time.

    fID = (hit->Channel() % 200000) * 10000 + (int(std::abs(t)) % 10000);
  }

  /// Constructor.
  ///
  /// Arguments:
  ///
  /// wireid  - Wire id.
  /// x       - X coordinate.
  /// xerr    - X error.
  ///
  KHitWireLine::KHitWireLine(const geo::WireID& wireid, double x, double xerr)
    : KHit(std::shared_ptr<const Surface>(new SurfWireLine(wireid, x)))
  {
    // Get plane number.

    setMeasPlane(wireid.Plane);

    // Update measurement vector and error matrix.
    // The measured value (aka impact parameter) is always zero.

    trkf::KVector<1>::type mvec(1, 0.);
    setMeasVector(mvec);

    trkf::KSymMatrix<1>::type merr(1);
    merr(0, 0) = xerr * xerr;
    setMeasError(merr);
  }

  bool KHitWireLine::subpredict(const KETrack& tre,
                                KVector<1>::type& pvec,
                                KSymMatrix<1>::type& perr,
                                KHMatrix<1>::type& hmatrix) const
  {
    // Make sure that the track surface and the measurement surface are the
    // same. Throw an exception if they are not.

    if (!getMeasSurface()->isEqual(*tre.getSurface()))
      throw cet::exception("KHitWireLine")
        << "Track surface not the same as measurement surface.\n";

    // Prediction is the signed impact parameter (parameter 0).

    int size = tre.getVector().size();
    pvec.resize(1, /* preserve */ false);
    pvec.clear();
    pvec(0) = tre.getVector()(0);

    perr.resize(1, /* preserve */ false);
    perr.clear();
    perr(0, 0) = tre.getError()(0, 0);

    // Update prediction error to include contribution from track slope.

    double pitch = art::ServiceHandle<geo::WireReadout>()->Get().Plane({0, 0, 0}).WirePitch();
    double phi = tre.getVector()(2);
    double cosphi = std::cos(phi);
    double slopevar = pitch * pitch * cosphi * cosphi / 12.;
    perr(0, 0) += slopevar;

    // Hmatrix - dr/dr = 1., all others are zero.

    hmatrix.resize(1, size, /* preserve */ false);
    hmatrix.clear();
    hmatrix(0, 0) = 1.;

    return true;
  }
} // end namespace trkf
