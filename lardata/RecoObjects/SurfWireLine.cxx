///////////////////////////////////////////////////////////////////////
///
/// \file   SurfWireLine.cxx
///
/// \brief  Linear surface defined by wire id and x coordinate.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/SurfWireLine.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TMath.h"

namespace trkf {

  /// Constructor.
  ///
  /// Arguments:
  ///
  /// wireid - Wire id.
  /// x      - X coordinate.
  ///
  SurfWireLine::SurfWireLine(const geo::WireID& wireid, double x)
  {
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();
    geo::WireGeo const& wgeom = wireReadoutGeom.Wire(wireid);

    // Get wire center and angle from the wire geometry.
    // Put local origin at center of wire.

    auto const xyz = wgeom.GetCenter();
    double phi = TMath::PiOver2() - wgeom.ThetaZ();

    // Update base class.

    *static_cast<SurfYZLine*>(this) = SurfYZLine(x, xyz.Y(), xyz.Z(), phi);
  }

  /// Destructor.
  SurfWireLine::~SurfWireLine() = default;

} // end namespace trkf
