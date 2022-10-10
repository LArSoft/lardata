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
#include "TMath.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/WireGeo.h"

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
    // Get geometry service.

    art::ServiceHandle<geo::Geometry const> geom;

    // Get wire geometry.

    geo::WireGeo const& wgeom = geom->WireIDToWireGeo(wireid);

    // Get wire center and angle from the wire geometry.
    // Put local origin at center of wire.

    auto const xyz = wgeom.GetCenter();
    double phi = TMath::PiOver2() - wgeom.ThetaZ();

    // Update base class.

    *static_cast<SurfYZLine*>(this) = SurfYZLine(x, xyz[1], xyz[2], phi);
  }

  /// Destructor.
  SurfWireLine::~SurfWireLine() = default;

} // end namespace trkf
