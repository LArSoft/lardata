///////////////////////////////////////////////////////////////////////
///
/// \file   SurfWireX.cxx
///
/// \brief  Planar surface defined by wire id and x-axis.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "lardata/RecoObjects/SurfWireX.h"
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
  ///
  SurfWireX::SurfWireX(const geo::WireID& wireid)
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

    *static_cast<SurfYZPlane*>(this) = SurfYZPlane(0., xyz.Y(), xyz.Z(), phi);
  }

  /// Destructor.
  SurfWireX::~SurfWireX() = default;

} // end namespace trkf
