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
  ///
  SurfWireX::SurfWireX(const geo::WireID& wireid)
  {
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();
    geo::WireGeo const& wgeom = wireReadoutGeom.Wire(wireid);

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
