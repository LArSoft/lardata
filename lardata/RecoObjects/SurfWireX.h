////////////////////////////////////////////////////////////////////////
///
/// \file   SurfWireX.h
///
/// \brief  Planar surface defined by wire id and x-axis.
///
/// \author H. Greenlee
///
/// This class derives from SurfYZPlane.  This class does not add any
/// new members, but has a constructor that allows construction from
/// a wire id.
///
////////////////////////////////////////////////////////////////////////

#ifndef SURFWIREX_H
#define SURFWIREX_H

#include "lardata/RecoObjects/SurfYZPlane.h"

namespace geo {
  struct WireID;
}

namespace trkf {

  class SurfWireX : public SurfYZPlane {
  public:
    /// Constructor.
    SurfWireX(const geo::WireID& wireid);

    /// Destructor.
    virtual ~SurfWireX();
  };
}

#endif
