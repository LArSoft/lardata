////////////////////////////////////////////////////////////////////////
///
/// \file   SurfWireLine.h
///
/// \brief  Linear surface defined by wire id and x coordinate.
///
/// \author H. Greenlee
///
/// This class derives from SurfYZLine.  This class does not add any
/// new members, but has a constructor that allows construction from
/// a wire id and x coordinate.
///
////////////////////////////////////////////////////////////////////////

#ifndef SURFWIRELINE_H
#define SURFWIRELINE_H

#include "lardata/RecoObjects/SurfYZLine.h"

namespace geo {
  struct WireID;
}

namespace trkf {

  class SurfWireLine : public SurfYZLine {
  public:
    /// Constructor.
    SurfWireLine(const geo::WireID& wireid, double x);

    /// Destructor.
    virtual ~SurfWireLine();
  };
}

#endif
