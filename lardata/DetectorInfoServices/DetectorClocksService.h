////////////////////////////////////////////////////////////////////////
// DetectorClocksService.h
//
// Pure virtual service interface for DetectorClocks functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORCLOCKSSERVICE_H
#define DETECTORCLOCKSSERVICE_H

#include "art/Framework/Principal/fwd.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"

namespace detinfo {
  class DetectorClocksService {
  public:
    using provider_type = detinfo::DetectorClocks;

    virtual ~DetectorClocksService() = default;

    virtual DetectorClocksData DataForJob() const = 0;
    virtual DetectorClocksData DataFor(art::Event const& e) const = 0;
  };
}

DECLARE_ART_SERVICE_INTERFACE(detinfo::DetectorClocksService, SHARED)

#endif // DETECTORCLOCKSSERVICE_H
