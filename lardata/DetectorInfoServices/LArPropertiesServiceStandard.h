////////////////////////////////////////////////////////////////////////
// LArPropertiesServiceStandard.h
//
// Service interface for Utility LAr functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef LARPROPERTIESSERVICESTANDARD_H
#define LARPROPERTIESSERVICESTANDARD_H

#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardataalg/DetectorInfo/LArPropertiesStandard.h"

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace art {
  class Run;
}

namespace detinfo {
  class LArPropertiesServiceStandard : public LArPropertiesService {
  public:
    using Parameters = art::ServiceTable<detinfo::LArPropertiesStandard::ConfigurationParameters_t>;
    LArPropertiesServiceStandard(Parameters const& params, art::ActivityRegistry& reg);

  private:
    provider_type const* provider() const override { return &fProp; }

    void preBeginRun(art::Run const& run);

    detinfo::LArPropertiesStandard fProp;
  };
} //namespace detinfo

DECLARE_ART_SERVICE_INTERFACE_IMPL(detinfo::LArPropertiesServiceStandard,
                                   detinfo::LArPropertiesService,
                                   SHARED)

#endif // LARPROPERTIESSERVICESTANDARD_H
