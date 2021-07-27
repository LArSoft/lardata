////////////////////////////////////////////////////////////////////////
// DetectorPropertiesService.h
//
// Pure virtual service interface for DetectorProperties functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef LARPROPERTIESSERVICE_H
#define LARPROPERTIESSERVICE_H

#include "lardataalg/DetectorInfo/LArProperties.h"

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

namespace detinfo {
  class LArPropertiesService {
  public:
    using provider_type = detinfo::LArProperties;

    virtual ~LArPropertiesService() = default;
    virtual provider_type const* provider() const = 0;
  };
} //namespace detinfo

DECLARE_ART_SERVICE_INTERFACE(detinfo::LArPropertiesService, SHARED)

#endif // LARPROPERTIESSERVICE_H
