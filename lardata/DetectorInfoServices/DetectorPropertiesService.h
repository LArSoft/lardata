////////////////////////////////////////////////////////////////////////
// DetectorPropertiesService.h
//
// Pure virtual service interface for DetectorProperties functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORPROPERTIESSERVICE_H
#define DETECTORPROPERTIESSERVICE_H

#include "art/Framework/Principal/fwd.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

namespace detinfo {
  class DetectorClocksData;

  class DetectorPropertiesService {
  public:
    using provider_type = detinfo::DetectorProperties;

    virtual ~DetectorPropertiesService() = default;
    DetectorPropertiesData
    DataForJob() const
    {
      auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
      return DataForJob(clockData);
    }
    DetectorPropertiesData
    DataForJob(detinfo::DetectorClocksData const& clockData) const
    {
      return getDataForJob(clockData);
    }
    DetectorPropertiesData
    DataFor(art::Event const& e) const
    {
      auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
      return DataFor(e, clockData);
    }
    DetectorPropertiesData
    DataFor(art::Event const& e, detinfo::DetectorClocksData const& clockData) const
    {
      return getDataFor(e, clockData);
    }

  private:
    virtual DetectorPropertiesData getDataForJob(
      detinfo::DetectorClocksData const& clockData) const = 0;
    virtual DetectorPropertiesData getDataFor(
      art::Event const& e,
      detinfo::DetectorClocksData const& clockData) const = 0;
  };
} //namespace detinfo

DECLARE_ART_SERVICE_INTERFACE(detinfo::DetectorPropertiesService, SHARED)

#endif // DETECTORPROPERTIESSERVICE_H
