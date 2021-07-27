////////////////////////////////////////////////////////////////////////
// \file DetectorPropertiesServiceArgoNeuT.h
//
// \brief service to contain information about detector electronics, etc
//
// \author brebel@fnal.gov
//
// From the original DetectorProperties.h ; this one preserves the dependency on
// DatabaseUtil service and the ability to read information from a database
// with direct DB connection.
// For new experiments, an indirect connection should be used instead.
//
// PLEASE DO NOT take this as a model to develop a service:
// this is just a backward-compatible hack.
//
////////////////////////////////////////////////////////////////////////
#ifndef UTIL_DETECTORPROPERTIESSERVICEARGONEUT_H
#define UTIL_DETECTORPROPERTIESSERVICEARGONEUT_H

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/DetectorPropertiesArgoNeuT.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataalg/DetectorInfo/ElecClock.h"

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace detinfo {
  class DetectorClocks;
}

/// General LArSoft Utilities
namespace util {
  class DetectorPropertiesServiceArgoNeuT : public detinfo::DetectorPropertiesService {
  public:
    DetectorPropertiesServiceArgoNeuT(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

  private:
    detinfo::DetectorPropertiesData
    getDataForJob(detinfo::DetectorClocksData const& clockData) const override
    {
      return fDetProp.DataFor(clockData);
    }

    detinfo::DetectorPropertiesData
    getDataFor(art::Event const&, detinfo::DetectorClocksData const& clockData) const override
    {
      return fDetProp.DataFor(clockData);
    }

    void postOpenFile(std::string const& filename);

    static bool isDetectorPropertiesServiceArgoNeuT(const fhicl::ParameterSet& ps);

    util::DetectorPropertiesArgoNeuT fDetProp;

    fhicl::ParameterSet fPS;        ///< Original parameter set.
    bool fInheritNumberTimeSamples; ///< Flag saying whether to inherit
                                    ///< NumberTimeSamples
  };                                // class DetectorPropertiesServiceArgoNeuT
} // namespace util
DECLARE_ART_SERVICE_INTERFACE_IMPL(util::DetectorPropertiesServiceArgoNeuT,
                                   detinfo::DetectorPropertiesService,
                                   SHARED)
#endif // UTIL_DETECTORPROPERTIESSERVICEARGONEUT_H
