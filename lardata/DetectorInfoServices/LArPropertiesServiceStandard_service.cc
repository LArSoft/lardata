#include "lardata/DetectorInfoServices/LArPropertiesServiceStandard.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

//------------------------------------------------
/// \todo these values should eventually come from a database
//-----------------------------------------------
detinfo::LArPropertiesServiceStandard::LArPropertiesServiceStandard(Parameters const& config,
                                                                    art::ActivityRegistry& reg)
  : fProp{config.get_PSet()}
{
  reg.sPreBeginRun.watch(this, &LArPropertiesServiceStandard::preBeginRun);
}

//----------------------------------------------
void
detinfo::LArPropertiesServiceStandard::preBeginRun(art::Run const& run)
{
  fProp.Update(run.run());
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(detinfo::LArPropertiesServiceStandard,
                                  detinfo::LArPropertiesService)
