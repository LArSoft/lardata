#include "lardata/DetectorInfoServices/LArPropertiesServiceStandard.h"

#include "art/Framework/Principal/Run.h"

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

