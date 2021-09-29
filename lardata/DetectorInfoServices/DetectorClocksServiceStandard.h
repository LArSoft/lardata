////////////////////////////////////////////////////////////////////////
// DetectorClocksServiceStandard.h
//
// Service interface for Detector Clock functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORCLOCKSSERVICESTANDARD_H
#define DETECTORCLOCKSSERVICESTANDARD_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include "fhiclcpp/fwd.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"

namespace detinfo {

  /**
   * @brief _art_ service managing `detinfo::DetectorClocksStandard`.
   * @see detinfo::DetectorClocksStandard, detinfo::DetectorClocks
   *
   * This _art_ service manages LArSoft's service provider
   * `detinfo::DetectorClocksStandard`, which implements
   * `detinfo::DetectorClocks` interface.
   *
   * For information about functionality of the service, see the documentation
   * of its interface, `detinfo::DetectorClocks`.
   * For information of the configuration, see also
   * `detinfo::DetectorClocksStandard`.
   *
   *
   * Configuration
   * ==============
   *
   * The configuration parameters are documented in the service provider
   * implementation: `detinfo::DetectorClocksStandard`.
   *
   *
   * Consistency check
   * ==================
   *
   * This service manager honors the `InheritClockConfig` configuration option
   * in the following way:
   * -# if the past jobs (explicitly excluding the current job) had inconsistent
   *     configuration, an exception is thrown claiming an "historical
   *     disagreement"
   * -# after the verification that the past configuration is consistent, the
   *     values from that configurations override the ones in the configuration
   *     of the current job; a value from the configuration of the current job
   *     is retained only if it was not present in the past (i.e. it is a new
   *     configuration parameter added since the input file was produced).
   *
   * The "past jobs" are the jobs that have produced the input file, and whose
   * configuration is stored by _art_ in the input file itself. The check and
   * reconfiguration is performed on each new input file.
   *
   *
   * Timing specifics
   * =================
   *
   * The trigger and beam gate times are set by this service before each event
   * is processed.
   * The logic is the following:
   * -# if the event contains a raw trigger (`raw::Trigger`) data product with
   *     input tag `TriggerName()` (from the configuration), that data product
   *     is read and the trigger and beam gate times stored in it are imported
   *     in the current service provider configuration; if there are more than
   *     one `raw::Trigger` objects in the data product, an exception is thrown
   * -# if no raw trigger is found with the specified label, the configuration
   *     of the service provider is updated using the default values of trigger
   *     and beam times specified in the service configuration
   *
   * The first set up happens on opening the first run in the first input file.
   * Accessing this service before (e.g. during `beginJob()` phase) yields
   * undefined behaviour.
   *
   */
  class DetectorClocksServiceStandard : public DetectorClocksService {
  public:
    DetectorClocksServiceStandard(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

  private:
    void preBeginRun(art::Run const& run);
    void postOpenFile(std::string const& filename);

    DetectorClocksData
    DataForJob() const override
    {
      return fClocks.DataForJob();
    }

    DetectorClocksData
    DataFor(art::Event const& e) const override;

    DetectorClocksStandard fClocks;
    bool fInheritClockConfig;
  };
} // namespace detinfo

DECLARE_ART_SERVICE_INTERFACE_IMPL(detinfo::DetectorClocksServiceStandard,
                                   detinfo::DetectorClocksService,
                                   SHARED)

#endif // DETECTORCLOCKSSERVICESTANDARD_H

