/**
 * @file   ServicePackTest_module.cc
 * @brief  Tests utilities in ServicePack.h
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   November 23, 2015
 * @see    ServicePack.h
 */

// LArSoft libraries
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/ServicePack.h" // lar::providerFrom<>()
#include "lardataalg/DetectorInfo/LArProperties.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <sstream>
#include <string>

namespace fhicl {
  class ParameterSet; //
} // namespace fhicl

namespace { // local

  template <typename T>
  inline std::string to_string(T const& v)
  {
    return std::to_string(v);
  }

  template <typename T>
  std::string to_string(T const* ptr)
  {
    std::ostringstream sstr;
    sstr << '<' << ((void*)ptr) << '>';
    return sstr.str();
  } // to_string()

} // local namespace

namespace lar {

  /**
   * @brief Test module for ServicePack.h utilities depending on art farmework
   *
   * Currently exercises:
   * - `lar::extractProviders()`
   *
   * Throws an exception on failure.
   *
   * Service requirements
   * =====================
   *
   * This module requires the following services to be configured:
   * - Geometry
   * - LArPropertiesService
   * - DetectorClocksService
   * - DetectorPropertiesService
   *
   * Configuration parameters
   * =========================
   *
   * Currently none.
   *
   */
  class ServicePackTest : public art::EDAnalyzer {
  public:
    explicit ServicePackTest(fhicl::ParameterSet const&);

  private:
    /// Run event-independent tests
    void beginJob() override;

    /// Run event-dependent tests (none so far)
    void analyze(const art::Event& /* evt */) override {}

    /// Throws if errors have been accumulated
    void endJob() override;

    /// @{
    /// @name Test functions

    /// Tests lar::extractProviders()
    void extractProviders_test_plain();

    /// Tests lar::extractProviders() and permuted constructor
    void extractProviders_test_permuted();

    /// Tests lar::extractProviders() and assignment to reduced pack
    void extractProviders_test_reduced();

    /// All tests on lar::extractProviders()
    void extractProviders_tests();

    /// @}

    std::vector<std::string> errors; ///< list of collected errors

  }; // ServicePackTest

  DEFINE_ART_MODULE(ServicePackTest)

} // namespace lar

//------------------------------------------------------------------------------
//--- implementation
//---
namespace lar {

  //----------------------------------------------------------------------------
  ServicePackTest::ServicePackTest(const fhicl::ParameterSet& pset) : EDAnalyzer(pset) {}

  //----------------------------------------------------------------------------
  void ServicePackTest::beginJob()
  {
    extractProviders_tests();
  } // ServicePackTest::beginJob()

  //----------------------------------------------------------------------------
  void ServicePackTest::endJob()
  {
    if (errors.empty()) {
      mf::LogInfo("ServicePackTest") << "All tests were successful.";
      return;
    }

    mf::LogError log("ServicePackTest");
    log << errors.size() << " errors detected:";

    for (std::string const& error : errors)
      log << "\n - " << error;

    throw art::Exception(art::errors::LogicError) << errors.size() << " errors detected";
  }

  //----------------------------------------------------------------------------
  void ServicePackTest::extractProviders_tests()
  {
    extractProviders_test_plain();
    extractProviders_test_permuted();
    extractProviders_test_reduced();
  }

  //----------------------------------------------------------------------------
  void ServicePackTest::extractProviders_test_plain()
  {
    /*
     * The test creates a ProviderPack and checks that its element as as
     * expected.
     *
     * The expected value is extracted from the framework in the "traditional"
     * way.
     */

    // these are the "solutions":
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
    detinfo::LArProperties const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();

    auto providers = lar::extractProviders<geo::Geometry, detinfo::LArPropertiesService>();

    // check time
    if (providers.get<geo::GeometryCore>() != geom) {
      errors.push_back("wrong geometry provider (got " +
                       ::to_string(providers.get<geo::GeometryCore>()) + ", expected " +
                       ::to_string(geom) + ")");
    }
    if (providers.get<detinfo::LArProperties>() != larprop) {
      errors.push_back("wrong LAr properties provider (got " +
                       ::to_string(providers.get<detinfo::LArProperties>()) + ", expected " +
                       ::to_string(larprop) + ")");
    }
  } // ServicePackTest::extractProviders_test_plain()

  //----------------------------------------------------------------------------
  void ServicePackTest::extractProviders_test_permuted()
  {
    /*
     * The test creates a ProviderPack and checks that its element as as
     * expected.
     *
     * The expected value is extracted from the framework in the "traditional"
     * way.
     *
     * The order of the providers is different from the order of the services;
     * in this way a "wrong" ProviderPack will be (deliberately) created,
     * and the code will have to convert it to the right pack.
     */

    // these are the "solutions":
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
    detinfo::LArProperties const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();

    auto providers = lar::extractProviders<geo::Geometry, detinfo::LArPropertiesService>();

    // check time
    if (providers.get<geo::GeometryCore>() != geom) {
      errors.push_back("wrong geometry provider (got " +
                       ::to_string(providers.get<geo::GeometryCore>()) + ", expected " +
                       ::to_string(geom) + ") [permuted]");
    }
    if (providers.get<detinfo::LArProperties>() != larprop) {
      errors.push_back("wrong LAr properties provider (got " +
                       ::to_string(providers.get<detinfo::LArProperties>()) + ", expected " +
                       ::to_string(larprop) + ") [permuted]");
    }
  }

  //----------------------------------------------------------------------------
  void ServicePackTest::extractProviders_test_reduced()
  {
    /*
     * The test creates a ProviderPack and checks that its element as as
     * expected.
     *
     * The expected value is extracted from the framework in the "traditional"
     * way.
     *
     * We use a smaller provider pack to store the providers;
     * DetectorProperties will be dropped.
     */

    // these are the "solutions":
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
    detinfo::LArProperties const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();
    auto providers = lar::extractProviders<geo::Geometry, detinfo::LArPropertiesService>();

    // check time
    if (providers.get<geo::GeometryCore>() != geom) {
      errors.push_back("wrong geometry provider (got " +
                       ::to_string(providers.get<geo::GeometryCore>()) + ", expected " +
                       ::to_string(geom) + ") [reduced]");
    }
    if (providers.get<detinfo::LArProperties>() != larprop) {
      errors.push_back("wrong LAr properties provider (got " +
                       ::to_string(providers.get<detinfo::LArProperties>()) + ", expected " +
                       ::to_string(larprop) + ") [reduced]");
    }
  }

} // namespace lar
