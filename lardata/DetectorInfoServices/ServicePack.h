/**
 * @file   ServicePack.h
 * @brief  Utilities to manage ProviderPack objects with art
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   November 22, 2015
 * @see    ProviderPack.h
 */

#ifndef DETECTORINFOSERVICE_SERVICEPACK_H
#define DETECTORINFOSERVICE_SERVICEPACK_H 1

// LArSoft libraries
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/CoreUtils/ProviderPack.h"

namespace lar {
  /*
  namespace details {

  } // namespace details
  */

  /// Type of provider pack with providers from all specified Services
  template <typename... Services>
  using ProviderPackFromServices = lar::ProviderPack<typename Services::provider_type...>;

  /**
   * @brief Returns a provider pack with providers from specified services
   * @tparam Services the services to extract the providers from
   * @return a ProviderPack containing the current service providers
   *
   * This convenience function automatically extracts all the service providers
   * from a list of services.
   * This is convenient if an algorithm or service accepts a provider pack
   * for setup:
   *
   *     algo->Setup(extractProviders<
   *       detinfo::DetectorPropertiesService, detinfo::LArPropertiesService
   *       >());
   *
   * Also note that the provider packs can rearrange their elements, so the call
   * above should work just the same as:
   *
   *     algo->Setup(extractProviders<
   *       detinfo::LArPropertiesService, detinfo::DetectorPropertiesService
   *       >());
   *
   * If a provider is needed in the setup argument that is not provided by any
   * of the specified services, a compilation error will occur.
   */
  template <typename... Services>
  ProviderPackFromServices<Services...> extractProviders()
  {
    return {lar::providerFrom<Services>()...};
  }

} // namespace lar

//==============================================================================
//=== template implementation
//===
/*
namespace lar {
  namespace details {
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
  } // namespace details
} // namespace lar
*/

#endif // DETECTORINFOSERVICE_SERVICEPACK_H
