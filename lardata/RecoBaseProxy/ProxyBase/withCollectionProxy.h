/**
 * @file   lardata/RecoBaseProxy/ProxyBase/withCollectionProxy.h
 * @brief  Creation of a collection proxy as auxiliary data for another proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_WITHCOLLECTIONPROXY_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_WITHCOLLECTIONPROXY_H

// LArSoft libraries
#include "lardata/RecoBaseProxy/ProxyBase/WithProxyAsAuxStructBase.h"

// framework libraries
#include "canvas/Utilities/InputTag.h"

// C/C++ standard
#include <tuple>       // also std::tuple_element_t<>
#include <type_traits> // std::is_convertible<>, std::decay_t<>, ...
#include <utility>     // std::forward(), std::move()

namespace proxy {

  // --- BEGIN Collection proxy infrastructure ---------------------------------
  /// @addtogroup LArSoftProxyCollections
  /// @{

  //----------------------------------------------------------------------------
  /// The same as `withCollectionProxy()`, but it also specified a tag.
  /// @bug Broken in many ways. Do not use.
  template <typename AuxProxy, typename AuxTag, typename... Args>
  auto withCollectionProxyAs(Args&&... args)
  {
    using ArgTuple_t = std::tuple<Args&&...>;
    static_assert(
      std::is_convertible<std::decay_t<std::tuple_element_t<0U, ArgTuple_t>>, art::InputTag>(),
      "The first argument of withCollectionProxyAs() must be art::InputTag.");
    ArgTuple_t argsTuple(std::forward<Args>(args)...);
    return details::WithProxyAsAuxStructBase<AuxProxy, ArgTuple_t, AuxTag>(std::move(argsTuple));
  } // withCollectionProxyAs()

  //----------------------------------------------------------------------------
  /**
   * @brief Helper function to merge an auxiliary proxy into the proxy.
   * @tparam AuxProxy type (proxy tag) of auxiliary collection proxy requested
   * @tparam Args types of constructor arguments for parallel data proxy
   * @param args constructor arguments for the parallel data collection proxy
   * @return a temporary object that `getCollection()` knows to handle
   * @bug Broken in many ways. Do not use.
   *
   * This function is meant to convey to `getCollection()` function the request
   * for merging a collection proxy to carry auxiliary data structured as
   * another collection proxy, parallel to the main collection.
   * The function also bridges the information required to create a proxy to
   * that auxiliary data.
   *
   * This data will be tagged with the type `AuxProxy`. To use a different type
   * as tag, use `withCollectionProxyAs()` instead, specifying the tag as second
   * template argument.
   *
   *
   * Customization of the auxiliary collection proxy
   * ================================================
   *
   * The customization of auxiliary collection proxy happens in a fashion
   * similar to the customization presented in `withParallelData()`.
   * The customization point here is `ProxyAsAuxProxyMaker`.
   */
  template <typename AuxProxy, typename... Args>
  auto withCollectionProxy(Args&&... args)
  {
    return withCollectionProxyAs<AuxProxy, AuxProxy>(std::forward<Args>(args)...);
  }

  /// @}
  // --- END Collection proxy infrastructure -----------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_WITHCOLLECTIONPROXY_H
