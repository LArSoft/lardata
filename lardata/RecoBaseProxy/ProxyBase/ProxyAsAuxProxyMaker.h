/**
 * @file   lardata/RecoBaseProxy/ProxyBase/ProxyAsAuxProxyMaker.h
 * @brief  Infrastructure for a collection proxy as auxiliary data for a proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase/withCollectionProxy.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_PROXYASAUXPROXYMAKER_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_PROXYASAUXPROXYMAKER_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...
#include "lardata/RecoBaseProxy/ProxyBase/ProxyAsParallelData.h"
#include "lardata/RecoBaseProxy/ProxyBase/getCollection.h"

// framework libraries
#include "canvas/Utilities/InputTag.h"

// C/C++ standard
#include <utility> // std::forward(), std::move()

namespace proxy {

  /// --- BEGIN LArSoftProxiesAuxProxy -----------------------------------------
  /// @addtogroup LArSoftProxiesAuxProxy
  /// @{

  /**
   * @brief Creates a proxy wrapper for merging into another proxy ("main").
   * @tparam Main type of main datum (element) of the main proxy
   * @tparam AuxProxy type ("proxy name") of the proxy being wrapped
   * @tparam AuxTag tag of the auxiliary proxy in the context of the main one
   *
   * By default, `AuxTag` is the same as the proxy name.
   *
   * This class works as a base class for `ProxyAsAuxProxyMaker` so that
   * the specializations of the latter can still inherit from this one if they
   * its facilities.
   */
  template <typename Main, typename AuxProxy, typename AuxTag = AuxProxy>
  struct ProxyAsAuxProxyMakerBase {

    /// Tag labelling the associated data we are going to produce.
    using data_tag = AuxTag;

    /// Type of the main datum.
    using main_element_t = Main;

    /// Tag-type of the auxiliary proxy (not the type of the proxy!).
    using aux_proxy_t = AuxProxy;

    /**
     * @brief Create a parallel data proxy collection using the specified tag.
     * @tparam Event type of the event to read data from
     * @tparam Handle (_unused_) type of handle to the main data product
     * @tparam MainArgs (_unused_) any type convertible to `art::InputTag`
     * @tparam AuxArgs type of arguments for the creation of the auxiliary proxy
     * @param event event to create the proxy from
     * @param auxProxyTag tag for the creation of the auxiliary collection proxy
     * @param args other arguments for the creation of the auxiliary proxy
     * @return a auxiliary proxy data object
     *
     * The returned object exposes a random access container interface, with
     * data indexed by the index of the corresponding object in the main
     * collection.
     *
     * The tag of the main collection proxy is ignored even if present, and
     * the caller must specify it.
     */
    template <typename Event, typename Handle, typename MainArgs, typename... AuxArgs>
    static auto make(Event const& event,
                     Handle&&,
                     MainArgs const&,
                     art::InputTag const& auxProxyTag,
                     AuxArgs&&... args)
    {
      auto auxProxy = makeAuxiliaryProxy(event, auxProxyTag, std::forward<AuxArgs>(args)...);
      return makeProxyAsParallelData<data_tag, util::collection_value_t<decltype(auxProxy)>>(
        std::move(auxProxy));
    }

  private:
    /// Creates the proxy to be used as parallel data.
    template <typename Event, typename... AuxArgs>
    static auto makeAuxiliaryProxy(Event const& event,
                                   art::InputTag const& auxProxyTag,
                                   AuxArgs&&... args)
    {
      return getCollection<aux_proxy_t>(event, auxProxyTag, std::forward<AuxArgs>(args)...);
    }

  }; // struct ProxyAsAuxProxyMakerBase<>

  //--------------------------------------------------------------------------
  /**
   * @brief Creates an auxiliary proxy wrapper for the specified proxy.
   * @tparam Main type of main datum (element) to associate from ("left")
   * @tparam AuxProxy type of proxy collection to be associated
   * @tparam CollProxy type of proxy this associated data works for
   * @tparam Tag tag for the association proxy to be created
   * @see `withCollectionProxy()`
   *
   * This class is (indirectly) called when using `proxy::withCollectionProxy()`
   * in `getCollection()`.
   * Its task is to supervise the creation of the collection proxy that is used
   * as auxiliary data for the main data type.
   * The interface required by `withCollectionProxy()` includes:
   * * a static `make()` method creating and returning the auxiliary data
   *   proxy with arguments an event, the main data product handle, a template
   *   argument representing the main collection information, and all the
   *   arguments required for the creation of the auxiliary collection proxy
   *   (coming from `withCollectionProxy()`); equivalent to the signature:
   *   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *   template <typename Event, typename Handle, typename MainArg, typename... Args>
   *   auto make(Event const&, Handle&&, MainArg const&, Args&&...);
   *   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * This class can be specialized.
   * The default implementation just uses `getCollection()` to create the
   * auxiliary proxy, and merges it to the main collection proxy in a fashion
   * similar to parallel data.
   *
   * The template argument `CollProxy` is designed for specialization of
   * auxiliary data in the context of a specific proxy type.
   */
  template <typename Main,
            typename AuxProxy,
            typename CollProxy,
            typename Tag = util::collection_value_t<AuxProxy>>
  class ProxyAsAuxProxyMaker : public ProxyAsAuxProxyMakerBase<Main, AuxProxy, Tag> {};

  /// @}
  /// --- END LArSoftProxiesAuxProxy -------------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_PROXYASAUXPROXYMAKER_H
