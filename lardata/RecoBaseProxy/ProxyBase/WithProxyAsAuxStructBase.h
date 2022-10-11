/**
 * @file   lardata/RecoBaseProxy/ProxyBase/WithProxyAsAuxStructBase.h
 * @brief  Infrastructure for a collection proxy as auxiliary data for a proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase/withCollectionProxy.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_WITHPROXYASAUXSTRUCTBASE_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_WITHPROXYASAUXSTRUCTBASE_H

// LArSoft libraries
#include "lardata/RecoBaseProxy/ProxyBase/ProxyAsAuxProxyMaker.h"

// framework libraries
#include "canvas/Utilities/InputTag.h"

// C/C++ standard
#include <cstdlib>     // std::size_t
#include <tuple>       // std::tuple_element_t<>, std::get(), ...
#include <type_traits> // std::is_convertible<>, std::decay_t<>, ...
#include <utility>     // std::forward(), std::move(), std::make_index_sequence()...

namespace proxy {

  namespace details {

    /// --- BEGIN LArSoftProxiesAuxProxy ---------------------------------------
    /// @addtogroup LArSoftProxiesAuxProxy
    /// @{

    //--------------------------------------------------------------------------
    /**
     * @brief Helper to create a proxy as auxiliary data for another proxy.
     * @tparam AuxProxy type of collection proxy associated to the main one
     * @tparam ArgTuple type of arguments required for the creation of proxy
     * @tparam AuxTag tag for the associated data (default: as `Aux`)
     *
     * This class stores user arguments for the construction of a collection
     * proxy to be used as auxiliary data for another proxy.
     * It can use those arguments plus some additional one to create the
     * collection proxy data itself. This additional information is provided by
     * `getCollection()`.
     *
     * The auxiliary data will be identified by type `AuxTag`.
     *
     * This is not a customization point: to have a custom associated data
     * produced, choose and then specialize the `ProxyAsAuxProxyMaker` class.
     */
    template <typename AuxProxy, typename ArgTuple, typename AuxTag = AuxProxy>
    class WithProxyAsAuxStructBase {

      static_assert(
        std::is_convertible<std::decay_t<std::tuple_element_t<0U, ArgTuple>>, art::InputTag>(),
        "The first argument of WithProxyAsAuxStructBase must be art::InputTag.");

      /// Type of main data product element from a proxy of type `CollProxy`.
      template <typename CollProxy>
      using main_t = typename CollProxy::main_element_t;

      /// Type of auxiliary proxy.
      using aux_proxy_t = AuxProxy;

      /// Tag for the associated data (same as the data type itself).
      using tag = AuxTag;

      /// Class to create the data proxy associated to a `CollProxy`.
      template <typename CollProxy>
      using proxy_maker_t = ProxyAsAuxProxyMaker<main_t<CollProxy>, aux_proxy_t, CollProxy, tag>;

    public:
      /// Constructor: steals the arguments, to be used by
      /// `createAuxProxyMaker()`.
      WithProxyAsAuxStructBase(ArgTuple&& args) : args(std::move(args)) {}

      /// Creates the associated data proxy by means of `ProxyAsAuxProxyMaker`.
      template <typename CollProxy, typename Event, typename Handle, typename MainArgs>
      auto createAuxProxyMaker(Event const& event, Handle&& mainHandle, MainArgs const& mainArgs)
      {
        return createAuxProxyImpl<CollProxy>(
          event, std::forward<Handle>(mainHandle), mainArgs, std::make_index_sequence<NArgs>());
      } // construct()

    protected:
      ArgTuple args; ///< Argument construction storage as tuple.

      /// Number of arguments stored.
      static constexpr std::size_t NArgs = std::tuple_size<ArgTuple>();

      // this method allows unpacking the arguments from the tuple
      template <typename CollProxy,
                typename Event,
                typename Handle,
                typename MainArgs,
                std::size_t... I>
      auto createAuxProxyImpl(Event const& event,
                              Handle&& mainHandle,
                              MainArgs const& mainArgs,
                              std::index_sequence<I...>)
      {
        return proxy_maker_t<CollProxy>::make(
          event, mainHandle, mainArgs, std::get<I>(std::forward<ArgTuple>(args))...);
      }

    }; // struct WithProxyAsAuxStructBase

    /// @}
    /// --- END LArSoftProxiesAuxProxy -----------------------------------------

  } // namespace details

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_WITHPROXYASAUXSTRUCTBASE_H
