/**
 * @file   lardata/RecoBaseProxy/ProxyBase/WithAssociatedStructBase.h
 * @brief  Template class to declare addition of associated data to a proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_WITHASSOCIATEDSTRUCTBASE_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_WITHASSOCIATEDSTRUCTBASE_H

// C/C++ standard
#include <cstdlib> // std::size_t
#include <tuple>   // std::tuple_size(), std::get()
#include <utility> // std::forward(), std::move(), std::index_sequence<>...

/**
 * @defgroup LArSoftProxiesAssociatedData Associated data infrastructure
 * @ingroup  LArSoftProxyCustom
 * @brief Infrastructure for support of associated data.
 *
 * Associated data is auxiliary data connected to the main data via _art_
 * associations.
 * The following associated data are currently supported:
 * * @ref LArSoftProxyDefinitionOneToManySeqAssn "one-to-many sequential association"
 *     implicitly supporting also one-to-any (one-to-one, one-to-zero/or/one) in
 *     a non-optimized way
 * * @ref LArSoftProxyDefinitionOneToZeroOrOneSeqAssn "one-to-(zero-or-one) sequential association"
 * * @ref LArSoftProxyDefinitionParallelData "parallel data product", implicit
 *     one-to-one associations
 *
 */

namespace proxy {

  //----------------------------------------------------------------------------
  namespace details {

    //--------------------------------------------------------------------------
    /**
     * @brief Helper to create associated data proxy.
     * @tparam Aux type of data associated to the main one
     * @tparam Metadata type of metadata of the association
     * @tparam ArgTuple type of arguments required for the creation of proxy
     * @tparam ProxyMaker template type of the proxy maker class
     * @tparam AuxTag tag for the associated data (default: as `Aux`)
     *
     * This class stores user arguments for the construction of a proxy to
     * associated data of type `Aux`.
     * It can use that information plus some additional one to create the
     * associated data itself. This additional information is provided by
     * `getCollection()`.
     *
     * The association will be identified by type `AuxTag`.
     *
     * This is not a customization point: to have a custom associated data
     * produced, choose and then specialize the `ProxyMaker` class.
     */
    template <typename Aux,
              typename Metadata,
              typename ArgTuple,
              template <typename CollProxy>
              class ProxyMaker,
              typename AuxTag /* = Aux */
              >
    class WithAssociatedStructBase {

      /// Type of main data product element from a proxy of type `CollProxy`.
      template <typename CollProxy>
      using main_t = typename CollProxy::main_element_t;

      /// Type of associated data.
      using aux_t = Aux;

      /// Type of associated metadata.
      using metadata_t = Metadata;

      /// Tag for the associated data (same as the data type itself).
      using tag = AuxTag;

      /// Class to create the data proxy associated to a `CollProxy`.
      template <typename CollProxy>
      using proxy_maker_t = ProxyMaker<CollProxy>;

    public:
      /// Type of association proxy created for the specified `CollProxy`.
      template <typename CollProxy>
      using aux_collection_proxy_t = typename proxy_maker_t<CollProxy>::aux_collection_proxy_t;

      /// Constructor: steals the arguments, to be used by
      /// `createAuxProxyMaker()`.
      WithAssociatedStructBase(ArgTuple&& args) : args(std::move(args)) {}

      /// Creates the associated data proxy by means of `ProxyMaker`.
      template <typename CollProxy, typename Event, typename Handle, typename MainArgs>
      auto createAuxProxyMaker(Event const& event, Handle&& mainHandle, MainArgs const& mainArgs)
      {
        return createAssnProxyMaker<CollProxy>(
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
      auto createAssnProxyMaker(Event const& event,
                                Handle&& mainHandle,
                                MainArgs const& mainArgs,
                                std::index_sequence<I...>)
      {
        return proxy_maker_t<CollProxy>::make(
          event, mainHandle, mainArgs, std::get<I>(std::forward<ArgTuple>(args))...);
      }

    }; // struct WithAssociatedStructBase<>

    //--------------------------------------------------------------------------

  } // namespace details

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_WITHASSOCIATEDSTRUCTBASE_H
