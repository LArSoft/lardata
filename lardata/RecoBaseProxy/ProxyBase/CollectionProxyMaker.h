/**
 * @file   lardata/RecoBaseProxy/ProxyBase/CollectionProxyMaker.h
 * @brief  Infrastructure for the creation of a collection proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase/getCollection.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_COLLECTIONPROXYMAKER_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_COLLECTIONPROXYMAKER_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...
#include "larcorealg/CoreUtils/MetaUtils.h"     // util::always_true_type
#include "lardata/RecoBaseProxy/ProxyBase/CollectionProxy.h"

// framework libraries
#include "canvas/Utilities/InputTag.h"

// C/C++ standard
#include <vector>

namespace proxy {

  // --- BEGIN Collection proxy infrastructure ---------------------------------
  /// @addtogroup LArSoftProxyCollections
  /// @{

  //----------------------------------------------------------------------------
  /**
   * @brief Collection of data type definitions for collection proxies.
   * @tparam Proxy type of proxy the traits refer to
   * @tparam Selector specialization helper type
   *
   * Expected traits:
   * * `main_collection_t`: type of the main data product collection
   * * `main_element_t`: type contained in the main data product collection
   * * `main_collection_proxy_t`: type wrapping the main data product collection
   * * `collection_proxy_impl_t`: type of implementation object used for the
   *     proxy; to use `proxy::CollectionProxy`, define this as:
   *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *     template <typename... Args>
   *     using collection_proxy_impl_t = proxy::CollectionProxyFromArgs<Args...>;
   *     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *     This trait is treated as optional by the default `CollectionProxyMaker`
   *     implementation, which uses `proxy::CollectionProxy` if the trait is
   *     missing.
   *
   * Note that the `Proxy` type is expected to be the same type as used in
   * `getCollection()` calls, and does not need to match the actual type of a
   * proxy collection.
   */
  template <typename Proxy, typename Selector = Proxy>
  struct CollectionProxyMakerTraits {
    static_assert(util::always_true_type<Proxy>(), "This class requires specialization.");
  };

  //----------------------------------------------------------------------------
  /**
   * @brief Class to assemble the required proxy.
   * @tparam CollProxy a type characterizing the produced proxy
   *
   * This is a base class suitable for deriving `CollectionProxyMaker`,
   * which is used by `getCollection()` to create the requested proxy.
   * Deriving a specialization of `CollectionProxyMaker` from
   * `CollectionProxyMakerBase` provides some basic definitions and
   * functionality that might be tedious to copy.
   *
   * In this, which can be considered as the default implementation of
   * `CollectionProxyMaker`, the template argument `CollProxy` is expected to
   * have the interface of `CollectionProxy`.
   * The arguments required are documented in the implementation of the `make()`
   * function. Note that the type of proxy returned does not need to be
   * `CollProxy`, but is in fact an instance of `proxy::CollectionProxy`.
   *
   * @attention This class is not meant to be specialized.
   */
  template <typename CollProxy>
  struct CollectionProxyMakerBase {

    /// Traits of the collection proxy for the collection proxy maker.
    using traits_t = CollectionProxyMakerTraits<CollProxy>;

    /// Type of main collection proxy.
    using main_collection_proxy_t = typename traits_t::main_collection_proxy_t;

    /// Type returned by the main collection indexing operator.
    using main_element_t = typename traits_t::main_element_t;

    /// Type of element of the main collection.
    using main_collection_t = typename traits_t::main_collection_t;

    /**
     * @brief Creates and returns a collection proxy based on `CollProxy` and
     *        with the requested associated data.
     * @tparam Event type of the event to read the information from
     * @tparam WithArgs type of arguments for associated data
     * @param event event to read the information from
     * @param tag input tag of the main data product
     * @param withArgs optional associated objects to be included
     * @return a collection proxy from data as retrieved via `tag`
     *
     * For each argument in `withArgs`, an action is taken. Usually that is to
     * add an association to the proxy.
     *
     * Only a few associated data collections are supported:
     * * `withAssociated<Aux>()` (optional argument: hit-track association
     *     tag, as track by default): adds to the proxy an association to the
     *     `Aux` data product (see `withAssociated()`)
     *
     */
    template <typename Event, typename... WithArgs>
    static auto make(Event const& event, art::InputTag const& tag, WithArgs&&... withArgs);

  }; // struct CollectionProxyMakerBase<>

  //----------------------------------------------------------------------------
  /**
   * @brief Class to assemble the required proxy.
   * @tparam CollProxy a type characterizing the produced proxy
   *
   * This class is used by `getCollection()` to create the requested proxy.
   * The required interface for this class is:
   * * `make()`: a static method returning the collection proxy, matching the
   *   signature
   *   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *   template <typename Event, typename... Args>
   *   auto make(Event const&, Args&&...);
   *   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *   where the first argument will be the event to read the information from.
   *
   * A default implementation is provided as `proxy::CollectionProxyMakerBase`.
   * Specializations of `proxy::CollectionProxyMaker` may choose to derive from
   * `proxy::CollectionProxyMakerBase` as well, for convenience.
   */
  template <typename CollProxy>
  struct CollectionProxyMaker : CollectionProxyMakerBase<CollProxy> {};

  /// @}
  // --- END Collection proxy infrastructure -----------------------------------

  //----------------------------------------------------------------------------
  //--- specializations of CollectionProxyMakerTraits
  //----------------------------------------------------------------------------

  template <typename T>
  struct CollectionProxyMakerTraits<std::vector<T>> {

    /// Type of element of the main collection.
    using main_collection_t = std::vector<T>;

    /// Type returned by the main collection indexing operator.
    using main_element_t = util::collection_value_t<main_collection_t>;

    /// Type of main collection proxy.
    using main_collection_proxy_t = details::MainCollectionProxy<main_collection_t>;

  }; // CollectionProxyMakerTraits<std::vector<T>>

  template <typename MainColl>
  struct CollectionProxyMakerTraits<CollectionProxy<MainColl>> {
    /// Type of main collection proxy.
    using main_collection_proxy_t = CollectionProxy<MainColl>;

    /// Type returned by the main collection indexing operator.
    using main_element_t = typename main_collection_proxy_t::main_element_t;

    /// Type of element of the main collection.
    using main_collection_t = typename main_collection_proxy_t::main_collection_t;

  }; // CollectionProxyMakerTraits<CollectionProxy>

  // ---------------------------------------------------------------------------

} // namespace proxy

// -----------------------------------------------------------------------------
// ---  Template implementation
// -----------------------------------------------------------------------------
namespace proxy {

  namespace details {

    // -------------------------------------------------------------------------
    //
    // We need to discover whether the trait class contains a
    // `collection_proxy_impl_t` type. Unfortunately, that is not a plain type
    // but rather a template type (it is expected to take the same template
    // arguments as `proxy::CollectionProxy`).
    //
    // Since Gianluca is not able to detect the presence of a templated type in
    // a class, but only of fully specified type, the following machinery is in
    // place. At the time we need to know the `collection_proxy_impl_t` type, we
    // are trying to create one and therefore we know the full argument set for
    // it already (let's call it `Args...`, that is `MainColl` plus
    // `AuxColls...`). Then, we test the presence of
    // `collection_proxy_impl_t<Args...>` instead of just
    // `collection_proxy_impl_t`.
    // An additional complication is the presence of a parameter pack, that
    // impedes the usual pattern of having `std::enable_if_t` as the last,
    // optional template parameter (`Args` must be the last). So we give up the
    // elegant feature of having a default value, and we force the caller to
    // explicitly spell it: `void` (always that one).
    // We absorb this ugliness and the need for the usual and cumbersome
    // `typename ...::type` in the definition of
    // `CollectionProxyImplFromTraits_t`, which should be the only access point
    // to this machinery.
    // About the naming: a "collection proxy" is the final proxy class.
    // `CollectionProxy` is also the name, already taken, of the class we
    // provide as reference for that role. So I have called the role an
    // "implementation" of the "collection proxy" idea, and `CollectionProxy` is
    // in fact an implementation (the "standard" one) of that idea.
    // The name of that idea in the proxy traits is `collection_proxy_impl_t`.
    // `CollectionProxyImpl` is supposed to be a generic name for that concept,
    // and the utility `CollectionProxyImplFromTraits` is in charge of
    // extracting it from the traits. This utility also needs an implementation,
    // hence the second "Impl": `CollectionProxyImplFromTraitsImpl`. Urgh.
    //
    // I don't know whether to be worried that the explanation of this code is
    // longer than the code itself...
    //
    template <typename Traits, typename, typename... Args>
    struct CollectionProxyImplFromTraitsImpl {
      using type = CollectionProxyFromArgs<Args...>;
    };

    template <typename Traits, typename... Args>
    struct CollectionProxyImplFromTraitsImpl<
      Traits,
      std::enable_if_t<
        util::always_true_v<typename Traits::template collection_proxy_impl_t<Args...>>>,
      Args...> {
      using type = typename Traits::template collection_proxy_impl_t<Args...>;
    };

    /// `Traits::collection_proxy_impl_t` if that (template) type exists,
    /// `proxy::CollectionProxy` otherwise.
    template <typename Traits, typename... Args>
    using CollectionProxyImplFromTraits_t =
      typename CollectionProxyImplFromTraitsImpl<Traits, void, Args...>::type;

    // -------------------------------------------------------------------------
    template <typename Traits, typename... Args>
    auto createCollectionProxyFromTraits(Args&&... args)
    {
      using collection_proxy_impl_t =
        CollectionProxyImplFromTraits_t<Traits, std::decay_t<Args>...>;
      return collection_proxy_impl_t(std::forward<Args>(args)...);
    } // createCollectionProxyFromTraits()

    // -------------------------------------------------------------------------

  } // namespace details

  // ---------------------------------------------------------------------------
  // ---  CollectionProxyMakerBase implementation
  // ---------------------------------------------------------------------------
  template <typename CollProxy>
  template <typename Event, typename... WithArgs>
  auto CollectionProxyMakerBase<CollProxy>::make(Event const& event,
                                                 art::InputTag const& tag,
                                                 WithArgs&&... withArgs)
  {
    auto mainHandle = event.template getValidHandle<main_collection_t>(tag);

    // The actual type of collection proxy implementation is extracted from
    // the traits of the proxy (`collection_proxy_impl_t`), but if that is not
    // provided a default implementation, `proxy::CollectionProxy`, is used:
    return details::createCollectionProxyFromTraits<traits_t>(
      *mainHandle,
      withArgs.template createAuxProxyMaker<main_collection_proxy_t>(event, mainHandle, tag)...);
  } // CollectionProxyMakerBase<>::make<>()

  // ---------------------------------------------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_COLLECTIONPROXYMAKER_H
