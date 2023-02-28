/**
 * @file   lardata/RecoBaseProxy/ProxyBase/ParallelDataProxyMaker.h
 * @brief  Infrastructure for the addition of parallel data to a proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_PARALLELDATAPROXYMAKER_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_PARALLELDATAPROXYMAKER_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...
#include "lardata/RecoBaseProxy/ProxyBase/ParallelData.h"
#include "lardata/RecoBaseProxy/ProxyBase/makeParallelDataFrom.h"

// framework libraries
#include "canvas/Utilities/InputTag.h"

// C/C++ standard
#include <utility> // std::forward()
#include <vector>

namespace proxy {

  // -- BEGIN Parallel data infrastructure -------------------------------------
  /**
   * @defgroup LArSoftProxiesParallelData Parallel data infrastructure
   * @ingroup LArSoftProxyCustom
   * @brief Infrastructure for support of parallel data structures.
   *
   * This infrastructure provides support for merging to a proxy data products
   * fulfilling the
   * @ref LArSoftProxyDefinitionParallelData "parallel data product"
   * requirements.
   *
   * The typical pattern is to read and merge a parallel data product with
   * `withParallelData()`. In alternative, a proxy can be augmented with any
   * existing collection, using `wrapParallelData()`.
   *
   * @{
   */

  /**
   * @brief Creates an parallel data wrapper for the specified types.
   * @tparam Main type of main datum (element)
   * @tparam AuxColl type of the parallel data collection
   * @tparam Aux type of the parallel data element
   * @tparam AuxTag tag labelling the parallel data collection
   *
   * Usually, `AuxTag` is also the type of datum (element) in the parallel
   * collection.
   *
   * This class works as a base class for `ParallelDataProxyMaker` so that
   * the specializations of the latter can still inherit from this one if they
   * its facilities.
   */
  template <typename Main,
            typename AuxColl,
            typename Aux,
            typename AuxTag = util::collection_value_t<AuxColl>>
  struct ParallelDataProxyMakerBase {

    /// Tag labelling the associated data we are going to produce.
    using data_tag = AuxTag;

    /// Type of the main datum.
    using main_element_t = Main;

    /// Type of the auxiliary data product.
    using aux_collection_t = AuxColl;

    /// Type of the auxiliary datum.
    using aux_element_t = Aux;

    /// Type of associated data proxy being created.
    using aux_collection_proxy_t = details::ParallelData<aux_collection_t, aux_element_t, data_tag>;

    /**
     * @brief Create a parallel data proxy collection using main collection tag.
     * @tparam Event type of the event to read data from
     * @tparam Handle type of handle to the main data product
     * @tparam MainArgs _(unused)_ any type convertible to `art::InputTag`
     * @param event the event to read data from
     * @param mainHandle _(unused)_ handle to the main collection data product
     * @param mainArgs an object describing the main data product
     * @return a parallel data proxy object
     *
     * The returned object exposes a random access container interface, with
     * data indexed by the index of the corresponding object in the main
     * collection.
     *
     * The `mainArgs` object is of an arbitrary type that must be convertible
     * by explicit type cast into a `art::InputTag`; that input tag will be
     * used to fetch the parallel data collection.
     */
    template <typename Event, typename Handle, typename MainArgs>
    static auto make(Event const& event, Handle&& mainHandle, MainArgs const& mainArgs)
    {
      return createFromTag(event, std::forward<Handle>(mainHandle), art::InputTag(mainArgs));
    }

    /**
     * @brief Create a parallel data proxy collection using the specified tag.
     * @tparam Event type of the event to read data from
     * @tparam Handle type of handle to the main data product
     * @tparam MainArgs (_unused_) any type convertible to `art::InputTag`
     * @param event the event to read associations from
     * @param mainHandle (_unused_) handle to the main collection data product
     * @param auxInputTag the tag of the data to be read
     * @return a parallel data proxy object
     *
     * The returned object exposes a random access container interface, with
     * data indexed by the index of the corresponding object in the main
     * collection.
     */
    template <typename Event, typename Handle, typename MainArgs>
    static auto make(Event const& event,
                     Handle&& mainHandle,
                     MainArgs const&,
                     art::InputTag const& auxInputTag)
    {
      return createFromTag(event, std::forward<Handle>(mainHandle), auxInputTag);
    }

    /**
     * @brief Create a parallel data proxy collection using the specified tag.
     * @tparam Event (_unused_) type of the event to read data from
     * @tparam Handle (_unused_) type of handle to the main data product
     * @tparam MainArgs (_unused_) any type convertible to `art::InputTag`
     * @param auxColl the collection to be wrapped
     * @return a parallel data proxy object
     *
     * The returned object exposes a random access container interface, with
     * data indexed by the index of the corresponding object in the main
     * collection.
     */
    template <typename Event, typename Handle, typename MainArgs>
    static auto make(Event const&, Handle&&, MainArgs const&, aux_collection_t const& auxColl)
    {
      return makeParallelDataFrom<aux_collection_t, aux_element_t, data_tag>(auxColl);
    }

  private:
    template <typename Event, typename Handle>
    static auto createFromTag(Event const& event, Handle&&, art::InputTag const& auxInputTag)
    {
      return makeParallelDataFrom<aux_collection_t, aux_element_t, data_tag>(event, auxInputTag);
    }

  }; // struct ParallelDataProxyMakerBase<>

  //--------------------------------------------------------------------------
  /**
   * @brief Creates an associated data wrapper for the specified types.
   * @tparam Main type of main datum (element) to associate from ("left")
   * @tparam Aux type of datum (element) to associate to ("right")
   * @tparam CollProxy type of proxy this associated data works for
   * @tparam Tag tag for the association proxy to be created
   * @tparam AuxColl type of the auxiliary data (default: `std::vector<Aux>`)
   * @see `withParallelDataAs()`, `wrapParallelDataAs()`
   *
   * This class is (indirectly) called when using `proxy::withParallelData()`
   * in `getCollection()`.
   * Its task is to supervise the creation of the proxy to the auxiliary data
   * parallel to the main data type.
   * The interface required by `withParallelData()` includes:
   * * a static `make()` method creating and returning the auxiliary data
   *   proxy with arguments an event, the main data product handle, a template
   *   argument representing the main collection information, and all the
   *   arguments required for the creation of the associated proxy (coming from
   *   `withParallelData()`); equivalent to the signature:
   *   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *   template <typename Event, typename Handle, typename MainArg, typename... Args>
   *   auto make(Event const&, Handle&&, MainArg const&, Args&&...);
   *   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * This class can be specialized (see `withParallelData()` for an example).
   * The default implementation just wraps a one-to-many
   * `std::vector<Aux>` data product fulfilling "parallel data product"
   * requirement (see the "Definitions" section in `ProxyBase.h` documentation).
   *
   * The last template argument is designed for specialization of auxiliary data
   * in the context of a specific proxy type.
   */
  template <typename Main,
            typename Aux,
            typename CollProxy,
            typename Tag = Aux,
            typename AuxColl = std::vector<Aux>>
  class ParallelDataProxyMaker : public ParallelDataProxyMakerBase<Main, AuxColl, Aux, Tag> {
    //
    // Note that this implementation is here only to document how to derive
    // a ParallelDataProxyMaker (specialization) from
    // ParallelDataProxyMakerBase. It's just mirroring the base class.
    //
    using base_t = ParallelDataProxyMakerBase<Main, AuxColl, Aux, Tag>;

  public:
    /// Type of the main datum.
    using typename base_t::main_element_t;

    /// Type of the auxiliary data product.
    using typename base_t::aux_collection_t;

    /// Type of the auxiliary datum.
    using typename base_t::aux_element_t;

    /// Type of collection data proxy being created.
    using typename base_t::aux_collection_proxy_t;

    /**
     * @brief Create a association proxy collection using main collection tag.
     * @tparam Event type of the event to read associations from
     * @tparam Handle type of data product handle
     * @tparam MainArgs any type convertible to `art::InputTag`
     * @tparam Args optional single type (`art::InputTag` required)
     * @param event the event to read associations from
     * @param mainHandle handle of the main collection data product
     * @param margs an object describing the main data product
     * @param args input tag for parallel data, if different from main
     * @return an auxiliary data proxy object
     *
     * The returned object exposes a random access container interface, with
     * data indexed by the index of the corresponding object in the main
     * collection.
     *
     * This implementation requires `margs` object to be convertible
     * by explicit type cast into a `art::InputTag`; that input tag will be
     * used to fetch the association.
     */
    template <typename Event, typename Handle, typename MainArgs, typename... Args>
    static auto make(Event const& event, Handle&& mainHandle, MainArgs const& margs, Args&&... args)
    {
      return base_t::make(
        event, std::forward<Handle>(mainHandle), margs, std::forward<Args>(args)...);
    }

  }; // struct ParallelDataProxyMaker<>

  /// @}
  // -- END Parallel data infrastructure ---------------------------------------

  //----------------------------------------------------------------------------
  namespace details {

    //--------------------------------------------------------------------------
    //--- stuff for parallel data collection (a form of auxiliary data)
    //--------------------------------------------------------------------------
    template <typename Aux, typename AuxTag, typename AuxColl = void>
    struct ParallelDataProxyMakerWrapper {
      template <typename CollProxy>
      using maker_t =
        ParallelDataProxyMaker<typename CollProxy::main_element_t, Aux, CollProxy, AuxTag, AuxColl>;
    }; // struct ParallelDataProxyMakerWrapper<Aux, AuxTag, AuxColl>

    template <typename Aux, typename AuxTag>
    struct ParallelDataProxyMakerWrapper<Aux, AuxTag, void> {
      template <typename CollProxy>
      using maker_t =
        ParallelDataProxyMaker<typename CollProxy::main_element_t, Aux, CollProxy, AuxTag>;
    }; // struct ParallelDataProxyMakerWrapper<Aux, AuxTag>

  } // namespace details

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_PARALLELDATAPROXYMAKER_H
