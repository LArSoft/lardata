/**
 * @file   lardata/RecoBaseProxy/ProxyBase/OneTo01DataProxyMaker.h
 * @brief  Infrastructure for merging optional associated data to a proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase/withZeroOrOne.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_ONETO01DATAPROXYMAKER_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_ONETO01DATAPROXYMAKER_H

// LArSoft libraries
#include "lardata/RecoBaseProxy/ProxyBase/OneTo01Data.h"
#include "lardata/RecoBaseProxy/ProxyBase/makeOneTo01dataFrom.h"

// framework libraries
#include "canvas/Utilities/InputTag.h"

// C/C++ standard libraries
#include <type_traits> // std::is_convertible<>
#include <utility>     // std::forward()

namespace proxy {

  // --- BEGIN LArSoftProxiesAssociatedData ------------------------------------
  /// @addtogroup LArSoftProxiesAssociatedData
  /// @{

  //----------------------------------------------------------------------------
  //--- one-to-(zero/one) associations
  //----------------------------------------------------------------------------
  /**
   * @brief Creates an one-to-(zero-or-one) wrapper for the specified types.
   * @tparam Main type of main datum (element) to associate from ("left")
   * @tparam Aux type of datum (element) to associate to ("right")
   * @tparam Metadata type of metadata coming with the association
   * @tparam AuxTag tag labelling this association
   *
   * Usually, `AuxTag` is also the type of datum (element) to associate to
   * ("right").
   *
   * This class works as a base class for `OneTo01DataProxyMaker` so that
   * the specializations of the latter can still inherit from this one if they
   * its facilities.
   */
  template <typename Main, typename Aux, typename Metadata = void, typename AuxTag = Aux>
  struct OneTo01DataProxyMakerBase {

    /// Tag labelling the associated data we are going to produce.
    using data_tag = AuxTag;

    /// Type of the main datum ("left").
    using main_element_t = Main;

    /// Type of the auxiliary associated datum ("right").
    using aux_element_t = Aux;

    /// Type of associated metadata.
    using metadata_t = Metadata;

    /// Type of associated data proxy being created.
    using aux_collection_proxy_t =
      details::OneTo01Data<main_element_t, aux_element_t, metadata_t, data_tag>;

    /// Type of _art_ association being used as input.
    using assns_t = typename aux_collection_proxy_t::assns_t;

    /**
     * @brief Create a association proxy collection using main collection tag.
     * @tparam Event type of the event to read associations from
     * @tparam Handle type of handle to the main data product
     * @tparam MainArgs any type convertible to `art::InputTag`
     * @param event the event to read associations from
     * @param mainHandle handle to the main collection data product
     * @param mainArgs an object describing the main data product
     * @return an auxiliary data proxy object
     *
     * The returned object exposes a random access container interface, with
     * data indexed by the index of the corresponding object in the main
     * collection.
     *
     * The `mainArgs` object is of an arbitrary type that must be convertible
     * by explicit type cast into a `art::InputTag`; that input tag will be
     * used to fetch the association.
     */
    template <typename Event, typename Handle, typename MainArgs>
    static auto make(Event const& event, Handle&& mainHandle, MainArgs const& mainArgs)
    {
      return createFromTag(event, std::forward<Handle>(mainHandle), art::InputTag(mainArgs));
    }

    /**
     * @brief Create a association proxy collection using the specified tag.
     * @tparam Event type of the event to read associations from
     * @tparam Handle type of handle to the main data product
     * @tparam MainArgs any type convertible to `art::InputTag` (unused)
     * @param event the event to read associations from
     * @param mainHandle handle to the main collection data product
     * @param auxInputTag the tag of the association to be read
     * @return a auxiliary data proxy object
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
     * @brief Create a association proxy collection using the specified tag.
     * @tparam Event type of the event to read associations from (unused)
     * @tparam Handle type of handle to the main data product
     * @tparam MainArgs any type convertible to `art::InputTag` (unused)
     * @param handle handle to the main data collection
     * @param assns the associations to be wrapped
     * @return a auxiliary data proxy object
     *
     * The returned object exposes a random access container interface, with
     * data indexed by the index of the corresponding object in the main
     * collection.
     */
    template <typename Event, typename Handle, typename MainArgs, typename Assns>
    static auto make(Event const&, Handle&& handle, MainArgs const&, Assns const& assns)
    {
      static_assert(std::is_convertible<typename Assns::right_t, aux_element_t>(),
                    "Improper right type for one-to-(zero-or-one) association.");
      return makeOneTo01dataFrom<data_tag>(assns, handle->size());
    }

  private:
    template <typename Event, typename Handle>
    static auto createFromTag(Event const& event,
                              Handle&& mainHandle,
                              art::InputTag const& auxInputTag)
    {
      return makeOneTo01dataFrom<main_element_t, aux_element_t, metadata_t, data_tag>(
        event, auxInputTag, mainHandle->size());
    }

  }; // struct OneTo01DataProxyMakerBase<>

  //--------------------------------------------------------------------------
  /**
   * @brief Creates an one-to-(zero-or-one) wrapper for the specified types.
   * @tparam Main type of main datum (element) to associate from ("left")
   * @tparam Aux type of datum (element) to associate to ("right")
   * @tparam Metadata type of metadata in the association
   * @tparam CollProxy type of proxy this associated data works for
   * @tparam Tag tag for the association proxy to be created
   * @see `withZeroOrOne()`
   * This class is (indirectly) called when using `proxy::withZeroOrOne()`
   * in `getCollection()`.
   * Its task is to supervise the creation of the proxy to the data
   * association between the main data type and an auxiliary one.
   * The interface required by `withZeroOrOne()` includes:
   * * a static `make()` method creating and returning the associated data
   *   proxy with arguments an event, the main data product handle, a template
   *   argument representing the main collection information, and all the
   *   arguments required for the creation of the associated proxy (coming from
   *   `withZeroOrOne()`); equivalent to the signature:
   *   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *   template <typename Event, typename Handle, typename MainArg, typename... Args>
   *   auto make(Event const&, Handle&&, MainArg const&, Args&&...);
   *   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * This class can be specialized (see `withAssociated()` for an example
   * of the general procedure).
   * The default implementation just wraps a one-to-(zero-or-one)
   * `art::Assns<Main, Aux>` data product fulfilling
   * @ref LArSoftProxyDefinitionOneToZeroOrOneSeqAssn "one-to-(zero or one) sequential association"
   * requirement.
   *
   * The last template argument is designed for specialization of associations
   * in the context of a specific proxy type.
   */
  template <typename Main, typename Aux, typename Metadata, typename CollProxy, typename Tag = Aux>
  class OneTo01DataProxyMaker : public OneTo01DataProxyMakerBase<Main, Aux, Metadata, Tag> {
    //
    // Note that this implementation is here only to document how to derive
    // a OneTo01DataProxyMaker (specialization) from
    // OneTo01DataProxyMakerBase. It's just mirroring the base class.
    //
    using base_t = OneTo01DataProxyMakerBase<Main, Aux, Metadata, Tag>;

  public:
    /// Type of the main datum ("left").
    using typename base_t::main_element_t;

    /// Type of the auxiliary associated datum ("right").
    using typename base_t::aux_element_t;

    /// Type of metadata in the association.
    using typename base_t::metadata_t;

    /// Type of associated data proxy being created.
    using typename base_t::aux_collection_proxy_t;

    /// Type of _art_ association being used as input.
    using typename base_t::assns_t;

    /**
     * @brief Create a association proxy collection using main collection tag.
     * @tparam Event type of the event to read associations from
     * @tparam Handle type of data product handle
     * @tparam MainArgs any type convertible to `art::InputTag`
     * @tparam Args optional single type (`art::InputTag` required)
     * @param event the event to read associations from
     * @param mainHandle handle of the main collection data product
     * @param margs an object describing the main data product
     * @param args input tag for associated data, if different from main
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

  }; // struct OneTo01DataProxyMaker<>

  /// @}
  // --- END LArSoftProxiesAssociatedData --------------------------------------

  namespace details {

    //--------------------------------------------------------------------------
    //---  Stuff for one-to-(zero or one) associated data
    //--------------------------------------------------------------------------
    template <typename Aux, typename Metadata /* = void */, typename AuxTag /* = Aux */
              >
    struct OneTo01DataProxyMakerWrapper {
      template <typename CollProxy>
      using maker_t =
        OneTo01DataProxyMaker<typename CollProxy::main_element_t, Aux, Metadata, CollProxy, AuxTag>;
    };

  } // namespace details

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_ONETO01DATAPROXYMAKER_H
