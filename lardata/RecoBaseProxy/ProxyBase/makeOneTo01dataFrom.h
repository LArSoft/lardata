/**
 * @file   lardata/RecoBaseProxy/ProxyBase/makeOneTo01dataFrom.h
 * @brief  Helper functions to create data structures associated to a proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase/withZeroOrOne.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_MAKEONETO01DATAFROM_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_MAKEONETO01DATAFROM_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t
#include "lardata/RecoBaseProxy/ProxyBase/OneTo01Data.h"

// framework libraries
#include "canvas/Utilities/InputTag.h"

// C/C++ standard libraries
#include <cstdlib> // std::size_t
#include <utility> // std::forward()

namespace proxy {

  // --- BEGIN LArSoftProxiesAssociatedData ------------------------------------
  /// @addtogroup LArSoftProxiesAssociatedData
  /// @{

  //----------------------------------------------------------------------------
  //--- one-to-(zero/one) associations
  //----------------------------------------------------------------------------
  /**
   * @brief Processes and returns an one-to-(zero/one) associated data object.
   * @tparam Tag the tag labelling this associated data
   *             (if omitted: second type of the association: `right_t`)
   * @tparam Assns type of association to be processed
   * @param assns association object to be processed
   * @param minSize minimum number of entries in the produced association data
   * @return a new `OneTo01Data` filled with associations from `tag`
   *
   * The content of the association object must fulfill the requirements of
   * @ref LArSoftProxyDefinitionOneToZeroOrOneSeqAssn "one-to-(zero or one) sequential association".
   * The `Assns` type is expected to be a `art::Assns` instance. At least,
   * the `Assns` type is required to have `left_t` and `right_t` definitions
   * representing respectively the main data type and the associated one, and
   * to respond to `begin()` and `end()` functions. The iterated object must
   * also respond to `std::get<0>()` with a `art::Ptr<left_t>` and to
   * `std::get<1>()` with a `art::Ptr<right_t>`.
   *
   * Elements in the main collection not associated with any object will present
   * an invalid _art_ pointer (`isNull()` true). If there is information for
   * less than `minSize` main objects, more records will be added to mark the
   * missing objects as not associated to anything.
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * art::Assns<recob::Track, recob::Vertex> trackVertexAssns;
   * // ...
   * auto assData = makeOneTo01dataFrom(assns);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will have `assData` tagged as `recob::Vertex`.
   */
  template <typename Tag, typename Assns>
  auto makeOneTo01dataFrom(Assns const& assns, std::size_t minSize = 0)
  {
    return proxy::makeOneTo01data<Tag>(assns, minSize);
  }

  template <typename Assns>
  auto makeOneTo01dataFrom(Assns const& assns, std::size_t minSize = 0)
  {
    return proxy::makeOneTo01data(assns, minSize);
  }

  /**
   * @brief Creates and returns an one-to-(zero/one) associated data object.
   * @tparam Main type of main object to be associated
   * @tparam Aux type of data to be associated to the main objects
   * @tparam Metadata type of metadata in the association
   * @tparam Tag the tag labelling this associated data (if omitted: `Aux`)
   * @tparam Event type of event to read associations from
   * @param event event to read associations from
   * @param tag input tag of the association object
   * @param minSize minimum number of entries in the produced association data
   * @return a new `OneTo01Data` filled with associations from `tag`
   * @see `makeOneTo01dataFrom(Assns, std::size_t)`
   *
   * The association being retrieved must fulfill the requirements of
   * @ref LArSoftProxyDefinitionOneToZeroOrOneSeqAssn "one-to-(zero or one) sequential association".
   *
   * Two template types must be explicitly specified, e.g.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto assData = makeOneTo01dataFrom<recob::Track, recob::Vertex>(event, tag);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename Main, typename Aux, typename Metadata, typename Tag, typename Event>
  auto makeOneTo01dataFrom(Event const& event, art::InputTag const& tag, std::size_t minSize = 0);

  template <typename Main, typename Aux, typename Metadata, typename Event>
  auto makeOneTo01dataFrom(Event const& event, art::InputTag const& tag, std::size_t minSize = 0)
  {
    return makeOneTo01dataFrom<Main, Aux, Metadata, Aux, Event>(event, tag, minSize);
  }

  /**
   * @brief Creates and returns an one-to-(zero/one) associated data object.
   * @tparam Aux type of data to be associated to the main objects
   * @tparam Tag the tag labelling this associated data (if omitted: `Aux`)
   * @tparam Handle type of handle to the main collection object
   * @tparam Event type of event to read associations from
   * @param handle handle to the main collection object
   * @param event event to read associations from
   * @param tag input tag of the association object
   * @return a new `OneTo01Data` wrapping the information in `assns`
   * @see `makeAssociatedDataFrom(Event const&, art::InputTag, std::size_t)`
   *
   * This function operates like
   * `makeOneTo01dataFrom(Event const&, art::InputTag, std::size_t)`, but it
   * extracts the information about the type of main object and the minimum
   * number of them from a handle.
   * The handle object is expected to behave as a smart pointer to a
   * collection of elements of the associated type.
   *
   * One template type must be explicitly specified, e.g.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto assData = makeOneTo01dataFrom<recob::Vertex>(handle, event, tag);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename Aux, typename Metadata, typename Tag, typename Handle, typename Event>
  auto makeOneTo01dataFrom(Handle&& handle, Event const& event, art::InputTag const& tag);

  template <typename Aux, typename Metadata, typename Handle, typename Event>
  auto makeOneTo01dataFrom(Handle&& handle, Event const& event, art::InputTag const& tag)
  {
    return makeOneTo01dataFrom<Aux, Metadata, Aux, Handle, Event>(
      std::forward<Handle>(handle), event, tag);
  }

  /**
   * @brief Creates and returns an one-to-(zero/one) associated data object.
   * @tparam Tag the tag labelling this associated data (if omitted: `Aux`)
   * @tparam MainColl type of the main collection object
   * @tparam Assns type of the association object
   * @param mainColl the main collection object
   * @param assns association data object
   * @return a new `OneTo01Data` wrapping the information in `assns`
   * @see `makeOneTo01dataFrom(Assns const&, std::size_t)`
   *
   * This function operates like
   * `makeOneTo01dataFrom(Assns const&, std::size_t)`, where the size is
   * extracted from the main data collection.
   */
  template <typename Tag, typename MainColl, typename Assns>
  auto makeOneTo01dataFrom(MainColl const& mainColl, Assns const& assns)
  {
    return proxy::makeOneTo01data<Tag>(assns, mainColl.size());
  }

  template <typename MainColl, typename Assns>
  auto makeOneTo01dataFrom(MainColl const& mainColl, Assns const& assns)
  {
    return proxy::makeOneTo01data<typename Assns::right_t>(mainColl, assns);
  }

  /// @}
  // --- END LArSoftProxiesAssociatedData --------------------------------------

} // namespace proxy

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace proxy {

  //----------------------------------------------------------------------------
  //--- makeOneTo01dataFrom() implementations
  //----------------------------------------------------------------------------
  template <typename Main, typename Aux, typename Metadata, typename Tag, typename Event>
  auto makeOneTo01dataFrom(Event const& event,
                           art::InputTag const& tag,
                           std::size_t minSize /* = 0 */
  )
  {
    using Main_t = Main;
    using Aux_t = Aux;
    using Metadata_t = Metadata;
    using AssociatedData_t = details::OneTo01Data<Main_t, Aux_t, Metadata_t, Tag>;
    using Assns_t = typename AssociatedData_t::assns_t;

    return makeOneTo01dataFrom<Tag>(*(event.template getValidHandle<Assns_t>(tag)), minSize);

  } // makeOneTo01dataFrom(tag)

  //----------------------------------------------------------------------------
  template <typename Aux, typename Metadata, typename Tag, typename Handle, typename Event>
  auto makeOneTo01dataFrom(Handle&& handle, Event const& event, art::InputTag const& tag)
  {
    // Handle::value_type is the main data product type (a collection)
    using Main_t = util::collection_value_t<typename Handle::value_type>;
    using Aux_t = Aux;
    using Metadata_t = Metadata;
    return makeOneTo01dataFrom<Main_t, Aux_t, Metadata_t, Tag>(event, tag, handle->size());
  } // makeOneTo01dataFrom(handle)

  //----------------------------------------------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_MAKEONETO01DATAFROM_H
