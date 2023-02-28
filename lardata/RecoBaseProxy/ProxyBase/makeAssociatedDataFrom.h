/**
 * @file   lardata/RecoBaseProxy/ProxyBase/makeAssociatedDataFrom.h
 * @brief  Helper functions to create an `AssociatedData` object.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase/withAssociated.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_MAKEASSOCIATEDDATAFROM_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_MAKEASSOCIATEDDATAFROM_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...
#include "lardata/RecoBaseProxy/ProxyBase/AssociatedData.h"

// framework libraries
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"

// C/C++ standard
#include <cstdlib> // std::size_t
#include <utility> // std::forward()

namespace proxy {

  // --- BEGIN Associated data infrastructure ----------------------------------
  /// @addtogroup LArSoftProxiesAssociatedData
  /// @{

  template <typename Tag, typename Assns>
  auto makeAssociatedDataFrom(Assns const& assns, std::size_t minSize = 0)
  {
    return proxy::makeAssociatedData<Tag>(assns, minSize);
  }

  template <typename Assns>
  auto makeAssociatedDataFrom(Assns const& assns, std::size_t minSize = 0)
  {
    return makeAssociatedDataFrom<typename Assns::right_t>(assns, minSize);
  }

  /**
   * @brief Creates and returns an associated data object.
   * @tparam Main type of main object to be associated
   * @tparam Aux type of data to be associated to the main objects
   * @tparam Metadata type of metadata in the association (if omitted: `void`)
   * @tparam Tag the tag labelling this associated data (if omitted: `Aux`)
   * @tparam Event type of event to read associations from
   * @param event event to read associations from
   * @param tag input tag of the association object
   * @param minSize minimum number of entries in the produced association data
   * @return a new `AssociatedData` filled with associations from `tag`
   *
   * The association being retrieved must fulfill the requirements of
   * @ref LArSoftProxyDefinitionOneToManySeqAssn "one-to-many sequential association".
   *
   * Elements in the main collection not associated with any object will be
   * recorded as such. If there is information for less than `minSize` main
   * objects, more records will be added to mark the missing objects as not
   * associated to anything.
   *
   * Two template types must be explicitly specified, e.g.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto assData = makeAssociatedDataFrom<recob::Track, recob::Hit>(event, tag);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename Main, typename Aux, typename Metadata, typename Tag, typename Event>
  auto makeAssociatedDataFrom(Event const& event,
                              art::InputTag const& tag,
                              std::size_t minSize = 0);

  template <typename Main, typename Aux, typename Metadata, typename Event>
  auto makeAssociatedDataFrom(Event const& event, art::InputTag const& tag, std::size_t minSize = 0)
  {
    return makeAssociatedDataFrom<Main, Aux, Metadata, Aux, Event>(event, tag, minSize);
  }

  template <typename Main, typename Aux, typename Event>
  auto makeAssociatedDataFrom(Event const& event, art::InputTag const& tag, std::size_t minSize = 0)
  {
    return makeAssociatedDataFrom<Main, Aux, void, Aux, Event>(event, tag, minSize);
  }

  /**
   * @brief Creates and returns an associated data object.
   * @tparam Aux type of data to be associated to the main objects
   * @tparam Metadata type of metadata in the association (if omitted: `void`)
   * @tparam Tag the tag labelling this associated data (if omitted: `Aux`)
   * @tparam Handle type of handle to the main collection object
   * @tparam Event type of event to read associations from
   * @param handle handle to the main collection object
   * @param event event to read associations from
   * @param tag input tag of the association object
   * @return a new `AssociatedData` filled with associations from `tag`
   * @see `makeAssociatedDataFrom(Event const&, art::InputTag, std::size_t)`
   *
   * This function operates like
   * `makeAssociatedDataFrom(Event const&, art::InputTag, std::size_t)`, but it
   * extracts the information about the type of main object and the minimum
   * number of them from a handle.
   * The handle object is expected to behave as a smart pointer to a
   * collection of elements of the associated type.
   *
   * One template type must be explicitly specified, e.g.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto assData = makeAssociatedDataFrom<recob::Hit>(handle, event, tag);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename Aux, typename Metadata, typename Tag, typename Handle, typename Event>
  auto makeAssociatedDataFrom(Handle&& handle, Event const& event, art::InputTag const& tag);

  template <typename Aux, typename Metadata, typename Handle, typename Event>
  auto makeAssociatedDataFrom(Handle&& handle, Event const& event, art::InputTag const& tag)
  {
    return makeAssociatedDataFrom<Aux, Metadata, Aux, Handle, Event>(
      std::forward<Handle>(handle), event, tag);
  }

  template <typename Aux, typename Handle, typename Event>
  auto makeAssociatedDataFrom(Handle&& handle, Event const& event, art::InputTag const& tag)
  {
    return makeAssociatedDataFrom<Aux, void, Aux, Handle, Event>(
      std::forward<Handle>(handle), event, tag);
  }

  /**
   * @brief Creates and returns an associated data object.
   * @tparam Tag the tag labelling this associated data (if omitted: `Aux`)
   * @tparam MainColl type of the main collection object
   * @tparam Assns type of the association object
   * @param mainColl the main collection object
   * @param assns association data object
   * @return a new `AssociatedData` wrapping the information in `assns`
   * @see `makeAssociatedDataFrom(Assns const&, std::size_t)`
   *
   * This function operates like
   * `makeAssociatedDataFrom(Assns const&, std::size_t)`, where the size is
   * extracted from the main data collection.
   */
  template <typename Tag, typename MainColl, typename Assns>
  auto makeAssociatedDataFrom(MainColl const& mainColl, Assns const& assns)
  {
    return proxy::makeAssociatedData<Tag>(assns, mainColl.size());
  }

  template <typename MainColl, typename Assns>
  auto makeAssociatedDataFrom(MainColl const& mainColl, Assns const& assns)
  {
    return makeAssociatedDataFrom<typename Assns::right_t>(mainColl, assns);
  }

  //----------------------------------------------------------------------------

  /// @}
  // --- END Associated data infrastructure ------------------------------------

} // namespace proxy

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace proxy {

  //----------------------------------------------------------------------------
  //--- makeAssociatedDataFrom() implementations
  //----------------------------------------------------------------------------
  template <typename Main, typename Aux, typename Metadata, typename Tag, typename Event>
  auto makeAssociatedDataFrom(Event const& event,
                              art::InputTag const& tag,
                              std::size_t minSize /* = 0 */
  )
  {
    using Main_t = Main;
    using Aux_t = Aux;
    using Metadata_t = Metadata;
    using AssociatedData_t = details::AssociatedData<Main_t, Aux_t, Metadata_t, Tag>;
    using Assns_t = typename AssociatedData_t::assns_t;

    return makeAssociatedDataFrom<Tag>(*(event.template getValidHandle<Assns_t>(tag)));

  } // makeAssociatedDataFrom(tag)

  //----------------------------------------------------------------------------
  template <typename Aux, typename Metadata, typename Tag, typename Handle, typename Event>
  auto makeAssociatedDataFrom(Handle&& handle, Event const& event, art::InputTag const& tag)
  {
    // Handle::value_type is the main data product type (a collection)
    using Main_t = util::collection_value_t<typename Handle::value_type>;
    using Aux_t = Aux;
    using Metadata_t = Metadata;
    return makeAssociatedDataFrom<Main_t, Aux_t, Metadata_t, Tag>(event, tag, handle->size());
  } // makeAssociatedDataFrom(handle)

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_MAKEASSOCIATEDDATAFROM_H
