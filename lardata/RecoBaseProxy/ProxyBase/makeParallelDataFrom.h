/**
 * @file   lardata/RecoBaseProxy/ProxyBase/makeParallelDataFrom.h
 * @brief  Helper functions to create `proxy::ParallelData` objects.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase/withParallelData.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_MAKEPARALLELDATAFROM_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_MAKEPARALLELDATAFROM_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...
#include "lardata/RecoBaseProxy/ProxyBase/ParallelData.h"

// framework libraries
#include "canvas/Utilities/InputTag.h"

namespace proxy {

  // -- BEGIN Parallel data infrastructure -------------------------------------
  /// @addtogroup LArSoftProxiesParallelData
  /// @{

  /**
   * @brief Wraps a collection into a parallel data collection object.
   * @tparam AuxColl type of parallel data data product container
   * @tparam Aux type of parallel data to be associated to the main objects
   *             (if omitted: `AuxColl::value_type`)
   * @tparam Tag the tag labelling this associated data (if omitted: as `Aux`)
   * @param data data collection to be wrapped
   * @return a new `ParallelData` wrapping the information in `data`
   *
   * The data collection must be non-temporary and it is treated as fulfilling
   * @ref LArSoftProxyDefinitionParallelData "parallel data product"
   * requirements.
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<recob::TrackFitHitInfo> trackData;
   * // ...
   * auto auxData = makeParallelDataFrom(trackData);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * where the tag of the parallel data is now `recob::TrackFitHitInfo` and
   * `auxData` behaviour becomes undefined as soon as `trackData` falls out of
   * scope.
   */
  template <typename AuxColl, typename Aux = util::collection_value_t<AuxColl>, typename Tag = Aux>
  auto makeParallelDataFrom(AuxColl const& data)
  {
    return proxy::makeParallelData<AuxColl, Aux, Tag>(data);
  }

  /**
   * @brief Creates and returns a parallel data collection object.
   * @tparam AuxColl type of parallel data data product container
   * @tparam Aux type of parallel data to be associated to the main objects
   *             (if omitted: `AuxColl::value_type`)
   * @tparam Tag the tag labelling this associated data (if omitted: as `Aux`)
   * @tparam Event type of event to read the data product from
   * @param event event to read the data product from
   * @param tag input tag of the parallel data product
   * @return a new `ParallelData` filled with data from `tag`
   *
   * The data product being retrieved must fulfill the requirements of
   * @ref LArSoftProxyDefinitionParallelData "parallel data product".
   *
   * At least one template type must be explicitly specified, e.g.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto auxData
   *   = makeParallelDataFrom<std::vector<recob::TrackFitHitInfo>>(event, tag);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this case, the `Aux` type is defined as `recob::TrackFitHitInfo`, as is
   * the tag.
   */
  template <typename AuxColl, typename Aux, typename Tag, typename Event>
  auto makeParallelDataFrom(Event const& event, art::InputTag const& tag);

  template <typename AuxColl, typename Aux, typename Event>
  auto makeParallelDataFrom(Event const& event, art::InputTag const& tag)
  {
    return makeParallelDataFrom<AuxColl, Aux, Aux, Event>(event, tag);
  }

  template <typename AuxColl, typename Event>
  auto makeParallelDataFrom(Event const& event, art::InputTag const& tag)
  {
    return makeParallelDataFrom<AuxColl, util::collection_value_t<AuxColl>, Event>(event, tag);
  }

  /// @}
  // -- END Parallel data infrastructure ---------------------------------------

} // namespace proxy

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace proxy {

  //----------------------------------------------------------------------------
  //--- makeParallelDataFrom() implementations
  //----------------------------------------------------------------------------
  template <typename AuxColl, typename Aux, typename Tag, typename Event>
  auto makeParallelDataFrom(Event const& event, art::InputTag const& tag)
  {
    return makeParallelDataFrom<AuxColl, Aux, Tag>(*(event.template getValidHandle<AuxColl>(tag)));
  } // makeParallelDataFrom()

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_MAKEPARALLELDATAFROM_H
