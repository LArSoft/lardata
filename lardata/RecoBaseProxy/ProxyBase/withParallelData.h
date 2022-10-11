/**
 * @file   lardata/RecoBaseProxy/ProxyBase/withParallelData.h
 * @brief  Interface to add auxiliary data from parallel collections to a proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_WITHPARALLELDATA_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_WITHPARALLELDATA_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...
#include "lardata/RecoBaseProxy/ProxyBase/ParallelDataProxyMaker.h"
#include "lardata/RecoBaseProxy/ProxyBase/WithAssociatedStructBase.h"

// C/C++ standard
#include <tuple>
#include <utility> // std::forward(), std::move()

namespace proxy {

  //----------------------------------------------------------------------------
  namespace details {

    template <typename Aux, typename ArgTuple, typename AuxTag = Aux>
    using WithParallelCollectionStruct =
      WithAssociatedStructBase<Aux,
                               void, // no metadata concept for parallel collections
                               ArgTuple,
                               ParallelDataProxyMakerWrapper<Aux, AuxTag>::template maker_t,
                               AuxTag>;

    template <typename Aux, typename ArgTuple, typename AuxColl, typename AuxTag = Aux>
    using WithWrappedParallelCollectionStruct = WithAssociatedStructBase<
      Aux,
      void, // no metadata concept for parallel collections
      ArgTuple,
      ParallelDataProxyMakerWrapper<Aux, AuxTag, AuxColl>::template maker_t,
      AuxTag>;

  } // namespace details

  // --- BEGIN Parallel data collections ---------------------------------------
  /**
   * @name Parallel data collections
   *
   * These functions allow to merge into a data collection proxy some auxiliary
   * data from other collections fulfilling the
   * @ref LArSoftProxyDefinitionParallelData "parallel data product requirement".
   *
   * Two categories of functions are available depending on the data source:
   *  * `proxy::withParallelData()` reads the relevant data product from an
   *      event
   *  * `proxy::wrapParallelData()` uses an existing collection
   *
   * Also, variants are available to customize the tag class.
   *
   * The implementation of this feature is documented in
   * @ref LArSoftProxiesParallelData "its own doxygen module".
   *
   * @{
   */

  //----------------------------------------------------------------------------
  /**
   * @brief Helper function to merge an auxiliary data product into the proxy.
   * @tparam Aux type of parallel data product requested
   * @tparam AuxTag the tag type to refer this auxiliary data as
   * @tparam Args types of constructor arguments for parallel data proxy
   * @param args constructor arguments for the parallel data collection proxy
   * @return a temporary object that `getCollection()` knows to handle
   * @ingroup LArSoftProxyBase
   * @see `proxy::withParallelData()`
   *
   * This function is meant to convey to `getCollection()` function the request
   * for merging a auxiliary data structured as a
   * @ref LArSoftProxyDefinitionParallelData "collection parallel" into the
   * collection proxy.
   *
   * It is functionally equivalent to `withParallelData()`, with the difference
   * that here the auxiliary data tag must be specified. `withParallelData()`
   * documentation also contains examples on how to use this function and the
   * proxy resulting from that.
   */
  template <typename Aux, typename AuxTag, typename... Args>
  auto withParallelDataAs(Args&&... args)
  {
    using ArgTuple_t = std::tuple<Args&&...>;
    ArgTuple_t argsTuple(std::forward<Args>(args)...);
    return details::WithParallelCollectionStruct<Aux, ArgTuple_t, AuxTag>(std::move(argsTuple));
  } // withParallelDataAs()

  //----------------------------------------------------------------------------
  /**
   * @brief Helper function to merge an auxiliary data product into the proxy.
   * @tparam Aux type of parallel data product requested
   * @tparam Args types of constructor arguments for parallel data proxy
   * @param args constructor arguments for the parallel data collection proxy
   * @return a temporary object that `getCollection()` knows to handle
   * @ingroup LArSoftProxyBase
   *
   * This function is meant to convey to `getCollection()` function the request
   * to merge auxiliary data structured as a
   * @ref LArSoftProxyDefinitionParallelData "collection parallel" into the
   * collection proxy.
   *
   * This data will be tagged with the type `Aux`. To use a different type as
   * tag, use `withParallelDataAs()` instead, specifying the tag as second
   * template argument, e.g.:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct MCS {};
   * auto tracks = proxy::getCollection<proxy::Tracks>(event, trackTag,
   *   withParallelData<recob::TrackMomentum>(defaultMomTag),
   *   withParallelDataAs<recob::TrackMomentum, MCS>(MCSmomTag)
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The first momentum association (`"defaultMomTag"`) will be accessed by
   * using the type `recob::TrackMomentum` as tag, while the second one will
   * be accessed by the `MCS` tag (which is better not be defined in a local
   * scope):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * for (auto&& track: tracks) {
   *   decltype(auto) trackMom = track.get<recob::TrackMomentum>();
   *   decltype(auto) trackMCSmom = track.get<MCS>();
   *   // ...
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * The default implementation of parallel data proxy returns for each element
   * query an object with the same interface as the element of the parallel data
   * collection. In the previous examples, that would be a constant reference to
   * an object with `recob::TrackMomentum` interface.
   *
   *
   * Customization of the parallel data proxy
   * =========================================
   *
   * To have a call like:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto tracks = getCollection<SpecialTracks>
   *   (event, tag, withParallelData<recob::TrackMomentum>(momTag, "special"));
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * to create something different than the standard parallel data proxy,
   * one needs to specialize `proxy::ParallelDataProxyMaker`, e.g.:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * namespace proxy {
   *   template <>
   *   struct ParallelDataProxyMaker<recob::Track, recob::TrackMomentum, SpecialTracks>
   *     : public ParallelDataProxyMakerBase<recob::Track, std::vector<recob::TrackMomentum>, recob::TrackMomentum>
   *   {
   *
   *     template<typename Event, typename MainArgs>
   *     static auto make(
   *       Event const& event,
   *       MainArgs const&,
   *       art::InputTag assnTag,
   *       std::string quality = "default"
   *       )
   *       {
   *         ::SpecialTrackHitsProxy myAuxProxy;
   *         // ... make it, and make it right
   *         return myAuxProxy;
   *       }
   *
   *   }; // struct ParallelDataProxyMaker<recob::Track, recob::TrackMomentum, SpecialTracks>
   *
   * } // namespace proxy
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   */
  template <typename Aux, typename... Args>
  auto withParallelData(Args&&... args)
  {
    return withParallelDataAs<Aux, Aux>(std::forward<Args>(args)...);
  }

  //----------------------------------------------------------------------------
  /**
   * @brief Uses a collection as auxiliary data for a collection proxy.
   * @tparam AuxTag the tag type to refer this auxiliary data as
   * @tparam AuxColl type of the auxiliary data collection
   * @param auxColl the data collection to be used as auxiliary data
   * @return an object making `getCollection()` add such auxiliary data
   * @ingroup LArSoftProxyBase
   * @see `withParallelDataAs()`, `wrapParallelData()`
   *
   * The specified collection is used directly as auxiliary data in a collection
   * proxy. It is required to fulfil the
   * @ref LArSoftProxyDefinitionParallelData "parallel data product"
   * requirements, but it does not have to actually be a data product (that is,
   * it does not have to be a collection read from _art_).
   *
   * The usage of the resulting proxy is the same as the ones created using
   * `proxy::withParallelDataAs()`, but the object `auxColl` must remain valid
   * as long as that proxy is being used.
   * Example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct MCS {};
   *
   * void checkMomenta(
   *   std::vector<recob::TrackMomentum> const& mom,
   *   std::vector<recob::TrackMomentum> const& MCSmom
   * ) {
   *   auto tracks = proxy::getCollection<proxy::Tracks>(event, trackTag,
   *     wrapParallelData<recob::TrackMomentum>(defaultMomTag),
   *     wrapParallelDataAs<recob::TrackMomentum, MCS>(MCSmomTag)
   *     );
   *   for (auto const& track: tracks) {
   *     auto const& trackMom = track.get<recob::TrackMomentum>();
   *     auto const& trackMCSmom = track.get<MCS>();
   *
   *     // ...
   *
   *   } // for tracks
   *
   * } // checkMomenta()
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The first momentum (`mom`) will be accessed by using the type
   * `recob::TrackMomentum` as tag, while the second one (`MCSmom`) will be
   * accessed by the `MCS` tag (which is better not be defined in a local
   * scope).
   */
  template <typename AuxTag, typename AuxColl>
  auto wrapParallelDataAs(AuxColl const& auxColl)
  {
    std::tuple<AuxColl const&> args = {auxColl};
    return details::WithWrappedParallelCollectionStruct<util::collection_value_t<AuxColl>,
                                                        decltype(args),
                                                        AuxColl,
                                                        AuxTag>(std::move(args));
  } // wrapParallelDataAs()

  /**
   * @brief Uses a collection as auxiliary data for a collection proxy.
   * @tparam AuxColl type of the auxiliary data collection
   * @param auxColl the data collection to be used as auxiliary data
   * @return an object making `getCollection()` add such auxiliary data
   * @ingroup LArSoftProxyBase
   * @see `proxy::withParallelData()`, `proxy::wrapParallelDataAs()`
   *
   * This function is meant to convey to `getCollection()` function the request
   * for merging an existing auxiliary data collection structured as a
   * @ref LArSoftProxyDefinitionParallelData "collection parallel" into the
   * collection proxy.
   *
   * It is functionally equivalent to `wrapParallelDataAs()`, with the
   * difference that here the auxiliary data tag is automatically defined after
   * the type of the data in the container. `withParallelDataAs()`
   * documentation also contains examples on how to use this function and the
   * proxy resulting from that.
   */
  template <typename AuxColl>
  auto wrapParallelData(AuxColl const& auxColl)
  {
    return wrapParallelDataAs<util::collection_value_t<AuxColl>>(auxColl);
  }

  /// @}
  // --- END Parallel data collections -----------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_WITHPARALLELDATA_H
