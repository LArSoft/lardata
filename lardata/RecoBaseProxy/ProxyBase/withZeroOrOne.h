/**
 * @file   lardata/RecoBaseProxy/ProxyBase/withZeroOrOne.h
 * @brief  Interface to add optional associated data to a collection proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_WITHZEROORONE_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_WITHZEROORONE_H

// LArSoft libraries
#include "lardata/RecoBaseProxy/ProxyBase/OneTo01DataProxyMaker.h"
#include "lardata/RecoBaseProxy/ProxyBase/WithAssociatedStructBase.h"

// C/C++ standard libraries
#include <tuple>
#include <utility> // std::forward(), std::move()

namespace proxy {

  namespace details {

    template <typename Aux, typename Metadata, typename ArgTuple, typename AuxTag = Aux>
    using WithOneTo01AssociatedStruct = WithAssociatedStructBase<
      Aux,
      Metadata,
      ArgTuple,
      OneTo01DataProxyMakerWrapper<Aux, Metadata, AuxTag>::template maker_t,
      AuxTag>;

  } // namespace details

  // --- BEGIN One-to-one (optional) associations ------------------------------
  /**
   * @name One-to-one (optional) associations
   *
   * These functions allow to merge into a data collection proxy some auxiliary
   * data via an _art_ association fulfilling the
   * @ref LArSoftProxyDefinitionOneToZeroOrOneSeqAssn "one-to-(zero-or-one) sequential association requirement".
   *
   * One category of functions is currently available:
   *  * `proxy::withZeroOrOne()` reads the relevant association from an event
   *
   * Variants of `proxy::withZeroOrOne()` called `proxy::withZeroOrOneMeta()`
   * will allow merging the metadata of an association too. This feature is not
   * supported yet, though.
   *
   * Also, variants are available to customize the tag class.
   *
   * The implementation of this feature is documented in
   * @ref LArSoftProxiesAssociatedData "its own doxygen module".
   *
   * @{
   */

  //----------------------------------------------------------------------------
  /// The same as `withZeroOrOneMeta()`, but it also specified a tag.
  /// @ingroup LArSoftProxyBase
  /// @todo Metadata is not supported yet.
  template <typename Aux, typename Metadata, typename AuxTag, typename... Args>
  auto withZeroOrOneMetaAs(Args&&... args)
  {
    using ArgTuple_t = std::tuple<Args&&...>;
    ArgTuple_t argsTuple(std::forward<Args>(args)...);
    return details::WithOneTo01AssociatedStruct<Aux, Metadata, ArgTuple_t, AuxTag>(
      std::move(argsTuple));
  } // withZeroOrOneAs()

  /// The same as `withZeroOrOne()`, but it also specified a tag for the data.
  /// @ingroup LArSoftProxyBase
  template <typename Aux, typename AuxTag, typename... Args>
  auto withZeroOrOneAs(Args&&... args)
  {
    return withZeroOrOneMetaAs<Aux, void, AuxTag>(std::forward<Args>(args)...);
  }

  /**
   * @brief Helper function to merge one-to-(zero-or-one) associated data.
   * @tparam Aux type of associated data requested
   * @tparam Metadata type of metadata coming with the associated data
   * @tparam Args types of constructor arguments for associated data collection
   * @param args constructor arguments for the associated data collection
   * @return a temporary object that `getCollection()` knows to handle
   * @ingroup LArSoftProxyBase
   * @see withZeroOrOneMetaAs(), withZeroOrOne()
   *
   * This function is meant to convey to `getCollection()` function the request
   * for the delivered collection proxy to carry auxiliary data from an
   * association fulfilling the
   * @ref LArSoftProxyDefinitionOneToManySeqAssn "one-to-many sequential association"
   * requirements.
   *
   * This data will be tagged with the type `Aux`. To use a different type as
   * tag, use `withZeroOrOneAs()` instead, specifying the tag as second
   * template argument, e.g.:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct QuestionableVertex {};
   * auto tracks = proxy::getCollection<proxy::Tracks>(event, trackTag,
   *   withZeroOrOneMeta<recob::Vertex, void>(defaultVertexTag),
   *   withZeroOrOneMetaAs<recob::Vertex, void, QuestionableVertex>
   *     (stinkyVertexTag)
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * and, since we are not requesting any metadata, this is equivalent to
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct QuestionableVertex {};
   * auto tracks = proxy::getCollection<proxy::Tracks>(event, trackTag,
   *   withZeroOrOne<recob::Vertex>(defaultVertexTag),
   *   withZeroOrOneAs<recob::Vertex, QuestionableVertex>(stinkyVertexTag)
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * The first vertex association (`"defaultVertexTag"`) will be accessed by
   * using the type `recob::Vertex` as tag, while the second one will be
   * accessed by the `QuestionableVertex` tag (which is better not be defined
   * in a local scope):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * for (auto&& track: tracks) {
   *   decltype(auto) vertex = track.get<recob::Vertex>();
   *   decltype(auto) maybeVertex = track.get<QuestionableVertex>();
   *   // ...
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   *
   * Customization of the association proxy
   * =======================================
   *
   * See the technical details about `withAssociated()`, which hold for this
   * function and related classes too.
   *
   *
   * Technical details
   * ==================
   *
   * See the technical details about `withAssociated()`, which hold for this
   * function and related classes too.
   *
   * @todo Metadata is not supported yet (the interface is apparently there though).
   */
  template <typename Aux, typename Metadata, typename... Args>
  auto withZeroOrOneMeta(Args&&... args)
  {
    return withZeroOrOneMetaAs<Aux, Metadata, Aux>(std::forward<Args>(args)...);
  }

  /// Works like `withZeroOrOneMeta()`, but for associations with no metadata.
  /// @ingroup LArSoftProxyBase
  /// @see withZeroOrOneAs(), withZeroOrOneMeta()
  template <typename Aux, typename... Args>
  auto withZeroOrOne(Args&&... args)
  {
    return withZeroOrOneMeta<Aux, void>(std::forward<Args>(args)...);
  }

  /// @}
  // --- END One-to-one (optional) associations --------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_WITHZEROORONE_H
