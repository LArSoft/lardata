/**
 * @file   lardata/RecoBaseProxy/ProxyBase/getCollection.h
 * @brief  Creation of a collection proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_GETCOLLECTION_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_GETCOLLECTION_H

// LArSoft libraries
#include "lardata/RecoBaseProxy/ProxyBase/CollectionProxyMaker.h"

// C/C++ standard
#include <utility> // std::forward()

namespace proxy {

  // ---------------------------------------------------------------------------
  /**
   * @brief Creates a proxy to a data product collection.
   * @tparam CollProxy type of target main collection proxy
   * @tparam Event type of event to read data from
   * @tparam OptionalArgs type of optional arguments
   * @param event event to read data from
   * @param optionalArgs optional arguments for construction of the proxy
   * @return a collection proxy object
   * @ingroup LArSoftProxyBase
   * @see @ref LArSoftProxyBase "ways to merge more data into a proxy"
   *
   * This function delivers a collection proxy related to `CollProxy`.
   *
   * The @ref LArSoftProxyQuirks "type of proxy delivered is arbitrary"
   * and usually not `CollProxy`.
   * The type of the collection proxy must be explicitly specified, e.g.:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto tracks = proxy::getCollection<proxy::Tracks>
   *   (event, tag, withAssociated<recob::Hits>());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this case, two optional arguments are passed: the input tag to the main
   * collection, and then `withAssociated<recob::Hits>()`. The meaning of both
   * is decided depending on the collection proxy being created, but it's common
   * to have the first argument be the input tag to the main collection, as in
   * this example.
   * `withAssociated()` is one of the ways for a proxy to have
   * @ref LArSoftProxyDefinitionAuxiliaryData "auxiliary data" "merged" into.
   * The options to @ref LArSoftProxyDefinitionMerging "merge" this data are
   * collected in the @ref LArSoftProxyBase "proxy interface documentation".
   *
   * The collection proxy name is arbitrary, but it's custom to have it live in
   * `proxy` namespace and have the same name as the base object, made plural:
   * a proxy to a `recob::Track` collection data product will have a proxy
   * called `proxy::Tracks`.
   *
   * Note that a proxy need to be explicitly supported in order to be available.
   * Nevertheless, a generic implementation is supported to create a proxy of a
   * data product which is a C++ vector, so that:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto tracks = proxy::getCollection<std::vector<recob::Track>>
   *   (event, tag, withAssociated<recob::Hits>());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will have an outcome similar to the previous example. In this case, though,
   * all the specific track interface that went into `proxy::Tracks` proxy will
   * not be available.
   *
   * The implementation of this feature is documented in
   * @ref LArSoftProxyCollections "its own doxygen module".
   *
   *
   * Customization
   * ==============
   *
   * To control which type of collection proxy is produced for the type
   * `CollProxy`, the class `CollectionProxyMaker` may be specialised.
   *
   */
  template <typename CollProxy, typename Event, typename... OptionalArgs>
  auto getCollection(Event const& event, OptionalArgs&&... optionalArgs)
  {
    return CollectionProxyMaker<CollProxy>::make(event,
                                                 std::forward<OptionalArgs>(optionalArgs)...);
  }

  // ---------------------------------------------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_GETCOLLECTION_H
