/**
 * @file   lardata/RecoBaseProxy/ProxyBase/ProxyAsParallelData.h
 * @brief  Data encapsulating a collection proxy as auxiliary data.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase/withCollectionProxy.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_PROXYASPARALLELDATA_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_PROXYASPARALLELDATA_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...
#include "lardata/RecoBaseProxy/ProxyBase/ParallelData.h"

// C/C++ standard
#include <utility> // std::move()

namespace proxy {

  namespace details {

    template <typename AuxProxyColl,
              typename Aux = util::collection_value_t<AuxProxyColl>,
              typename Tag = Aux>
    struct ProxyAsParallelData;

  } // namespace details

  // --- BEGIN Infrastructure for proxies as auxiliary data --------------------
  /**
   * @defgroup LArSoftProxiesAuxProxy Infrastructure for proxies as auxiliary
   *           data.
   * @ingroup LArSoftProxyCustom
   * @brief Infrastructure to use a collection proxy as auxiliary data for
   *        another proxy.
   * @bug Broken in many ways. Do not use.
   *
   * @{
   */

  template <typename Tag /* = Aux */,
            typename Aux /* = util::collection_value_t<AuxProxyColl>*/,
            typename AuxProxyColl>
  auto makeProxyAsParallelData(AuxProxyColl&& auxProxy)
  {
    return details::ProxyAsParallelData<AuxProxyColl, Aux, Tag>(std::move(auxProxy));
  } // makeProxyAsParallelData()

  /// @}
  /// --- END Infrastructure for proxies as auxiliary data ---------------------

  //----------------------------------------------------------------------------
  namespace details {

    //--------------------------------------------------------------------------
    /**
     * @brief Object presenting a proxy as parallel data for another one.
     * @tparam AuxProxyColl type of the parallel data collection
     * @tparam Aux type of the associated object
     * @tparam Tag tag this data is labeled with
     * @ingroup LArSoftProxiesAuxProxy
     *
     * This object inherits its interface from `proxy::ParallelData`.
     * In addition, it owns the proxy it wraps.
     */
    template <typename AuxProxyColl,
              typename Aux /* = util::collection_value_t<AuxProxyColl> */,
              typename Tag /* = Aux */
              >
    struct ProxyAsParallelData : private AuxProxyColl, public ParallelData<AuxProxyColl, Aux, Tag> {
      /// Steals and wraps collection `proxy`.
      ProxyAsParallelData(AuxProxyColl&& proxy)
        : AuxProxyColl(std::move(proxy))
        , ParallelData<AuxProxyColl, Aux, Tag>(static_cast<AuxProxyColl const*>(this))
      {}

      // explicitly select the tag from the parallel data (same as Tag)
      using typename ParallelData<AuxProxyColl, Aux, Tag>::tag;

    }; // class ProxyAsParallelData<>

    //--------------------------------------------------------------------------

  } // namespace details

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_PROXYASPARALLELDATA_H
