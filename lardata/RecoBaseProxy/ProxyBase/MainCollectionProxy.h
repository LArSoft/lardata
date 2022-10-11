/**
 * @file   lardata/RecoBaseProxy/ProxyBase/MainCollectionProxy.h
 * @brief  Utilities for the main collection of a collection proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_MAINCOLLECTIONPROXY_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_MAINCOLLECTIONPROXY_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...

// C/C++ standard
#include <cstdlib> // std::size_t

namespace proxy {

  //----------------------------------------------------------------------------
  namespace details {

    //--------------------------------------------------------------------------
    //--- stuff for the main collection
    //--------------------------------------------------------------------------
    /**
     * @brief Wrapper for the main collection of a proxy.
     * @tparam MainColl type of the collection being wrapped
     *
     * The wrapper contains a pointer to the original collection, which must
     * persist. The original collection is not modified.
     *
     * The `MainColl` type must expose a random access container interface.
     */
    template <typename MainColl>
    struct MainCollectionProxy {

      /// Type of the original collection.
      using main_collection_t = MainColl;

      /// Type of the elements in the original collection.
      using main_element_t = util::collection_value_t<MainColl>;

      /// Constructor: wraps the specified collection.
      MainCollectionProxy(main_collection_t const& main) : fMain(&main) {}

      /// Returns the wrapped collection.
      main_collection_t const& main() const { return mainRef(); }

      /// Returns a reference to the wrapped collection.
      main_collection_t const& mainRef() const { return *fMain; }

      /// Returns a pointer to the wrapped collection.
      main_collection_t const* mainPtr() const { return fMain; }

    protected:
      /// This type.
      using this_t = MainCollectionProxy<main_collection_t>;

      /// Return this object as main collection proxy.
      this_t& mainProxy() { return *this; }

      /// Return this object as main collection proxy.
      this_t const& mainProxy() const { return *this; }

      /// Returns the specified item in the original collection.
      auto getMainAt(std::size_t i) const -> decltype(auto) { return fMain->operator[](i); }

    private:
      main_collection_t const* fMain; /// Pointer to the original collection.

    }; // struct MainCollectionProxy<>

    //--------------------------------------------------------------------------

  } // namespace details

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_MAINCOLLECTIONPROXY_H
