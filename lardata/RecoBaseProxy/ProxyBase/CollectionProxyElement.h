/**
 * @file   lardata/RecoBaseProxy/ProxyBase/CollectionProxyElement.h
 * @brief  Utilities for a single element of a collection proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_COLLECTIONPROXYELEMENT_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_COLLECTIONPROXYELEMENT_H

// LArSoft libraries
#include "larcorealg/CoreUtils/DebugUtils.h"    // lar::debug::demangle()
#include "larcorealg/CoreUtils/MetaUtils.h"     // util::always_true_type, ...
#include "lardata/Utilities/TupleLookupByTag.h" // util::index_of_tag_v, ...

// C/C++ standard
#include <cstdlib>     // std::size_t
#include <stdexcept>   // std::logic_error
#include <tuple>       // also std::tuple_element_t<>, std::get()
#include <type_traits> // std::integral_constant<>
#include <utility>     // std::move()

namespace proxy {

  namespace details {

    template <typename AuxCollTuple>
    struct SubstituteWithAuxList;

  } // namespace details

  //--- BEGIN Proxy element infrastructure -------------------------------------
  /**
   * @defgroup LArSoftProxiesElement Proxy element infrastructure
   * @ingroup LArSoftProxyCustom
   * @brief Infrastructure to describe the element of a proxy
   *
   * A collection proxy element is a complex entity in that it has to figure out
   * all the data it can return, and their type.
   *
   * The default type of a collection proxy element is
   * `proxy::CollectionProxyElement`.
   *
   * @{
   */

  //----------------------------------------------------------------------------
  /**
   * @brief An element of a collection proxy.
   * @tparam CollProxy the type of collection proxy this is the element of
   *
   * The collection proxy element represents, unsurprisingly, a single element
   * of the collection proxy.
   * It exposes all the connections between the components merged into the
   * collection proxy, for a specific element.
   *
   * The interface of a proxy element allows it to access the main and auxiliary
   * data as follows:
   * * main data element (always present!):
   *     * the whole object, by dereference method (`*track`)
   *     * methods and fields, by member access (`track->Length()`)
   * * auxiliary data: accessed by tag (usually it's the default, that is the
   *   same as the data type):
   *     * check if available: `has()` method, static
   *       (`track.has<recob::Hit>()`)
   *     * get the data:
   *         * `get()` method; the data *must* be available or this
   *           call will not compile (`track.get<recob::Hit>()`)
   *         * `getIf()` method; if the data is not available, an exception is
   *           thrown
   * * _index_ of the element in the collection (a bonus!): `index()` method
   *
   * The for loop block illustrates how to use some of them:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto tracks = proxy::getCollection<std::vector<recob::Track>>
   *   (event, tracksTag, proxy::withAssociated<recob::Hit>());
   *
   * for (auto track: tracks) {
   *
   *   recob::Track const& trackObj = *track; // access to the track
   *
   *   double const length = track->Length(); // access to track members
   *
   *   // access to associated data
   *   // (returns random-access collection-like object)
   *   decltype(auto) hits = track.get<recob::Hit>();
   *   double charge = 0.0;
   *   for (auto const& hitPtr: hits) charge += hitPtr->Integral();
   *
   *   std::vector<recob::TrackFitHitInfo> const* fitInfo
   *     = track.has<std::vector<recob::TrackFitHitInfo>>()
   *     ? &(track.getIf<std::vector<recob::TrackFitHitInfo>>())
   *     : nullptr
   *     ;
   *
   *   mf::LogVerbatim("Info")
   *     << "[#" << trackInfo.index() << "] track ID=" << trackObj.ID()
   *     << " (" << length << " cm) deposited charge=" << charge
   *     << " with " << hits.size() << " hits and"
   *     << (fitInfo? "": " no") << " fit information";
   *
   * } // for tracks
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Please note that `getIf()` method is _not easy to use correctly_.
   * In this example, we rely on the knowledge that `track.getIf()` would return
   * a reference to a (non-temporary) object, and we take the address of that
   * object: this is not necessarily the case (e.g., it is not the case for
   * `getIf<recob::Hit>()`).
   *
   *
   * Customization
   * ==============
   *
   * A proxy element class can (and should) be derived from
   * `proxy::CollectionProxyElement`.
   * On top of it, the derived class can extend or completely rewrite the
   * interface:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename CollProxy>
   * class TrackProxy: public proxy::CollectionProxyElement<CollProxy> {
   *   using base_t = proxy::CollectionProxyElement<CollProxy>;
   *     public:
   *   unsigned int nHits() const
   *     { return base_t::template get<recob::Hit>().size(); }
   * }; // class TrackProxy
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Specialization is also possible, but discouraged.
   *
   * A `proxy::CollectionProxyElement` has knowledge of the type of collection
   * proxy it's an element of, but it has no connection to a collection proxy
   * object and in fact it never interacts with any, not even during
   * construction.
   *
   * Once the element is customized, a `proxy::CollectionProxy` can use it by
   * having it as first template argument:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename MainColl, typename... AuxColl>
   * using Proxy_t
   *   = proxy::CollectionProxyBase<TrackProxy, MainColl, AuxColl...>;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * A `CollectionProxyMaker` specialization will have to `make()` such objects.
   *
   *
   * @note There is an indirect dependency of this class on the whole
   *       collection proxy class it is an element of. The type of that
   *       collection proxy is the template argument `CollProxy`.
   *       This class uses `CollProxy` to discovers some types it needs, but it
   *       does not contain or link to `CollProxy` itself.
   *       It should be possible therefore to have `CollProxy` types and
   *       `CollectionProxyElement` depend on a third, "trait" type delivering
   *       the relevant data types to both objects.
   *
   */
  template <typename CollProxy>
  struct CollectionProxyElement {

  public:
    using collection_proxy_t = CollProxy;
    using main_element_t = typename collection_proxy_t::main_element_t;

    /// Tuple of elements (expected to be tagged types).
    using aux_elements_t =
      typename details::SubstituteWithAuxList<typename collection_proxy_t::aux_collections_t>::type;

    /// Constructor: sets the element index, the main element and steals
    /// auxiliary data.
    CollectionProxyElement(std::size_t index, main_element_t const& main, aux_elements_t&& auxData)
      : fIndex(index), fMain(&main), fAuxData(std::move(auxData))
    {}

    /// Returns a pointer to the main element.
    main_element_t const* operator->() const { return fMain; }

    /// Returns a reference to the main element.
    main_element_t const& operator*() const { return *fMain; }

    /// Returns the index of this element in the collection.
    std::size_t index() const { return fIndex; };

    /// Returns the auxiliary data specified by type (`Tag`).
    template <typename Tag>
    auto get() const -> decltype(auto)
    {
      return std::get<util::index_of_tag_v<Tag, aux_elements_t>>(fAuxData);
    }

    /**
     * @brief Returns the auxiliary data specified by type (`Tag`).
     * @tparam Tag tag of the data to fetch (usually, its type)
     * @tparam T type to return (by default, a constant reference to `Tag`)
     * @return the auxiliary data specified by type (`Tag`).
     * @throw std::logic_error if the tag is not available.
     * @see get(), has()
     * 
     * @deprecated C++17 `if constexpr` should be used instead
     *             (see the example below)
     *             
     * This method is a `get()` which forgives when the requested type is not
     * available (because this proxy was configured not to hold it).
     *
     * The difference with `get()` is the following:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * decltype(auto) elem = tracks[0];
     * if (elem.has<recob::Hit>()) {
     *   auto hits = elem.get<recob::Hit>();
     *   // ...
     * }
     * if (elem.has<recob::SpacePoint>()) {
     *   auto spacepoints = elem.getIf<recob::SpacePoint>();
     *   // ...
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * If the proxy `tracks` has _not_ been coded with `recob::Hit` data,
     * the code snippet will _not_ compile, because `get()` will not compile
     * for tags that were not coded in. On the other end, if `recob::Hit`
     * is coded in `tracks` but `recob::SpacePoint` is not, the snippet _will_
     * compile. In that case, if the `has()` check had been omitted, `getIt()`
     * would throw a `std::logic_error` exception when executed.
     * With C++17, this will not be necessary any more by using "constexpr if":
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * decltype(auto) elem = tracks[0];
     * if constexpr (elem.has<recob::Hit>()) {
     *   auto hits = elem.get<recob::Hit>();
     *   // ...
     * }
     * if constexpr (elem.has<recob::SpacePoint>()) {
     *   auto spacepoints = elem.get<recob::SpacePoint>();
     *   // ...
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @note The second template argument contains the *exact* type returned by
     *       the function. That information is needed because this method needs
     *       to return the same type (or compatible types) whether the required
     *       tag is present or not, in order for user code like
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     *       if (elem.has<recob::SpacePoint>()) {
     *         recob::SpacePoint const* spacepoints
     *           = elem.getIf<recob::SpacePoint>(); // won't do
     *         // ...
     *       }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *       to work; it becomes:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     *       if (elem.has<recob::SpacePoint>()) {
     *         recob::SpacePoint const* spacepoints
     *           = elem.getIf<recob::SpacePoint, recob::SpacePoint const*>();
     *         // ...
     *       }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *       This is necessary because when the tag (in the example,
     *       `recob::SpacePoint`) is not registered in the proxy, `getIf()` has
     *       no clue of what the return value should be ("which is the type of
     *       the data that does not exist?").
     *
     */
    template <typename Tag, typename T = Tag const&>
    [[deprecated("Use C++17 constexpr if instead and get() instead")]] auto getIf() const
      -> decltype(auto);

    /// Returns whether this class knowns about the specified type (`Tag`).
    template <typename Tag>
    static constexpr bool has()
    {
      return util::has_tag_v<Tag, aux_elements_t>;
    }

  private:
    std::size_t fIndex;          ///< Index of this element in the proxy.
    main_element_t const* fMain; ///< Pointer to the main object of the element.

    // note that the auxiliary data is not tagged, we need to learn which
    // the tags are from the collection.
    aux_elements_t fAuxData; ///< Data associated to the main object.

    template <typename Tag, typename>
    auto getIfHas(std::bool_constant<true>) const -> decltype(auto);

    template <typename Tag, typename T>
    [[noreturn]] auto getIfHas(std::bool_constant<false>) const -> T;

  }; // CollectionProxyElement<>

  /// @}
  //--- END Proxy element infrastructure ---------------------------------------

  //----------------------------------------------------------------------------
  namespace details {

    //--------------------------------------------------------------------------
    //---  Stuff for the whole collection proxy
    //--------------------------------------------------------------------------
    /**
     * @brief Creates a collection proxy element object from data structures.
     * @tparam ProxyElement type of proxy element to be created
     * @tparam AuxData types of auxiliary data structures being included
     * @param index index in main collection of the element being represented
     * @param main main collection proxy data
     * @param auxData auxiliary data collections
     * @return a `ProxyElement` object bound to the specified data element
     */
    template <typename ProxyElement, typename... AuxData>
    auto makeCollectionProxyElement(std::size_t index,
                                    typename ProxyElement::main_element_t const& main,
                                    AuxData&&... auxData)
    {
      return ProxyElement(
        index, main, typename ProxyElement::aux_elements_t(std::forward<AuxData>(auxData)...));
    } // makeCollectionProxyElement()

    //--------------------------------------------------------------------------

  } // namespace details

} // namespace proxy

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace proxy {

  namespace details {

    //--------------------------------------------------------------------------
    //--- stuff for auxiliary data
    //--------------------------------------------------------------------------
    // Trait replacing each element of the specified tuple with its
    // `auxiliary_data_t`
    template <typename Tuple>
    struct SubstituteWithAuxList {
      static_assert(util::always_true_type<Tuple>(), "Template argument must be a tuple");
    }; // SubstituteWithAuxList<>

    template <typename... T>
    struct SubstituteWithAuxList<std::tuple<T...>> {
      using type = std::tuple<typename T::auxiliary_data_t...>;
    }; // SubstituteWithAuxList<tuple>

    //--------------------------------------------------------------------------

  } // namespace details

  //----------------------------------------------------------------------------
  //---  CollectionProxyElement
  //----------------------------------------------------------------------------
  template <typename CollProxy>
  template <typename Tag, typename T>
  auto CollectionProxyElement<CollProxy>::getIf() const -> decltype(auto)
  {
    return getIfHas<Tag, T>(std::bool_constant<has<Tag>()>{});
  }

  //----------------------------------------------------------------------------
  template <typename CollProxy>
  template <typename Tag, typename>
  auto CollectionProxyElement<CollProxy>::getIfHas(std::bool_constant<true>) const -> decltype(auto)
  {
    return get<Tag>();
  }

  template <typename CollProxy>
  template <typename Tag, typename T>
  auto CollectionProxyElement<CollProxy>::getIfHas(std::bool_constant<false>) const -> T
  {
    throw std::logic_error("Tag '" + lar::debug::demangle<Tag>() + "' not available.");
  }

  //----------------------------------------------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_COLLECTIONPROXYELEMENT_H
