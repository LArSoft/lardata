/**
 * @file   lardata/RecoBaseProxy/ProxyBase/CollectionProxy.h
 * @brief  Utilities for the collection proxy object.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_COLLECTIONPROXY_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_COLLECTIONPROXY_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...
#include "lardata/RecoBaseProxy/ProxyBase/CollectionProxyElement.h"
#include "lardata/RecoBaseProxy/ProxyBase/MainCollectionProxy.h"
#include "lardata/Utilities/TupleLookupByTag.h" // util::type_with_tag_t, ...

// C/C++ standard
#include <cstdlib> // std::size_t
#include <limits>  // std::numeric_limits<>
#include <tuple>
#include <utility> // std::move()
#include <vector>

namespace proxy {

  namespace details {

    template <typename Cont>
    class IndexBasedIterator;

    template <template <typename, typename...> class F, typename...>
    struct TemplateAdaptorOnePlus;

  } // namespace details

  // --- BEGIN Collection proxy infrastructure ---------------------------------
  /**
   * @defgroup LArSoftProxyCollections Collection proxy infrastructure
   * @ingroup LArSoftProxyCustom
   * @brief Infrastructure to define a proxy of collection data product.
   *
   * A data collection proxy connects a main collection data product with other
   * data products whose elements have relation with one of the main product
   * elements.
   *
   * A collection proxy is created via a call to `getCollection()`, in a fashion
   * non dissimilar from calling, e.g., `art::Event::getValidHandle()`.
   *
   * A proxy is characterized by a proxy tag, which is the type used in the
   * `getCollection()` call, but that is not necessarily in direct relation with
   * the type of the collection proxy itself. In other words, calling
   * `proxy::getCollection<proxy::Tracks>(event, trackTag)` will return a proxy
   * object whose type is likely not `proxy::Tracks`.
   *
   * Proxy collections are created by _merging_ auxiliary data into a main
   * collection data product. See `getCollection()` for more details.
   *
   * Proxies can be customized, meaning that collection proxies may be written
   * to have a specific interface. There are different levels of customization,
   * of the collection proxy itself, of its element type, and of the auxiliary
   * data merged. Customization may require quite a bit of coding, and the
   * easiest approach is to start from an already customized proxy implementing
   * features similar to the desired ones.
   *
   * It's worth stressing @ref LArSoftProxyQuirks "again" that collection
   * proxies composed by merging different auxiliary data have different C++
   * types, and that if such proxies need to be propagated as function arguments
   * those arguments need to be of template type.
   *
   * @{
   */

  //----------------------------------------------------------------------------
  /**
   * @brief Base representation of a collection of proxied objects.
   * @tparam Element type of element of the collection proxy
   * @tparam MainColl type of the collection of the main data product
   * @tparam AuxColls type of all included auxiliary data proxies
   * @see proxy::CollectionProxyElement
   *
   * This object exposes a collection interface.
   * The collection proxy is driven by a data product containing the main
   * objects. The size of the collection proxy is the same as the one of this
   * main data product, and all associated data is referring to its elements.
   *
   * Thus, the elements of this collection proxy are objects that collect the
   * information of a single element in the main data product and all the data
   * associated with it.
   *
   * The `AuxColls` types are tagged types: all must define a `tag` type.
   * Their data is accessed specifying that tag, i.e. via `get<Tag>()`.
   * Therefore, tags must be unique.
   *
   * The type `Element` is expected to expose the same interface of
   * `CollectionProxyElement`, from which it can derive. It is a template
   * that needs to take as only argument the type of collection proxy it is the
   * element of. This is a way to customize the interface of access to single
   * element of proxy.
   *
   * @note This class depends on an `Element` type, which indirectly depends on
   *       this class for discovering some relevant data types. This is a
   *       circular dependency that might be solved by introducing a third class
   *       with the definition of the types that both classes need.
   *
   */
  template <template <typename CollProxy> class Element, typename MainColl, typename... AuxColls>
  class CollectionProxyBase : public details::MainCollectionProxy<MainColl>, public AuxColls... {
    /// This type.
    using collection_proxy_t = CollectionProxyBase<Element, MainColl, AuxColls...>;

    /// Type of wrapper used for the main data product.
    using main_collection_proxy_t = details::MainCollectionProxy<MainColl>;

  public:
    /// Type of element of the main data product.
    using typename main_collection_proxy_t::main_element_t;

    /// Type of collection in the main data product.
    using typename main_collection_proxy_t::main_collection_t;

    /// Type of element of this collection proxy.
    using element_proxy_t = Element<collection_proxy_t>;

    /// Tuple of all auxiliary data collections (wrappers).
    using aux_collections_t = std::tuple<AuxColls...>;

    /// Type of element of this collection proxy.
    using value_type = element_proxy_t;

    /// Type of iterator to this collection (constant).
    using const_iterator = details::IndexBasedIterator<collection_proxy_t>;

    /// Type of iterator to this collection (still constant).
    using iterator = const_iterator;

    /**
     * @brief Constructor: uses the specified data.
     * @param main the original main data product collection
     * @param aux all auxiliary data collections and structures
     *
     * The auxiliary data structures are stolen (moved) from the arguments.
     * They are expected to be wrappers around the original associated data,
     * not owning the auxiliary data itself.
     */
    CollectionProxyBase(main_collection_t const& main, AuxColls&&... aux)
      : main_collection_proxy_t(main), AuxColls(std::move(aux))...
    {}

    /**
     * @brief Returns the element of this collection with the specified index.
     * @param i the index in the collection
     * @return a value representing an element of the collection
     *
     * The returned value is an object created on the spot, not a reference to
     * an existing structure.
     * The structure exposes the `i`-th element in the main collection, plus all
     * objects that are associated to it.
     */
    element_proxy_t const operator[](std::size_t i) const
    {
      return details::makeCollectionProxyElement<element_proxy_t>(
        i, getMainAt(i), aux<AuxColls>().operator[](i)...);
    }

    /// Returns an iterator to the first element of the collection.
    const_iterator begin() const { return makeIterator(0U); }

    /// Returns an iterator past the last element of the collection.
    const_iterator end() const { return makeIterator(size()); }

    /// Returns whether this collection is empty.
    bool empty() const { return main().empty(); }

    /// Returns the size of this collection.
    std::size_t size() const { return main().size(); }

    /// Returns the associated data proxy specified by `AuxTag`.
    template <typename AuxTag>
    auto get() const -> decltype(auto)
    {
      return auxByTag<AuxTag>();
    }

    /**
     * @brief Returns the auxiliary data specified by type (`Tag`).
     * @tparam Tag tag of the data to fetch (usually, its type)
     * @tparam T exact type returned by the method (by default a vector of tags)
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
     * if (tracks.has<recob::Hit>()) {
     *   auto hits = tracks.get<recob::Hit>();
     *   // ...
     * }
     * if (tracks.has<recob::SpacePoint>()) {
     *   auto spacepoints = tracks.getIf<recob::SpacePoint>();
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
     * if constexpr (tracks.has<recob::Hit>()) {
     *   auto hits = tracks.get<recob::Hit>();
     *   // ...
     * }
     * if constexpr (tracks.has<recob::SpacePoint>()) {
     *   auto spacepoints = tracks.get<recob::SpacePoint>();
     *   // ...
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @note If the wrapped data product is something different than a vector
     *       of space points (which in the example is likely, if space points
     *       are associated to tracks), the almost-correct type of return value
     *       needs to be specified as second template parameter `T`.
     *
     * @warning This functionality is not trivial to use!
     *          It's mostly meant for implementation of higher level wrappers.
     *
     * @deprecated Use C++17 constexpr if instead.
     */
    template <typename Tag, typename T = std::vector<Tag> const&>
    [[deprecated("Use C++17 constexpr if instead and get() instead")]] auto getIf() const
      -> decltype(auto);

    /// Returns whether this class knowns about the specified type (`Tag`).
    template <typename Tag>
    static constexpr bool has()
    {
      return util::has_tag_v<Tag, aux_collections_t>;
    }

  protected:
    using main_collection_proxy_t::getMainAt;
    using main_collection_proxy_t::main;
    using main_collection_proxy_t::mainProxy;

    /// Returns the auxiliary data specified by type.
    template <typename AuxColl>
    AuxColl const& aux() const
    {
      return static_cast<AuxColl const&>(*this);
    }

    /// Returns the auxiliary data specified by type.
    template <typename AuxTag>
    auto auxByTag() const -> decltype(auto)
    {
      return aux<util::type_with_tag_t<AuxTag, aux_collections_t>>();
    }

    template <typename Tag, typename>
    auto getIfHas(std::bool_constant<true>) const -> decltype(auto);
    template <typename Tag, typename T>
    [[noreturn]] auto getIfHas(std::bool_constant<false>) const -> T;

    /// Returns an iterator pointing to the specified index of this collection.
    const_iterator makeIterator(std::size_t i) const { return {*this, i}; }

    static_assert(!util::has_duplicate_tags<aux_collections_t>(),
                  "Some auxiliary data collections share the same tag. They should not.");

  }; // struct CollectionProxyBase

  //----------------------------------------------------------------------------
  /**
   * @brief Base representation of a collection of proxied objects.
   * @tparam MainColl type of the collection of the main data product
   * @tparam AuxColls type of all included auxiliary data proxies
   * @see proxy::CollectionProxyElement
   *
   * This object is a "specialization" of `proxy::CollectionProxyBase` using
   * `proxy::CollectionProxyElement` as element type.
   */
  template <typename MainColl, typename... AuxColls>
  using CollectionProxy = CollectionProxyBase<CollectionProxyElement, MainColl, AuxColls...>;

  // this joke is necessary because expanding directly CollectionProxy<Args...>
  // into CollectionProxy<Main, Aux...>  template arguments does not work
  template <typename... Args>
  using CollectionProxyFromArgs =
    typename details::TemplateAdaptorOnePlus<CollectionProxy, Args...>::type;

  /// @}
  // --- END Collection proxy infrastructure -----------------------------------

  //----------------------------------------------------------------------------
  namespace details {

    //--------------------------------------------------------------------------
    /// Creates a collection proxy of a specified type with the given arguments.
    template <template <typename...> class CollProxy, typename MainColl, typename... AuxColl>
    auto createCollectionProxy(MainColl const& main, AuxColl&&... aux)
    {
      return CollProxy<MainColl, AuxColl...>(main, std::forward<AuxColl>(aux)...);
    }

    //--------------------------------------------------------------------------
    /// Creates a `CollectionProxy` object with the given arguments.
    template <typename MainColl, typename... AuxColl>
    auto makeCollectionProxy(MainColl const& main, AuxColl&&... aux)
    {
      return createCollectionProxy<CollectionProxy>(main, std::forward<AuxColl>(aux)...);
    }

    //--------------------------------------------------------------------------
    /**
     * @brief Iterator to random access collection storing a current index.
     * @tparam Cont type of random-access container to iterate
     *
     * `Cont` is a type providing a public `operator[](std::size_t)` method.
     */
    template <typename Cont>
    class IndexBasedIterator {

    public:
      using container_t = Cont;

      using value_type = util::collection_value_t<container_t>;
      using const_iterator = IndexBasedIterator;

      /// Default constructor (required by iterator protocol): an unusable iterator.
      IndexBasedIterator() = default;

      /// Constructor: initializes from an iterator of the proxy main collection.
      IndexBasedIterator(container_t const& cont, std::size_t index = 0)
        : fCont(&cont), fIndex(index)
      {}

      /// Returns the value pointed by this iterator.
      auto operator*() const -> decltype(auto) { return fCont->operator[](fIndex); }

      /// Returns the value pointed by this iterator.
      const_iterator& operator++()
      {
        ++fIndex;
        return *this;
      }

      /// Returns whether the iterators point to the same element.
      bool operator!=(const_iterator const& other) const
      {
        return (other.fIndex != fIndex) || (other.fCont != fCont);
      }

    protected:
      container_t const* fCont = nullptr; ///< Pointer to the original container.

      /// Current index in the main collection.
      std::size_t fIndex = std::numeric_limits<std::size_t>::max();

    }; // IndexBasedIterator<>

  } // namespace details

} // namespace proxy

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace proxy {

  namespace details {

    //--------------------------------------------------------------------------
    template <template <typename, typename...> class F, typename First, typename... Others>
    struct TemplateAdaptorOnePlus<F, First, Others...> {
      using type = F<First, Others...>;
    };

  } // namespace details

  //----------------------------------------------------------------------------
  //---  CollectionProxyBase
  //----------------------------------------------------------------------------
  template <template <typename CollProxy> class Element, typename MainColl, typename... AuxColls>
  template <typename Tag, typename T>
  auto CollectionProxyBase<Element, MainColl, AuxColls...>::getIf() const -> decltype(auto)
  {
    return getIfHas<Tag, T>(std::bool_constant<has<Tag>()>{});
  }

  template <template <typename CollProxy> class Element, typename MainColl, typename... AuxColls>
  template <typename Tag, typename>
  auto CollectionProxyBase<Element, MainColl, AuxColls...>::getIfHas(std::bool_constant<true>) const
    -> decltype(auto)
  {
    return get<Tag>();
  }

  template <template <typename CollProxy> class Element, typename MainColl, typename... AuxColls>
  template <typename Tag, typename T>
  auto CollectionProxyBase<Element, MainColl, AuxColls...>::getIfHas(
    std::bool_constant<false>) const -> T
  {
    throw std::logic_error("Tag '" + lar::debug::demangle<Tag>() + "' not available.");
  }

  //----------------------------------------------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_COLLECTIONPROXY_H
