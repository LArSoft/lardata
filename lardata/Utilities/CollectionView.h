/**
 * @file   lardata/Utilities/CollectionView.h
 * @brief  Provides the features of a collections, from begin and end iterators.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   August 3rd, 2017
 *
 * This library is header-only.
 *
 * @note It is likely that a skilful use of Nibbler's range library will
 *       provide the same functionality (and loads more).
 *
 *
 *
 * @section InterfaceSubstitution Interface substitution technique
 *
 * A technique that is used in this implementation is to replace (or extend) the
 * interface of an existing object.
 * A key requirement is that the new interface object must not have any
 * additional state.
 *
 * The interface class is superimposed to the _existing_ data without
 * replication by _reinterpreting_ its content. This is achieved deriving the
 * new interface class from the data class:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 *
 * class Data {
 *   double chiSq;
 *   double NDF;
 *     public:
 *   Data(double chiSq, double NDF): chiSq(chiSq), NDF(NDF) {}
 *
 *   double chiSquare() const { return chiSq; }
 *   double DegreesOfFreedom() const { return NDF; }
 *
 * }; // class Data
 *
 *
 * class DataInterface: private Data {
 *   Data const& asData() const { return static_cast<Data const&>(*this); }
 *     public:
 *
 *   double normalizedChiSquare() const
 *     { return asData().chiSquare() / asData().DegreesOfFreedom(); }
 *
 *     protected:
 *   DataInterface(Data const& from): Data(data) {}
 *
 *   friend DataInterface const& makeDataInterface(Data const&);
 *
 * }; // class DataInterface
 *
 *
 * DataInterface const& makeDataInterface(Data const& data)
 *   { return static_const<DataInterface const&>(data); }
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * With this pattern, an interface object can be obtained only by calling
 * `makeDataInterface()` on the base object, and in this way it will be returned
 * only as a reference (in this case, constant). The interface object can't
 * be copied, and it must be passed around as reference. It's not possible to
 * convert it back to `Data`, because the base class is private.
 * There is a single protected constructor. This choice, compared to deleting
 * all constructors, allows for a derived class to acquire the same interface:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * struct DataWithInterface: public DataInterface {
 *   DataWithInterface(Data const& from): DataInterface(from) {}
 * }; // class DataWithInterface
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * This simple class provides the storage for `Data` in addition to exposing
 * `DataInterface`.
 * There are other ways to achieve the same goal (e.g., multiple inheritance).
 *
 * The presence of a friend function should raise a warning. Friendship is
 * required because the function is attempting a downcast from a private base
 * class. If it is intended that the full `Data` interface is exposed, then
 * the inheritance can be public and no special friendship will be needed.
 * Another way is to replace the `static_cast` with a `reinterpret_cast`.
 *
 *
 */

#ifndef LARDATA_UTILITIES_COLLECTIONVIEW_H
#define LARDATA_UTILITIES_COLLECTIONVIEW_H

// C/C++ standard libraries
#include <cstddef>     // std::size_t
#include <iterator>    // std::iterator_traits, std::reverse_iterator
#include <stdexcept>   // std::out_of_range
#include <string>      // std::to_string()
#include <type_traits> // std::declval(), ...

namespace lar {

  // forward declarations
  template <typename Range>
  class CollectionView;

  namespace details {

    //--------------------------------------------------------------------------
    template <typename Range>
    struct RangeTraits {
      using range_t = Range;
      using traits_t = RangeTraits<range_t>;
      static auto getCBegin(range_t const& range)
      {
        using std::cbegin;
        return cbegin(range);
      }
      static auto getCEnd(range_t const& range)
      {
        using std::cend;
        return cend(range);
      }
      using begin_iterator_t = std::decay_t<decltype(traits_t::getCBegin(std::declval<range_t>()))>;
      using end_iterator_t = std::decay_t<decltype(traits_t::getCEnd(std::declval<range_t>()))>;
    }; // RangeTraits<>

    //--------------------------------------------------------------------------
    /// Class storing a begin and a end iterator.
    template <typename BeginIter, typename EndIter = BeginIter>
    class CollectionExtremes {
    public:
      using begin_iterator_t = BeginIter;
      using end_iterator_t = EndIter;

      struct FromContainerTag {};
      static constexpr FromContainerTag fromContainer{};

      /// Constructor: copies the specified iterators.
      CollectionExtremes(BeginIter b, EndIter e) : b(b), e(e) {}

      /// Returns the stored begin iterator.
      begin_iterator_t const& begin() const { return b; }

      /// Returns the stored end iterator.
      end_iterator_t const& end() const { return e; }

    private:
      begin_iterator_t b; ///< Stored copy of begin iterator.
      end_iterator_t e;   ///< Stored copy of end iterator.
    };                    // class CollectionExtremes<>

    //--------------------------------------------------------------------------
    /// Helper to create a CollectionExtremes object from two iterators.
    template <typename BeginIter, typename EndIter>
    auto makeCollectionExtremes(BeginIter const& b, EndIter const& e)
    {
      return CollectionExtremes<BeginIter, EndIter>(b, e);
    }

    //--------------------------------------------------------------------------
    /// Helper to create a CollectionExtremes object from a range object.
    template <typename Range>
    auto makeCollectionExtremes(Range const& range)
    {
      using std::cbegin;
      using std::cend;
      return makeCollectionExtremes(cbegin(range), cend(range));
    } // makeCollectionExtremes(range)

    //--------------------------------------------------------------------------
    // forward declaration
    template <typename Range>
    CollectionView<Range> makeCollectionView(Range&&);
    //--------------------------------------------------------------------------

  } // namespace details

  //----------------------------------------------------------------------------
  /**
   * @brief Provides features of a collections, from begin and end iterators.
   * @tparam Range the type of collection to be wrapped
   *
   * A collection view is a class that offers a collection-like interface,
   * mostly like the C+ standard library containers, based on two iterators.
   *
   * The base, wrapped collection is required to respond to `begin()` and a
   * `end()` global functions. If the desired view is not described by such
   * an object, a temporary one must be created (see `makeCollectionView()`).
   *
   * There are two ways to use this class:
   * 1. to wrap an existing `Range`-like container
   * 2. to turn two iterators into a container
   *
   * The two use cases are both addressed by this class, but helper functions
   * are provided to make it easier to create them as needed.
   *
   * @note While the object is currently copiable and moveable, this is not
   *       guaranteed for the future.
   *
   *
   * Wrap an existing `Range`-like container
   * ----------------------------------------
   *
   * Here we assume there is somewhere an instance of the object `range` which
   * fulfills the requirement of the `Range` type above, that is it responds
   * to the requests of a `begin()` and a `end()` iterator.
   *
   * To create a collection view of `range`, the easiest way is to use
   * `wrapCollectionIntoView()`. In the following example the `range` object
   * is a STL vector (which does not really need any added interface...):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> range(5);
   * std::iota(range.begin(), range.end(), 1); // { 1, 2, 3, 4, 5 }
   *
   * for (int d: lar::wrapCollectionIntoView(range)) {
   *   std::cout << d << " ";
   * }
   * std::cout << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which will print "1 2 3 4 5 ".
   * Here the wrapped collection object, returned by `wrapCollectionIntoView()`,
   * is insubstantial. It can be saved with
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * decltype(auto) view = lar::wrapCollectionIntoView(range);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * or even:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const& view = lar::wrapCollectionIntoView(range);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * but it will be just a (constant) reference to something else.
   *
   *
   * Turn two iterators into a container
   * ------------------------------------
   *
   * In this approach, we have two iterators to an existing collection, and we
   * want to use them as extremes of a "virtual" collection.
   * Again we use a STL vector as a base container for the example. Here we want
   * to see a subrange of it as a new collection. We use `makeCollectionView()`.
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> v(10);
   * std::iota(v.begin(), v.end(), 0); // { 0, 1, ..., 9 }
   *
   * for (int d: lar::makeCollectionView(v.begin() + 4, v.begin() + 7)) {
   *   std::cout << d << " ";
   * }
   * std::cout << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will print "0 1 2 3 4 5 6 7 8 9 ".
   *
   *
   * Declaring a wrapper class
   * --------------------------
   *
   * The function `lar::makeCollectionView()` creates a view owning the
   * information the view requires. Similarly, a new class can be defined which
   * does the same, by simply deriving it from `lar::CollectionView`:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * class IntVector {
   *    using vector_t = std::vector<int>;
   *
   *    vector_t data;
   *
   *      public:
   *    IntVector(vector_t&& data): data(std::move(data)) {}
   *
   *    auto begin() const -> decltype(auto) { return data.cbegin(); }
   *    auto end() const -> decltype(auto) { return data.cend(); }
   *
   * }; // struct IntVector
   *
   * using IntViewBase_t = lar::CollectionView<IntVector>;
   *
   * struct MyCollection: public IntViewBase_t {
   *   MyCollection(std::vector<int>&& data) : IntViewBase_t(std::move(data)) {}
   * }; // class MyCollection
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * aftter which `MyCollection` interface can be enriched as needed.
   * The `IntViewBase_t` alias is a way to overcome the fact that the name
   * `IntVector` can't be used inside `MyCollection` because it's actually a
   * private base class (the base class, even if not direct, will hide
   * `IntVector`, even if, being private, it can't even be used). Another way
   * is to fully qualify its name (e.g. `::IntVector`).
   *
   * Note that to avoid accidental copies, `lar::CollectionView` objects can't
   * be directly instantiated: using directly `IntViewBase_t` will _not_ be
   * allowed.
   *
   */
  template <typename Range>
  class CollectionView : private Range {
    using this_t = CollectionView<Range>;           ///< This type.
    using range_t = Range;                          ///< Type of the range being wrapped.
    using traits_t = details::RangeTraits<range_t>; ///< Range traits.

    /// Type of the begin iterator.
    using begin_iter_t = typename traits_t::begin_iterator_t;
    /// Type of the end iterator.
    using end_iter_t = typename traits_t::end_iterator_t;

    /// Type of traits of iterator.
    using iter_traits_t = std::iterator_traits<begin_iter_t>;

  public:
    using collection_type = range_t; ///< Type of collection being wrapped.

    using value_type = typename iter_traits_t::value_type;
    // reference not provided
    using const_reference = typename iter_traits_t::reference;
    // pointer not provided
    using const_pointer = typename iter_traits_t::pointer;
    // iterator not provided
    using const_iterator = begin_iter_t;
    // reverse_iterator not provided
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using difference_type = typename iter_traits_t::difference_type;
    using size_type = std::size_t;

    /// @{
    /// @name Forward access.

    /// Returns whether the collection is empty.
    bool empty() const noexcept { return cbegin() == cend(); }

    /// Returns the size of the collection.
    size_type size() const noexcept { return std::distance(cbegin(), cend()); }

    /// Returns an iterator to the begin of the collection.
    const_iterator cbegin() const noexcept
    {
      using std::cbegin;
      return cbegin(asRange());
    }

    /// Returns an iterator past the end of the collection.
    end_iter_t cend() const noexcept
    {
      using std::cend;
      return cend(asRange());
    }

    /// Returns an iterator to the begin of the collection.
    const_iterator begin() const noexcept { return cbegin(); }

    /// Returns an iterator past the end of the collection.
    end_iter_t end() const noexcept { return cend(); }

    /// Returns the first element in the collection.
    auto front() const -> decltype(auto) { return *cbegin(); }

    /// @}

    /// @{
    /// @name Backward access.

    /// Returns a reverse iterator to the begin of the collection.
    const_reverse_iterator rbegin() const noexcept { return crbegin(); }

    /// Returns a reverse iterator past the end of the collection.
    const_reverse_iterator rend() const noexcept { return crend(); }

    /// Returns a reverse iterator to the begin of the collection.
    const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(cend()); }

    /// Returns a reverse iterator past the end of the collection.
    const_reverse_iterator crend() const noexcept { return const_reverse_iterator(cbegin()); }

    /// Returns the last element in the collection.
    auto back() const -> decltype(auto) { return *crbegin(); }

    /// @}

    /// @{
    /// @name Random access.

    /// Returns the content of the `i`-th element.
    auto operator[](size_type i) const -> decltype(auto) { return cbegin()[i]; }

    /// Returns the content of the `i`-th element.
    auto at(size_type i) const -> decltype(auto)
    {
      if (i >= size()) {
        throw std::out_of_range("CollectionView index out of range: " + std::to_string(i) +
                                " >= " + std::to_string(size()));
      }
      return operator[](i);
    }

    /// @}

    /// @{
    /// @name Contiguous access.

    const_pointer data() const { return &front(); }

    /// @}

  protected:
    /*
    /// @{
    /// @name This object can't be directly constructed.
    ///
    /// A class deriving from this will be used directly as the base collection
    /// to extract the iterators from.
    CollectionViewWrapper() = default;
    CollectionViewWrapper(this_t const&) = default;
    CollectionViewWrapper(this_t&&) = default;
    CollectionViewWrapper& operator=(this_t const&) = default;
    CollectionViewWrapper& operator=(this_t&&) = default;
    /// @}
    */

    /// Constructor: steals the data from the specified range.
    explicit CollectionView(range_t&& range) : range_t(std::move(range)) {}

    // we license the creation of this class from ranges to a single function
    // (and to derived classes)
    friend this_t details::makeCollectionView<range_t>(range_t&&);

  private:
    /// Returns this very object, cast back to `range_t`.
    range_t const& asRange() const { return static_cast<range_t const&>(*this); }

  }; // class CollectionView<>

  //----------------------------------------------------------------------------
  /// Returns the specified container, wrapped in the view.
  template <typename Range>
  CollectionView<Range> const& wrapCollectionIntoView(Range const& c)
  {
    return reinterpret_cast<CollectionView<Range> const&>(c);
  }

  //----------------------------------------------------------------------------
  /// Creates a `CollectionView` from the specified iterators.
  template <typename BeginIter, typename EndIter>
  auto makeCollectionView(BeginIter const& b, EndIter const& e)
  {
    return details::makeCollectionView(details::makeCollectionExtremes(b, e));
  }

  /// Type of collection view owning the two range boundary iterators.
  template <typename BeginIter, typename EndIter = BeginIter>
  using RangeAsCollection_t =
    decltype(makeCollectionView(std::declval<BeginIter>(), std::declval<EndIter>()));

  //----------------------------------------------------------------------------
  namespace details {
    template <typename Range>
    CollectionView<Range> makeCollectionView(Range&& range)
    {
      return CollectionView<Range>(std::move(range));
    }
  } // namespace details
  //----------------------------------------------------------------------------

} // namespace lar

#endif // LARDATA_UTILITIES_COLLECTIONVIEW_H
