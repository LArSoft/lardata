/**
 * @file   lardata/RecoBaseProxy/ProxyBase/AssociatedData.h
 * @brief  Auxiliary data from one-to-many sequential association.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_ASSOCIATEDDATA_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_ASSOCIATEDDATA_H

// LArSoft libraries
#include "larcorealg/CoreUtils/MetaUtils.h" // util::is_not_same<>
#include "lardata/RecoBaseProxy/ProxyBase/AssnsNodeAsTuple.h"
#include "lardata/RecoBaseProxy/ProxyBase/AssnsTraits.h"
#include "lardata/Utilities/CollectionView.h"
#include "lardata/Utilities/TupleLookupByTag.h" // util::add_tag_t, ...

// framework libraries
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"

// C/C++ standard libraries
#include <vector>
// #include <tuple> // std::tuple_element_t<>, std::get()
#include <algorithm> // std::min()
#include <cassert>
#include <cstdlib>     // std::size_t
#include <iterator>    // std::distance(), std::forward_iterator_tag, ...
#include <memory>      // std::addressof()
#include <type_traits> // std::is_same<>, std::enable_if_t<>, ...
#include <utility>     // std::forward(), std::declval(), ...

namespace proxy {

  //----------------------------------------------------------------------------
  /**
   * @ingroup LArSoftProxiesAssociatedData
   * @{
   */

  //----------------------------------------------------------------------------
  namespace details {

    //--------------------------------------------------------------------------
    //---  general infrastructure
    //--------------------------------------------------------------------------

    /**
     * @brief Simple iterator wrapper for manipulation of dereferenced result.
     * @tparam Iter final iterator
     * @tparam BaseIter base iterator
     * @tparam ValueType type returned by the new dereference operator
     *
     * This class is designed to be used as base class for an iterator that
     * redefines the dereference operations without changing anything else.
     * An example of such iterator is:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * using int_iterator_t = std::vector<int>::const_iterator;
     *
     * class twice_iterator
     *   : public IteratorWrapperBase<twice_iterator, int_iterator_t, float>
     * {
     *   using base_iterator_t
     *     = IteratorWrapperBase<twice_iterator, int_iterator_t, float>;
     *
     *     public:
     *
     *   using base_iterator_t::base_iterator_t; // inherit constructors
     *
     *   float operator*() const -> decltype(auto)
     *     { return 2.0F * base_iterator_t::operator*(); }
     *
     * }; // twice_iterator
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * This new iterator `twice_iterator` effectively inherits `int_iterator_t`
     * constructors (via `base_iterator_t`). Its `operator++()` returns a
     * `twice_iterator`, so that operations like `*++it` are correctly handled.
     *
     * Note that the derived classes need to redefine all the dereferencing
     * operators that are supported:
     *
     * * operator*()
     *
     * The `operator->()` is not supported. Typically, `operator*()` returns a
     * reference and `operator->()` a pointer to an existing structure.
     * The `operator*()` is not required to return a reference when the iterator
     * does not allow to modify the pointed data, and the wrapper can manipulate
     * values and return temporary results. To have a similar feature for
     * `operator->()` requires more work, as it should return a smart pointer
     * that would perform the transformation on dereference.
     *
     */
    template <typename Iter, typename DataIter, typename ValueType = typename DataIter::value_type>
    class IteratorWrapperBase : private DataIter {
    protected:
      using data_iterator_t = DataIter;

    public:
      using iterator = Iter; ///!< The type of base iterator wrapper.

      /// @{
      /// @name Iterator traits
      using typename data_iterator_t::difference_type;
      using value_type = ValueType;
      using pointer = std::add_pointer_t<value_type>;
      using reference = std::add_lvalue_reference_t<value_type>;
      using iterator_category = std::forward_iterator_tag;
      /// @}

      /// Default constructor: default-constructs the underlying iterator.
      IteratorWrapperBase() = default;

      /// Copy-from-base constructor.
      IteratorWrapperBase(data_iterator_t const& from) : data_iterator_t(from) {}

      /// Prefix increment operator.
      iterator& operator++()
      {
        data_iterator_t::operator++();
        return asIterator();
      }

      /// Comparison with a data iterator (makes unnecessary to wrap end iterators).
      bool operator!=(data_iterator_t const& other) const { return other != asDataIterator(); }

      /// Comparison with another iterator.
      bool operator!=(iterator const& other) const { return operator!=(other.asDataIterator()); }

      auto operator[](std::size_t index) const -> decltype(auto)
      {
        return asIterator().transform(asDataIterator() + index);
      }

      /// Dereference operator; need to be redefined by derived classes.
      auto operator*() const -> decltype(auto) { return asIterator().transform(asDataIterator()); }

      /// Dereference operator; need to be redefined by derived classes.
      auto operator-> () const -> decltype(auto) { return makeValuePointer(operator*()); }

    protected:
      /// Transforms and returns the value at the specified data iterator.
      static auto transform(data_iterator_t const&) -> decltype(auto)
      {
        return data_iterator_t::operator*();
      }

      data_iterator_t const& asDataIterator() const
      {
        return static_cast<data_iterator_t const&>(*this);
      }

    private:
      /// Value box for use with pointer dereference `operator->()`.
      template <typename Value>
      class ValuePtr {
        Value value; ///< Value to return the address of (may be reference).
      public:
        ValuePtr(Value const& value) : value(value) {}
        /// Access the contained value via its pointer.
        auto operator-> () const -> decltype(auto) { return std::addressof(value); }
      }; // class ValuePtr<>
      template <typename Value>
      static ValuePtr<Value> makeValuePointer(Value&& value)
      {
        return {std::forward<Value>(value)};
      }

      iterator const& asIterator() const { return static_cast<iterator const&>(*this); }
      iterator& asIterator() { return static_cast<iterator&>(*this); }

    }; // IteratorWrapperBase<>

#if 0

    /// Modified iterator returning the `N`-th element out of the pointed tuple.
    template <std::size_t N, typename TupleIter>
    class tuple_element_iterator:
      public IteratorWrapperBase<
        tuple_element_iterator<N, TupleIter>,
        TupleIter,
        std::tuple_element_t<N, typename TupleIter::value_type>
        >
    {
      using base_iterator_t = IteratorWrapperBase<
        tuple_element_iterator<N, TupleIter>,
        TupleIter,
        std::tuple_element_t<N, typename TupleIter::value_type>
        >;

        public:
      using base_iterator_t::base_iterator_t;

      /// Constructor from a base iterator (explicitly allowed).
      tuple_element_iterator(base_iterator_t const& from)
        : base_iterator_t(from) {}

      static auto transform(TupleIter const& v) -> decltype(auto)
        {return std::get<N>(*v); }

    }; // tuple_element_iterator

#endif // 0
    //--------------------------------------------------------------------------
    //--- BEGIN iterators for art::Assns
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------

    /// This type extends the interface of the art pointer to Assns right side.
    template <typename ArtAssnsIterValue>
    class AssnsNode : private ArtAssnsIterValue {

      using base_t = ArtAssnsIterValue;            ///< Base class type.
      using this_t = AssnsNode<ArtAssnsIterValue>; ///< This class.
      using node_t = ArtAssnsIterValue;            ///< Type of the wrapped node.

      /// Set of traits of the node.
      using assns_node_traits_t = lar::util::assns_traits<node_t>;

    public:
      /// Type of the main (left) object in the association.
      using main_t = typename assns_node_traits_t::left_t;

      /// Type of the associated (right) object.
      using value_t = typename assns_node_traits_t::right_t;

      /// Type of the associated additional data (`void` if none).
      using data_t = typename assns_node_traits_t::data_t;

      /// Type of _art_ pointer to the main (left) object in the association.
      using mainptr_t = typename assns_node_traits_t::leftptr_t;

      /// Type of _art_ pointer to the associated (right) object.
      using valueptr_t = typename assns_node_traits_t::rightptr_t;

      /// Type of the pointer to associated additional data.
      using dataptr_t = typename assns_node_traits_t::dataptr_t;

      /// @{
      /// @name Access to the associated (right) value

      /// Returns the _art_ pointer to the associated value.
      valueptr_t const& valuePtr() const { return base_t::second; }

      /// Returns a reference to the associated value.
      value_t const& value() const { return *valuePtr(); }

      /// @}

      /// @{
      /// @name Access to the key (left) value

      /// Returns the _art_ pointer to the main value, key of the association.
      mainptr_t const& mainPtr() const { return base_t::first; }

      /// Returns the main value, key of the association.
      main_t const& main() const { return *mainPtr(); }

      /// @}

      /// @{
      /**
       * @name Metadata access
       *
       * The complete interface is available only if the association has
       * metadata. Otherwise, only the static member `hasMetadata()` is
       * available.
       */

      // the templates are needed to activate "SFINAE" on std::enable_if
      /// Returns whether this node type supports metadata.
      template <typename Node = node_t>
      static constexpr bool hasMetadata()
      {
        return lar::util::assns_has_metadata_v<Node>;
      }

      /// Returns the pointer to the metadata on this association node.
      template <typename Node = node_t>
      std::enable_if_t<hasMetadata<Node>(), dataptr_t> dataPtr() const
      {
        return base_t::data;
      }

      // this is even more complicate, since if `data_t` is void we can't write
      // `data_t const&` as type of enable_if, because it does not depend on
      // templates and therefore it may be unconditionally evaluated;
      // and C++ does not like references to `void`...
      /// Returns a reference to the metadata on this association node.
      template <typename Node = node_t>
      std::enable_if_t<hasMetadata<Node>(), typename lar::util::assns_traits<Node>::data_t const&>
      data() const
      {
        return *dataPtr();
      }

      /// @}

      /// @{
      /// @name Interface to the art pointer to the associated (right) value

      /// Implicit conversion to _art_ pointer of the associated object.
      operator valueptr_t const&() const& { return valuePtr(); }

      /// Implicit conversion to _art_ pointer of the associated object.
      operator valueptr_t() const&& { return valuePtr(); }

      /// Returns a reference to the associated value (alias of `value()`).
      value_t const& operator*() const { return value(); }

      /// Returns the associated value (alias of `valuePtr()`).
      valueptr_t operator->() const { return valuePtr(); }

      /// Returns the key of the _art_ pointer to the value.
      auto key() const -> decltype(auto) { return valuePtr().key(); }

      /// Returns the product ID of the _art_ pointer to the value.
      auto id() const -> decltype(auto) { return valuePtr().id(); }
      /// @}

      /// Reinterprets the specified association node as a `AssnsNode`.
      static this_t const& makeFrom(node_t const& from) { return static_cast<this_t const&>(from); }

    }; // class AssnsNode<>

    template <typename ArtAssnsIterValue>
    bool operator==(AssnsNode<ArtAssnsIterValue> const& A,
                    typename AssnsNode<ArtAssnsIterValue>::valueptr_t const& B)
    {
      return A.valuePtr() == B;
    }
    template <typename ArtAssnsIterValue>
    bool operator==(typename AssnsNode<ArtAssnsIterValue>::valueptr_t const& A,
                    AssnsNode<ArtAssnsIterValue> const& B)
    {
      return A == B.valuePtr();
    }
    template <typename ArtAssnsIterValue>
    bool operator!=(AssnsNode<ArtAssnsIterValue> const& A,
                    typename AssnsNode<ArtAssnsIterValue>::valueptr_t const& B)
    {
      return A.valuePtr() != B;
    }
    template <typename ArtAssnsIterValue>
    bool operator!=(typename AssnsNode<ArtAssnsIterValue>::valueptr_t const& A,
                    AssnsNode<ArtAssnsIterValue> const& B)
    {
      return A != B.valuePtr();
    }

    /// Reinterprets the specified association node as a `AssnsNode`.
    template <typename ArtAssnsIterValue>
    AssnsNode<ArtAssnsIterValue> const& makeAssnsNode(ArtAssnsIterValue const& from)
    {
      return AssnsNode<ArtAssnsIterValue>::makeFrom(from);
    }

  } // namespace details
} // namespace proxy

// we interrupt this namespace for an urgent specialization...
namespace lar {
  namespace util {

    // specialization for the art node wrapper
    template <typename ArtAssnsIterValue>
    struct assns_metadata_type<proxy::details::AssnsNode<ArtAssnsIterValue>>
      : assns_metadata_type<ArtAssnsIterValue> {};

  } // namespace util
} // namespace lar

// back to what we were doing:
namespace proxy {
  namespace details {

    //--------------------------------------------------------------------------
    /// Traits for a association iterator.
    template <typename ArtAssnsIter>
    struct AssnsIterTraits : public lar::util::assns_traits<typename ArtAssnsIter::value_type> {
      using art_node_t = typename ArtAssnsIter::value_type;
      using node_t = AssnsNode<art_node_t>;
    }; // struct AssnsIterTraits

    /// Modified iterator returning a association node interface.
    /// The basic iterator interface is to the associated (right) _art_ pointer.
    template <typename ArtAssnsIter>
    class assns_node_iterator
      : public IteratorWrapperBase<assns_node_iterator<ArtAssnsIter>,
                                   ArtAssnsIter,
                                   typename AssnsIterTraits<ArtAssnsIter>::node_t> {
      using base_iterator_t = IteratorWrapperBase<assns_node_iterator<ArtAssnsIter>,
                                                  ArtAssnsIter,
                                                  typename AssnsIterTraits<ArtAssnsIter>::node_t>;

      using art_assns_iter_t = ArtAssnsIter;
      using traits_t = AssnsIterTraits<art_assns_iter_t>;

      /// Type of node for this association iterator.
      using AssnsNode_t = typename traits_t::node_t;
      using ArtAssnsNode_t = typename traits_t::art_node_t;

    public:
      using base_iterator_t::base_iterator_t;

      /// Constructor from a base iterator (explicitly allowed).
      assns_node_iterator(base_iterator_t const& from) : base_iterator_t(from) {}

      /// Returns the full information the iterator points to.
      AssnsNode_t const& info() const { return base_iterator_t::operator*(); }

      /// Returns the full information the iterator points to.
      AssnsNode_t const& operator()() const { return info(); }

      //--- BEGIN Access to the full association information -------------------
      /// @name Access to the full association information
      /// This interface is a replica of the one of `AssnsNode_t`.
      /// @{

      using main_t = typename AssnsNode_t::main_t;
      using value_t = typename AssnsNode_t::value_t;
      using data_t = typename AssnsNode_t::data_t;
      using mainptr_t = typename AssnsNode_t::mainptr_t;
      using valueptr_t = typename AssnsNode_t::valueptr_t;
      using dataptr_t = typename AssnsNode_t::dataptr_t;

      /// Returns the _art_ pointer to the associated value.
      valueptr_t valuePtr() const { return info().valuePtr(); }

      /// Returns the _art_ pointer to the associated value.
      value_t const& value() const { return info().value(); }

      /// Returns the _art_ pointer to the main value, key of the association.
      mainptr_t mainPtr() const { return info().mainPtr(); }

      /// Returns the main value, key of the association.
      main_t const& main() const { return info().main(); }

      // see the comments on AssnsNode for the need of all these decorations
      /// Returns whether this node type supports metadata.
      template <typename Node = AssnsNode_t>
      static constexpr bool hasMetadata()
      {
        return lar::util::assns_has_metadata_v<Node>;
      }

      /// Returns the pointer to the metadata on this association node.
      template <typename ArtNode = ArtAssnsNode_t>
      std::enable_if_t<hasMetadata<ArtNode>(), dataptr_t> dataPtr() const
      {
        return info().dataPtr();
      }

      /// Returns a reference to the metadata on this association node.
      template <typename ArtNode = ArtAssnsNode_t>
      std::enable_if_t<hasMetadata<ArtNode>(),
                       typename lar::util::assns_traits<ArtNode>::data_t const&>
      data() const
      {
        return info().data();
      }

      /// @}
      //--- END Access to the full association information ---------------------

      /*
       * Associations with metadata have an iterator with value type
       * art::AssnsNode, while the value for the ones without have just a
       * std::pair.
       * The std::pair returned by the one without metadata is a reference to an
       * element of the original collection.
       * Instead, the art::AssnsNode returned by the art::Assns iterator is a
       * temporary put together copying information from the original pair
       * collection and from the parallel metadata collection.
       * Therefore, in the first case we can return references to the existing
       * data, while in the latter we can't and we have to return the results by
       * value.
       * In this implementation we compromise and return always the data by
       * value; the values are art pointers, or plain pointers, so the copy
       * should not be extremely taxing. It is possible to change this, at the
       * cost of additional complexity of the implementation.
       */
      static AssnsNode_t const& transform(art_assns_iter_t const& v) { return makeAssnsNode(*v); }

    }; // class assns_node_iterator<>

    //--- END iterators for art::Assns -----------------------------------------

    //--------------------------------------------------------------------------
    //--- stuff for associated data (a form of auxiliary data)
    //--------------------------------------------------------------------------

    /// Interface providing begin and end iterator of a range.
    /// @tparam BoundaryIter iterator to the first of the range iterators.
    template <typename BoundaryIter>
    class BoundaryListRangeBase : private BoundaryIter {
      using boundary_iterator_t = BoundaryIter;

      /// Returns the iterator to the begin iterator.
      auto boundaryIter() const { return static_cast<BoundaryIter const&>(*this); }

    public:
      /// Constructor: copies the specified base iterator.
      BoundaryListRangeBase(boundary_iterator_t const& it) : boundary_iterator_t(it) {}

      /// Returns the begin iterator of the range.
      auto begin() const -> decltype(auto) { return *(boundaryIter()); }

      /// Returns the end iterator of the range (next to the begin iterator).
      auto end() const -> decltype(auto) { return *(std::next(boundaryIter())); }

    }; // BoundaryListRangeBase<>

    /// A `BoundaryListRangeBase` with a full container interface.
    template <typename BoundaryIter>
    class BoundaryListRange : public lar::CollectionView<BoundaryListRangeBase<BoundaryIter>> {
      // A CollectionView can't be constructed except from deriver classes;
      // we define here such a class.

      using rangebase_t = BoundaryListRangeBase<BoundaryIter>;
      using base_t = lar::CollectionView<rangebase_t>;

    public:
      using boundary_iterator_t = BoundaryIter;

      /// Constructor: from an iterator to the begin iterator.
      BoundaryListRange(boundary_iterator_t const& iBegin) : base_t(rangebase_t(iBegin)) {}

    }; // class BoundaryListRange<>

    /**
     * @brief Reinterprets a iterator to boundaries list as a range collection.
     * @tparam BoundaryIter type of iterator to boundary collection
     * @param iBegin iterator to the begin iterator of a range
     * @return a collection view (`lar::CollectionView`) of the range
     *
     * A range is conceptually defined as a sequence of data between a begin and
     * and end iterator.
     * The argument of this function is an iterator to the begin iterator of the
     * range. The begin iterator itself is obtained by dereferencing the
     * argument: `*iBegin`. The end iterator of the range is required to be
     * immediately after the begin iterator (@code *std::next(iBegin) @endcode).
     * This pair of iterators is exposed via the `lar::CollectionView` view,
     * that presents a vector-like interface.
     * For this to fully work, the data iterators (e.g., `*iBegin`) must comply
     * with the random-access iterator requirements.
     *
     * An example of `BoundaryIter` is the iterator to the list of boundaries
     * in `BoundaryList`: `BoundaryList::boundaries_t::const_iterator`.
     */
    template <typename BoundaryIter>
    BoundaryListRange<BoundaryIter> makeBoundaryListRange(BoundaryIter const& iBegin)
    {
      return {iBegin};
    }

    /**
     * @brief Iterator exposing elements of a boundary list as ranges.
     * @tparam BoundaryIter type of iterator to boundary iterators
     *
     * This iterator wraps an iterator (`BoundaryIter`) to a sequence of
     * iterators representing ranges of data.
     * As an example, let the data be a collection of real numbers:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * std::vector<float> data;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * and we have a list of iterators that define subranges within that data:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * using data_iterator_t = std::vector<data>::const_iterator;
     * std::vector<data_iterator_t> rangeStart;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Here `rangeStart` holds the begin iterator of each of the subranges,
     * plus a "global" end iterator, as in the convention described in
     * `lar::makeBoundaryListRange()` (which takes an iterator to `rangeStart`
     * as argument).
     * The `BoundaryIter` type in this example would be
     * @code std::vector<data_iterator_t>::const_iterator @endcode.
     *
     * When dereferenced, this iterator returns a view of the range currently
     * pointed. This view has a STL-vector-like interface (again, see
     * `lar::makeBoundaryListRange()`).
     */
    template <typename BoundaryIter>
    class BoundaryListRangeIterator
      : public IteratorWrapperBase<BoundaryListRangeIterator<BoundaryIter>, BoundaryIter> {
      using base_t = IteratorWrapperBase<BoundaryListRangeIterator<BoundaryIter>, BoundaryIter>;

    public:
      using boundary_iterator_t = BoundaryIter; ///< Type of boundary iterator.

      /// Type of range returned when dereferencing.
      using rangeview_t = decltype(makeBoundaryListRange(std::declval<boundary_iterator_t>()));

      using value_type = rangeview_t;
      using pointer = std::add_pointer_t<std::decay_t<value_type>>;
      using reference = std::add_lvalue_reference_t<std::decay_t<value_type>>;

      using base_t::base_t; // import constructors (explicitly allowed)

      /**
       * @brief Returns the pointed range.
       * @return the pointed range, as a view
       *
       * The returned value may be a view object, that is a pure interface
       * overlaid on a different data structure.
       * As such, it may be not copiable and need to be propagated by reference.
       */
      static auto transform(BoundaryIter const& iter) { return makeBoundaryListRange(iter); }

    }; // class BoundaryListRangeIterator<>

    /**
     * @brief Builds and keeps track of internal boundaries in a sequence.
     * @tparam Iter type of iterators to the original sequence
     *
     * This class manages a sequence of boundary iterators defining the
     * beginning of contiguous subsequences. Each iterator marks the begin of a
     * subsequence, whose end is marked by the beginning of the next one.
     * The last iterator in the boundary list marks the end of the last
     * subsequence, but it does not mark the beginning of a following one.
     * Therefore, for a list of _N_ subsequences there will be _N + 1_ boundary
     * iterators in the list: _N_ marking the beginning of the respective
     * subsequences, plus another marking the end of the last subsequence.
     *
     * It is likely that the first iterator in the boundary list is the begin
     * iterator of the underlying sequence (usually a container) being
     * partitioned, while the last one is the end iterator of that sequence.
     *
     * This is a data class which does not contain any logic to define the
     * subsequences, but rather acquires the result of an algorithm which is
     * expected to have established which the boundaries are.
     *
     * The underlying representation of the class is a random-access sequence
     * of boundaries. The exposed value, `range_t`, is a range of data elements
     * (a view with the interface of a random access container), which is
     * internally represented as a single iterator pointing to the begin
     * iterator of the range.
     */
    template <typename Iter>
    class BoundaryList {
      using boundarylist_t = BoundaryList<Iter>;

    public:
      using data_iterator_t = Iter;
      using boundaries_t = std::vector<data_iterator_t>;

      /// Iterator on the ranges contained in the collection.
      using range_iterator_t = BoundaryListRangeIterator<typename boundaries_t::const_iterator>;

      /// Structure holding begin and end iterator for a single range.
      // BoundaryListRange<data_iterator_t> const&
      using range_ref_t = typename range_iterator_t::value_type;

      /// Range object directly containing the boundary iterators.
      using range_t = lar::RangeAsCollection_t<data_iterator_t>;

      /// Constructor: steals the specified boundary list.
      explicit BoundaryList(boundaries_t&& boundaries) : boundaries(std::move(boundaries))
      {
        assert(this->boundaries.size() >= 1);
      }

      /// Returns the number of ranges contained in the list.
      std::size_t nRanges() const { return boundaries.size() - 1; }
      /// Returns the begin iterator of the `i`-th range (end if overflow).
      data_iterator_t const& rangeBegin(std::size_t i) const
      {
        return boundaries[std::min(i, nRanges())];
      }
      /// Returns the end iterator of the `i`-th range (end if overflow).
      data_iterator_t const& rangeEnd(std::size_t i) const { return rangeBegin(i + 1); }

      /// Returns the number of ranges contained in the list.
      std::size_t size() const { return nRanges(); }
      /// Returns the begin iterator of the first range.
      range_iterator_t begin() const { return {boundaries.begin()}; }
      /// Returns the end iterator of the last range.
      range_iterator_t end() const { return {std::prev(boundaries.end())}; }
      /**
       * @brief Returns the specified range.
       * @param i index of the range to be returned
       * @return a proxy object with container interface
       * @see `range()`
       *
       * The returned object exposes the range as a random access container
       * interface.
       *
       * Internally, it refers to the relevant information from this
       * `BoundaryList` object, and therefore it becomes invalid when this
       * `BoundaryList` object is destroyed.
       * If this is not acceptable, use `range()` instead.
       */
      range_ref_t rangeRef(std::size_t i) const { return {std::next(boundaries.begin(), i)}; }
      /**
       * @brief Returns the specified range in an object holding the iterators.
       * @param i index of the range to be returned
       * @return a new object with container interface
       * @see `rangeRef()`
       *
       * The returned object contains copies of the begin and end iterators of
       * the range. This object is self-contained and valid even after this
       * BoundaryList object is destroyed.
       *
       * Note the content of the range itself is _not_ copied: just the boundary
       * iterators of the range are.
       */
      range_t range(std::size_t i) const
      {
        return lar::makeCollectionView(rangeBegin(i), rangeEnd(i));
      }

      /// Returns the begin iterator of the `i`-th range (unchecked).
      /// @see `range()`
      auto operator[](std::size_t i) const -> decltype(auto) { return range(i); }

    private:
      /// Begin iterator of each range, plus end iterator of whole sequence.
      boundaries_t boundaries;

    }; // class BoundaryList

    /**
     * @brief Object to draft associated data interface.
     * @tparam Main type of the main associated object (one)
     * @tparam Aux type of the additional associated objects (many)
     * @tparam Metadata type of metadata in the association (default: `void`)
     * @tparam Tag tag this data is labeled with
     *
     * Allows:
     *  * random access (no index check guarantee)
     *  * forward iteration
     *
     * Construction is not part of the interface.
     *
     * The `AssociatedData` object, on creation, finds the borders surrounding
     * the associated `Aux` objects for each `Main` one, and keep a record of
     * them (this is actually delegated to `BoundaryList` class).
     * The `AssociatedData` object also provides a container-like view of this
     * information, where each element in the container is associated to a
     * single `Main` and it is a container (actually, another view) of `Right`.
     * Both levels of containers are random access, so that the set of `Right`
     * objects associated to a `Left` can be accessed by index, and the `Right`
     * objects within can be accessed with `Right` index in the `Left`.
     */
    template <typename Main, typename Aux, typename Metadata /* = void */, typename Tag /* = Aux */
              >
    class AssociatedData {
      using This_t = AssociatedData<Main, Aux, Metadata, Tag>; ///< This type.

    public:
      /// Type of _art_ association.
      using assns_t = art::Assns<Main, Aux, Metadata>;

    private:
      using associated_data_iterator_t = assns_node_iterator<lar::util::assns_iterator_t<assns_t>>;
      //  = tuple_element_iterator<1U, lar::util::assns_iterator_t<assns_t>>;

    public:
      using tag = Tag; ///< Tag of this association proxy.

      using group_ranges_t = BoundaryList<associated_data_iterator_t>;

      /// Type of collection of auxiliary data associated with a main item.
      using auxiliary_data_t = util::add_tag_t<typename group_ranges_t::range_t, tag>;

      // constructor is not part of the interface
      AssociatedData(group_ranges_t&& groups) : fGroups(std::move(groups)) {}

      /// Returns an iterator pointing to the first associated data range.
      auto begin() const -> decltype(auto) { return fGroups.begin(); }

      /// Returns an iterator pointing past the last associated data range.
      auto end() const -> decltype(auto) { return fGroups.end(); }

      /// Returns the range with the specified index (no check performed).
      auto getRange(std::size_t i) const -> decltype(auto)
      {
        return util::makeTagged<tag>(fGroups.range(i));
      }

      /// Returns the range with the specified index (no check performed).
      auto operator[](std::size_t index) const -> decltype(auto)
      {
        static_assert(std::is_convertible<decltype(getRange(index)), auxiliary_data_t>(),
                      "Inconsistent data types.");
        return getRange(index);
      }

      /// Returns whether this data is labeled with the specified tag.
      template <typename TestTag>
      static constexpr bool hasTag()
      {
        return std::is_same<TestTag, tag>();
      }

    private:
      group_ranges_t fGroups;

    }; // class AssociatedData<>

    //--------------------------------------------------------------------------

  } // namespace details

  //----------------------------------------------------------------------------
  //@{
  /**
   * @brief Processes and returns an associated data object.
   * @tparam Tag the tag labelling this associated data
   *             (if omitted: second type of the association: `right_t`)
   * @tparam Assns type of association to be processed
   * @param assns association object to be processed
   * @param minSize minimum number of entries in the produced association data
   * @return a new `AssociatedData` filled with associations from `tag`
   *
   * The content of the association object must fulfill the requirements of
   * @ref LArSoftProxyDefinitionOneToManySeqAssn "one-to-many sequential association".
   * The `Assns` type is expected to be a `art::Assns` instance. At least,
   * the `Assns` type is required to have `left_t` and `right_t` definitions
   * representing respectively the main data type and the associated one, and
   * to respond to `begin()` and `end()` functions. The iterated object must
   * also respond to `std::get<0>()` with a `art::Ptr<left_t>` and to
   * `std::get<1>()` with a `art::Ptr<right_t>`.
   *
   * Elements in the main collection not associated with any object will be
   * recorded as such. If there is information for less than `minSize` main
   * objects, more records will be added to mark the missing objects as not
   * associated to anything.
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * art::Assns<recob::Track, recob::Hit> trackHitAssns;
   * // ...
   * auto assData = proxy::makeAssociatedData(assns);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will have `assData` tagged as `recob::Hit`.
   */
  template <typename Tag, typename Assns>
  auto makeAssociatedData(Assns const& assns, std::size_t minSize = 0);

  template <typename Assns>
  auto makeAssociatedData(Assns const& assns, std::size_t minSize = 0)
  {
    return makeAssociatedData<typename Assns::right_t>(assns, minSize);
  }
  //@}

  //@{
  /**
   * @brief Creates and returns an associated data object.
   * @tparam Tag the tag labelling this associated data
   *             (if omitted: second type of the association: `right_t`)
   * @tparam MainColl type of the main collection object
   * @tparam Assns type of the association object
   * @param mainColl the main collection object
   * @param assns association data object
   * @return a new `AssociatedData` wrapping the information in `assns`
   * @see `makeAssociatedData(Assns const&, std::size_t)`
   *
   * This function operates like
   * `makeAssociatedData(Assns const&, std::size_t)`, where the size is
   * extracted from the main data collection.
   */
  template <typename Tag, typename MainColl, typename Assns>
  auto makeAssociatedData(MainColl const& mainColl, Assns const& assns)
  {
    return makeAssociatedData<Tag>(assns, mainColl.size());
  }

  template <typename MainColl, typename Assns>
  auto makeAssociatedData(MainColl const& mainColl, Assns const& assns)
  {
    return makeAssociatedData<typename Assns::right_t>(mainColl, assns);
  }
  //@}

  //----------------------------------------------------------------------------

} // namespace proxy

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace proxy {

  namespace details {

    //--------------------------------------------------------------------------
    //--- associationRangeBoundaries() implementation
    //--------------------------------------------------------------------------
    template <std::size_t GroupKey, typename Iter>
    typename BoundaryList<Iter>::boundaries_t
    associationRangesImpl(Iter begin, Iter end, std::size_t expectedSize /* = 0 */)
    {
      constexpr auto KeyIndex = GroupKey;

      auto extractKey = [](auto const& assn) { return std::get<KeyIndex>(assn).key(); };

      typename BoundaryList<Iter>::boundaries_t boundaries;
      boundaries.reserve(expectedSize + 1);
      boundaries.push_back(begin);
      std::size_t current = 0;
      for (auto it = begin; it != end; ++it) {
        auto const key = extractKey(*it);
        if (key == current) continue;
        if (key < current) {
          auto index = std::distance(begin, it);
          throw std::runtime_error(
            "associationRanges() got input element #" + std::to_string(index - 1) + " with key " +
            std::to_string(current) + " and the next with key " + std::to_string(key) + "!");
        }
        boundaries.insert(boundaries.end(), key - current, it);
        current = key;
      } // for
      boundaries.push_back(end);
      return boundaries;
    } // associationRangesImpl()

    //--------------------------------------------------------------------------
    template <std::size_t GroupKey, typename Iter>
    auto associationRangeBoundaries(Iter begin, Iter end)
    {
      return associationRangesImpl<GroupKey, Iter>(begin, end);
    }

    //--------------------------------------------------------------------------
    template <std::size_t GroupKey, typename Iter>
    auto associationRangeBoundaries(Iter begin, Iter end, std::size_t n)
    {
      auto boundaries = associationRangesImpl<GroupKey, Iter>(begin, end, n);
      if (boundaries.size() <= n) {
        boundaries.insert(boundaries.end(), n + 1 - boundaries.size(), boundaries.back());
        assert(boundaries.size() == (n + 1));
      }
      return boundaries;
    } // associationRangeBoundaries(Iter, Iter, std::size_t)

    //--------------------------------------------------------------------------
    /**
     * @brief Groups associations by the first key.
     * @tparam GroupKey index of the key in the tuple pointed by the iterator
     * @tparam Iter type of iterators delimiting the data (same type required)
     * @param begin iterator to the first association in the list
     * @param end iterator past the last association in the list
     * @return a list of range boundaries marking the different groups.
     * @throw std::runtime_error if input key is not monotonic
     *
     * The input iterators are expected to point to a tuple-like structure whose
     * key element can be accessed as `std::get<GroupKey>()` and is an _art_
     * pointer of some sort.
     *
     * The index of the grouping key is expected to be monotonically increasing.
     * Gaps are supported except that at the end: if e.g. an association of 5
     * keys associates objects to only elements #0, #1 and #3, the resulting
     * list will cover 4 ranges for elements #0 to #3 included, but excluding
     * the end elements, the existence of which can't be inferred from the
     * association list in input. In this example, the range #2 will exist and
     * be empty. To enforce a minimum number of elements, use
     * `associationRanges(Iter, Iter, std::size_t)`.
     */
    template <std::size_t GroupKey, typename Iter>
    BoundaryList<Iter> associationRanges(Iter begin, Iter end)
    {
      return BoundaryList<Iter>(associationRangeBoundaries<GroupKey>(begin, end));
    }

    /**
     * @brief Groups associations by the first key.
     * @tparam GroupKey index of the key in the tuple pointed by the iterator
     * @tparam Iter type of iterators delimiting the data (same type required)
     * @param begin iterator to the first association in the list
     * @param end iterator past the last association in the list
     * @param n minimum number of ranges to be produced.
     * @return a list of range boundaries marking the different groups.
     * @throw std::runtime_error if input key is not monotonic
     * @see `associationRanges(Iter, Iter)`
     *
     * This function operates almost like `associationRanges(Iter, Iter)`.
     * The only difference is that at least `n` ranges are guaranteed to be
     * produced: if the input defines less than `n`, the missing ones will be
     * added at the end, empty.
     * This allows to work around the possibility of empty ranges at the end,
     * which the `associationRanges(Iter, Iter)` algorithm does not support.
     */
    template <std::size_t GroupKey, typename Iter>
    BoundaryList<Iter> associationRanges(Iter begin, Iter end, std::size_t n)
    {
      return BoundaryList<Iter>(associationRangeBoundaries<GroupKey>(begin, end, n));
    }

    //--------------------------------------------------------------------------

  } // namespace details

  //----------------------------------------------------------------------------
  template <typename Tag, typename Assns>
  auto makeAssociatedData(Assns const& assns, std::size_t minSize /* = 0 */)
  {
    using Main_t = typename Assns::left_t;
    using Aux_t = typename Assns::right_t;
    using Metadata_t = lar::util::assns_metadata_t<Assns>;
    using AssociatedData_t = details::AssociatedData<Main_t, Aux_t, Metadata_t, Tag>;

    // associationRangeBoundaries() produces iterators to association elements,
    // (i.e. tuples)
    using std::begin;
    using std::end;
    auto ranges = details::associationRangeBoundaries<0U>(begin(assns), end(assns), minSize);
    // we convert those iterators into iterators to the right associated item
    // (it takes a few steps)
    using group_ranges_t = typename AssociatedData_t::group_ranges_t;
    return AssociatedData_t(
      group_ranges_t(typename group_ranges_t::boundaries_t(ranges.begin(), ranges.end())));
  } // makeAssociatedDataFrom(assns)

  //----------------------------------------------------------------------------

} // namespace proxy
//------------------------------------------------------------------------------

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_ASSOCIATEDDATA_H
