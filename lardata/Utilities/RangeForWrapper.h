/**
 * @file    lardata/Utilities/RangeForWrapper.h
 * @brief   Utility function to enable range-for on different type iterators
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    December 12, 2016
 *
 * This header will be needed likely until C++20, since range library is not
 * providing some definitions that are needed for range for loops.
 */

#ifndef LARDATA_UTILITIES_RANGEFORWRAPPER_H
#define LARDATA_UTILITIES_RANGEFORWRAPPER_H

// Boost libraries
#include "boost/variant.hpp"

// C/C++ standard libraries
#include <stdexcept> // std::logic_error
#include <utility> // std::move(), std::declval(), ...
#include <iterator> // std::iterator_traits
#include <type_traits> // std::is_same<>, std::enable_if_t<>, ...


namespace util {

  namespace details {

    /// Iterator wrapping one of two types of iterators.
    /// @see RangeForWrapperBox
    template <typename BeginIter, typename EndIter>
    class RangeForWrapperIterator {

      using traits_t = std::iterator_traits<BeginIter>;

        public:
      using begin_t = BeginIter; ///< Type of begin iterator we can store.
      using end_t = EndIter; ///< Type of end iterator we can store.
      using this_t = RangeForWrapperIterator<begin_t, end_t>; ///< This class.

      /// @{
      /// @brief Iterator traits, imported from the wrapped begin iterator.
      using difference_type   = typename traits_t::difference_type;
      using value_type        = typename traits_t::value_type;
      using pointer           = typename traits_t::pointer;
      using reference         = typename traits_t::reference;
      //
      // This wrapper fully supports up to bidirectional iterators;
      // if the wrapped iterator is a random or contiguous iterator,
      // the wrapper will still expose only a bidirectional iterator interface.
      // Supporting random access is possible, but writing a proper unit test
      // is tedious... open a feature request if needed.
      //
      using iterator_category = std::conditional_t<
        std::is_base_of<std::bidirectional_iterator_tag, typename traits_t::iterator_category>::value,
        std::bidirectional_iterator_tag,
        typename traits_t::iterator_category
        >;
      /// @}

      /// Constructor: initializes with a end-type default-constructed iterator.
      explicit RangeForWrapperIterator()
        : fIter(end_t{})
        {}

      /// Constructor: initializes with a begin-type iterator.
      explicit RangeForWrapperIterator(begin_t&& begin)
        : fIter(std::move(begin))
        {}

      /// Constructor: initializes with a end-type iterator.
      explicit RangeForWrapperIterator(end_t&& end)
        : fIter(std::move(end))
        {}

      /// Returns the pointed value (just like the original iterator).
      reference operator*() const
        { return boost::apply_visitor(Dereferencer(), fIter); }

      /// Returns the pointed value (just like the original iterator).
      pointer operator->() const
        { return boost::apply_visitor(MemberAccessor(), fIter); }

      /// Increments the iterator (prefix operator).
      this_t& operator++()
        { boost::apply_visitor(Incrementer(), fIter); return *this; }

      /// Decrements the iterator (prefix operator).
      this_t& operator--()
        { boost::apply_visitor(Decrementer(), fIter); return *this; }

      /// Increments the iterator (postfix operator).
      this_t operator++(int)
        { auto old = *this; this_t::operator++(); return old; }

      /// Decrements the iterator (postfix operator).
      this_t operator--(int)
        { auto old = *this; this_t::operator--(); return old; }

      /// Returns whether the other iterator is not equal to this one.
      bool operator!=(this_t const& other) const
        { return boost::apply_visitor(Comparer(), fIter, other.fIter); }

      /// Returns whether the other iterator is equal to this one.
      bool operator==(this_t const& other) const
        { return !(this->operator!=(other)); }


      reference operator[] (difference_type offset) const
        { return boost::apply_visitor(IndexAccessor(offset), fIter); }

      difference_type operator- (this_t const& other) const
        { return boost::apply_visitor(Difference(), fIter, other.fIter); }

        private:
      static_assert(!std::is_same<begin_t, end_t>::value,
        "RangeForWrapperIterator requires two different iterator types."
        );

      boost::variant<begin_t, end_t> fIter; ///< The actual iterator we store.

      //
      // We opt for allowing all the operations if the underlying operators do.
      // While it is true that, for example, an end iterator should not be
      // dereferenced, if it's bidirectional, its operator--() may make it
      // useful, but its type will still be the one of the end iterator.
      // Therefore we don't judge by the type, but by the action.
      //

      /// Visitor to dereference an iterator.
      struct Dereferencer: public boost::static_visitor<reference> {

        template <typename Iter>
        auto operator() (Iter& iter) const -> decltype(auto)
          { return DereferencerImpl<reference, Iter>::dereference(iter); }

          private:
        template <typename Result, typename Iter, typename = void>
        struct DereferencerImpl;

      }; // Dereferencer

      /// Visitor to access a data member of the pointed class.
      struct MemberAccessor: public boost::static_visitor<pointer> {

        template <typename Iter>
        auto operator() (Iter& iter) const -> decltype(auto)
          { return MemberAccessorImpl<pointer, Iter>::access(iter); }

          private:
        template <typename Result, typename Iter, typename = void>
        struct MemberAccessorImpl;

      }; // MemberAccessor

      /// Visitor to increment an iterator.
      struct Incrementer: public boost::static_visitor<> {

        template <typename Iter>
        void operator() (Iter& iter) const
          { IncrementerImpl<Iter>::increment(iter); }

          private:
        template <typename Iter, typename = void>
        struct IncrementerImpl;

      }; // Incrementer

      /// Visitor to decrement an iterator.
      struct Decrementer: public boost::static_visitor<> {

        template <typename Iter>
        void operator() (Iter& iter) const
          { DecrementerImpl<Iter>::decrement(iter); }

          private:
        template <typename Iter, typename = void>
        struct DecrementerImpl;

      }; // Decrementer


      /// Visitor to compare iterators (returns whether they differ).
      struct Comparer: public boost::static_visitor<bool> {

        template <typename A, typename B>
        bool operator() (A const& left, B const& right) const
          { return ComparerImpl<A, B>::compare(left, right); }

          private:
        template <typename A, typename B, typename = void>
        struct ComparerImpl;

      }; // Comparer

      /// Visitor to access element by index.
      struct IndexAccessor: public boost::static_visitor<reference> {

        difference_type offset;

        IndexAccessor(difference_type offset): offset(offset) {}

        template <typename Iter>
        bool operator() (Iter& iter) const
          { return IndexAccessorImpl<reference, Iter>(offset).access(iter); }

          private:
        template <typename Result, typename Iter, typename = void>
        struct IndexAccessorImpl;

      }; // IndexAccessor

      /// Visitor to compare iterators (returns whether they differ).
      struct Difference: public boost::static_visitor<difference_type > {

        template <typename A, typename B>
        difference_type operator() (A const& minuend, B const& subtrahend) const
          { return DifferenceImpl<A, B>::subtract(minuend, subtrahend); }

          private:
        template <typename A, typename B, typename = void>
        struct DifferenceImpl;

      }; // Difference

    }; // class RangeForWrapperIterator<>



    /// Class defining types and traits for RangeForWrapperBox
    template <typename RangeRef>
    struct RangeForWrapperTraits {

      using RangeRef_t = RangeRef; ///< Type of the stored reference.

      ///< Type of the stored range (constantness is preserved).
      using Range_t = std::remove_reference_t<RangeRef_t>;

      /// Extractor of the begin iterator from a range.
      static auto extractBegin(RangeRef_t range)
        { using namespace std; return begin(range); }

      /// Extracts the end iterator from a range object.
      static auto extractEnd(RangeRef_t range)
        { using namespace std; return end(range); }

      /// Type of wrapped begin iterator.
      using BeginIter_t = decltype(extractBegin(std::declval<RangeRef_t>()));

      /// Type of wrapped end iterator.
      using EndIter_t = decltype(extractEnd(std::declval<RangeRef_t>()));

      /// True if the range has iterators of the same type.
      static constexpr bool sameIteratorTypes
        = std::is_same<BeginIter_t, EndIter_t>();

      /// Type of wrapper iterators (same for begin and end iterators).
      using Iterator_t = RangeForWrapperIterator<BeginIter_t, EndIter_t>;

      using value_type = typename BeginIter_t::value_type;
      using size_type = std::size_t;
      using difference_type = typename BeginIter_t::difference_type;
      using reference = typename BeginIter_t::value_type;
      using pointer = typename BeginIter_t::pointer;

    }; // class RangeForWrapperTraits<>


    /**
     * @brief Class offering begin/end iterators of the same type out of a range
     *        of iterators of different types.
     * @tparam RangeRef type of reference to be stored (constantness embedded)
     *
     * The class steals (moves) the value if `RangeRef` is a rvalue reference
     * type, while it just references the original one otherwise.
     */
    template <typename RangeRef>
    class RangeForWrapperBox {

      static_assert(std::is_reference<RangeRef>::value,
        "RangeForWrapperBox requires a reference type.");

      using Traits_t = RangeForWrapperTraits<RangeRef>;

        public:

      // Import traits
      using RangeRef_t = typename Traits_t::RangeRef_t;
      using Range_t = typename Traits_t::Range_t;

      /// Type of wrapper iterators (same for begin and end iterators).
      using Iterator_t = typename Traits_t::Iterator_t;

      /// Type of number of stored elements.
      using size_type = typename Traits_t::size_type;

      /// Type of difference between element positions.
      using difference_type = typename Traits_t::difference_type;

      /// Type of value contained in the wrapped sequence.
      using value_type = typename Traits_t::value_type;


      /// Constructor: references the specified range (lvalue reference).
      RangeForWrapperBox(Range_t& range)
        : fRange(range)
        {}

      /// Constructor: references the specified range (rvalue reference).
      RangeForWrapperBox(Range_t&& range)
        : fRange(std::move(range))
        {}

      /// Returns a begin-of-range iterator.
      Iterator_t begin() const
        { return Iterator_t(wrappedBegin()); }

      /// Returns a end-of-range iterator.
      Iterator_t end() const
        { return Iterator_t(wrappedEnd()); }

      /// @{
      /// @name Reduced container interface.

      auto size() const { return std::distance(begin(), end()); }

      bool empty() const { return !(wrappedBegin() != wrappedEnd()); }

      auto operator[] (difference_type index) const -> decltype(auto)
        { return wrappedBegin()[index]; }

      /// @}


        private:

      struct DataBox {

        using Stored_t = std::conditional_t<
          std::is_rvalue_reference<RangeRef_t>::value,
          std::remove_reference_t<RangeRef_t>,
          RangeRef_t
          >;
        using Data_t = std::remove_reference_t<Stored_t>;

        Stored_t data;

        // only one of these is valid...
        DataBox(Data_t& data): data(data) {}
        DataBox(Data_t&& data): data(std::move(data)) {}

        operator RangeRef_t() const { return RangeRef_t(data); }
        operator RangeRef_t() { return RangeRef_t(data); }

      }; // DataBox

      DataBox fRange; ///< A reference to the original range.

      auto wrappedBegin() const -> decltype(auto)
        { return Traits_t::extractBegin(static_cast<RangeRef_t>(fRange)); }
      auto wrappedEnd() const -> decltype(auto)
        { return Traits_t::extractEnd(static_cast<RangeRef_t>(fRange)); }

    }; // class RangeForWrapperBox<>


    /// Tag for internal use.
    struct SameIterTag {};

    /// Tag for internal use.
    struct DiffIterTag {};

    /// Wraps an object for use in a range-for loop
    /// (same iterator types: pass through)
    // the return type decltype(auto) is necessary to preserve the forwarded
    // referenceness
    template <
      typename BaseRange,
      bool SameIteratorsType
        = details::RangeForWrapperTraits<std::decay_t<BaseRange>>::sameIteratorTypes
      >
    struct WrapRangeForDispatcher;

    // Template specialization for same iterator types
    template <typename BaseRange>
    struct WrapRangeForDispatcher<BaseRange, true> {

      using BaseRange_t = std::decay_t<BaseRange>;

      static BaseRange_t wrap(BaseRange_t&& range) { return std::move(range); }
      static BaseRange_t& wrap(BaseRange_t& range) { return range; }
      static BaseRange_t const& wrap(BaseRange_t const& range) { return range; }
    }; // WrapRangeForDispatcher<BaseRange, true>


    // Template specialization for different-iterator types
    template <typename BaseRange>
    struct WrapRangeForDispatcher<BaseRange, false> {
      template <typename Range>
      static auto wrap(Range&& range)
        {
          return RangeForWrapperBox<decltype(range)>
            (static_cast<decltype(range)>(range));
        }
    }; // WrapRangeForDispatcher<BaseRange, false>

  } // namespace details


  /**
   * @brief Wraps an object for use in a range-for loop
   * @tparam Range type of range object (anything with begin() and end())
   * @param range instance of the range object to be wrapped
   *
   * This is necessary only when the argument provides different types for
   * the begin-of-range and end-of-range iterators.
   * This is also superfluous for compilers adhering to C++ 2017 standard,
   * which accepts iterators of different types by requirement.
   * Example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * Range data; // initialization
   * for (auto&& value: util::wrapRangeFor(data)) // ...
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * where data is supposed to gave begin and end iterators of different types.
   */
  template <typename Range>
  auto wrapRangeFor(Range&& range) -> decltype(auto)
    {
      return details::WrapRangeForDispatcher<Range>::wrap
        (std::forward<Range>(range));
    }


  /// Tag marking the use of RangeForWrapperBox
  struct RangeForWrapperTag {};

  /// Constant to be used with
  /// `operator|(Range&&, details::RangeForWrapperTag)`.
  constexpr RangeForWrapperTag range_for;

  /**
   * @brief Transforms a range so that it can be used in a range-for loop
   * @tparam Range the type of range to be transformed
   * @param range the range to be transformed
   * @return an equivalent range object to be used in a range-for loop
   *
   * This is necessary only when the argument provides different types for
   * the begin-of-range and end-of-range iterators.
   * This is also superfluous for compilers adhering to C++ 2017 standard,
   * which accepts iterators of different types by requirement.
   * Example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * Range data; // initialization
   * for (auto&& value: data | util::range_for) // ...
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * where data is supposed to gave begin and end iterators of different types.
   */
  template <typename Range>
  auto operator| (Range&& range, RangeForWrapperTag) -> decltype(auto)
    { return wrapRangeFor(std::forward<Range>(range)); }


} // namespace util


//------------------------------------------------------------------------------
//---  template implementation
//------------------------------------------------------------------------------
namespace util {

  namespace details {

    //--------------------------------------------------------------------------
    template <typename T>
    struct is_type: public std::true_type {};

    template <typename T>
    constexpr bool is_type_v = is_type<T>();


    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename A, typename B, typename /* = void */>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Comparer::ComparerImpl {
      // this would be worth a static_assert(), but apparently boost::variant
      // visitor instantiates it even when it's not called
      static bool compare(A const&, B const&)
        { throw std::logic_error("These iterators can't be compared!"); }
    }; //

    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename A, typename B>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Comparer::ComparerImpl<
      A, B, std::enable_if_t<
        std::is_convertible
          <decltype(std::declval<A>() != std::declval<B>()), bool>::value
        >
      >
    {
      static bool compare(A const& left, B const& right)
        { return left != right; }
    }; //


    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Result, typename Iter, typename /* = void */>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Dereferencer::DereferencerImpl {
      // this would be worth a static_assert(), but apparently boost::variant
      // visitor instantiates it even when it's not called
      [[noreturn]] static Result dereference(Iter const&)
        { throw std::logic_error("This iterator can't be dereferenced!"); }
    }; //

    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Result, typename Iter>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Dereferencer::DereferencerImpl<
      Result, Iter, std::enable_if_t<is_type_v<decltype(*(std::declval<Iter>()))>>
      >
    {
      static Result dereference(Iter const& iter)
        { return *iter; }
    }; //


    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Result, typename Iter, typename /* = void */>
    struct RangeForWrapperIterator<BeginIter, EndIter>::MemberAccessor::MemberAccessorImpl {
      // this would be worth a static_assert(), but apparently boost::variant
      // visitor instantiates it even when it's not called
      [[noreturn]] static Result access(Iter const&)
        { throw std::logic_error("This iterator can't be dereferenced!"); }
    }; //

    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Result, typename Iter>
    struct RangeForWrapperIterator<BeginIter, EndIter>::MemberAccessor::MemberAccessorImpl<
      Result, Iter, std::enable_if_t<is_type_v<decltype(std::declval<Iter>().operator->())>>
      >
    {
      static Result access(Iter const& iter)
        { return iter.operator->(); }
    }; //


    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Iter, typename /* = void */>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Incrementer::IncrementerImpl {
      // this would be worth a static_assert(), but apparently boost::variant
      // visitor instantiates it even when it's not called
      [[noreturn]] static void increment(Iter&)
        { throw std::logic_error("This iterator can't be incremented!"); }
    }; //

    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Iter>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Incrementer::IncrementerImpl<
      Iter, std::enable_if_t<is_type_v<decltype(++(std::declval<Iter>()))>>
      >
    {
      static void increment(Iter& iter)
        { ++iter; }
    }; //


    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Iter, typename /* = void */>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Decrementer::DecrementerImpl {
      // this would be worth a static_assert(), but apparently boost::variant
      // visitor instantiates it even when it's not called
      [[noreturn]] static void decrement(Iter&)
        { throw std::logic_error("This iterator can't be decremented!"); }
    }; //

    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Iter>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Decrementer::DecrementerImpl<
      Iter, std::enable_if_t<is_type_v<decltype(--(std::declval<Iter>()))>>
      >
    {
      static void decrement(Iter& iter)
        { --iter; }
    }; //


    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Result, typename Iter, typename /* = void */>
    struct RangeForWrapperIterator<BeginIter, EndIter>::IndexAccessor::IndexAccessorImpl {
      // this would be worth a static_assert(), but apparently boost::variant
      // visitor instantiates it even when it's not called

      IndexAccessorImpl(difference_type) {}

      [[noreturn]] Result access(Iter const&) const
        { throw std::logic_error("This iterator can't be indexed!"); }
    }; //

    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename Result, typename Iter>
    struct RangeForWrapperIterator<BeginIter, EndIter>::IndexAccessor::IndexAccessorImpl<
      Result, Iter, std::enable_if_t<is_type_v<decltype((std::declval<Iter>())[0])>>
      >
    {
      difference_type offset;

      IndexAccessorImpl(difference_type offset): offset(offset) {}

      Result dereference(Iter const& iter) const
        { return iter[offset]; }
    }; //


    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename A, typename B, typename /* = void */>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Difference::DifferenceImpl {
      // this would be worth a static_assert(), but apparently boost::variant
      // visitor instantiates it even when it's not called
      static difference_type subtract(A const&, B const&)
        { throw std::logic_error("These iterators can't be subtracted!"); }
    }; //

    //--------------------------------------------------------------------------
    template <typename BeginIter, typename EndIter>
    template <typename A, typename B>
    struct RangeForWrapperIterator<BeginIter, EndIter>::Difference::DifferenceImpl<
      A, B, std::enable_if_t<
        std::is_convertible<
          decltype(std::declval<A>() - std::declval<B>()),
          typename RangeForWrapperIterator<BeginIter, EndIter>::difference_type>::value
        >
      >
    {
      static difference_type subtract(A const& minuend, B const& subtrahend)
        { return minuend - subtrahend; }
    }; //


    //--------------------------------------------------------------------------


  } // namespace details
} // namespace util


//------------------------------------------------------------------------------


#endif // LARDATA_UTILITIES_RANGEFORWRAPPER_H
