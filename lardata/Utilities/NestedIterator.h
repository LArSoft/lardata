/**
 * @file   NestedIterator.h
 * @brief  Iterators recursing though nested collections
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   August 22th, 2014
 *
 * This header is currently a draft with reduced (essential) functionality.
 * Its ambition is to become general enough to be used transparently.
 */

#ifndef NESTEDITERATOR_H
#define NESTEDITERATOR_H

// interface include
#include <iterator>    // std::forward_iterator_tag
#include <type_traits> // std::true_type, std::false_type, std::conditional

namespace lar {

  namespace details {
    namespace type_traits {

      //------------------------------------------------------------------------
      // declarations of helpers for has_const_iterator() function
      template <typename T, bool>
      struct has_const_iterator_struct;

      template <typename T>
      constexpr auto has_const_iterator_helper(T* = nullptr)
        -> decltype(typename T::const_iterator(), bool());

      // Used as fallback when SFINAE culls the template method
      constexpr bool has_const_iterator_helper(...);
      //------------------------------------------------------------------------

    } // namespace type_traits
  }   // namespace details

  template <typename T>
  struct has_const_iterator : public details::type_traits::has_const_iterator_struct<
                                T,
                                details::type_traits::has_const_iterator_helper((T*)(nullptr))> {};

  /// Functor returning the object specified as argument
  template <typename T>
  class Identity;

  template <typename T>
  class PairSecond;

  /// Internal helper class: actual implementation of nested iterator
  template <typename ITER, typename INNERCONTEXTRACT = Identity<typename ITER::value_type>>
  class deep_const_fwd_iterator_nested;

  /// Deep iterator
  /// @todo write documentation about how to use it
  /*  template <typename CITER>
  using deep_fwd_const_iterator = std::conditional<
    has_const_iterator<typename CITER::value_type>(),
    deep_const_fwd_iterator_nested<CITER>,
    CITER
    >;
*/

  template <typename CITER, typename INNERCONTEXTRACT = Identity<typename CITER::value_type>>
  using double_fwd_const_iterator = deep_const_fwd_iterator_nested<CITER, INNERCONTEXTRACT>;

} // namespace lar

namespace std {
  template <typename CITER, typename INNERCONTEXTRACT>
  void swap(lar::deep_const_fwd_iterator_nested<CITER, INNERCONTEXTRACT>& a,
            lar::deep_const_fwd_iterator_nested<CITER, INNERCONTEXTRACT>& b);

} // namespace std

//------------------------------------------------------------------------------

namespace lar {

  namespace details {
    namespace type_traits {

      //------------------------------------------------------------------------
      // helpers for has_const_iterator() function
      template <typename T, bool>
      struct has_const_iterator_struct : public std::false_type {};

      template <typename T>
      struct has_const_iterator_struct<T, true> : public std::true_type {};

      // Culled by SFINAE if T::const_iterator does not exist
      // or is not accessible or not default-constructable
      template <typename T>
      constexpr auto has_const_iterator_helper(T * /* = nullptr */)
        -> decltype(typename T::const_iterator(), bool())
      {
        return true;
      }

      // Used as fallback when SFINAE culls the template method
      constexpr bool has_const_iterator_helper(...) { return false; }
      //------------------------------------------------------------------------

    } // namespace type_traits
  }   // namespace details

  //---
  //--- Identity declaration
  //---
  template <class T>
  class Identity {
  public:
    typedef T argument_type;
    typedef T result_type;

    result_type& operator()(argument_type& x) const { return x; }
    const result_type& operator()(const argument_type& x) const { return x; }
  }; // class Identity<>

  //---
  //--- PairSecond declaration
  //---
  template <class T>
  class PairSecond {
  public:
    typedef T argument_type;
    typedef typename T::second_type result_type;

    result_type& operator()(argument_type& p) const { return p.second; }
    const result_type& operator()(const argument_type& p) const { return p.second; }
  }; // class PairSecond<>

  //---
  //--- deep_const_fwd_iterator_nested declaration
  //---
  template <typename ITER, typename INNERCONTEXTRACT>
  class deep_const_fwd_iterator_nested : public std::forward_iterator_tag {
  public:
    using OuterIterator_t = ITER;
    using InnerContainerExtractor_t = INNERCONTEXTRACT;
    using InnerContainer_t = typename InnerContainerExtractor_t::result_type;
    using InnerIterator_t
      //    = deep_fwd_const_iterator<typename ITER::value_type::const_iterator>;
      = typename InnerContainer_t::const_iterator;

    //  using iterator_type = deep_fwd_const_iterator<OuterIterator_t>;
    using iterator_type =
      deep_const_fwd_iterator_nested<OuterIterator_t, InnerContainerExtractor_t>;

    /// Type of the value pointed by the iterator
    using value_type = typename InnerIterator_t::value_type;

    struct BeginPositionTag {};
    struct EndPositionTag {};

    static constexpr BeginPositionTag begin = {};
    static constexpr EndPositionTag end = {};

    /**
     * @brief Default constructor: invalid iterator
     *
     * This constructor sets the iterator in an invalid, end-of-container state.
     */
    deep_const_fwd_iterator_nested() = default;

    /**
     * @brief Constructor: starts from the container at the specified iterator
     * @param src the starting point of the iterator
     * @param end the end point of the iterator
     *
     * This constructor does not set and end. Due to how the class works,
     * if the outer container has an "end", reaching it with this iterator has
     * a undefined behaviour. You most likely want to use the constructor where
     * you can also specify the end of the container.
     */
    deep_const_fwd_iterator_nested(OuterIterator_t src, OuterIterator_t end);

    /**
     * @brief Constructor: starts from the beginning of the specified container
     * @param cont the container
     * @param [anonymous] tag to select the begin-of-container behaviour
     *
     * The second parameter is used just to choose which constructor to use.
     * Two constants are provided, begin and end, defined in the iterator itself
     * (no explicit namespace is required).
     * Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * using Data_t = std::vector<std::vector<float>>;
     * Data_t data(5, { 1., 3., 5. });
     * deep_const_fwd_iterator_nested<Data_t::const_iterator> iData(data, begin),
     *   data_end(data, end);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <class CONT>
    deep_const_fwd_iterator_nested(const CONT& cont, BeginPositionTag)
      : deep_const_fwd_iterator_nested(std::begin(cont), std::end(cont))
    {
      skip_empty();
    }

    /**
     * @brief Constructor: starts from the end of the specified container
     * @param cont the container
     * @param [anonymous] tag to select the end-of-container behaviour
     *
     * The second parameter is used just to choose which constructor to use.
     * Two constants are provided, begin and end, defined in the iterator itself
     * (no explicit namespace is required).
     * Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * using Data_t = std::vector<std::vector<float>>;
     * Data_t data(5, { 1., 3., 5. });
     * deep_const_fwd_iterator_nested<Data_t::const_iterator> iData(data, begin),
     *   data_end(data, end);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <class CONT>
    deep_const_fwd_iterator_nested(const CONT& cont, EndPositionTag)
      : deep_const_fwd_iterator_nested(std::end(cont))
    {}

    /**
     * @brief Prefix increment operator: points to the next element
     * @return this iterator, already incremented
     *
     * The behaviour of this method on a past-the-end iterator is undefined
     * (currently, chances are it will access invalid memory).
     */
    iterator_type& operator++();

    /**
     * @brief Postfix increment operator: points to the next element
     * @return a copy of this iterator before the increment
     *
     * The behaviour of this method on a past-the-end iterator is undefined
     * (currently, chances are it will access invalid memory).
     */
    iterator_type operator++(int)
    {
      iterator_type old(*this);
      this->operator++();
      return old;
    }

    ///@{
    /// @name Dereference operators
    const value_type& operator*() const { return *inner_iter; }
    const value_type& operator->() const { return *inner_iter; }
    ///@}

    ///@{
    /// @name Comparison operators
    /// Returns true if the two iterators are equivalent
    bool operator==(const iterator_type& as) const
    {
      return (as.outer_iter == outer_iter) &&
             ((as.inner_iter == inner_iter) || (is_end() && as.is_end()));
    }
    /// Returns true if the two iterators are not equivalent
    bool operator!=(const iterator_type& as) const
    {
      return (as.outer_iter != outer_iter) ||
             ((as.inner_iter != inner_iter) && (!is_end() || !as.is_end()));
    }
    ///@}

    ///@{
    /**
     * @brief Bonus: convert to bool to find out if we are at the end
     * @return whether this operator is past-the-end
     */
    operator bool() const { return !is_end(); }
    bool operator!() const { return is_end(); }
    ///@}

    /// Swaps this with the specified iterator
    void swap(iterator_type& with);

  protected:
    OuterIterator_t outer_iter; ///< points to current inner container
    OuterIterator_t outer_end;  ///< points to past-the-end inner container

    InnerIterator_t inner_iter; ///< points to the current element
    InnerIterator_t inner_end;  ///< stores the end of current inner container

    /// Internal constructor: past-the-end iterator pointing to specified place
    deep_const_fwd_iterator_nested(OuterIterator_t end) : deep_const_fwd_iterator_nested(end, end)
    {}

  private:
    void init_inner();
    void reset_inner();
    void skip_empty(); ///< points to the next item

    bool is_end() const { return outer_iter == outer_end; }

    /// Extracts the value out of the inner iterator
    const typename InnerContainerExtractor_t::result_type& extract_container(
      const typename OuterIterator_t::value_type& v)
    {
      return InnerContainerExtractor_t()(v);
    }

  }; // class deep_const_fwd_iterator_nested<>

  //---
  //--- deep_const_fwd_iterator_nested implementation
  //---

  template <typename ITER, typename INNERCONTEXTRACT>
  deep_const_fwd_iterator_nested<ITER, INNERCONTEXTRACT>::deep_const_fwd_iterator_nested(
    OuterIterator_t src,
    OuterIterator_t end)
    : outer_iter(src), outer_end(end)
  {
    if (is_end()) return;
    init_inner();
    skip_empty();
  } // deep_const_fwd_iterator_nested(OuterIterator_t, OuterIterator_t)

  template <typename ITER, typename INNERCONTEXTRACT>
  deep_const_fwd_iterator_nested<ITER, INNERCONTEXTRACT>&
  deep_const_fwd_iterator_nested<ITER, INNERCONTEXTRACT>::operator++()
  {
    ++inner_iter;
    skip_empty();
    return *this;
  } // deep_const_fwd_iterator_nested<>::operator++()

  template <typename ITER, typename INNERCONTEXTRACT>
  void deep_const_fwd_iterator_nested<ITER, INNERCONTEXTRACT>::swap(iterator_type& with)
  {
    std::swap(outer_iter, with.outer_iter);
    std::swap(outer_end, with.outer_end);
    std::swap(inner_iter, with.inner_iter);
    std::swap(inner_end, with.inner_end);
  } // deep_const_fwd_iterator_nested<>::swap()

  template <typename ITER, typename INNERCONTEXTRACT>
  void deep_const_fwd_iterator_nested<ITER, INNERCONTEXTRACT>::init_inner()
  {
    inner_iter = std::begin(extract_container(*outer_iter));
    inner_end = std::end(extract_container(*outer_iter));
  } // deep_const_fwd_iterator_nested<>::init_inner()

  template <typename ITER, typename INNERCONTEXTRACT>
  void deep_const_fwd_iterator_nested<ITER, INNERCONTEXTRACT>::reset_inner()
  {
    inner_end = inner_iter = {};
  }

  template <typename ITER, typename INNERCONTEXTRACT>
  void deep_const_fwd_iterator_nested<ITER, INNERCONTEXTRACT>::skip_empty()
  {
    while (inner_iter == inner_end) {
      ++outer_iter;
      if (is_end()) {
        reset_inner();
        return;
      } // if
      init_inner();
    } // while inner iterator ended
  }   // skip_empty()

} // namespace lar

namespace std {
  template <typename CITER, typename INNERCONTEXTRACT>
  inline void swap(lar::deep_const_fwd_iterator_nested<CITER, INNERCONTEXTRACT>& a,
                   lar::deep_const_fwd_iterator_nested<CITER, INNERCONTEXTRACT>& b)
  {
    a.swap(b);
  }
} // namespace std

#endif // NESTEDITERATOR_H
