/**
 * @file   lardata/Utilities/filterRangeFor.h
 * @brief  Utilities to manipulate range for loops.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   May 1, 2018
 *
 * This is a header-only library.
 */

#ifndef LARDATA_UTILITIES_FILTERRANGEFOR_H
#define LARDATA_UTILITIES_FILTERRANGEFOR_H

// Boost libraries
#include <boost/iterator/filter_iterator.hpp>

namespace util {

  /**
   * @brief Provides iteration only through elements passing a condition.
   * @tparam Range the data to be iterated
   * @tparam Pred the type of the predicate to be fulfilled
   * @param range the data to be iterated through
   * @param pred the predicate to be tested
   * @return an object suitable to be used in a range-for loop
   *
   * This adapter makes the range for loop iterate only through the elements of
   * `range` which fulfil the predicate `pred`.
   *
   * This example will print: "0 3 6 9 ":
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *
   * std::vector<int> data = { 0, 1, 2, 3, 4, 5, 6 ,7, 8, 9 };
   * for (int v: util::filterRangeFor(data, [](int v){ return v % 3 == 0; })) {
   *
   *   std::cout << v << " ";
   *
   * } // for
   * std::cout << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * Note that `pred` may be copied (`range` will not be).
   *
   *
   * Requirements
   * -------------
   *
   * * `Range` is an object which can itself go through a range-for:
   *
   *     for (auto&& v: range);
   *
   *   is valid
   * * `Pred` is a copiable unary function type, whose single argument can be
   *     converted from the value type of `Range`, and whose return value can be
   *     converted into a `bool` vaule
   *
   */
  template <typename Range, typename Pred>
  auto filterRangeFor(Range&& range, Pred&& pred) -> decltype(auto);

} // namespace util

//------------------------------------------------------------------------------
//--- implementation
//------------------------------------------------------------------------------
namespace util {

  namespace details {

    //--------------------------------------------------------------------------
    template <typename Range, typename Pred>
    class FilterRangeForStruct {

      /// Extract the begin iterator from a range.
      static auto getBegin(Range&& range) -> decltype(auto)
      {
        using std::begin;
        return begin(range);
      }

      /// Extract the end iterator from a range.
      static auto getEnd(Range&& range) -> decltype(auto)
      {
        using std::end;
        return end(range);
      }

      /// Create a Boost filter iterator pointing to the beginning of data.
      static auto makeBeginIterator(Range&& range, Pred&& pred)
      {
        auto begin = getBegin(std::forward<Range>(range));
        auto end = getEnd(std::forward<Range>(range));
        return boost::make_filter_iterator(std::forward<Pred>(pred), begin, end);
      }

      /// Create a Boost filter iterator pointing to the end of data.
      static auto makeEndIterator(Range&& range, Pred&& pred)
      {
        auto end = getEnd(std::forward<Range>(range));
        return boost::make_filter_iterator(std::forward<Pred>(pred), end, end);
      }

      using begin_iterator_t =
        decltype(makeBeginIterator(std::declval<Range&&>(), std::declval<Pred&&>()));
      using end_iterator_t =
        decltype(makeEndIterator(std::declval<Range&&>(), std::declval<Pred&&>()));

      begin_iterator_t fBegin;
      end_iterator_t fEnd;

    public:
      /// Extracts the iterators from the specified range.
      FilterRangeForStruct(Range&& range, Pred&& pred)
        : fBegin(makeBeginIterator(std::forward<Range>(range), std::forward<Pred>(pred)))
        , fEnd(makeEndIterator(std::forward<Range>(range), std::forward<Pred>(pred)))
      {}

      auto begin() const { return fBegin; }
      auto end() const { return fEnd; }

    }; // class FilterRangeForStruct<>

    //--------------------------------------------------------------------------

  } // namespace details

  //----------------------------------------------------------------------------
  template <typename Range, typename Pred>
  auto filterRangeFor(Range&& range, Pred&& pred) -> decltype(auto)
  {
    return details::FilterRangeForStruct<Range, Pred>(std::forward<Range>(range),
                                                      std::forward<Pred>(pred));
  }

  //----------------------------------------------------------------------------

} // namespace util

//------------------------------------------------------------------------------

#endif // LARDATA_UTILITIES_FILTERRANGEFOR_H
