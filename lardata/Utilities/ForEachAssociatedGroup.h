/**
 * @file lardata/Utilities/ForEachAssociatedGroup.h
 * @brief  Helper functions to access associations in order.
 *
 * No additional linking is required to use these functions.
 *
 * Provided functions:
 *
 * * `util::associated_groups()` providing a sequence of objects associated to
 *   the same object, for each object
 *
 */


#ifndef LARDATA_UTILITIES_FOR_EACH_ASSOCIATED_GROUP_H
#define LARDATA_UTILITIES_FOR_EACH_ASSOCIATED_GROUP_H

// LArSoft libraries
#include "lardata/Utilities/RangeForWrapper.h"

// framework libraries
#include "canvas/Persistency/Common/AssnsAlgorithms.h" // art::for_each_group()

// range library
#include "range/v3/algorithm/for_each.hpp"
#include "range/v3/view/group_by.hpp"
#include "range/v3/view/transform.hpp"
#include "range/v3/view/map.hpp" // range::views::values
#include "range/v3/view/all.hpp"

// C/C++ standard libraries
#include <utility> // std::make_pair()
#include <iterator> // std::next()


namespace util {
  /**
   * @brief  Helper functions to access associations in order.
   * @tparam A type of association being read
   * @tparam F type of functor to be called on each associated group
   * @param assns the association being read
   * @param func functor to be called on each associated group
   * @see associated_groups() art::for_each_group()
   *
   * @deprecated Moved into _canvas_: `art::for_each_group()`.
   */
  template <class A, class F>
  [[deprecated("Use art::for_each_group() instead")]]
  void for_each_associated_group(A const & assns, F & func)
     { art::for_each_group(assns, func); }


  /**
   * @brief  Helper functions to access associations in order.
   * @tparam A type of association being read
   * @param assns the association being read
   * @see for_each_associated_group()
   *
   * This function provides a functionality equivalent to
   * `art::for_each_group()`, but it grants the caller additional control on the
   * external loop and on the function.
   *
   * Example: assuming that a module with input tag stored in `fTrackTag` has
   * created associations of each track to its hits, the total charge for each
   * track can be extracted by:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto assns = art::getValidHandle<art::Assns<recob::Track, recob::Hit>>
   *   (fTrackTag);
   *
   * std::vector<double> totalCharge;
   *
   * for (auto const& hits: util::associated_groups(*assns)) {
   *   double total = 0.;
   *   for (art::Ptr<recob::Hit> const& hit: hits)
   *     total += hit->Integral();
   *   totalCharge.push_back(total);
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * A number of important points need to be realised about this example:
   *
   *  * the requirements of this function on its input association are the same
   *    as for `art::for_each_group()`
   *  * we can code the action on each group of hits directly in a loop, if
   *    like in this case the code is succinct
   *  * again, there is one outer loop iteration for every track;
   *  * the value of `hits` is an object representing a range of _art_ pointers
   *    (`art::Ptr<recob::Hit>`) which can be navigated with the
   *    `begin()`/`end()` free functions, or in a range-for loop;
   *  * on each iteration, the information of which track the hits are
   *    associated to is not available; if that is also needed, use
   *    `util::associated_groups_with_left()` instead.
   */
  template <class A>
  auto associated_groups(A const & assns) {
     return assns |
            ranges::views::all |
            ranges::views::group_by([](auto a1, auto a2) { return a1.first == a2.first;}) |
            ranges::views::transform([] (auto pairs) {return pairs | ranges::views::values | util::range_for;}) |
            util::range_for
            ;
  } // associated_groups()


  /**
   * @brief  Helper functions to access associations in order, also with key.
   * @tparam A type of association being read
   * @param assns the association being read
   * @see for_each_associated_group()
   *
   * This function provides a functionality equivalent to
   * `art::for_each_group_with_left()`, but it grants the caller additional
   * control on the external loop and on the function.
   *
   * Example: assuming that a module with input tag stored in `fTrackTag` has
   * created associations of each track to its hits, the total charge for each
   * track can be extracted by:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto assns = art::getValidHandle<art::Assns<recob::Track, recob::Hit>>
   *   (fTrackTag);
   *
   * std::map<int, double> totalCharge;
   *
   * for (auto const& trackWithHits: util::associated_groups_with_left(*assns))
   * {
   *   art::Ptr<recob::Track> const& track = trackWithHits.first;
   *   auto const& hits = trackWithHits.second;
   *
   *   if (totalCharge.count(track->ID()) > 0) {
   *     throw art::Exception(art::errors::LogicError)
   *       << "Multiple tracks have ID " << track->ID() << "!\n";
   *   }
   *
   *   double& total = totalCharge[track->ID()];
   *   total = 0.0;
   *   for (art::Ptr<recob::Hit> const& hit: hits)
   *     total += hit->Integral();
   *
   * } // for
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * A number of important points need to be realised about this example:
   *
   *  * the requirements of this function on its input association are the same
   *    as for `art::for_each_group_with_left()`
   *  * we can code the action on each group of hits directly in a loop, if
   *    like in this case the code is succinct
   *  * again, there is one outer loop iteration for every track;
   *  * the value of `hits` is an object representing a range of _art_ pointers
   *    (`art::Ptr<recob::Hit>`) which can be navigated with the
   *    `begin()`/`end()` free functions, or in a range-for loop.
   */
  template <class A>
  auto associated_groups_with_left(A const & assns) {
     return assns
        | ranges::views::all
        | ranges::views::group_by([](auto a1, auto a2) { return a1.first == a2.first;})
        | ranges::views::transform([] (auto pairs)
           {
             return std::make_pair(
               pairs.front().first, // assuming they're all the same, pick first
               pairs | ranges::views::values | util::range_for
               );
           })
        | util::range_for
        ;
  } // associated_groups_with_left()


  /**
   * @brief  Returns the group within `groups` with the specified index.
   * @tparam Groups the type of collection of groups
   * @param groups the collection of all groups
   * @param index the index of the group to be accessed
   * @return the group with specified index (may be a reference)
   * @see `associated_groups()`
   *
   * The `groups` argument is expected to be the one returned by
   * `associated_groups`.
   */
  template <typename Groups>
  auto groupByIndex(Groups&& groups, std::size_t index) -> decltype(auto)
    { return *(std::next(groups.begin(), index)); }

} // namespace util

#endif // LARDATA_UTILITIES_FOR_EACH_ASSOCIATED_GROUP_H
