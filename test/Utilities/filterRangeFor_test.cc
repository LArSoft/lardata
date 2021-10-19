/**
 * @file   filterRangeFor_test.cc
 * @brief  Test for `util::filterRangeFor()`.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   May 1, 2018
 * @see    `lardata/Utilities/filterRangeFor.h`
 *
 * The test is run with no arguments.
 *
 */

// LArSoft libraries
#include "lardata/Utilities/filterRangeFor.h"

// Boost libraries
#define BOOST_TEST_MODULE ( filterRangeFor_test )
#include "boost/test/unit_test.hpp"

// C/C++ standard libraries
#include <algorithm> // std::count_if()
#include <numeric> // std::iota()


//-----------------------------------------------------------------------------
template <typename Cont, typename Pred>
void testPredicate(Cont& data, Pred pred) {

  auto const nPass = std::count_if(data.begin(), data.end(), pred);

  unsigned int n = 0;
  for (auto const& v: util::filterRangeFor(data, pred)) {
    ++n;
    BOOST_TEST_CHECKPOINT("  testing value: " << v);
    BOOST_TEST(pred(v));
  } // for

  BOOST_TEST(n == nPass);
} // testPredicate()


//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(filterRangeFor_testCase) {

  std::vector<int> data(20);
  std::iota(data.begin(), data.end(), 0);

  BOOST_TEST_MESSAGE("Selecting multiples of 3");
  testPredicate<std::vector<int> const>
    (data, [](int v){ return (v % 3) == 0; });
  testPredicate(data, [](int v){ return (v % 3) == 0; });

  BOOST_TEST_MESSAGE("Selecting values that are not 9");
  testPredicate<std::vector<int> const>(data, [](int v){ return v != 9; });
  testPredicate(data, [](int v){ return v != 9; });

  BOOST_TEST_MESSAGE("Selecting values that are 50");
  testPredicate<std::vector<int> const>(data, [](int v){ return v == 50; });
  testPredicate(data, [](int v){ return v == 50; });

} // BOOST_AUTO_TEST_CASE(filterRangeFor_testCase)
