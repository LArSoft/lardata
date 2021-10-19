/**
 * @file    NestedIterator_test.cc
 * @brief   Tests nested iterators
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    20140822
 * @version 1.0
 *
 * See http://www.boost.org/libs/test for the Boost test library home page.
 *
 * Timing:
 * version 1.0 takes negligible time on a 3 GHz machine
 */

// C/C++ standard libraries
#include <map>
#include <random>
#include <iostream>

// Boost libraries
/*
 * Boost Magic: define the name of the module;
 * and do that before the inclusion of Boost unit test headers
 * because it will change what they provide.
 * Among the those, there is a main() function and some wrapping catching
 * unhandled exceptions and considering them test failures, and probably more.
 * This also makes fairly complicate to receive parameters from the command line
 * (for example, a random seed).
 */
#define BOOST_TEST_MODULE ( NestedIterator_test )
#include "boost/test/unit_test.hpp"

// LArSoft libraries
#include "lardata/Utilities/NestedIterator.h"


/// The seed for the default random engine
constexpr unsigned int RandomSeed = 12345;


//------------------------------------------------------------------------------
//--- Test code
//

/**
 * @brief Tests bulk allocator with a vector of vectors
 *
 * The test consists in filling a sequence of integers in a two-level structure,
 * and then iterating to recover the sequence.
 *
 * The test fails if the extracted sequence is not correct.
 */
void RunVectorVectorTest() {

  // fill the double tier structure
  using DoubleVectorI_t = std::vector<std::vector<int>>;
  DoubleVectorI_t data(1); // get the first vector started
  constexpr size_t NElements = 10000;
  constexpr float SwitchProbability = 0.1; // expect about 1000 containers

  static std::default_random_engine random_engine(RandomSeed);
  std::uniform_real_distribution<float> uniform(0., 1.);

  unsigned int nEmpty = 0;
  for (size_t i = 0; i < NElements; ++i) {
    // add a new vector (some times)
    if (uniform(random_engine) < SwitchProbability) {
      if (data.back().empty()) ++nEmpty;
      data.emplace_back();
    }
    // add the element i to the last vector
    data.back().push_back(i);
  } // for
  std::cout << "Working with " << NElements << " elements in " << data.size()
    << " vectors (" << nEmpty << " empty) in a vector" << std::endl;

  unsigned int nMismatches = 0;

  int expected = 0;
  lar::double_fwd_const_iterator<DoubleVectorI_t::const_iterator>
    iElem(data, lar::double_fwd_const_iterator<DoubleVectorI_t::const_iterator>::begin),
    eend(data, lar::double_fwd_const_iterator<DoubleVectorI_t::const_iterator>::end);
  while (iElem != eend) {
    int elem = *iElem;
    if (elem != expected) ++nMismatches;
    ++expected;
    ++iElem;
  } // while

  BOOST_TEST((unsigned int) expected == NElements);
  BOOST_TEST(nMismatches == 0U);
} // RunVectorVectorTest()


/**
 * @brief Tests bulk allocator with a map of vectors
 *
 * The test consists in filling a sequence of integers in a two-level structure,
 * and then iterating to recover the sequence.
 *
 * The test fails if the extracted sequence is not correct.
 */
void RunVectorMapTest() {

  // fill the double tier structure
  using VectorMapI_t = std::map<int, std::vector<int>>;
  VectorMapI_t data; // get the first vector started
  //data.emplace(0, std::vector<int>{});
  data[0] = {};
  constexpr size_t NElements = 10000;
  constexpr float SwitchProbability = 0.1; // expect about 1000 containers

  static std::default_random_engine random_engine(RandomSeed);
  std::uniform_real_distribution<float> uniform(0., 1.);

  unsigned int nEmpty = 0;
  for (size_t i = 0; i < NElements; ++i) {
    // add a new map (some times)
    if (uniform(random_engine) < SwitchProbability) {
      if (data.rbegin()->second.empty()) ++nEmpty;
      data.insert({ data.size(), {} });
    }
    // add the element i to the last vector
    data.rbegin()->second.push_back(i);
  } // for
  std::cout << "Working with " << NElements << " elements in " << data.size()
    << " vectors (" << nEmpty << " empty) in a map" << std::endl;

  unsigned int nMismatches = 0;

  int expected = 0;
  using ConstIterator_t = lar::double_fwd_const_iterator
    <VectorMapI_t::const_iterator, lar::PairSecond<VectorMapI_t::value_type>>;
  ConstIterator_t iElem(data, ConstIterator_t::begin),
    eend(data, ConstIterator_t::end);
  while (iElem != eend) {
    int elem = *iElem;
    if (elem != expected) ++nMismatches;
    ++expected;
    ++iElem;
  } // while

  BOOST_TEST((unsigned int) expected == NElements);
  BOOST_TEST(nMismatches == 0U);
} // RunVectorMapTest()


//------------------------------------------------------------------------------
//--- registration of tests
//
// Boost needs now to know which tests we want to run.
// Tests are "automatically" registered, hence the BOOST_AUTO_TEST_CASE()
// macro name. The argument is a name for the test; each test will have a
// number of checks and it will fail if any of them does.
//

BOOST_AUTO_TEST_CASE(RunVectorVector) {
  RunVectorVectorTest();
}

BOOST_AUTO_TEST_CASE(RunVectorMap) {
  RunVectorMapTest();
}
