/**
 * @file    CountersMap_test.cc
 * @brief   Tests the counter map
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    20140822
 * @version 1.0
 *
 * See http://www.boost.org/libs/test for the Boost test library home page.
 *
 * Timing:
 * version 1.0 takes about 30" on a 3 GHz machine.
 */

// C/C++ standard libraries
#include <iostream>
#include <map>
#include <random>

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
#define BOOST_TEST_MODULE (CountersMap_test)
#include "boost/test/unit_test.hpp"

// LArSoft libraries
#include "lardata/Utilities/CountersMap.h"

/// The seed for the default random engine
constexpr unsigned int RandomSeed = 12345;

//------------------------------------------------------------------------------
//--- Test code
//

/**
 * @brief Tests with a vector of counter maps (Hough transfort use case)
 *
 * The test consists in filling a lot of points into a 2D sparse "image" (or
 * histogram).
 * Two structures are maintained, one with the standard STL map, another
 * with the CountersMap.
 * The test fails if the two images do not match.
 */
void RunHoughTransformTreeTest()
{

  // the structure we are testing is a 2D "image" of integers;
  // image is mostly empty (zero), but each abscissa has roughly the same
  // number of non-empty pixels (NPoints), and at least one of them.

  constexpr unsigned int NPoints = 1000;
  constexpr unsigned int NAngles = 10800;
  constexpr unsigned int NDist = 2500; // half distance

  typedef std::map<int, int> BaseMap_t;

  // STL container
  typedef std::vector<BaseMap_t> MapVectorI_t;
  MapVectorI_t stl_image(NAngles);

  // CountersMap; uses chunks of 8 counters per block
  std::vector<lar::CountersMap<int, int, 8>> cm_image(NAngles);
  // the following should fail compilation
  //  std::vector<lar::CountersMap<int, int, 9>> cm_image_broken(NAngles);

  static std::default_random_engine random_engine(RandomSeed);
  std::uniform_real_distribution<float> uniform(-1., 1.);

  for (unsigned int iPoint = 0; iPoint != NPoints; ++iPoint) {
    // we add here some simple image, not to strain the test;
    // this is a straight line
    const float offset = uniform(random_engine) * NDist;
    const float slope = uniform(random_engine);
    float d = offset;
    for (size_t iAngle = 0; iAngle < NAngles; ++iAngle) {
      // add one entry on the (angle ; distance) plane
      ++(stl_image[iAngle][int(d)]);
      cm_image[iAngle].increment(int(d)); // different interface than usual map
      // prepare for the next point; wrap in the [-NDist, NDist[ range
      d += slope;
      while (d >= (float)NDist)
        d -= 2 * NDist;
      while (d < 0)
        d += 2 * NDist;
    } // for iAngle
  }   // for iPoint

  std::cout << "Filling complete, now checking." << std::endl;

  // we have to provide a comparison between two "different" structures
  // (having different allocators is enough to make them unrelated)
  unsigned int nExtraKeys = 0, nMismatchValue = 0, nMissingKeys = 0;
  auto stl_begin = stl_image.cbegin();
  unsigned int iMap = 0;
  for (const auto& cm_map : cm_image) {

    const MapVectorI_t::value_type& stl_map = *(stl_begin++);

    std::cout << "Map #" << iMap << " (" << cm_map.n_counters() << " counters, " << stl_map.size()
              << " real)" << std::endl;

    // compare the two maps; the CountersMap one has more elements,
    // since the counters are allocated in blocks;
    // if a key is in STL map, it must be also in the CountersMap;
    // if a key is not in STL map, counter in CountersMap must be missing or 0
    MapVectorI_t::value_type::const_iterator stl_iter = stl_map.begin(), stl_end = stl_map.end();
    for (auto p : cm_map) { // this should be a pair (index, counter)

      if (stl_iter != stl_end) { // we have still counters to find
        // if counter is already beyond the next non-empty one froml STL map,
        // then we are missing some
        while (p.first > stl_iter->first) {
          ++nMissingKeys;
          std::cout << "ERROR missing key " << stl_iter->first << std::endl;
          if (++stl_iter == stl_end) break;
        }
      } // if

      if (stl_iter != stl_end) { // we have still counters to find
        if (p.first == stl_iter->first) {
          // if the counter is in SLT map, the two counts must match
          //  std::cout << "  " << p.first << " " << p.second << std::endl;
          if (stl_iter->second != p.second) {
            std::cout << "ERROR wrong counter value " << p.second << ", expected "
                      << stl_iter->second << std::endl;
            ++nMismatchValue;
          }
          ++stl_iter; // done with it
        }
        else if (p.first < stl_iter->first) {
          // if the counter is not in STL map, then it must be 0
          if (p.second != 0) {
            ++nExtraKeys;
            std::cout << "ERROR extra key " << p.first << " (" << p.second << ")" << std::endl;
          }
          //  else {
          //    std::cout << "  " << p.first << " " << p.second << " (not in STL)"
          //      << std::endl;
          //  }
        }
      }
      else {
        // no more keys in STL map
        if (p.second != 0) {
          ++nExtraKeys;
          std::cout << "ERROR extra key " << p.first << " (" << p.second << ")" << std::endl;
        }
      }
    } // for element in map

    BOOST_TEST(cm_map.is_equal(stl_map));

    // if they were the same, make sure that now they differ
    const_cast<MapVectorI_t::value_type&>(stl_map)[NDist / 2]++;
    BOOST_TEST(!cm_map.is_equal(stl_map));

    ++iMap;
  } // for map

  BOOST_TEST(nMismatchValue == 0U);
  BOOST_TEST(nMissingKeys == 0U);
  BOOST_TEST(nExtraKeys == 0U);

} // RunHoughTransformTreeTest()

//------------------------------------------------------------------------------
//--- registration of tests
//
// Boost needs now to know which tests we want to run.
// Tests are "automatically" registered, hence the BOOST_AUTO_TEST_CASE()
// macro name. The argument is a name for the test; each test will have a
// number of checks and it will fail if any of them does.
//

BOOST_AUTO_TEST_CASE(RunHoughTransformTree)
{
  RunHoughTransformTreeTest();
  std::cout << "Done." << std::endl;
}
