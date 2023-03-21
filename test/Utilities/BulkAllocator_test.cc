/**
 * @file    BulkAllocator_test.cc
 * @brief   Tests the bulk allocator
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    20140817
 * @version 1.0
 *
 * See http://www.boost.org/libs/test for the Boost test library home page.
 *
 * Timing:
 * version 1.0 takes less than 3" on a 3 GHz machine FIXME
 */

// C/C++ standard libraries
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
#define BOOST_TEST_MODULE (BulkAllocator_test)
#include "boost/test/unit_test.hpp"

// LArSoft libraries
#include "lardata/Utilities/BulkAllocator.h"

/// The seed for the default random engine
constexpr unsigned int RandomSeed = 12345;

//------------------------------------------------------------------------------
//--- Test code
//

/**
 * @brief Tests bulk allocator with a vector of maps (Hough transfort use case)
 *
 * The test consists in filling a lot of points into a 2D sparse "image" (or
 * histogram).
 * Two structures are maintained, one with the standard STL allocator, another
 * with the BulkAllocator.
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

  // STL allocator
  std::vector<std::map<int, int, BaseMap_t::key_compare, std::allocator<BaseMap_t::value_type>>>
    stl_image(NAngles);

  // BulkAllocator
  // Normally, I would declare an allocator like
  // lar::BulkAllocator<BaseMap_t::value_type> MyAllocator(500000);
  // and use it to initialize the STL container; then the STL would use the same
  // allocator or a copy of it or a default-constructed version of the same
  // allocator class, and that does not matter because all of them refer
  // internally to the same base allocator for BaseMap_t::value_type.
  // But map rebinds the allocator to one based not on BaseMap_t::value_type,
  // but on the tree node used in the map implementation; that allocator is
  // of a different type, and does not refer to
  // lar::BulkAllocator<BaseMap_t::value_type> any more.
  // Therefore I let map do its tricks...
  std::vector<std::map<int, int, BaseMap_t::key_compare, lar::BulkAllocator<BaseMap_t::value_type>>>
    bulk_image(NAngles);
  // ... and then I change the underlying allocator directly;
  // Unfortunately this trick is implementation-dependent, in that there is no
  // requirement for map to expose the allocator it really uses.
  // Here the non-standard part is the name of the node type, std::_Rb_tree_node
  lar::BulkAllocator<std::_Rb_tree_node<BaseMap_t::value_type>>::SetChunkSize(500000);

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
      ++(bulk_image[iAngle][int(d)]);
      // prepare for the next point; wrap in the [-NDist, NDist[ range
      d += slope;
      while (d >= (float)NDist)
        d -= 2 * NDist;
      while (d < 0)
        d += 2 * NDist;
    } // for iAngle
  }   // for iPoint

  // we have to provide a comparison between two "different" structures
  // (having different allocators is enough to make them unrelated)
  bool bSame = true;
  auto stl_begin = stl_image.cbegin();
  for (const auto& bulk_map : bulk_image) {
    // the std::equal() call compares pairs (int, int) of each map
    if (std::equal(bulk_map.begin(), bulk_map.end(), (stl_begin++)->begin())) continue;
    bSame = false;
    break;
  } // for

  BOOST_TEST(bSame);

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
