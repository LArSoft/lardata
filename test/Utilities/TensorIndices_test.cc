/**
 * @file   TensorIndices_test.cc
 * @brief  Test for TensorIndices class
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 28, 2016
 * @see    TensorIndices.h
 *
 * This test instantiates some TensorIndices and verifies the reaction to some
 * hard-coded settings and queries.
 *
 * The test is run with no arguments.
 *
 * Tests are run:
 *
 * * `VectorTest`: one-dimension tensor test
 * * `MatrixTest`: two-dimension tensor test
 * * `TensorRank3Test`: test rank 3 tensor
 *
 * See the documentation of the three functions for more information.
 *
 */

// LArSoft libraries
#include "lardata/Utilities/TensorIndices.h"

// Boost libraries
#define BOOST_TEST_MODULE (PointIsolationAlg_test)
#include "boost/test/unit_test.hpp"

// C/C++ standard libraries
#include <array>
#include <stdexcept> // std::out_of_range

//------------------------------------------------------------------------------
//--- Test code
//---
/**
 * @brief Test for a rank 1 tensor (vector)
 */
void VectorTest()
{

  // indices for a vector of size 4
  util::TensorIndices<1> const indices{4};

  //
  // reflection
  //
  BOOST_TEST(indices.dim<0>() == 4U);

  BOOST_TEST(indices.size<0>() == 4U);

  BOOST_TEST(indices.size() == 4U);

  //
  // indexing
  //
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<char, 1> ii{{1}}; // char is for test; it should be size_t
                               /* TODO
  BOOST_TEST(indices[0] == 0U);
  BOOST_TEST(indices[1] == 1U);
  BOOST_TEST(indices[2] == 2U);
  BOOST_TEST(indices[3] == 3U);
  BOOST_TEST(indices[4] == 4U); // no bound check
*/
  BOOST_TEST(indices(0) == 0U);
  BOOST_TEST(indices(1) == 1U);
  BOOST_TEST(indices(ii.begin()) == 1U);
  BOOST_TEST(indices(2) == 2U);
  BOOST_TEST(indices(3) == 3U);
  BOOST_CHECK_NO_THROW(indices(4)); // no bound check

  BOOST_TEST(indices.at(0) == 0U);
  BOOST_TEST(indices.at(1) == 1U);
  BOOST_TEST(indices.at(ii.begin()) == 1U);
  BOOST_TEST(indices.at(2) == 2U);
  BOOST_TEST(indices.at(3) == 3U);
  BOOST_CHECK_THROW(indices.at(4), std::out_of_range);

  BOOST_TEST(indices.has(0));
  BOOST_TEST(indices.has(1));
  BOOST_TEST(indices.has(ii.begin()));
  BOOST_TEST(indices.has(2));
  BOOST_TEST(indices.has(3));
  BOOST_TEST(!indices.has(4));

  BOOST_TEST(indices.hasIndex<0>(0));
  BOOST_TEST(indices.hasIndex<0>(1));
  BOOST_TEST(indices.hasIndex<0>(2));
  BOOST_TEST(indices.hasIndex<0>(3));
  BOOST_TEST(!indices.hasIndex<0>(4));
  BOOST_TEST(!indices.hasIndex<0>(5));

  BOOST_TEST(indices.hasLinIndex(0U));
  BOOST_TEST(indices.hasLinIndex(indices.size() - 1));
  BOOST_TEST(!indices.hasLinIndex(indices.size()));

  BOOST_TEST(indices.size() == 4U);

  // check that the function also works
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<int, 1> ia{{4}};
  util::TensorIndices<1> indicesAgain(ia.begin());
  util::TensorIndices<1> indicesOther = util::makeTensorIndices(3);
  util::TensorIndices<2> indicesRank = util::makeTensorIndices(4, 3);

  // check the comparison operators
  BOOST_TEST((indicesAgain == indices));
  BOOST_TEST(!(indicesAgain != indices));
  BOOST_TEST((indicesOther != indices));
  BOOST_TEST(!(indicesOther == indices));
  BOOST_TEST((indicesRank != indices));
  BOOST_TEST(!(indicesRank == indices));

} // VectorTest()

/**
 * @brief Test for a rank 2 tensor (matrix)
 */
void MatrixTest()
{

  // indices for a matrix 4 x 3
  util::TensorIndices<2> const indices{4, 3};

  //
  // reflection
  //
  BOOST_TEST(indices.dim<0>() == 4U);
  BOOST_TEST(indices.dim<1>() == 3U);

  BOOST_TEST(indices.size<0>() == 4U * 3U);
  BOOST_TEST(indices.size<1>() == 3U);

  BOOST_TEST(indices.size() == 4U * 3U);

  BOOST_TEST(indices.minorTensor().rank() == 1U);
  BOOST_TEST(indices.minorTensor().size() == 3U);

  //
  // indexing
  //
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<char, 2> ii{{1, 2}}; // char is for test; it should be size_t
                                  /* TODO
  BOOST_TEST(indices[0] == 0U);
  BOOST_TEST(indices[1] == 1U);
  BOOST_TEST(indices[2] == 2U);
  BOOST_TEST(indices[3] == 3U);
  BOOST_TEST(indices[4] == 4U); // no bound check
*/
  BOOST_TEST(indices(0, 0) == 0U);
  BOOST_TEST(indices(1, 2) == 5U);
  BOOST_TEST(indices(ii.begin()) == 5U);
  BOOST_CHECK_NO_THROW(indices(1, 3));
  BOOST_CHECK_NO_THROW(indices(4, 2));
  BOOST_CHECK_NO_THROW(indices(7, 6));

  BOOST_TEST(indices.at(0, 0) == 0U);
  BOOST_TEST(indices.at(1, 2) == 5U);
  BOOST_TEST(indices.at(ii.begin()) == 5U);
  BOOST_CHECK_THROW(indices.at(1, 3), std::out_of_range);
  BOOST_CHECK_THROW(indices.at(4, 2), std::out_of_range);
  BOOST_CHECK_THROW(indices.at(7, 6), std::out_of_range);

  BOOST_TEST(indices.has(0, 0));
  BOOST_TEST(indices.has(ii.begin()));
  BOOST_TEST(indices.has(2, 2));
  BOOST_TEST(!indices.has(1, 3));
  BOOST_TEST(!indices.has(4, 2));
  BOOST_TEST(!indices.has(7, 6));

  BOOST_TEST(indices.hasIndex<0>(0));
  BOOST_TEST(indices.hasIndex<0>(3));
  BOOST_TEST(!indices.hasIndex<0>(4));
  BOOST_TEST(indices.hasIndex<1>(0));
  BOOST_TEST(indices.hasIndex<1>(2));
  BOOST_TEST(!indices.hasIndex<1>(3));

  BOOST_TEST(indices.hasLinIndex(0U));
  BOOST_TEST(indices.hasLinIndex(indices.size() - 1));
  BOOST_TEST(!indices.hasLinIndex(indices.size()));

  //
  // comparisons
  //
  // check that the function and other constructors also works
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<int, 2> ia{{4, 3}};
  util::TensorIndices<2> indicesAgain(ia.begin());
  util::TensorIndices<2> indicesOther = util::makeTensorIndices(4, 4);
  util::TensorIndices<3> indicesRank = util::makeTensorIndices(5, 4, 3);

  // check the comparison operators
  BOOST_TEST((indicesAgain == indices));
  BOOST_TEST(!(indicesAgain != indices));
  BOOST_TEST((indicesOther != indices));
  BOOST_TEST(!(indicesOther == indices));
  BOOST_TEST((indicesRank != indices));
  BOOST_TEST(!(indicesRank == indices));

} // MatrixTest()

/**
 * @brief Test for a rank 3 tensor
 */
void TensorRank3Test()
{

  // indices for a matrix 2 x 3 x 4
  util::TensorIndices<3> const indices{2, 3, 4};

  //
  // reflection
  //
  BOOST_TEST(indices.dim<0>() == 2U);
  BOOST_TEST(indices.dim<1>() == 3U);
  BOOST_TEST(indices.dim<2>() == 4U);

  BOOST_TEST(indices.size<0>() == 2U * 3U * 4U);
  BOOST_TEST(indices.size<1>() == 3U * 4U);
  BOOST_TEST(indices.size<2>() == 4U);

  BOOST_TEST(indices.size() == 2U * 3U * 4U);

  BOOST_TEST(indices.minorTensor().rank() == 2U);
  BOOST_TEST(indices.minorTensor().size() == 3U * 4U);

  //
  // indexing
  //
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<char, 3> ii{{1, 2, 3}}; // char is for test; it should be size_t
                                     /* TODO
  BOOST_TEST(indices[0] == 0U);
  BOOST_TEST(indices[1] == 1U);
  BOOST_TEST(indices[2] == 2U);
  BOOST_TEST(indices[3] == 3U);
  BOOST_TEST(indices[4] == 4U); // no bound check
*/
  BOOST_TEST(indices(0, 0, 0) == 0U);
  BOOST_TEST(indices(1, 2, 3) == 23U);
  BOOST_TEST(indices(ii.begin()) == 23U);
  BOOST_CHECK_NO_THROW(indices(1, 3, 1));
  BOOST_CHECK_NO_THROW(indices(4, 2, 1));
  BOOST_CHECK_NO_THROW(indices(1, 2, 6));
  BOOST_CHECK_NO_THROW(indices(7, 6, 6));

  BOOST_TEST(indices.at(0, 0, 0) == 0U);
  BOOST_TEST(indices.at(1, 2, 3) == 23U);
  BOOST_TEST(indices.at(ii.begin()) == 23U);
  BOOST_CHECK_THROW(indices.at(1, 3, 1), std::out_of_range);
  BOOST_CHECK_THROW(indices.at(4, 2, 1), std::out_of_range);
  BOOST_CHECK_THROW(indices.at(1, 2, 6), std::out_of_range);
  BOOST_CHECK_THROW(indices.at(7, 6, 6), std::out_of_range);

  BOOST_TEST(indices.has(0, 0, 0));
  BOOST_TEST(indices.has(1, 2, 3));
  BOOST_TEST(indices.has(ii.begin()));
  BOOST_TEST(!indices.has(1, 3, 1));
  BOOST_TEST(!indices.has(4, 2, 1));
  BOOST_TEST(!indices.has(1, 2, 6));
  BOOST_TEST(!indices.has(7, 6, 6));

  BOOST_TEST(indices.hasIndex<0>(0));
  BOOST_TEST(indices.hasIndex<0>(1));
  BOOST_TEST(!indices.hasIndex<0>(2));
  BOOST_TEST(indices.hasIndex<1>(0));
  BOOST_TEST(indices.hasIndex<1>(2));
  BOOST_TEST(!indices.hasIndex<1>(3));
  BOOST_TEST(indices.hasIndex<2>(0));
  BOOST_TEST(indices.hasIndex<2>(3));
  BOOST_TEST(!indices.hasIndex<2>(4));

  BOOST_TEST(indices.hasLinIndex(0U));
  BOOST_TEST(indices.hasLinIndex(indices.size() - 1));
  BOOST_TEST(!indices.hasLinIndex(indices.size()));

  //
  // comparisons
  //
  // check that the function also works
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<int, 3> ia{{2, 3, 4}};
  util::TensorIndices<3> indicesAgain(ia.begin());
  util::TensorIndices<3> indicesOther = util::makeTensorIndices(2, 3, 5);
  util::TensorIndices<2> indicesRank = util::makeTensorIndices(2, 3);

  // check the comparison operators
  BOOST_TEST((indicesAgain == indices));
  BOOST_TEST(!(indicesAgain != indices));
  BOOST_TEST((indicesOther != indices));
  BOOST_TEST(!(indicesOther == indices));
  BOOST_TEST((indicesRank != indices));
  BOOST_TEST(!(indicesRank == indices));

} // TensorRank3Test()

//------------------------------------------------------------------------------
//--- tests
//
BOOST_AUTO_TEST_CASE(VectorTestCase)
{
  VectorTest();
} // VectorTestCase

BOOST_AUTO_TEST_CASE(MatrixTestCase)
{
  MatrixTest();
} // MatrixTestCase

BOOST_AUTO_TEST_CASE(TensorRank3TestCase)
{
  TensorRank3Test();
} // TensorRank3TestCase
