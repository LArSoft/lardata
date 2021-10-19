/**
 * @file    ChiSquareAccumulator_test.cc
 * @brief   Tests the classes in `ChiSquareAccumulator.h`
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    July 26, 2018
 * @version 1.0
 * @see     `lardata/Utilities/ChiSquareAccumulator.h`
 *
 * See http://www.boost.org/libs/test for the Boost test library home page.
 *
 * Timing:
 * not given yet
 */


// Boost libraries
#define BOOST_TEST_MODULE ( ChiSquareAccumulator_test )
#include "boost/test/unit_test.hpp"

// LArSoft libraries
#include "lardata/Utilities/ChiSquareAccumulator.h"

// C/C++ standard libraries
#include <type_traits> // std::is_same<>

using boost::test_tools::tolerance;

//------------------------------------------------------------------------------
void testChiSquareAccumulator() {

  auto one = [](double){ return 1.0; };
  auto chiSquare = lar::util::makeChiSquareAccumulator(one);

  BOOST_TEST(chiSquare.expected(1.0) == 1.0);
  BOOST_TEST(chiSquare.expected(2.0) == 1.0);
  BOOST_TEST(chiSquare.expected(3.0) == 1.0);

  BOOST_TEST(chiSquare.N() == 0U);
  BOOST_TEST(chiSquare() == 0.0);
  BOOST_TEST(double(chiSquare) == 0.0);
  BOOST_TEST(chiSquare.chiSquare() == 0.0);

  chiSquare.add(1.0, 1.0); // uncertainty: 1
  BOOST_TEST(chiSquare.N() == 1U);
  BOOST_TEST(chiSquare() == 0, 1e-5% tolerance());
  BOOST_TEST(double(chiSquare) == 0, 1e-5% tolerance());
  BOOST_TEST(chiSquare.chiSquare() == 0, 1e-5% tolerance());

  chiSquare.add(2.0, 0.5); // uncertainty: 1
  BOOST_TEST(chiSquare.N() == 2U);
  BOOST_TEST(chiSquare() == 0.25, 1e-4% tolerance());
  BOOST_TEST(double(chiSquare) == 0.25, 1e-4% tolerance());
  BOOST_TEST(chiSquare.chiSquare() == 0.25, 1e-4% tolerance());

  chiSquare.add(3.0, 2.0, 0.5);
  BOOST_TEST(chiSquare.N() == 3U);
  BOOST_TEST(chiSquare() == 4.25, 1e-4% tolerance());
  BOOST_TEST(double(chiSquare) == 4.25, 1e-4% tolerance());
  BOOST_TEST(chiSquare.chiSquare() == 4.25, 1e-4% tolerance());

} // testChiSquareAccumulator()


//------------------------------------------------------------------------------
void testChiSquareAccumulator_documentation() {
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * double const a =  2.0;
   * double const b = -1.0;
   * auto f = [a,b](double x){ return a + b * x; };
   * lar::util::ChiSquareAccumulator<decltype(f)> chiSquare;
   *
   * chiSquare.add(0.0, 1.0, 0.5); // add ( 0 ; 1.0 +/- 0.5 )
   * chiSquare.add(1.0, 1.0, 0.5); // add ( 1 ; 1.0 +/- 0.5 )
   * chiSquare.add(2.0, 1.0, 0.5); // add ( 2 ; 1.0 +/- 0.5 )
   *
   * double const chi2value = chiSquare();
   * int degreesOfFreedom = int(chiSquare.N()) - 3;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * promised `chi2value` `8.0` and `degreeOfFreedom` `0`.
   */
  double const a =  2.0;
  double const b = -1.0;
  auto f = [a,b](double x){ return a + b * x; };
  lar::util::ChiSquareAccumulator<decltype(f)> chiSquare(f);

  chiSquare.add(0.0, 1.0, 0.5); // add ( 0 ; 1.0 +/- 0.5 )
  chiSquare.add(1.0, 1.0, 0.5); // add ( 1 ; 1.0 +/- 0.5 )
  chiSquare.add(2.0, 1.0, 0.5); // add ( 2 ; 1.0 +/- 0.5 )

  double const chi2value = chiSquare();
  int degreesOfFreedom = chiSquare.N() - 3;

  BOOST_TEST(chi2value == 8.0, 0.001% tolerance());
  BOOST_TEST(degreesOfFreedom == 0U);

} // testChiSquareAccumulator_documentation();


//------------------------------------------------------------------------------
void testMakeChiSquareAccumulator_documentation1() {

  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto zero = [](double){ return 0.0; }; // expectation function
   * auto chiSquare = lar::util::makeChiSquareAccumulator(zero);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * declare `chiSquare` in a way equivalent to:
   * `lar::util::ChiSquareAccumulator<decltype(zero)> chiSquare(zero)`.
   */
  auto zero = [](double){ return 0.0; }; // expectation function
  auto const& chiSquare = lar::util::makeChiSquareAccumulator(zero);

  BOOST_TEST(chiSquare.expected(-2.0) == 0.0);
  BOOST_TEST(chiSquare.expected(0.0) == 0.0);
  BOOST_TEST(chiSquare.expected(2.0) == 0.0);
  static_assert(std::is_same<decltype(chiSquare()), double>::value,
    "makeChiSquareAccumulator() returned an unexpected type!"
    );

} // testMakeChiSquareAccumulator_documentation1()


//------------------------------------------------------------------------------
void testMakeChiSquareAccumulator_documentation2() {

  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto zero = [](float){ return 0.0F; }; // expectation function
   * auto chiSquare = lar::util::makeChiSquareAccumulator<float>(zero);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * declare `chiSquare` in a way equivalent to:
   * `lar::util::ChiSquareAccumulator<decltype(zero), float> chiSquare(zero)`.
   */
  auto zero = [](float){ return 0.0F; }; // expectation function
  auto chiSquare = lar::util::makeChiSquareAccumulator<float>(zero);

  BOOST_TEST(chiSquare.expected(-2.0F) == 0.0F);
  BOOST_TEST(chiSquare.expected(0.0F) == 0.0F);
  BOOST_TEST(chiSquare.expected(2.0F) == 0.0F);
  static_assert(std::is_same<decltype(chiSquare()), float>::value,
    "makeChiSquareAccumulator<float>() returned an unexpected type!"
    );

} // testMakeChiSquareAccumulator_documentation2()


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(ChiSquareAccumulatorTestCase) {

  testChiSquareAccumulator();
  testChiSquareAccumulator_documentation();
  testMakeChiSquareAccumulator_documentation1();
  testMakeChiSquareAccumulator_documentation2();

} // ChiSquareAccumulatorTestCase


//------------------------------------------------------------------------------
