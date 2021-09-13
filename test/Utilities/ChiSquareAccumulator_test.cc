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
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()
#include <boost/test/tools/floating_point_comparison.hpp> // BOOST_CHECK_CLOSE()

// LArSoft libraries
#include "lardata/Utilities/ChiSquareAccumulator.h"

// C/C++ standard libraries
#include <type_traits> // std::is_same<>


//------------------------------------------------------------------------------
void testChiSquareAccumulator() {

  auto one = [](double){ return 1.0; };

  auto chiSquare = lar::util::makeChiSquareAccumulator(one);

  BOOST_CHECK_CLOSE(chiSquare.expected(1.0), 1.0, 1e-4);
  BOOST_CHECK_CLOSE(chiSquare.expected(2.0), 1.0, 1e-4);
  BOOST_CHECK_CLOSE(chiSquare.expected(3.0), 1.0, 1e-4);

  BOOST_CHECK_EQUAL(chiSquare.N(), 0U);
  BOOST_CHECK_EQUAL(chiSquare(), 0.0);
  BOOST_CHECK_EQUAL(double(chiSquare), 0.0);
  BOOST_CHECK_EQUAL(chiSquare.chiSquare(), 0.0);

  chiSquare.add(1.0, 1.0); // uncertainty: 1
  BOOST_CHECK_EQUAL(chiSquare.N(), 1U);
  BOOST_CHECK_SMALL(chiSquare(), 1e-5);
  BOOST_CHECK_SMALL(double(chiSquare), 1e-5);
  BOOST_CHECK_SMALL(chiSquare.chiSquare(), 1e-5);

  chiSquare.add(2.0, 0.5); // uncertainty: 1
  BOOST_CHECK_EQUAL(chiSquare.N(), 2U);
  BOOST_CHECK_CLOSE(chiSquare(), 0.25, 1e-4); // at 0.0001%
  BOOST_CHECK_CLOSE(double(chiSquare), 0.25, 1e-4);
  BOOST_CHECK_CLOSE(chiSquare.chiSquare(), 0.25, 1e-4);

  chiSquare.add(3.0, 2.0, 0.5);
  BOOST_CHECK_EQUAL(chiSquare.N(), 3U);
  BOOST_CHECK_CLOSE(chiSquare(), 4.25, 1e-4); // at 0.0001%
  BOOST_CHECK_CLOSE(double(chiSquare), 4.25, 1e-4);
  BOOST_CHECK_CLOSE(chiSquare.chiSquare(), 4.25, 1e-4);

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

  BOOST_CHECK_CLOSE(chi2value, 8.0, 0.001); // up to 10^-5
  BOOST_CHECK_EQUAL(degreesOfFreedom, 0U);

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

  BOOST_CHECK_EQUAL(chiSquare.expected(-2.0), 0.0);
  BOOST_CHECK_EQUAL(chiSquare.expected(0.0), 0.0);
  BOOST_CHECK_EQUAL(chiSquare.expected(2.0), 0.0);
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

  BOOST_CHECK_EQUAL(chiSquare.expected(-2.0F), 0.0F);
  BOOST_CHECK_EQUAL(chiSquare.expected(0.0F), 0.0F);
  BOOST_CHECK_EQUAL(chiSquare.expected(2.0F), 0.0F);
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

