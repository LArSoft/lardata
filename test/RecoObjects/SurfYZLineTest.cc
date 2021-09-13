#define BOOST_TEST_MODULE ( SurfYZTest )
#include "cetlib/quiet_unit_test.hpp"
#include "boost/test/tools/floating_point_comparison.hpp"

//
// File: SurfYZLineTest.cxx
//
// Purpose: Unit test for SurfYZLine.
//

#include <cmath>
#include "lardata/RecoObjects/SurfYZLine.h"
#include "lardata/RecoObjects/KalmanLinearAlgebra.h"
#include "cetlib_except/exception.h"

struct SurfYZLineTestFixture
{
  SurfYZLineTestFixture() :
    surf1(),
    surf2(0., 0., 0., 0.),
    surf3(1., 1., 1., 0.),
    surf4(2., 2., 2., 1.) {}
  trkf::SurfYZLine surf1;  // Default constructed.
  trkf::SurfYZLine surf2;  // Same as surf1.
  trkf::SurfYZLine surf3;  // Different origin, parallel to surf1 and surf2.
  trkf::SurfYZLine surf4;  // Not parallel.
};

BOOST_FIXTURE_TEST_SUITE(SurfYZLineTest, SurfYZLineTestFixture)

// Test equality comparisons.

BOOST_AUTO_TEST_CASE(Equality) {
  BOOST_CHECK(surf1.isEqual(surf2));
  BOOST_CHECK(!surf1.isEqual(surf3));
  BOOST_CHECK(!surf1.isEqual(surf4));
  BOOST_CHECK(!surf2.isEqual(surf3));
  BOOST_CHECK(!surf2.isEqual(surf4));
  BOOST_CHECK(!surf3.isEqual(surf4));
}

// Test parallel comparisions.

BOOST_AUTO_TEST_CASE(Parallel) {
  BOOST_CHECK(surf1.isParallel(surf2));
  BOOST_CHECK(surf1.isParallel(surf3));
  BOOST_CHECK(!surf1.isParallel(surf4));
  BOOST_CHECK(surf2.isParallel(surf3));
  BOOST_CHECK(!surf2.isParallel(surf4));
  BOOST_CHECK(!surf3.isParallel(surf4));
}

// Test coordinate transformations.

BOOST_AUTO_TEST_CASE(Transformation) {
  double xyz1[3] = {1., 2., 3.};
  double xyz2[3];
  double uvw[3];
  surf4.toLocal(xyz1, uvw);
  surf4.toGlobal(uvw, xyz2);
  for(int i=0; i<3; ++i)
    BOOST_CHECK_CLOSE(xyz1[i], xyz2[i], 1.e-6);
}

// Test separation.

BOOST_AUTO_TEST_CASE(Separation) {
  BOOST_CHECK(surf1.distanceTo(surf2) == 0.);
  BOOST_CHECK_CLOSE(surf1.distanceTo(surf3), sqrt(2.), 1.e-6);
  BOOST_CHECK_CLOSE(surf3.distanceTo(surf1), sqrt(2.), 1.e-6);
}

// Should throw exception (not parallel).

BOOST_AUTO_TEST_CASE(NotParallel) {
  BOOST_CHECK_EXCEPTION( surf1.distanceTo(surf4), cet::exception, \
                         [](cet::exception const & e)	          \
			 {				          \
			   return e.category() == "SurfYZLine";  \
			 } );
}

// Test track parameters.

BOOST_AUTO_TEST_CASE(TrackParameters) {
  trkf::TrackVector v(5);
  v(0) = 0.1;   // r.
  v(1) = 0.2;   // v.
  v(2) = 2.;    // phi.
  v(3) = 1.;    // eta.
  v(4) = 0.5;   // p = 2 GeV.

  double xyz[3];
  double mom[3];
  surf1.getPosition(v, xyz);
  BOOST_CHECK_CLOSE(xyz[0], -0.1*std::sin(2.), 1.e-6);
  BOOST_CHECK_CLOSE(xyz[1], 0.2, 1.e-6);
  BOOST_CHECK_CLOSE(xyz[2], 0.1*std::cos(2.), 1.e-6);
  surf3.getPosition(v, xyz);
  BOOST_CHECK_CLOSE(xyz[0], 1. - 0.1*std::sin(2.), 1.e-6);
  BOOST_CHECK_CLOSE(xyz[1], 1.2, 1.e-6);
  BOOST_CHECK_CLOSE(xyz[2], 1. + 0.1*std::cos(2.), 1.e-6);
  surf1.getMomentum(v, mom, trkf::Surface::FORWARD);
  BOOST_CHECK_CLOSE(mom[0], 2. * std::cos(2.) / std::cosh(1.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[1], 2. * std::tanh(1.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[2], 2. * std::sin(2.) / cosh(1.), 1.e-6);
  surf1.getMomentum(v, mom, trkf::Surface::BACKWARD);
  BOOST_CHECK_CLOSE(mom[0], 2. * std::cos(2.) / std::cosh(1.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[1], 2. * std::tanh(1.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[2], 2. * std::sin(2.) / cosh(1.), 1.e-6);
  surf1.getMomentum(v, mom);
  BOOST_CHECK_CLOSE(mom[0], 2. * std::cos(2.) / std::cosh(1.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[1], 2. * std::tanh(1.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[2], 2. * std::sin(2.) / cosh(1.), 1.e-6);
}

BOOST_AUTO_TEST_SUITE_END()
