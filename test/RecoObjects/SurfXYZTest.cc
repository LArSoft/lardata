#define BOOST_TEST_MODULE ( SurfXYZTest )
#include "cetlib/quiet_unit_test.hpp"
#include "boost/test/tools/floating_point_comparison.hpp"

//
// File: SurfXYZTest.cxx
//
// Purpose: Unit test for SurfXYZPlane.
//

#include <cmath>
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "lardata/RecoObjects/KalmanLinearAlgebra.h"
#include "cetlib_except/exception.h"

struct SurfXYZTestFixture
{
  SurfXYZTestFixture() :
    surf1(),
    surf2(0., 0., 0., 0., 0.),
    surf3(1., 1., 1., 0., 0.),
    surf4(1., 2., 2., 1., 0.1) {}
  trkf::SurfXYZPlane surf1;  // Default constructed.
  trkf::SurfXYZPlane surf2;  // Same as surf1.
  trkf::SurfXYZPlane surf3;  // Different origin, parallel to surf1 and surf2.
  trkf::SurfXYZPlane surf4;  // Not parallel.
};

BOOST_FIXTURE_TEST_SUITE(SurfXYZTest, SurfXYZTestFixture)

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
  BOOST_CHECK(surf1.distanceTo(surf3) == 1.);
  BOOST_CHECK(surf3.distanceTo(surf1) == -1.);
}

// Should throw exception (not parallel).

BOOST_AUTO_TEST_CASE(NotParallel) {
  BOOST_CHECK_EXCEPTION( surf1.distanceTo(surf4), cet::exception, \
                         [](cet::exception const & e)	          \
			 {				          \
			   return e.category() == "SurfXYZPlane";  \
			 } );
}

// Test track parameters.

BOOST_AUTO_TEST_CASE(TrackParameters) {
  trkf::TrackVector v(5);
  v(0) = 0.1;   // u.
  v(1) = 0.2;   // v.
  v(2) = 2.;    // du/dw.
  v(3) = 3.;    // dv/dw.
  v(4) = 0.5;   // p = 2 GeV.

  // For this vector, the direction cosines are.
  // du/ds = 2./sqrt(14.);
  // dv/ds = 3./sqrt(14.);
  // dw/ds = 1./sqrt(14.);

  double xyz[3];
  double mom[3];
  surf1.getPosition(v, xyz);
  BOOST_CHECK_CLOSE(xyz[0], 0.1, 1.e-6);
  BOOST_CHECK_CLOSE(xyz[1], 0.2, 1.e-6);
  BOOST_CHECK_CLOSE(xyz[2], 0.0, 1.e-6);
  surf3.getPosition(v, xyz);
  BOOST_CHECK_CLOSE(xyz[0], 1.1, 1.e-6);
  BOOST_CHECK_CLOSE(xyz[1], 1.2, 1.e-6);
  BOOST_CHECK_CLOSE(xyz[2], 1.0, 1.e-6);
  surf1.getMomentum(v, mom, trkf::Surface::FORWARD);
  BOOST_CHECK_CLOSE(mom[0], 4./std::sqrt(14.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[1], 6./std::sqrt(14.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[2], 2./std::sqrt(14.), 1.e-6);
  surf1.getMomentum(v, mom, trkf::Surface::BACKWARD);
  BOOST_CHECK_CLOSE(mom[0], -4./std::sqrt(14.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[1], -6./std::sqrt(14.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[2], -2./std::sqrt(14.), 1.e-6);
  surf4.getMomentum(v, mom, trkf::Surface::FORWARD);
  BOOST_CHECK_CLOSE(mom[0], (4.*std::cos(0.1) + 2.*std::sin(0.1))/std::sqrt(14.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[1], (4.*std::sin(0.1)*std::sin(1.) + 6.*std::cos(1.) - 2.*std::cos(0.1)*std::sin(1.))/std::sqrt(14.), 1.e-6);
  BOOST_CHECK_CLOSE(mom[2], (-4.*std::sin(0.1)*std::cos(1.) + 6.*std::sin(1.) + 2.*std::cos(0.1)*std::cos(1.))/std::sqrt(14.), 1.e-6);

  // Should throw exception (no direction).

  BOOST_CHECK_EXCEPTION( surf1.getMomentum(v, mom), cet::exception, \
                         [](cet::exception const & e)		    \
			 {					    \
			   return e.category() == "SurfXYZPlane";    \
			 } );
}

BOOST_AUTO_TEST_SUITE_END()
