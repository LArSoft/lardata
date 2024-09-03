//
// File: LATest.cxx
//
// Purpose: Test linear algebra.
//

#include "boost/numeric/ublas/io.hpp"
#include "lardata/RecoObjects/KalmanLinearAlgebra.h"

#include <cassert>
#include <cmath>
#include <iostream>

int main()
{
  // Make sure assert is enabled.

  bool assert_flag = false;
  assert((assert_flag = true, assert_flag));
  if (!assert_flag) {
    std::cerr << "Assert is disabled" << std::endl;
    return 1;
  }

  // Test matrix inversion.

  // 1x1 KSymMatrix.
  trkf::KSymMatrix<1>::type m1(1);
  for (unsigned int i = 0; i < m1.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j)
      m1(i, j) = i + j + 1.;
  }
  trkf::KSymMatrix<1>::type minv1(m1);
  bool ok = trkf::syminvert(minv1);
  assert(ok);
  trkf::ublas::matrix<double> unit1 = prod(m1, minv1);
  std::cout << m1 << std::endl;
  std::cout << minv1 << std::endl;
  std::cout << unit1 << std::endl;
  for (unsigned int i = 0; i < m1.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit1(i, j) - val) < 1.e-10);
    }
  }

  // 2x2 KSymMatrix.

  trkf::KSymMatrix<2>::type m2(2);
  for (unsigned int i = 0; i < m2.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j)
      m2(i, j) = i + j + 1.;
  }
  trkf::KSymMatrix<2>::type minv2(m2);
  ok = trkf::syminvert(minv2);
  assert(ok);
  trkf::ublas::matrix<double> unit2 = prod(m2, minv2);
  std::cout << m2 << std::endl;
  std::cout << minv2 << std::endl;
  std::cout << unit2 << std::endl;
  for (unsigned int i = 0; i < m2.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit2(i, j) - val) < 1.e-10);
    }
  }

  // 3x3 KSymMatrix.

  trkf::KSymMatrix<3>::type m3(3);
  for (unsigned int i = 0; i < m3.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      m3(i, j) = i + j + 1.;
      if (i == j) m3(i, j) += 1.;
    }
  }
  trkf::KSymMatrix<3>::type minv3(m3);
  ok = trkf::syminvert(minv3);
  assert(ok);
  trkf::ublas::matrix<double> unit3 = prod(m3, minv3);
  std::cout << m3 << std::endl;
  std::cout << minv3 << std::endl;
  std::cout << unit3 << std::endl;
  for (unsigned int i = 0; i < m3.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit3(i, j) - val) < 1.e-10);
    }
  }

  // 4x4 KSymMatrix.

  trkf::KSymMatrix<4>::type m4(4);
  for (unsigned int i = 0; i < m4.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      m4(i, j) = i + j + 1.;
      if (i == j) m4(i, j) += 1.;
    }
  }
  trkf::KSymMatrix<4>::type minv4(m4);
  ok = trkf::syminvert(minv4);
  assert(ok);
  trkf::ublas::matrix<double> unit4 = prod(m4, minv4);
  std::cout << m4 << std::endl;
  std::cout << minv4 << std::endl;
  std::cout << unit4 << std::endl;
  for (unsigned int i = 0; i < m4.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit4(i, j) - val) < 1.e-10);
    }
  }

  // 5x5 TrackError.

  trkf::TrackError m5;
  for (unsigned int i = 0; i < m5.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      m5(i, j) = i + j + 1.;
      if (i == j) m5(i, j) += 1.;
    }
  }
  trkf::TrackError minv5(m5);
  ok = trkf::syminvert(minv5);
  assert(ok);
  trkf::ublas::matrix<double> unit5 = prod(m5, minv5);
  std::cout << m5 << std::endl;
  std::cout << minv5 << std::endl;
  std::cout << unit5 << std::endl;
  for (unsigned int i = 0; i < m5.size1(); ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit5(i, j) - val) < 1.e-10);
    }
  }

  // 1x1 KMatrix.

  trkf::KMatrix<1, 1>::type m6(1, 1);
  for (unsigned int i = 0; i < m6.size1(); ++i) {
    for (unsigned int j = 0; j < m6.size2(); ++j)
      m6(i, j) = i + 2 * j + 1.;
  }
  trkf::KMatrix<1, 1>::type minv6(m6);
  ok = trkf::invert(minv6);
  assert(ok);
  trkf::ublas::matrix<double> unit6 = prod(m6, minv6);
  std::cout << m6 << std::endl;
  std::cout << minv6 << std::endl;
  std::cout << unit6 << std::endl;
  for (unsigned int i = 0; i < m6.size1(); ++i) {
    for (unsigned int j = 0; j < m6.size2(); ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit6(i, j) - val) < 1.e-10);
    }
  }

  // 2x2 KMatrix.

  trkf::KMatrix<2, 2>::type m7(2, 2);
  for (unsigned int i = 0; i < m7.size1(); ++i) {
    for (unsigned int j = 0; j < m7.size2(); ++j)
      m7(i, j) = i + 2 * j + 1.;
  }
  trkf::KMatrix<2, 2>::type minv7(m7);
  ok = trkf::invert(minv7);
  assert(ok);
  trkf::ublas::matrix<double> unit7 = prod(m7, minv7);
  std::cout << m7 << std::endl;
  std::cout << minv7 << std::endl;
  std::cout << unit7 << std::endl;
  for (unsigned int i = 0; i < m7.size1(); ++i) {
    for (unsigned int j = 0; j < m7.size2(); ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit7(i, j) - val) < 1.e-10);
    }
  }

  // 3x3 KMatrix.

  trkf::KMatrix<3, 3>::type m8(3, 3);
  for (unsigned int i = 0; i < m8.size1(); ++i) {
    for (unsigned int j = 0; j < m8.size2(); ++j)
      m8(i, j) = i + 2 * j + (i == j ? 1. : 0.);
  }
  trkf::KMatrix<3, 3>::type minv8(m8);
  ok = trkf::invert(minv8);
  assert(ok);
  trkf::ublas::matrix<double> unit8 = prod(m8, minv8);
  std::cout << m8 << std::endl;
  std::cout << minv8 << std::endl;
  std::cout << unit8 << std::endl;
  for (unsigned int i = 0; i < m8.size1(); ++i) {
    for (unsigned int j = 0; j < m8.size2(); ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit8(i, j) - val) < 1.e-10);
    }
  }

  // 4x4 KMatrix.

  trkf::KMatrix<4, 4>::type m9(4, 4);
  for (unsigned int i = 0; i < m9.size1(); ++i) {
    for (unsigned int j = 0; j < m9.size2(); ++j)
      m9(i, j) = i + 2 * j + (i == j ? 1. : 0.);
  }
  trkf::KMatrix<4, 4>::type minv9(m9);
  ok = trkf::invert(minv9);
  assert(ok);
  trkf::ublas::matrix<double> unit9 = prod(m9, minv9);
  std::cout << m9 << std::endl;
  std::cout << minv9 << std::endl;
  std::cout << unit9 << std::endl;
  for (unsigned int i = 0; i < m9.size1(); ++i) {
    for (unsigned int j = 0; j < m9.size2(); ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit9(i, j) - val) < 1.e-10);
    }
  }

  // 5x5 TrackMatrix.

  trkf::TrackMatrix m10(5, 5);
  for (unsigned int i = 0; i < m10.size1(); ++i) {
    for (unsigned int j = 0; j < m10.size2(); ++j)
      m10(i, j) = i + 2 * j + (i == j ? 1. : 0.);
  }
  trkf::TrackMatrix minv10(m10);
  ok = trkf::invert(minv10);
  assert(ok);
  trkf::ublas::matrix<double> unit10 = prod(m10, minv10);
  std::cout << m10 << std::endl;
  std::cout << minv10 << std::endl;
  std::cout << unit10 << std::endl;
  for (unsigned int i = 0; i < m10.size1(); ++i) {
    for (unsigned int j = 0; j < m10.size2(); ++j) {
      double val = (i == j ? 1. : 0.);
      assert(std::abs(unit10(i, j) - val) < 1.e-10);
    }
  }

  // Done (success).

  std::cout << "LATest: All tests passed." << std::endl;

  return 0;
}
