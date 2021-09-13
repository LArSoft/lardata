/**
 * @file    FastMatrixMath_test.cc
 * @brief   Tests the classes in SimpleFits.h
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    20141229
 * @version 1.0
 * @see     SimpleFits.h
 *
 * See http://www.boost.org/libs/test for the Boost test library home page.
 *
 * Timing:
 * not given yet
 */

// define the following symbol to print the matrices
// #define FASTMATRIXMATH_TEST_DEBUG

// C/C++ standard libraries
#include <cmath>
#include <limits> // std::numeric_limits<>
#include <array>
#include <algorithm> // std::generate()
#include <random>
#include <string>
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
#define BOOST_TEST_MODULE ( FastMatrixMath_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()
#include <boost/test/tools/floating_point_comparison.hpp> // BOOST_CHECK_CLOSE()

// LArSoft libraries
#include "lardata/Utilities/FastMatrixMathHelper.h"


//==============================================================================
//=== Test code
//===

constexpr unsigned int static_sqrt(unsigned int n) {
  // too lazy to copy sqrt from Stroustrup
  return
    ((n ==  1)? 1:
    ((n ==  4)? 2:
    ((n ==  9)? 3:
    ((n == 16)? 4:
    std::numeric_limits<unsigned int>::max()))));
} // static_sqrt()


template <typename Array>
void PrintMatrix
  (std::ostream& out, Array const& m, std::string name = "Matrix")
{
#ifdef FASTMATRIXMATH_TEST_DEBUG
  const size_t Dim = size_t(std::sqrt(m.size()));

  out << name << " " << Dim << "x" << Dim << ":";
  for (size_t r = 0; r < Dim; ++r) {
    out << "\n |";
    for (size_t c = 0; c < Dim; ++c)
      out << " " << m[r * Dim + c];
    out << " |";
  } // for
  out << std::endl;
#endif // FASTMATRIXMATH_TEST_DEBUG
} // PrintMatrix()


template <typename Array>
void CheckSymmetric(Array const& m) {
  const size_t Dim = size_t(std::sqrt(m.size()));

  for (size_t r = 1; r < Dim; ++r)
    for (size_t c = r + 1; c < Dim; ++c)
      BOOST_CHECK_CLOSE(m[r * Dim + c], m[r * Dim + c], 1e-3); // at 0.001%

} // CheckSymmetric()


template <typename Array>
void CheckInverse(Array const& a, Array const& a_inv) {

  const size_t Dim = size_t(std::sqrt(a.size()));
  using Data_t = typename Array::value_type;

  BOOST_CHECK_EQUAL(a.size(), a_inv.size());

  for (size_t r = 0; r < Dim; ++r) {
    for (size_t c = 0; c < Dim; ++c) {
      Data_t v = 0.;
      for (size_t k = 0; k < Dim; ++k) {
        v += a[r * Dim + k] * a_inv[k * Dim + c];
      } // for
      if (r == c) BOOST_CHECK_CLOSE(v, Data_t(1), 0.01); // 0.01%
      else        BOOST_CHECK_SMALL(v, 1e-5);
    } // for column
  } // for row

} // CheckInverse()


template <typename Array>
void MatrixTest(Array const& mat) {

  using Data_t = typename Array::value_type;

  constexpr unsigned int Dim = static_sqrt(std::tuple_size<Array>::value);
  static_assert((Dim >= 1) && (Dim <= 4), "Dimension not supported");

  using FastMatrixOperations
    = lar::util::details::FastMatrixOperations<Data_t, Dim>;

  Array mat_inv = FastMatrixOperations::InvertMatrix(mat);
  PrintMatrix(std::cout, mat_inv, "Alleged inverse matrix");
  CheckInverse(mat, mat_inv);

} // MatrixTest()


template <typename Array>
void MatrixTest(Array const& mat, typename Array::value_type det) {

  using Data_t = typename Array::value_type;

  constexpr unsigned int Dim = static_sqrt(std::tuple_size<Array>::value);
  static_assert((Dim >= 1) && (Dim <= 4), "Dimension not supported");

  using FastMatrixOperations
    = lar::util::details::FastMatrixOperations<Data_t, Dim>;

  const Data_t my_det = FastMatrixOperations::Determinant(mat);
  BOOST_CHECK_CLOSE(my_det, det, 1e-4); // 0.0001%

  if (std::isnormal(det)) {
    Array mat_inv = FastMatrixOperations::InvertMatrix(mat, det);
    PrintMatrix(std::cout, mat_inv, "Alleged inverse matrix");
    CheckInverse(mat, mat_inv);
  }
} // MatrixTest()


template <typename Array>
void SymmetricMatrixTest(Array const& mat) {

  using Data_t = typename Array::value_type;

  constexpr unsigned int Dim = static_sqrt(std::tuple_size<Array>::value);
  static_assert((Dim >= 1) && (Dim <= 4), "Dimension not supported");

  using FastMatrixOperations
    = lar::util::details::FastMatrixOperations<Data_t, Dim>;

  CheckSymmetric(mat);

  Array mat_inv = FastMatrixOperations::InvertSymmetricMatrix(mat);
  PrintMatrix(std::cout, mat_inv, "Alleged inverse matrix");
  CheckInverse(mat, mat_inv);
  CheckSymmetric(mat_inv);
} // SymmetricMatrixTest()


template <typename Array>
void SymmetricMatrixTest(Array const& mat, typename Array::value_type det) {

  using Data_t = typename Array::value_type;

  constexpr unsigned int Dim = static_sqrt(std::tuple_size<Array>::value);
  static_assert((Dim >= 1) && (Dim <= 4), "Dimension not supported");

  using FastMatrixOperations
    = lar::util::details::FastMatrixOperations<Data_t, Dim>;

  const Data_t my_det = FastMatrixOperations::Determinant(mat);
  BOOST_CHECK_CLOSE(my_det, det, 1e-4); // 0.0001%

  if (std::isnormal(det)) {
    Array mat_inv = FastMatrixOperations::InvertSymmetricMatrix(mat, det);
    PrintMatrix(std::cout, mat_inv, "Alleged inverse matrix");
    CheckInverse(mat, mat_inv);
    CheckSymmetric(mat_inv);
  }
} // SymmetricMatrixTest()



template <typename T>
void TestMatrix2x2() {
  using Data_t = T;
  constexpr unsigned int Dim = 2; // we are testing 2x2 matrices

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<Data_t, Dim*Dim> matrix = {{
    Data_t(2), Data_t(3),
    Data_t(4), Data_t(1)
  }}; // matrix
  const Data_t true_det = Data_t(-10);

  PrintMatrix(std::cout, matrix, "Matrix");
  MatrixTest(matrix, true_det);
} // TestMatrix2x2()


template <typename T>
void TestSymmetricMatrix2x2() {
  using Data_t = T;
  constexpr unsigned int Dim = 2; // we are testing 2x2 matrices

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<Data_t, Dim*Dim> matrix = {{
    Data_t(2), Data_t(3),
    Data_t(3), Data_t(1)
  }}; // matrix
  const Data_t true_det = Data_t(-7);

  PrintMatrix(std::cout, matrix, "Symmetric matrix");
  SymmetricMatrixTest(matrix, true_det);

} // TestSymmetricMatrix2x2()


template <typename T>
void TestMatrix3x3_1() {

  using Data_t = T;
  constexpr unsigned int Dim = 3; // we are testing 3x3 matrices

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<Data_t, Dim*Dim> matrix = {{
    Data_t(2), Data_t(0), Data_t(3),
    Data_t(0), Data_t(3), Data_t(0),
    Data_t(4), Data_t(0), Data_t(1),
  }}; // matrix
  const Data_t true_det = Data_t(-30);

  PrintMatrix(std::cout, matrix, "Matrix");
  MatrixTest(matrix, true_det);

} // TestMatrix3x3_1()


template <typename T>
void TestMatrix3x3_2() {

  using Data_t = T;
  constexpr unsigned int Dim = 3; // we are testing 3x3 matrices

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<Data_t, Dim*Dim> matrix = {{
    Data_t(2), Data_t(4), Data_t(3),
    Data_t(0), Data_t(3), Data_t(0),
    Data_t(4), Data_t(0), Data_t(1),
  }}; // matrix
  const Data_t true_det = Data_t(-30);

  PrintMatrix(std::cout, matrix, "Matrix");
  MatrixTest(matrix, true_det);

} // TestMatrix3x3_2()



template <typename T>
void TestSymmetricMatrix3x3() {

  using Data_t = T;
  constexpr unsigned int Dim = 3; // we are testing 3x3 matrices

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<Data_t, Dim*Dim> matrix = {{
    Data_t(2), Data_t(0), Data_t(3),
    Data_t(0), Data_t(3), Data_t(0),
    Data_t(3), Data_t(0), Data_t(1),
  }}; // matrix
  const Data_t true_det = Data_t(-21);

  PrintMatrix(std::cout, matrix, "Symmetric matrix");
  SymmetricMatrixTest(matrix, true_det);

} // TestSymmetricMatrix3x3()


template <typename T>
void TestMatrix4x4_1() {

  using Data_t = T;
  constexpr unsigned int Dim = 4; // we are testing 4x4 matrices

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<Data_t, Dim*Dim> matrix = {{
    Data_t(2), Data_t(0), Data_t(3), Data_t(0),
    Data_t(0), Data_t(3), Data_t(0), Data_t(6),
    Data_t(4), Data_t(0), Data_t(1), Data_t(0),
    Data_t(0), Data_t(2), Data_t(0), Data_t(7)
  }}; // matrix
  const Data_t true_det = Data_t(-90);

  PrintMatrix(std::cout, matrix, "Matrix");
  MatrixTest(matrix, true_det);

} // TestMatrix4x4_1()


template <typename T>
void TestMatrix4x4_2() {

  using Data_t = T;
  constexpr unsigned int Dim = 4; // we are testing 4x4 matrices

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<Data_t, Dim*Dim> matrix = {{
    Data_t(2), Data_t(0), Data_t(3), Data_t(0),
    Data_t(5), Data_t(3), Data_t(0), Data_t(6),
    Data_t(4), Data_t(0), Data_t(1), Data_t(0),
    Data_t(3), Data_t(2), Data_t(0), Data_t(7)
  }}; // matrix
  const Data_t true_det = Data_t(-90);

  PrintMatrix(std::cout, matrix, "Matrix");
  MatrixTest(matrix, true_det);

} // TestMatrix4x4_2()


template <typename T, unsigned int Dim>
void TestMatrix_N(unsigned int N = 100) {

  using Data_t = T;

  std::default_random_engine engine;
  std::uniform_real_distribution<Data_t> uniform(Data_t(-10.), Data_t(10.));

  std::array<Data_t, Dim*Dim> matrix;
  for (unsigned int i = 0; i < N; ++i) {

    std::generate(matrix.begin(), matrix.end(),
      [&engine, &uniform] { return uniform(engine); }
      );

    PrintMatrix(std::cout, matrix, "Matrix");
    MatrixTest(matrix);
  } // for

} // TestMatrix_N()


template <typename T, unsigned int Dim>
void TestNullMatrix() {

  using Data_t = T;

  std::array<Data_t, Dim*Dim> matrix;
  matrix.fill(Data_t(0));

  PrintMatrix(std::cout, matrix, "Empty matrix");
  MatrixTest(matrix, Data_t(0));
  PrintMatrix(std::cout, matrix, "Empty symmetric matrix");
  SymmetricMatrixTest(matrix, Data_t(0));

} // TestNullMatrix()



template <typename T>
void TestSymmetricMatrix4x4() {

  using Data_t = T;
  constexpr unsigned int Dim = 4; // we are testing 4x4 matrices

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  std::array<Data_t, Dim*Dim> matrix = {{
    Data_t(2), Data_t(0), Data_t(3), Data_t(0),
    Data_t(0), Data_t(3), Data_t(0), Data_t(2),
    Data_t(3), Data_t(0), Data_t(1), Data_t(0),
    Data_t(0), Data_t(2), Data_t(0), Data_t(7)
  }}; // matrix
  const Data_t true_det = Data_t(-119);

  PrintMatrix(std::cout, matrix, "Symmetric matrix");
  SymmetricMatrixTest(matrix, true_det);

} // TestSymmetricMatrix4x4()


//------------------------------------------------------------------------------
//--- registration of tests
//
// Boost needs now to know which tests we want to run.
// Tests are "automatically" registered, hence the BOOST_AUTO_TEST_CASE()
// macro name. The argument is a name for the test; each test will have a
// number of checks and it will fail if any of them does.
//

//
// Matrix tests
//
BOOST_AUTO_TEST_CASE(Matrix2x2RealTest) {
  TestMatrix2x2<double>();
  TestSymmetricMatrix2x2<double>();

  TestMatrix_N<double, 2>();
  TestNullMatrix<double, 2>();
}

BOOST_AUTO_TEST_CASE(Matrix3x3RealTest) {
  TestMatrix3x3_1<double>();
  TestMatrix3x3_2<double>();
  TestSymmetricMatrix3x3<double>();

  TestMatrix_N<double, 3>();
  TestNullMatrix<double, 3>();
}

BOOST_AUTO_TEST_CASE(Matrix4x4RealTest) {
  TestMatrix4x4_1<double>();
  TestMatrix4x4_2<double>();
  TestSymmetricMatrix4x4<double>();

  TestMatrix_N<double, 4>();
  TestNullMatrix<double, 4>();
}
