////////////////////////////////////////////////////////////////////////
///
/// \file   KalmanLinearAlgebra.h
///
/// \brief  Kalman filter linear algebra typedefs.
///
/// \author H. Greenlee
///
/// There are various linear algebra typedefs defined in this header:
///
/// 1. KVector<N>::type - Vector, dimension N.
/// 2. KSymMatrix<N>::type - Symmetric matrix, dimension NxN.
/// 3. KMatrix<N,M>::type - A matrix with dimension NxM.
/// 4. KHMatrix<N>::type - Matrix with dimension Nx5 (H-matrix).
/// 5. KGMatrix<N>::type - Matrix with dimension 5xN (gain matrix).
/// 6. TrackVector - Track state vector, fixed dimension 5
/// 7. TrackError - Track error matrix, fixed dimension 5x5.
/// 8. TrackMatrix - General matrix, fixed dimension 5x5.
///
/// The above typedefs refer to some ublas (boost/numeric/ublas)
/// linear algebra class.  The templated typedefs are defined as a
/// member of a template class because c++ doesn't have template
/// typedefs.
///
/// All linear algebra objects use the following storage model.
///
/// 1.  Matrices are stored in row major order (normal c/c++ convention).
/// 2.  Symmetric matrices are stored in lower triangular format.
/// 3.  Amount of preallocated stack memory is specified at compilation
///     time (the actual size of objects must still be specified at run
///     time).
///
/// Surprisingly, the ublas linear algebra package does not have any
/// built-in symmetric matrix inverse function.  We provide one here
/// as free function syminvert.
///
////////////////////////////////////////////////////////////////////////

#ifndef KALMANLINEARALGEBRA_H
#define KALMANLINEARALGEBRA_H

#include "boost/numeric/ublas/lu.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/symmetric.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/serialization/array_wrapper.hpp" // workaround for deficiency in boost 1.64
#include <cmath>

namespace trkf {

  /// Define a shortened alias for ublas namespace.
  namespace ublas = boost::numeric::ublas;

  /// Vector, dimension N.
  template <int N>
  struct KVector {
    typedef ublas::vector<double, ublas::bounded_array<double, N>> type;
  };

  /// Symmetric matrix, dimension NxN.
  template <int N>
  struct KSymMatrix {
    typedef ublas::symmetric_matrix<double,
                                    ublas::lower,
                                    ublas::row_major,
                                    ublas::bounded_array<double, N*(N + 1) / 2>>
      type;
  };

  /// General matrix, dimension NxM.
  template <int N, int M>
  struct KMatrix {
    typedef ublas::matrix<double, ublas::row_major, ublas::bounded_array<double, N * M>> type;
  };

  /// Kalman H-matrix, dimension Nx5.
  template <int N>
  struct KHMatrix {
    typedef typename KMatrix<N, 5>::type type;
  };

  /// Kalman gain matrix, dimension 5xN.
  template <int N>
  struct KGMatrix {
    typedef typename KMatrix<5, N>::type type;
  };

  /// Track state vector, dimension 5.
  typedef typename KVector<5>::type TrackVector;

  /// Track error matrix, dimension 5x5.
  typedef typename KSymMatrix<5>::type TrackError;

  /// General 5x5 matrix.
  typedef typename KMatrix<5, 5>::type TrackMatrix;

  /// Invert symmetric matrix (return false if singular).
  ///
  /// The method used is Cholesky decomposition.
  /// This method is efficient and stable for positive-definite matrices.
  /// In case the matrix is not positive-definite, this method will usually
  /// work, but there can be some numerical pathologies, including "false
  /// singular" failures, and numerical instability.
  /// In the Kalman filter, we expect that this method will be used
  /// exclusively for positive-definite matrices.
  ///
  template <class T, class TRI, class L, class A>
  bool syminvert(ublas::symmetric_matrix<T, TRI, L, A>& m)
  {
    typedef typename ublas::symmetric_matrix<T, TRI, L, A>::size_type size_type;
    typedef typename ublas::symmetric_matrix<T, TRI, L, A>::value_type value_type;

    // In situ Cholesky decomposition m = LDL^T.
    // D is diagonal matrix.
    // L is lower triangular with ones on the diagonal (ones not stored).

    for (size_type i = 0; i < m.size1(); ++i) {
      for (size_type j = 0; j <= i; ++j) {

        value_type ele = m(i, j);

        for (size_type k = 0; k < j; ++k)
          ele -= m(k, k) * m(i, k) * m(j, k);

        // Diagonal elements (can't have zeroes).

        if (i == j) {
          if (ele == 0.) return false;
        }

        // Off-diagonal elements.

        else
          ele = ele / m(j, j);

        // Replace element.

        m(i, j) = ele;
      }
    }

    // In situ inversion of D by simple division.
    // In situ inversion of L by back-substitution.

    for (size_type i = 0; i < m.size1(); ++i) {
      for (size_type j = 0; j <= i; ++j) {

        value_type ele = m(i, j);

        // Diagonal elements.

        if (i == j) m(i, i) = 1. / ele;

        // Off diagonal elements.

        else {
          value_type sum = -ele;
          for (size_type k = j + 1; k < i; ++k)
            sum -= m(i, k) * m(k, j);
          m(i, j) = sum;
        }
      }
    }

    // Recompose the inverse matrix in situ by matrix multiplication m = L^T DL.

    for (size_type i = 0; i < m.size1(); ++i) {
      for (size_type j = 0; j <= i; ++j) {

        value_type sum = m(i, i);
        if (i != j) sum *= m(i, j);

        for (size_type k = i + 1; k < m.size1(); ++k)
          sum += m(k, k) * m(k, i) * m(k, j);

        m(i, j) = sum;
      }
    }

    // Done (success).

    return true;
  }

  /// Invert general square matrix by LU decomposition with partial pivoting.
  /// Return false if singular or not square.
  ///
  template <class T, class L, class A>
  bool invert(ublas::matrix<T, L, A>& m)
  {
    // Make sure matrix is square.

    if (m.size1() != m.size2()) return false;

    // Create permutation matrix for pivoting.

    ublas::permutation_matrix<std::size_t> pm(m.size1());

    // Make temp copy of input matrix.

    ublas::matrix<T, L, A> mcopy(m);

    // Do LU factorization with partial pivoting.
    // This step will fail if matrix is singular.

    int res = lu_factorize(mcopy, pm);
    if (res != 0) return false;

    // Set original matrix to the identity matrix.

    m.assign(ublas::identity_matrix<T>(m.size1()));

    // Do backsubstitution.

    lu_substitute(mcopy, pm, m);

    // Done (success).

    return true;
  }

  /// Trace of matrix.
  ///
  template <class M>
  typename M::value_type trace(const M& m)
  {
    typename M::size_type n = std::min(m.size1(), m.size2());
    typename M::value_type result = 0.;

    for (typename M::size_type i = 0; i < n; ++i)
      result += m(i, i);

    return result;
  }
}

#endif
