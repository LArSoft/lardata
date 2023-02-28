/**
 * @file    FastMatrixMathHelper.h
 * @brief   Classes with hard-coded (hence "fast") matrix math
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 31st, 2015
 *
 * Currently includes:
 *  - determinant (2x2, 3x3, 4x4)
 *  - inversion (2x2, 3x3, 4x4)
 *
 */

#ifndef FASTMATRIXMATHHELPER_H
#define FASTMATRIXMATHHELPER_H 1

// C/C++ standard libraries
#include <algorithm> // std::for_each()
#include <array>
#include <cmath>     // std::sqrt()
#include <iterator>  // std::begin(), std::end()
#include <stdexcept> // std::range_error
#include <tuple>

#include "lardataalg/Utilities/StatCollector.h" // lar::util::identity

namespace lar {
  namespace util {

    namespace details {

      /// Base class with common definitions for "fast" matrix operations.
      template <typename T, unsigned int DIM>
      struct FastMatrixOperationsBase {
        using Data_t = T;

        static constexpr unsigned int Dim = DIM; ///< matrix dimensions

        using Matrix_t = std::array<Data_t, Dim * Dim>;
        using Vector_t = std::array<Data_t, Dim>;

        /// Returns the product of a square matrix times a column vector
        static Vector_t MatrixVectorProduct(Matrix_t const& mat, Vector_t const& vec);

        static constexpr Data_t sqr(Data_t v) { return v * v; }
      }; // struct FastMatrixOperationsBase<>

      /**
       * @brief Provides "fast" matrix operations.
       * @tparam T data type for the elements of the matrix
       * @tparam DIM the dimension of the (square) matrix
       *
       * Actually this class does nothing: specialize it!
       *
       * Once the specialization is in place, this class offers:
       *
       *     constexpr unsigned int Dim = 2;
       *     std::array<float, Dim*Dim> matrix;
       *
       *     float det = FastMatrixOperations<float, Dim>::Determinant();
       *
       *     std::array<float, Dim*Dim> inverse;
       *
       *     // generic inversion
       *     inverse = FastMatrixOperations<float, Dim>::InvertMatrix(matrix);
       *
       *     // faster inversion if we already have the determinant
       *     inverse = FastMatrixOperations<float, Dim>::InvertMatrix
       *       (matrix, det);
       *
       *     // faster inversion if we know the matrix is symmetric
       *     inverse = FastMatrixOperations<float, Dim>::InvertSymmetricMatrix
       *       (matrix);
       *
       *     // even faster inversion if we also know the determinant already
       *     inverse = FastMatrixOperations<float, Dim>::InvertSymmetricMatrix
       *       (matrix, det);
       *
       * Note that the inversion functions do not have a defined policy for
       * non-invertible matrices. If you need to check (and you usually do),
       * compute the determinant first, and invert only if `std::isnormal(det)`.
       */
      template <typename T, unsigned int DIM>
      struct FastMatrixOperations : public FastMatrixOperationsBase<T, DIM> {
        using Base_t = FastMatrixOperationsBase<T, DIM>;
        using Data_t = typename Base_t::Data_t;
        using Matrix_t = typename Base_t::Matrix_t;

        /// Computes the determinant of a matrix
        static Data_t Determinant(Matrix_t const& mat);

        /// Computes the determinant of a matrix, using the provided determinant
        static Matrix_t InvertMatrix(Matrix_t const& mat, Data_t det);

        /// Computes the determinant of a symmatric matrix,
        /// using the provided determinant
        static Matrix_t InvertSymmetricMatrix(Matrix_t const& mat, Data_t det);

        /// Computes the determinant of a matrix
        static Matrix_t InvertMatrix(Matrix_t const& mat)
        {
          return InvertMatrix(mat, Determinant(mat));
        }

        /// Computes the determinant of a matrix
        static Matrix_t InvertSymmetricMatrix(Matrix_t const& mat)
        {
          return InvertSymmetricMatrix(mat, Determinant(mat));
        }

      }; // struct FastMatrixOperations<>

    } // namespace details
  }   // namespace util
} // namespace lar

//******************************************************************************
//***  template implementation
//***
namespace lar {
  namespace util {

    namespace details {

      template <unsigned int NCols>
      constexpr size_t MatrixIndex(unsigned int row, unsigned int col)
      {
        return row * NCols + col;
      }

      template <unsigned int N>
      struct DeterminantHelperBase {
        static constexpr size_t index(unsigned int row, unsigned int col)
        {
          return MatrixIndex<N>(row, col);
        }
      }; // struct DeterminantHelperBase<>

      template <typename T,
                unsigned int N,
                unsigned int... RnC // all row indices, then all column indices
                >
      struct DeterminantHelper : public DeterminantHelperBase<N> {
        using DeterminantHelperBase<N>::index;

        static T compute(T const* data);
      }; // struct DeterminantHelper<>

      /// Determinant of a 1x1 submatrix
      template <typename T, unsigned int N, unsigned int R, unsigned int C>
      struct DeterminantHelper<T, N, R, C> : public DeterminantHelperBase<N> {
        using DeterminantHelperBase<N>::index;
        static_assert(R < N, "invalid row index specified");
        static_assert(C < N, "invalid column index specified");

        static T compute(T const* data) { return data[index(R, C)]; }
      }; // struct DeterminantHelper<>

      /// Determinant of a 2x2 submatrix
      template <typename T,
                unsigned int N,
                unsigned int R1,
                unsigned int R2,
                unsigned int C1,
                unsigned int C2>
      struct DeterminantHelper<T, N, R1, R2, C1, C2> : public DeterminantHelperBase<N> {
        using DeterminantHelperBase<N>::index;
        static_assert(R1 < N, "invalid row index specified");
        static_assert(R2 < N, "invalid row index specified");
        static_assert(C1 < N, "invalid column index specified");
        static_assert(C2 < N, "invalid column index specified");

        static constexpr T sqr(T v) { return v * v; }

        static T compute(T const* data)
        {
          return data[index(R1, C1)] * data[index(R2, C2)] -
                 data[index(R1, C2)] * data[index(R2, C1)];
          /* // also
                data[index(R1, C1)] * UpHelper<R2, C2>::compute(data)
              - data[index(R1, C2)] * UpHelper<R2, C1>::compute(data)
             */
        }
      }; // struct DeterminantHelper<>

      /// Determinant of a 3x3 submatrix
      template <typename T,
                unsigned int N,
                unsigned int R1,
                unsigned int R2,
                unsigned int R3,
                unsigned int C1,
                unsigned int C2,
                unsigned int C3>
      struct DeterminantHelper<T, N, R1, R2, R3, C1, C2, C3> : public DeterminantHelperBase<N> {
        using DeterminantHelperBase<N>::index;
        static_assert(R1 < N, "invalid row index specified");
        static_assert(R2 < N, "invalid row index specified");
        static_assert(R3 < N, "invalid row index specified");
        static_assert(C1 < N, "invalid column index specified");
        static_assert(C2 < N, "invalid column index specified");
        static_assert(C3 < N, "invalid column index specified");
        template <unsigned int sR1, unsigned int sR2, unsigned int sC1, unsigned int sC2>
        using UpHelper = DeterminantHelper<T, N, sR1, sR2, sC1, sC2>;
        static T compute(T const* data)
        {
          return data[index(R1, C1)] * UpHelper<R2, R3, C2, C3>::compute(data) -
                 data[index(R1, C2)] * UpHelper<R2, R3, C1, C3>::compute(data) +
                 data[index(R1, C3)] * UpHelper<R2, R3, C1, C2>::compute(data);
        }
      }; // struct DeterminantHelper<>

      /// Determinant of a 4x4 submatrix
      template <typename T,
                unsigned int N,
                unsigned int R1,
                unsigned int R2,
                unsigned int R3,
                unsigned int R4,
                unsigned int C1,
                unsigned int C2,
                unsigned int C3,
                unsigned int C4>
      struct DeterminantHelper<T, N, R1, R2, R3, R4, C1, C2, C3, C4>
        : public DeterminantHelperBase<N> {
        using DeterminantHelperBase<N>::index;
        static_assert(R1 < N, "invalid row index specified");
        static_assert(R2 < N, "invalid row index specified");
        static_assert(R3 < N, "invalid row index specified");
        static_assert(R4 < N, "invalid row index specified");
        static_assert(C1 < N, "invalid column index specified");
        static_assert(C2 < N, "invalid column index specified");
        static_assert(C3 < N, "invalid column index specified");
        static_assert(C4 < N, "invalid column index specified");
        template <unsigned int sR1,
                  unsigned int sR2,
                  unsigned int sR3,
                  unsigned int sC1,
                  unsigned int sC2,
                  unsigned int sC3>
        using UpHelper = DeterminantHelper<T, N, sR1, sR2, sR3, sC1, sC2, sC3>;
        static T compute(T const* data)
        {
          return data[index(R1, C1)] * UpHelper<R2, R3, R4, C2, C3, C4>::compute(data) -
                 data[index(R1, C2)] * UpHelper<R2, R3, R4, C1, C3, C4>::compute(data) +
                 data[index(R1, C3)] * UpHelper<R2, R3, R4, C1, C2, C4>::compute(data) -
                 data[index(R1, C4)] * UpHelper<R2, R3, R4, C1, C2, C3>::compute(data);
        }
      }; // struct DeterminantHelper<>

      /// Routines for 2x2 matrices
      template <typename T>
      struct FastMatrixOperations<T, 2> : public FastMatrixOperationsBase<T, 2> {
        using Base_t = FastMatrixOperationsBase<T, 2>;
        static constexpr unsigned int Dim = Base_t::Dim;
        using Data_t = typename Base_t::Data_t;
        using Matrix_t = typename Base_t::Matrix_t;

        /// Computes the determinant of a matrix
        static Data_t Determinant(Matrix_t const& mat)
        {
          return DeterminantHelper<T, Dim, 0, 1, 0, 1>::compute(mat.data());
        }

        /// Computes the determinant of a matrix, using the provided determinant
        static Matrix_t InvertMatrix(Matrix_t const& mat, Data_t det);

        /// Computes the determinant of a symmatric matrix,
        /// using the provided determinant
        static Matrix_t InvertSymmetricMatrix(Matrix_t const& mat, Data_t det);

        /// Computes the determinant of a matrix
        static Matrix_t InvertMatrix(Matrix_t const& mat)
        {
          return InvertMatrix(mat, Determinant(mat));
        }

        /// Computes the determinant of a matrix
        static Matrix_t InvertSymmetricMatrix(Matrix_t const& mat)
        {
          return InvertSymmetricMatrix(mat, Determinant(mat));
        }

      }; // struct FastMatrixOperations<T, 2>

      /// Routines for 3x3 matrices
      template <typename T>
      struct FastMatrixOperations<T, 3> : public FastMatrixOperationsBase<T, 3> {
        using Base_t = FastMatrixOperationsBase<T, 3>;
        static constexpr unsigned int Dim = Base_t::Dim;
        using Data_t = typename Base_t::Data_t;
        using Matrix_t = typename Base_t::Matrix_t;

        /// Computes the determinant of a matrix
        static Data_t Determinant(Matrix_t const& mat)
        {
          return DeterminantHelper<T, Dim, 0, 1, 2, 0, 1, 2>::compute(mat.data());
        }

        /// Computes the determinant of a matrix, using the provided determinant
        static Matrix_t InvertMatrix(Matrix_t const& mat, Data_t det);

        /// Computes the determinant of a symmatric matrix,
        /// using the provided determinant
        static Matrix_t InvertSymmetricMatrix(Matrix_t const& mat, Data_t det);

        /// Computes the determinant of a matrix
        static Matrix_t InvertMatrix(Matrix_t const& mat)
        {
          return InvertMatrix(mat, Determinant(mat));
        }

        /// Computes the determinant of a matrix
        static Matrix_t InvertSymmetricMatrix(Matrix_t const& mat)
        {
          return InvertSymmetricMatrix(mat, Determinant(mat));
        }

      }; // struct FastMatrixOperations<T, 3>

      /// Routines for 4x4 matrices
      template <typename T>
      struct FastMatrixOperations<T, 4> : public FastMatrixOperationsBase<T, 4> {
        using Base_t = FastMatrixOperationsBase<T, 4>;
        static constexpr unsigned int Dim = Base_t::Dim;
        using Data_t = typename Base_t::Data_t;
        using Matrix_t = typename Base_t::Matrix_t;

        /// Computes the determinant of a matrix
        static Data_t Determinant(Matrix_t const& mat)
        {
          return DeterminantHelper<T, Dim, 0, 1, 2, 3, 0, 1, 2, 3>::compute(mat.data());
        }

        /// Computes the determinant of a matrix, using the provided determinant
        static Matrix_t InvertMatrix(Matrix_t const& mat, Data_t det);

        /// Computes the determinant of a symmatric matrix,
        /// using the provided determinant
        static Matrix_t InvertSymmetricMatrix(Matrix_t const& mat, Data_t det);

        /// Computes the determinant of a matrix
        static Matrix_t InvertMatrix(Matrix_t const& mat)
        {
          return InvertMatrix(mat, Determinant(mat));
        }

        /// Computes the determinant of a matrix
        static Matrix_t InvertSymmetricMatrix(Matrix_t const& mat)
        {
          return InvertSymmetricMatrix(mat, Determinant(mat));
        }

      }; // struct FastMatrixOperations<T, 3>

    } // namespace details
  }   // namespace util
} // namespace lar

template <typename T, unsigned int DIM>
auto lar::util::details::FastMatrixOperationsBase<T, DIM>::MatrixVectorProduct(Matrix_t const& mat,
                                                                               Vector_t const& vec)
  -> Vector_t
{
  // not really fast, but there is probably not much to fasten...
  Vector_t res;
  Data_t const* mat_row = mat.data();
  for (size_t r = 0; r < Dim; ++r) {
    Data_t elem = Data_t(0);
    for (size_t c = 0; c < Dim; ++c)
      elem += *(mat_row++) * vec[c];
    res[r] = elem;
  } // for
  return res;
} // FastMatrixOperationsBase<>::MatrixVectorProduct()

template <typename T>
auto lar::util::details::FastMatrixOperations<T, 2>::InvertMatrix(Matrix_t const& mat, Data_t det)
  -> Matrix_t
{
  Matrix_t Inverse;
  Inverse[0 * Dim + 0] = mat[1 * Dim + 1] / det;
  Inverse[0 * Dim + 1] = -mat[0 * Dim + 1] / det;
  Inverse[1 * Dim + 0] = -mat[1 * Dim + 0] / det;
  Inverse[1 * Dim + 1] = mat[0 * Dim + 0] / det;
  return Inverse;
} // FastMatrixOperations<T, 2>::InvertMatrix()

template <typename T>
auto lar::util::details::FastMatrixOperations<T, 2>::InvertSymmetricMatrix(Matrix_t const& mat,
                                                                           Data_t det) -> Matrix_t
{
  Matrix_t Inverse;
  Inverse[0 * Dim + 0] = mat[1 * Dim + 1] / det;
  Inverse[0 * Dim + 1] = Inverse[1 * Dim + 0] = -mat[1 * Dim + 0] / det;
  Inverse[1 * Dim + 1] = mat[0 * Dim + 0] / det;
  return Inverse;
} // FastMatrixOperations<T, 2>::InvertMatrix()

template <typename T>
auto lar::util::details::FastMatrixOperations<T, 3>::InvertMatrix(Matrix_t const& mat, Data_t det)
  -> Matrix_t
{
  Data_t const* data = mat.data();
  Matrix_t Inverse;
  //
  // Basically using Cramer's rule;
  // each element [r,c] gets assigned the determinant of the submatrix
  //   after removing c from the rows and r from the columns
  //   (effectively assigning the transpose of the minor matrix)
  //   with the usual sign -1^(r+c)
  //
  //
  Inverse[0 * Dim + 0] = DeterminantHelper<T, 3, 1, 2, 1, 2>::compute(data) / det;
  Inverse[0 * Dim + 1] = -DeterminantHelper<T, 3, 0, 2, 1, 2>::compute(data) / det;
  Inverse[0 * Dim + 2] = DeterminantHelper<T, 3, 0, 1, 1, 2>::compute(data) / det;
  Inverse[1 * Dim + 0] = -DeterminantHelper<T, 3, 1, 2, 0, 2>::compute(data) / det;
  Inverse[1 * Dim + 1] = DeterminantHelper<T, 3, 0, 2, 0, 2>::compute(data) / det;
  Inverse[1 * Dim + 2] = -DeterminantHelper<T, 3, 0, 1, 0, 2>::compute(data) / det;
  Inverse[2 * Dim + 0] = DeterminantHelper<T, 3, 1, 2, 0, 1>::compute(data) / det;
  Inverse[2 * Dim + 1] = -DeterminantHelper<T, 3, 0, 2, 0, 1>::compute(data) / det;
  Inverse[2 * Dim + 2] = DeterminantHelper<T, 3, 0, 1, 0, 1>::compute(data) / det;
  return Inverse;
} // FastMatrixOperations<T, 3>::InvertMatrix()

template <typename T>
auto lar::util::details::FastMatrixOperations<T, 3>::InvertSymmetricMatrix(Matrix_t const& mat,
                                                                           Data_t det) -> Matrix_t
{
  //
  // Same algorithm as InvertMatrix(), but use the fact that the result is
  // also symmetric
  //
  Data_t const* data = mat.data();
  Matrix_t Inverse;
  Inverse[0 * Dim + 0] = DeterminantHelper<T, 3, 1, 2, 1, 2>::compute(data) / det;
  Inverse[0 * Dim + 1] = Inverse[1 * Dim + 0] =
    -DeterminantHelper<T, 3, 1, 2, 0, 2>::compute(data) / det;
  Inverse[0 * Dim + 2] = Inverse[2 * Dim + 0] =
    DeterminantHelper<T, 3, 1, 2, 0, 1>::compute(data) / det;
  Inverse[1 * Dim + 1] = DeterminantHelper<T, 3, 0, 2, 0, 2>::compute(data) / det;
  Inverse[1 * Dim + 2] = Inverse[2 * Dim + 1] =
    -DeterminantHelper<T, 3, 0, 2, 0, 1>::compute(data) / det;
  Inverse[2 * Dim + 2] = DeterminantHelper<T, 3, 0, 1, 0, 1>::compute(data) / det;
  return Inverse;
} // FastMatrixOperations<T, 3>::InvertSymmetricMatrix()

template <typename T>
auto lar::util::details::FastMatrixOperations<T, 4>::InvertMatrix(Matrix_t const& mat, Data_t det)
  -> Matrix_t
{
  //
  // Basically using Cramer's rule;
  // each element [r,c] gets assigned the determinant of the submatrix
  //   after removing c from the rows and r from the columns
  //   (effectively assigning the transpose of the minor matrix)
  //   with the usual sign -1^(r+c)
  //
  //
  Data_t const* data = mat.data();
  Matrix_t Inverse;
  Inverse[0 * Dim + 0] = DeterminantHelper<T, 4, 1, 2, 3, 1, 2, 3>::compute(data) / det;
  Inverse[0 * Dim + 1] = -DeterminantHelper<T, 4, 0, 2, 3, 1, 2, 3>::compute(data) / det;
  Inverse[0 * Dim + 2] = DeterminantHelper<T, 4, 0, 1, 3, 1, 2, 3>::compute(data) / det;
  Inverse[0 * Dim + 3] = -DeterminantHelper<T, 4, 0, 1, 2, 1, 2, 3>::compute(data) / det;
  Inverse[1 * Dim + 0] = -DeterminantHelper<T, 4, 1, 2, 3, 0, 2, 3>::compute(data) / det;
  Inverse[1 * Dim + 1] = DeterminantHelper<T, 4, 0, 2, 3, 0, 2, 3>::compute(data) / det;
  Inverse[1 * Dim + 2] = -DeterminantHelper<T, 4, 0, 1, 3, 0, 2, 3>::compute(data) / det;
  Inverse[1 * Dim + 3] = DeterminantHelper<T, 4, 0, 1, 2, 0, 2, 3>::compute(data) / det;
  Inverse[2 * Dim + 0] = DeterminantHelper<T, 4, 1, 2, 3, 0, 1, 3>::compute(data) / det;
  Inverse[2 * Dim + 1] = -DeterminantHelper<T, 4, 0, 2, 3, 0, 1, 3>::compute(data) / det;
  Inverse[2 * Dim + 2] = DeterminantHelper<T, 4, 0, 1, 3, 0, 1, 3>::compute(data) / det;
  Inverse[2 * Dim + 3] = -DeterminantHelper<T, 4, 0, 1, 2, 0, 1, 3>::compute(data) / det;
  Inverse[3 * Dim + 0] = -DeterminantHelper<T, 4, 1, 2, 3, 0, 1, 2>::compute(data) / det;
  Inverse[3 * Dim + 1] = DeterminantHelper<T, 4, 0, 2, 3, 0, 1, 2>::compute(data) / det;
  Inverse[3 * Dim + 2] = -DeterminantHelper<T, 4, 0, 1, 3, 0, 1, 2>::compute(data) / det;
  Inverse[3 * Dim + 3] = DeterminantHelper<T, 4, 0, 1, 2, 0, 1, 2>::compute(data) / det;
  return Inverse;
} // FastMatrixOperations<T, 4>::InvertMatrix()

template <typename T>
auto lar::util::details::FastMatrixOperations<T, 4>::InvertSymmetricMatrix(Matrix_t const& mat,
                                                                           Data_t det) -> Matrix_t
{
  //
  // Same algorithm as InvertMatrix(), but use the fact that the result is
  // also symmetric
  //
  Data_t const* data = mat.data();
  Matrix_t Inverse;
  Inverse[0 * Dim + 0] = DeterminantHelper<T, 4, 1, 2, 3, 1, 2, 3>::compute(data) / det;
  Inverse[0 * Dim + 1] = Inverse[1 * Dim + 0] =
    -DeterminantHelper<T, 4, 1, 2, 3, 0, 2, 3>::compute(data) / det;
  Inverse[0 * Dim + 2] = Inverse[2 * Dim + 0] =
    DeterminantHelper<T, 4, 1, 2, 3, 0, 1, 3>::compute(data) / det;
  Inverse[0 * Dim + 3] = Inverse[3 * Dim + 0] =
    -DeterminantHelper<T, 4, 1, 2, 3, 0, 1, 2>::compute(data) / det;
  Inverse[1 * Dim + 1] = DeterminantHelper<T, 4, 0, 2, 3, 0, 2, 3>::compute(data) / det;
  Inverse[1 * Dim + 2] = Inverse[2 * Dim + 1] =
    -DeterminantHelper<T, 4, 0, 2, 3, 3, 0, 1>::compute(data) / det;
  Inverse[1 * Dim + 3] = Inverse[3 * Dim + 1] =
    DeterminantHelper<T, 4, 0, 2, 3, 0, 1, 2>::compute(data) / det;
  Inverse[2 * Dim + 2] = DeterminantHelper<T, 4, 0, 1, 3, 0, 1, 3>::compute(data) / det;
  Inverse[2 * Dim + 3] = Inverse[3 * Dim + 2] =
    -DeterminantHelper<T, 4, 0, 1, 3, 0, 1, 2>::compute(data) / det;
  Inverse[3 * Dim + 3] = DeterminantHelper<T, 4, 0, 1, 2, 0, 1, 2>::compute(data) / det;
  return Inverse;
} // FastMatrixOperations<T, 4>::InvertSymmetricMatrix()

#endif // FASTMATRIXMATHHELPER_H
