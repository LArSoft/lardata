/**
 * @file   TensorIndices.h
 * @brief  TensorIndices class to flatten multi-dimension indices into linear
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 28, 2016
 *
 * This header provides:
 *
 * * util::TensorIndices: general template
 * * util::MatrixIndices: alias for 2D tensors
 *
 * This is a pure header that contains only template classes.
 */

#ifndef LARDATA_UTILITIES_TENSORINDICES_H
#define LARDATA_UTILITIES_TENSORINDICES_H

// C/C++ standard libraries
#include <cstddef>     // std::size_t
#include <stdexcept>   // std::out_of_range
#include <string>      // std::to_string()
#include <type_traits> // std::enable_if

namespace util {

  /// Types for TensorIndices class
  struct TensorIndicesBasicTypes {

    /// Type of a single index in the tensor
    using Index_t = std::ptrdiff_t;

    /// Type of size of a dimension in the tensor
    using DimSize_t = std::size_t;

    /// Type of the linear index
    using LinIndex_t = std::size_t;

  }; // TensorIndicesBasicTypes

  // This is the class declaration; it has a specialisation for rank 1 and
  // a generic implementation (valid for rank larger than 1).
  template <unsigned int RANK>
  class TensorIndices;

  namespace details {
    template <unsigned int RANK, unsigned int DIM>
    struct ExtractTensorDimension;
  } // namespace details

  /**
   * Specialisation for vector index (rank 1).
   *
   * Documentation is in the general template.
   * No public methods are added respect to it (but the `minor()` method is
   * removed).
   */
  template <>
  class TensorIndices<1U> {
  public:
    using Index_t = TensorIndicesBasicTypes::Index_t;

    using DimSize_t = TensorIndicesBasicTypes::DimSize_t;

    using LinIndex_t = TensorIndicesBasicTypes::LinIndex_t;

    static constexpr unsigned int rank() { return 1U; }

    TensorIndices(DimSize_t dim) : dimSize(dim) {}

    template <
      typename ITER,
      typename = std::enable_if_t<std::is_convertible<decltype(*(ITER())), DimSize_t>::value, void>>
    TensorIndices(ITER dimIter) : TensorIndices(*dimIter)
    {}

    LinIndex_t operator()(Index_t index) const { return index; }

    template <typename ITER>
    std::enable_if_t<std::is_convertible<decltype(*(ITER())), DimSize_t>::value, LinIndex_t>
    operator()(ITER indexIter) const
    {
      return *indexIter;
    }

    LinIndex_t at(Index_t index) const { return checkOuterIndex(index); }

    template <typename ITER>
    std::enable_if_t<std::is_convertible<decltype(*(ITER())), DimSize_t>::value, LinIndex_t> at(
      ITER indexIter) const
    {
      return checkOuterIndex(*indexIter);
    }

    bool has(Index_t index) const { return (index >= 0) && ((DimSize_t)index < size()); }

    template <typename ITER>
    std::enable_if_t<std::is_convertible<decltype(*(ITER())), DimSize_t>::value, bool> has(
      ITER indexIter) const
    {
      return (*indexIter >= 0) && ((DimSize_t)*indexIter < size());
    }

    template <unsigned int DIM = 0>
    DimSize_t size() const
    {
      return (DIM == 0) ? totalSize() : dim<DIM>();
    }

    template <unsigned int DIM>
    DimSize_t dim() const
    {
      static_assert(DIM == 0, "Invalid dimension requested");
      return dimSize;
    } // dim()

    template <unsigned int DIM>
    bool hasIndex(Index_t index) const
    {
      static_assert(DIM == 0, "Invalid dimension requested");
      return has(index);
    }

    bool hasLinIndex(LinIndex_t linIndex) const { return has((Index_t)linIndex); }

    bool operator==(TensorIndices<1U> const& t) const { return t.size() == size(); }
    bool operator!=(TensorIndices<1U> const& t) const { return t.size() != size(); }

  protected:
    Index_t checkOuterIndex(Index_t index) const
    {
      if (has(index)) return index; // good
      throw std::out_of_range("Requested index " + std::to_string(index) +
                              " for a dimension of size " + std::to_string(size()));
    }

  protected:
    /// Returns the size of the outer dimension
    DimSize_t dim0() const { return dimSize; }

  private:
    template <unsigned int R, unsigned int D>
    friend struct details::ExtractTensorDimension;

    DimSize_t dimSize; ///< size of the largest dimension

    /// Returns the total size of this tensor
    DimSize_t totalSize() const { return dim0(); }

  }; // class TensorIndices<1>

  /**
   * @brief Converts a tensor element specification into a linear index
   *
   * Example to use a 6 x 3 x 2 x 4 tensor of rank 4:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * TensorIndices<4> indices(6, 3, 2, 4);
   *
   * std::vector<double> v(indices.size(), 0.);
   * std::cout << v[indices(4, 1, 1, 0)] << std::endl;
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * This will map the content of the vector `v` as a tensor of rank 4.
   *
   */
  template <unsigned int RANK>
  class TensorIndices : private TensorIndices<1U> {
    static_assert(RANK > 1, "TensorIndices must have rank 1 or higher");

    using Base_t = TensorIndices<1U>; ///< type of base class

  public:
    /// Type of a single index in the tensor
    using Index_t = Base_t::Index_t;

    /// Type for the specification of a dimension size
    using DimSize_t = Base_t::DimSize_t;

    /// Type of the linear index
    using LinIndex_t = Base_t::LinIndex_t;

    /// Rank of this tensor
    static constexpr unsigned int rank() { return RANK; }

    /// Type of the tensor indices with rank smaller by one
    using MinorTensor_t = TensorIndices<rank() - 1>;

    /**
     * @brief Constructor: initialises the dimension of the tensor
     * @tparam OTHERINDICES types of the other index values
     * @param first the size of the first dimension
     * @param others size of each of the other dimensions
     *
     * Example to initialise indices for a 6 x 3 x 2 x 4 tensor of rank 4:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * TensorIndices<4> indices(6, 3, 2, 4);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     */
    template <typename... OTHERDIMS>
    TensorIndices(DimSize_t first, OTHERDIMS... others)
      : Base_t(first), m(others...), totSize(dim0() * m.size())
    {}

    /**
     * @brief Constructor: initialises the dimension of the tensor
     * @tparam ITER type of iterator to dimension sizes
     * @param dimIter iterator pointing to the first dimension size
     *
     * Dimensions are initialised from the values pointed by the specified
     * iterator. The iterator is required to be a forward iterator pointing to
     * a type convertible to a `DimSize_t` type.
     * The iterator must be valid when increased up to `rank() - 1` times, so
     * that all the `rank()` dimensions can be extracted.
     *
     * Example to initialise indices for a 6 x 3 x 2 x 4 tensor of rank 4:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * std::array<size_t, 4> dims { 6, 3, 2, 4 };
     * TensorIndices<4> indices(dims.begin());
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Note that no `end` iterator is required.
     *
     */
    // only if ITER can be dereferenced into something convertible to DimSize_t
    template <
      typename ITER,
      typename = std::enable_if_t<std::is_convertible<decltype(*(ITER())), DimSize_t>::value, void>>
    TensorIndices(ITER dimIter) : Base_t(dimIter), m(++dimIter), totSize(dim0() * m.size())
    {}

    /**
     * @brief Returns the linear index corresponding to the tensor indices
     * @tparam OTHERINDICES types of the other index values
     * @param first the index of the first dimension
     * @param others the indices of the other dimensions
     * @return the offset of the specified element with respect to the first one
     * @see at(), has()
     *
     * No check is performed on the validity of the indices (that is, whether
     * they fall within the dimension range).
     *
     * Use it as:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * TensorIndices<4> indices(6, 3, 2, 4);
     * auto validIndex = indices(4, 2, 0, 1);
     * auto invalidIndex = indices(4, 3, 2, 1); // invalid index is returned!
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * (in the example, note that the last line has two indices out of range and
     * the return value will not be valid!).
     */
    template <typename... OTHERINDICES>
    LinIndex_t operator()(Index_t first, OTHERINDICES... others) const
    {
      return first * minorTensor().size() + minorTensor()(others...);
    }

    /**
     * @brief Returns the linear index corresponding to the tensor indices
     * @tparam ITER type of iterator to dimension sizes
     * @param indexIter iterator pointing to the first dimension size
     * @return the offset of the specified element with respect to the first one
     * @see at(), has()
     *
     * No check is performed on the validity of the indices (that is, whether
     * they fall within the dimension range).
     *
     * Use it as:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * TensorIndices<4> indices(6, 3, 2, 4);
     * std::array<int, 4> ix { 4, 2, 0, 1 };
     * auto validIndex = indices(ix.begin());
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The behaviour is the same as in
     * `operator()(Index_t first, OTHERINDICES... others)`.
     */
    template <typename ITER>
    std::enable_if_t<std::is_convertible<decltype(*(ITER())), DimSize_t>::value, LinIndex_t>
    operator()(ITER indexIter) const
    {
      auto const baseSize = (*indexIter) * minorTensor().size();
      return baseSize + minorTensor()(++indexIter);
    }

    /**
     * @brief Returns the linear index corresponding to the tensor indices
     * @tparam OTHERINDICES types of the other index values
     * @param first the index of the first dimension
     * @param others the indices of the other dimensions
     * @throw std::out_of_range if any of the indices is not valid
     * @return the offset of the specified element with respect to the first one
     * @see operator()(), has()
     *
     * If any of the index values is not in the valid index range, an exception
     * is thrown.
     *
     * Use it as:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * TensorIndices<4> indices(6, 3, 2, 4);
     * auto validIndex = indices.at(4, 2, 0, 1);
     * auto invalidIndex = indices.at(4, 3, 2, 1); // throws std::out_of_range
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * (in the example, note that the last line has two indices out of range and
     * an exception will be thrown).
     *
     */
    template <typename... OTHERINDICES>
    LinIndex_t at(Index_t first, OTHERINDICES... others) const
    {
      return Base_t::checkOuterIndex(first) * minorTensor().size() + minorTensor().at(others...);
    }

    /**
     * @brief Returns the linear index corresponding to the tensor indices
     * @tparam ITER type of iterator to dimension sizes
     * @param indexIter iterator pointing to the first dimension size
     * @throw std::out_of_range if any of the indices is not valid
     * @return the offset of the specified element with respect to the first one
     * @see operator()(), has()
     *
     * If any of the index values is not in the valid index range, an exception
     * is thrown.
     *
     * Use it as:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * TensorIndices<4> indices(6, 3, 2, 4);
     * std::array<int, 4> ix { 4, 2, 0, 1 };
     * auto validIndex = indices.at(ix.begin());
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The behaviour is the same as in
     * `at()(Index_t first, OTHERINDICES... others)`.
     *
     */
    template <typename ITER>
    std::enable_if_t<std::is_convertible<decltype(*(ITER())), DimSize_t>::value, LinIndex_t> at(
      ITER indexIter) const
    {
      auto const baseSize = Base_t::checkOuterIndex(*indexIter) * minorTensor().size();
      return baseSize + minorTensor()(++indexIter);
    }

    /**
     * @brief Returns whether the specified set of indices is valid
     * @tparam OTHERINDICES types of the other index values
     * @param first the index of the first dimension
     * @param others the indices of the other dimensions
     * @return whether the specified set of indices is valid
     * @see operator(), at()
     *
     * If any of the index values is not in the valid index range, an exception
     * is thrown.
     *
     * Use it as:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * TensorIndices<4> indices(6, 3, 2, 4);
     * if (indices.has(4, 2, 0, 1)) {
     *   // ...
     * }
     * if (indices.has(4, 3, 2, 1)) {
     *   // ... (won't be executed)
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * (in the example, note that the last line has two indices out of range and
     * an exception will be thrown).
     */
    template <typename... OTHERINDICES>
    bool has(Index_t first, OTHERINDICES... others) const
    {
      return Base_t::has(first) && minorTensor().has(others...);
    }

    /**
     * @brief Returns whether the specified set of indices is valid
     * @tparam ITER type of iterator to dimension sizes
     * @param indexIter iterator pointing to the first dimension size
     * @return whether the specified set of indices is valid
     * @see operator(), at()
     *
     * If any of the index values is not in the valid index range, an exception
     * is thrown.
     *
     * Use it as:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * TensorIndices<4> indices(6, 3, 2, 4);
     * std::array<int, 4> ix { 4, 2, 0, 1 };
     * if (indices.has(ix.begin())) {
     *   // ...
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * The behaviour is the same as in
     * `has()(Index_t first, OTHERINDICES... others)`.

     */
    template <typename ITER>
    std::enable_if_t<std::is_convertible<decltype(*(ITER())), DimSize_t>::value, bool> has(
      ITER indexIter) const
    {
      return Base_t::has(*indexIter) ? minorTensor().has(++indexIter) : false;
    }

    /**
     * @brief Returns the size of the specified dimension
     * @tparam DIM the dimension (0 is the first one)
     * @return the size of the specified dimension
     *
     * Code won't compile for DIM larger than the tensor rank.
     */
    template <unsigned int DIM>
    DimSize_t dim() const
    {
      return details::ExtractTensorDimension<rank(), DIM>::dim(*this);
    }

    /**
     * @brief Returns whether a index is valid within a specified dimension
     * @tparam DIM the dimension (0 is the first one)
     * @param index the value of the index to be tested
     * @return whether a index is valid within a specified dimension
     *
     * Use it as:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * TensorIndices<4> indices(6, 3, 2, 4);
     * if (indices.hasIndex<0>(2)) { // dimension #0 has 6 indices, [ 0 ; 5 ]
     *   // ... (will be executed)
     * }
     * if (indices.hasIndex<0>(6)) {
     *   // ... (will not be executed)
     * }
     * if (indices.hasIndex<2>(0)) { // dimension #2 has 2 indices, { 0 , 1 }
     *   // ... (will be executed)
     * }
     * if (indices.hasIndex<3>(6)) {// dimension #3 has 4 indices, [ 0 ; 3 ]
     *   // ... (will not be executed)
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     */
    template <unsigned int DIM>
    bool hasIndex(Index_t index) const
    {
      return (index >= 0U) && ((DimSize_t)index < dim<DIM>());
    }

    /**
     * @brief Returns the size of the minor tensor
     * @tparam DIM the index of the outer dimension of the minor tensor
     * @return the size of the specified tensor
     *
     * The total size of the tensor obtained after stripping the DIM most outer
     * dimensions if returned. Therefore, `size<0>()` is the size of the full
     * tensor, while `size<1>()` is the size of the tensor obtained by stripping
     * the first dimension (that is, the tensor returned by `minor()`), and so
     * on.
     *
     * Code won't compile for DIM larger than the tensor rank.
     */
    template <unsigned int DIM = 0>
    DimSize_t size() const
    {
      return (DIM == 0) ? totalSize() : details::ExtractTensorDimension<rank(), DIM>::size(*this);
    }

    /// Returns whether the specified linear index is valid in this tensor
    bool hasLinIndex(LinIndex_t linIndex) const
    {
      return (linIndex >= 0U) && ((DimSize_t)linIndex < size());
    }

    /// Returns the tensor of rank `Rank-1` from stripping the first dimension
    /// Note that the minorTensor variable was previously named "minor".
    /// However, minor is also defined in sys/types.h to mean something completely different.
    MinorTensor_t const& minorTensor() const { return m; }

    /// Returns whether all sizes of the tensor t are the same as this one
    bool operator==(TensorIndices<RANK> const& t) const
    {
      return Base_t::operator==(t) && (minorTensor() == t.minorTensor());
    }

    /// Returns whether any size of the tensor t is different from this one
    bool operator!=(TensorIndices<RANK> const& t) const
    {
      return Base_t::operator!=(t) || (minorTensor() != t.minorTensor());
    }

  protected:
    // need to be friend of the dimension extractor
    template <unsigned int R, unsigned int D>
    friend struct details::ExtractTensorDimension;

    MinorTensor_t m;   ///< the rest of the tensor indices
    DimSize_t totSize; ///< size of this tensor

    /// Returns the total size of this tensor (the same as `size()`)
    DimSize_t totalSize() const { return totSize; }

  }; // class TensorIndices<>

  /// Comparison operator with tensors of different rank
  template <unsigned int RANK1,
            unsigned int RANK2,
            typename = std::enable_if_t<(RANK1 != RANK2), bool>>
  bool operator==(TensorIndices<RANK1> const& a, TensorIndices<RANK2> const& b)
  {
    return false;
  }

  /// Comparison operator with tensors of different rank
  template <unsigned int RANK1,
            unsigned int RANK2,
            typename = std::enable_if_t<(RANK1 != RANK2), bool>>
  bool operator!=(TensorIndices<RANK1> const& a, TensorIndices<RANK2> const& b)
  {
    return true;
  }

  /**
   * @brief Instantiates a TensorIndices class with the specified dimensions
   * @tparam DIMS types for each of the arguments
   * @param dims size of each of the dimensions
   * @return a TensorIndices object with properly initialised dimensions
   *
   * The rank of the tensor is determined by the number of arguments; example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * auto indices = util::makeTensorIndices(3, 4);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will initialise a `TensorIndices<2>` (that's matrix indices), for a 3 x 4
   * (3 rows, 4 columns) disposition.
   */
  template <typename... DIMS>
  auto makeTensorIndices(DIMS... dims)
  {
    return TensorIndices<sizeof...(DIMS)>{TensorIndicesBasicTypes::DimSize_t(dims)...};
  }

  /// Type for indexing a 2D-tensor (matrix)
  using MatrixIndices = TensorIndices<2U>;

} // namespace util

//------------------------------------------------------------------------------
//--- details implementation
//---
namespace util {
  namespace details {

    /**
     * @tparam RANK rank of the TensorIndices to be queried
     * @tparam DIM the dimension to ask about (0 is the most outer one)
     *
     * Helper to recourse tensor dimensions with C++ metaprogramming.
     * Currently supported:
     *
     * * `dim(t)`: if `DIM` is 0, the request is for the outer dimension of `t`,
     *   which is immediately returned; otherwise, another extractor is queried
     *   on the `minor()` of `t`, decreasing `DIM` by 1
     *
     * * `size(t)`: if `DIM` is 0, the request is for the outer dimension of `t`,
     *   which is immediately returned; otherwise, another extractor is queried
     *   on the `minor()` of `t`, decreasing `DIM` by 1
     *
     */
    template <unsigned int RANK, unsigned int DIM>
    struct ExtractTensorDimension {
      static_assert(RANK > DIM, "Invalid dimension requested");
      static TensorIndicesBasicTypes::DimSize_t dim(TensorIndices<RANK> const& t)
      {
        return ExtractTensorDimension<RANK - 1U, DIM - 1U>::dim(t.minorTensor());
      }

      static TensorIndicesBasicTypes::DimSize_t size(TensorIndices<RANK> const& t)
      {
        return ExtractTensorDimension<RANK - 1U, DIM - 1U>::size(t.minorTensor());
      }

    }; // ExtractTensorDimension<>()

    template <unsigned int RANK>
    struct ExtractTensorDimension<RANK, 0U> {
      static_assert(RANK > 0, "Invalid rank 0 for TensorIndices");
      static TensorIndicesBasicTypes::DimSize_t dim(TensorIndices<RANK> const& t)
      {
        return t.dim0();
      }

      static TensorIndicesBasicTypes::DimSize_t size(TensorIndices<RANK> const& t)
      {
        return t.size();
      }

    }; // ExtractTensorDimension<RANK, 0>()

  } // namespace details
} // namespace util

#endif // LARDATA_UTILITIES_TENSORINDICES_H
