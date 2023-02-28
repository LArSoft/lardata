/**
 * @file   GridContainerIndices.h
 * @brief  Classes to manage containers with indices in 1, 2 and 3 dimensions
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 29, 2016
 *
 * This header provides:
 *
 * * GridContainer2DIndices: index manager for object in a 2D space
 * * GridContainer3DIndices: index manager for object in a 3D space
 *
 * These classes have methods whose names reflect the idea of a physical space
 * ("x", "y", "z"). The functionality is provided by TensorIndices class.
 *
 * This is a pure header that contains only template classes.
 */

#ifndef LARDATA_UTILITIES_GRIDCONTAINERINDICES_H
#define LARDATA_UTILITIES_GRIDCONTAINERINDICES_H

// LArSoft libraries
#include "lardata/Utilities/TensorIndices.h"

// C/C++ standard libraries
#include <array>
#include <cstddef> // std::ptrdiff_t

namespace util {

  namespace details {

    /// Index manager for a container of data arranged on a DIMS-dimension grid
    template <unsigned int DIMS>
    class GridContainerIndicesBase {
      using IndexManager_t = util::TensorIndices<DIMS>;

    public:
      /// Returns the number of dimensions in this object
      static constexpr unsigned int dims() { return DIMS; }

      /// type of index for direct access to the cell
      using CellIndex_t = typename IndexManager_t::LinIndex_t;

      /// type of difference between indices
      using CellIndexOffset_t = std::ptrdiff_t;

      /// type of difference between indices along a dimension
      using CellDimIndex_t = CellIndexOffset_t;

      /// type of cell coordinate (x, y, z)
      using CellID_t = std::array<CellDimIndex_t, dims()>;

      /// Constructor: specifies the size of the container and allocates it
      GridContainerIndicesBase(std::array<size_t, dims()> const& new_dims)
        : indices(new_dims.begin())
      {}

      /// @{
      /// @name Grid structure

      /// Returns whether the specified index is valid
      bool has(CellIndexOffset_t index) const { return indices.hasLinIndex(index); }

      /// Returns the number of cells in the grid
      size_t size() const { return indices.size(); }

      /// @}

      /// @{
      /// @name Indexing

      /// Returns the index of the element from its cell coordinates (no check!)
      CellIndex_t operator[](CellID_t id) const { return index(id); }

      /// Returns the difference in index of cellID respect to origin
      CellIndexOffset_t offset(CellID_t const& origin, CellID_t const& cellID) const
      {
        return index(cellID) - index(origin);
      }

      /// @}

    protected:
      IndexManager_t indices; ///< the actual worker

      /// Returns the index of the element from its cell coordinates (no check!)
      CellIndex_t index(CellID_t id) const { return indices(id.begin()); }

    }; // GridContainerIndicesBase

  } // namespace details

  /// Index manager for a container of data arranged on a >=1-dim grid
  template <unsigned int DIMS = 1U>
  class GridContainerIndicesBase1D : public details::GridContainerIndicesBase<DIMS> {
    static_assert(DIMS >= 1U, "Dimensions for GridContainerIndicesBase1D must be at least 1");
    using Base_t = details::GridContainerIndicesBase<DIMS>;

  public:
    using Base_t::GridContainerIndicesBase;

    /// @{
    /// @name Grid structure

    /// Returns whether the specified x index is valid
    bool hasX(typename Base_t::CellDimIndex_t index) const
    {
      return Base_t::indices.template hasIndex<0>(index);
    }

    /// Returns the number of cells on the x axis of the grid
    size_t sizeX() const { return Base_t::indices.template dim<0>(); }

    /// @}

  }; // GridContainerIndicesBase1D

  /// Index manager for a container of data arranged on a >=2-dim grid
  template <unsigned int DIMS = 2U>
  class GridContainerIndicesBase2D : public GridContainerIndicesBase1D<DIMS> {
    static_assert(DIMS >= 2U, "Dimensions for GridContainerIndicesBase2D must be at least 2");

    using Base_t = GridContainerIndicesBase1D<DIMS>;

  public:
    using Base_t::GridContainerIndicesBase1D;

    /// @{
    /// @name Grid structure

    /// Returns whether the specified y index is valid
    bool hasY(typename Base_t::CellDimIndex_t index) const
    {
      return Base_t::indices.template hasIndex<1>(index);
    }

    /// Returns the number of cells on the y axis of the grid
    size_t sizeY() const { return Base_t::indices.template dim<1>(); }

    /// @}

  }; // GridContainerIndicesBase2D

  /// Index manager for a container of data arranged on a >=3-dim grid
  template <unsigned int DIMS = 3U>
  class GridContainerIndicesBase3D : public GridContainerIndicesBase2D<DIMS> {
    static_assert(DIMS >= 3U, "Dimensions for GridContainerIndicesBase3D must be at least 3");
    using Base_t = GridContainerIndicesBase2D<DIMS>;

  public:
    using Base_t::GridContainerIndicesBase2D;

    /// @{
    /// @name Grid structure

    /// Returns whether the specified z index is valid
    bool hasZ(typename Base_t::CellDimIndex_t index) const
    {
      return Base_t::indices.template hasIndex<2>(index);
    }

    /// Returns the number of cells on the z axis of the grid
    size_t sizeZ() const { return Base_t::indices.template dim<2>(); }

    /// @}

  }; // GridContainerIndicesBase3D

  /// Index manager for a container of data arranged on a 2D grid
  using GridContainer2DIndices = GridContainerIndicesBase2D<>;

  /// Index manager for a container of data arranged on a 3D grid
  using GridContainer3DIndices = GridContainerIndicesBase3D<>;

} // namespace util

#endif // LARDATA_UTILITIES_GRIDCONTAINERINDICES_H
