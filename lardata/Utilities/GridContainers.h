/**
 * @file   GridContainers.h
 * @brief  Containers with indices in 1, 2 and 3 dimensions
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 29, 2016
 *
 * This header provides:
 *
 * * GridContainer2D: container on data in 2D space
 * * GridContainer3D: container of data in 3D space
 * * GridContainerBase: base class for containers in a N-dimension space
 *
 * This is a pure header that contains only template classes.
 */

#ifndef LARDATA_UTILITIES_GRIDCONTAINERS_H
#define LARDATA_UTILITIES_GRIDCONTAINERS_H

// LArSoft libraries
#include "lardata/Utilities/GridContainerIndices.h"

// C/C++ standard libraries
#include <array>
#include <vector>

namespace util {

  namespace details {

    /**
     * @brief Base class for a container of data arranged on a grid
     * @tparam DATUM type of datum to be contained
     * @tparam IXMAN type of the grid index manager
     *
     * This is base class for `GridContainer#D` classes.
     * It provides the full functionality, to which the other classes add
     * some dimension-specific interface.
     *
     */
    template <typename DATUM, typename IXMAN>
    class GridContainerBase {

    public:
      using Datum_t = DATUM;   ///< type of contained datum
      using Indexer_t = IXMAN; /// type of index manager

      using Grid_t = GridContainerBase<Datum_t, Indexer_t>; ///< this type

      static constexpr unsigned int dims() { return IXMAN::dims(); }

      /// type of index for direct access to the cell
      using CellIndex_t = typename Indexer_t::CellIndex_t;

      /// type of difference between indices
      using CellIndexOffset_t = typename Indexer_t::CellIndexOffset_t;

      /// type of difference between indices
      using CellDimIndex_t = typename Indexer_t::CellDimIndex_t;

      /// type of cell coordinate (x, y, z)
      using CellID_t = typename Indexer_t::CellID_t;

      /// type of a single cell container
      using Cell_t = std::vector<Datum_t>;

      /// type of container holding all cells
      using Cells_t = std::vector<Cell_t>;

      /// type of iterator to all cells
      using const_iterator = typename Cells_t::const_iterator;

      /// Constructor: specifies the size of the container and allocates it
      GridContainerBase(std::array<size_t, dims()> const& dims)
        : indices(dims), data(indices.size())
      {}

      /// @{
      /// @name Data structure

      /// Returns the total size of the container
      size_t size() const { return indices.size(); }

      /// Returns whether the specified index is valid
      bool has(CellIndexOffset_t index) const { return indices.has(index); }

      /// @}

      /// @{
      /// @name Data access

      /// Return the index of the element from its cell coordinates (no check!)
      CellIndex_t index(CellID_t const& id) const { return indices[id]; }

      /// Returns the difference in index from two cells
      CellIndexOffset_t indexOffset(CellID_t const& origin, CellID_t const& cellID) const
      {
        return indices.offset(origin, cellID);
      }

      /// Returns a reference to the specified cell
      Cell_t& operator[](CellID_t const& id) { return cell(id); }

      /// Returns a constant reference to the specified cell
      Cell_t const& operator[](CellID_t const& id) const { return cell(id); }

      /// Returns a reference to to the cell with specified index
      Cell_t& operator[](CellIndex_t index) { return data[index]; }

      /// Returns a constant reference to the cell with specified index
      Cell_t const& operator[](CellIndex_t index) const { return data[index]; }

      ///@}

      /// @{
      /// @name Data insertion

      /// Copies an element into the specified cell
      void insert(CellID_t const& cellID, Datum_t const& elem) { cell(cellID).push_back(elem); }

      /// Moves an element into the specified cell
      void insert(CellID_t const& cellID, Datum_t&& elem)
      {
        cell(cellID).push_back(std::move(elem));
      }

      /// Copies an element into the cell with the specified index
      void insert(CellIndex_t index, Datum_t const& elem) { data[index].push_back(elem); }

      /// Moves an element into the cell with the specified index
      void insert(CellIndex_t index, Datum_t&& elem) { data[index].push_back(std::move(elem)); }

      /// @}

      /// Returns the index manager of the grid
      Indexer_t const& indexManager() const { return indices; }

    protected:
      Indexer_t indices; ///< manager of the indices of the container

      Cells_t data; ///< organised collection of points

      /// Returns a reference to the specified cell
      Cell_t& cell(CellID_t const& cellID) { return data[index(cellID)]; }

      /// Returns a constant reference to the specified cell
      Cell_t const& cell(CellID_t const& cellID) const { return data[index(cellID)]; }

    }; // GridContainerBase<>

  } // namespace details

  /**
   * @brief Base class for a container of data arranged on a 1D-grid
   * @tparam DATUM type of datum to be contained
   * @tparam IXMAN type of the grid index manager
   *
   *
   */
  template <typename DATUM, typename IXMAN>
  class GridContainerBase1D : public details::GridContainerBase<DATUM, IXMAN> {
    using Base_t = details::GridContainerBase<DATUM, IXMAN>;
    static_assert(Base_t::dims() >= 1, "GridContainerBase1D must have dimensions 1 or larger.");

  public:
    using Base_t::GridContainerBase;

    /// @{
    /// @name Data structure

    /// Returns whether the specified x index is valid
    bool hasX(typename Base_t::CellDimIndex_t index) const { return Base_t::indices.hasX(index); }

    /// Returns the size of the container in the first dimension (x)
    size_t sizeX() const { return Base_t::indices.sizeX(); }

    /// @}

  protected:
  }; // GridContainerBase1D<>

  /**
   * @brief Base class for a container of data arranged on a 2D-grid
   * @tparam DATUM type of datum to be contained
   * @tparam IXMAN type of the grid index manager
   *
   *
   */
  template <typename DATUM, typename IXMAN>
  class GridContainerBase2D : public GridContainerBase1D<DATUM, IXMAN> {
    using Base_t = GridContainerBase1D<DATUM, IXMAN>;
    static_assert(Base_t::dims() >= 2, "GridContainerBase2D must have dimensions 2 or larger.");

  public:
    using Base_t::GridContainerBase1D;

    /// @{
    /// @name Data structure

    /// Returns the size of the container in the second dimension (y)
    size_t sizeY() const { return Base_t::indices.sizeY(); }

    /// Returns whether the specified x index is valid
    bool hasY(typename Base_t::CellDimIndex_t index) const { return Base_t::indices.hasY(index); }

    /// @}

  protected:
  }; // GridContainerBase2D<>

  /**
   * @brief Base class for a container of data arranged on a 3D-grid
   * @tparam DATUM type of datum to be contained
   * @tparam IXMAN type of the grid index manager
   *
   *
   */
  template <typename DATUM, typename IXMAN>
  class GridContainerBase3D : public GridContainerBase2D<DATUM, IXMAN> {
    using Base_t = GridContainerBase2D<DATUM, IXMAN>;
    static_assert(Base_t::dims() >= 3, "GridContainerBase3D must have dimensions 3 or larger.");

  public:
    using Base_t::GridContainerBase2D;

    /// @{
    /// @name Data structure

    /// Returns the size of the container in the third dimension (z)
    size_t sizeZ() const { return Base_t::indices.sizeZ(); }

    /// Returns whether the specified x index is valid
    bool hasZ(typename Base_t::CellDimIndex_t index) const { return Base_t::indices.hasZ(index); }

    /// @}

  protected:
  }; // GridContainerBase3D<>

  /**
   * @brief Container allowing 2D indexing
   * @tparam DATUM type of contained data
   * @see GridContainer2DIndices
   *
   * This is an alias for GridContainerBase2D, with a proper index manager.
   * See the documentation of GridContainerBase2D.
   */
  template <typename DATUM>
  using GridContainer2D = GridContainerBase2D<DATUM, GridContainer2DIndices>;

  /**
   * @brief Container allowing 3D indexing
   * @tparam DATUM type of contained data
   * @see GridContainer3DIndices
   *
   * This is an alias for GridContainerBase3D, with a proper index manager.
   * See the documentation of GridContainerBase3D.
   */
  template <typename DATUM>
  using GridContainer3D = GridContainerBase3D<DATUM, GridContainer3DIndices>;

} // namespace util

#endif // LARDATA_UTILITIES_GRIDCONTAINERS_H
