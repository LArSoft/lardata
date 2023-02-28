/**
 * @file   lardata/RecoBaseProxy/ProxyBase/ParallelData.h
 * @brief  Auxiliary data from parallel data products.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_PARALLELDATA_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_PARALLELDATA_H

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h" // util::collection_value_t, ...
#include "lardata/Utilities/TupleLookupByTag.h" // util::makeTagged(), ...

// C/C++ standard
#include <cstdlib>     // std::size_t
#include <type_traits> // std::is_convertible<>, ...
#include <utility>     // std::declval()

namespace proxy {

  // --- BEGIN LArSoftProxiesParallelData --------------------------------------
  /// @addtogroup LArSoftProxiesParallelData
  /// @{

  /**
   * @brief Wraps a collection into a parallel data collection object.
   * @tparam AuxColl type of parallel data data product container
   * @tparam Aux type of parallel data to be associated to the main objects
   *             (if omitted: `AuxColl::value_type`)
   * @tparam Tag the tag labelling this associated data (if omitted: as `Aux`)
   * @param data data collection to be wrapped
   * @return a new `ParallelData` wrapping the information in `data`
   *
   * The data collection must be non-temporary and it is treated as fulfilling
   * @ref LArSoftProxyDefinitionParallelData "parallel data product"
   * requirements.
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<recob::TrackFitHitInfo> trackData;
   * // ...
   * auto auxData = makeParallelData(trackData);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * where the tag of the parallel data is now `recob::TrackFitHitInfo` and
   * `auxData` behaviour becomes undefined as soon as `trackData` falls out of
   * scope.
   */
  template <typename AuxColl, typename Aux = util::collection_value_t<AuxColl>, typename Tag = Aux>
  auto makeParallelData(AuxColl const& data);

  //----------------------------------------------------------------------------
  namespace details {

    /**
     * @brief Object to draft parallel data interface.
     * @tparam AuxColl type of the parallel data collection
     * @tparam Aux type of the associated object
     * @tparam Tag tag this data is labeled with
     *
     * Allows:
     *  * random access (no index check guarantee)
     *  * forward iteration
     *
     * Construction is not part of the interface.
     */
    template <typename AuxColl,
              typename Aux /* = util::collection_value_t<AuxColl> */,
              typename Tag /* = Aux */
              >
    class ParallelData {
      using This_t = ParallelData<AuxColl, Aux, Tag>; ///< This type.

      /// Type of auxiliary collection.
      using parallel_data_t = AuxColl;

      /// Type of the value of auxiliary collection element.
      using aux_t = Aux; // unused

      /// Type returned when accessing an auxiliary collection element.
      using aux_element_t = util::collection_value_constant_access_t<AuxColl>;

      using parallel_data_iterator_t = typename parallel_data_t::const_iterator;

    public:
      using tag = Tag; ///< Tag of this association proxy.

      /// Type returned when accessing auxiliary data.
      using auxiliary_data_t = decltype(util::makeTagged<tag>(std::declval<aux_element_t>()));

      /// Constructor: points to the specified data collection.
      ParallelData(parallel_data_t const& data) : fData(&data) {}

      /// Returns an iterator pointing to the first data element.
      auto begin() const -> decltype(auto) { return fData->begin(); }

      /// Returns an iterator pointing past the last data element.
      auto end() const -> decltype(auto) { return fData->end(); }

      /// Returns the element with the specified index (no check performed).
      auto operator[](std::size_t index) const -> decltype(auto)
      {
        static_assert(std::is_convertible<decltype(getElement(index)), auxiliary_data_t>(),
                      "Inconsistent data types.");
        return getElement(index);
      }

      /// Returns whether this data is labeled with the specified tag.
      template <typename TestTag>
      static constexpr bool hasTag()
      {
        return std::is_same<TestTag, tag>();
      }

      /// Returns a pointer to the whole data collection.
      parallel_data_t const* data() const { return fData; }

      /// Returns a reference to the whole data collection.
      parallel_data_t const& dataRef() const { return *(data()); }

    private:
      parallel_data_t const* fData; ///< Reference to the original data product.

      auto getElement(std::size_t index) const -> decltype(auto)
      {
        return util::makeTagged<tag>(fData->operator[](index));
      }

    }; // class ParallelData<>

    //--------------------------------------------------------------------------

  } // namespace details

  /// @}
  // --- END LArSoftProxiesParallelData ----------------------------------------

} // namespace proxy

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace proxy {

  //----------------------------------------------------------------------------
  //--- makeParallelData() implementations
  //----------------------------------------------------------------------------
  template <typename AuxColl,
            typename Aux /* = util::collection_value_t<AuxColl>*/,
            typename Tag /* = Aux */
            >
  auto makeParallelData(AuxColl const& data)
  {

    // Ahh, simplicity.
    return details::ParallelData<AuxColl, Aux, Tag>(data);

  } // makeParallelData(AuxColl)

  //----------------------------------------------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_PARALLELDATA_H
