/**
 * @file   lardata/RecoBaseProxy/ProxyBase/OneTo01Data.h
 * @brief  Auxiliary data from one-to-(zero-or-one) sequential association.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_ONETO01DATA_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_ONETO01DATA_H

// LArSoft libraries
#include "lardata/RecoBaseProxy/ProxyBase/AssnsTraits.h"
#include "lardata/Utilities/TupleLookupByTag.h" // util::add_tag_t(), ...

// framework libraries
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"

// C/C++ standard
#include <cstdlib>     // std::size_t
#include <iterator>    // std::cbegin(), std::cend()
#include <tuple>       // std::tuple_element_t<>, std::get()
#include <type_traits> // std::is_convertible<>
#include <utility>     // std::move()
#include <vector>

namespace proxy {

  namespace details {

    /**
     * @brief Object for one-to-zero/or/one associated data interface.
     * @tparam Main type of the main associated object (one)
     * @tparam Aux type of the additional associated objects (zero or one)
     * @tparam Metadata type of associated metadata
     * @tparam Tag tag this data is labeled with
     *
     * Allows:
     *  * random access (no index check guarantee)
     *  * forward iteration
     *
     * Construction is not part of the interface.
     *
     * The `OneTo01Data` object acquires a vector of _art_ pointers, one for
     * each element in the main collection.
     * It is an implementation detail for associations fulfilling the
     * @ref LArSoftProxyDefinitionOneToZeroOrOneSeqAssn "one-to-(zero-or-one) sequential association"
     * requirement.
     *
     * @note This data structure marks the main elements which have no
     *       associated data with an invalid _art_ pointer (default-constructed)
     *       and it does not distinguish that from the element being actually
     *       associated to a default-constructed _art_ pointer.
     *
     * The `OneTo01Data` object also provides a container-like view of this
     * information, where each element in the container is associated to a
     * single `Main` and it is an _art_ pointer to the `Right` element.
     *
     * Association metadata is not accessible from this object.
     *
     * @todo Metadata for `proxy::details::OneTo01Data` is not supported yet.
     */
    template <typename Main, typename Aux, typename Metadata /* = void */, typename Tag /* = Aux */
              >
    class OneTo01Data {
      using This_t = OneTo01Data<Main, Aux, Metadata, Tag>; ///< This type.

    public:
      /// Type of associated datum.
      using aux_t = Aux;

      /// Type of associated metadata.
      using metadata_t = Metadata;

      /// Type of tag.
      using tag = Tag;

      /// Type of main datum.
      using main_t = Main;

      /// Type of _art_ pointer to associated datum.
      using aux_ptr_t = art::Ptr<aux_t>;

      /// Type of auxiliary data associated with a main item.
      using auxiliary_data_t = util::add_tag_t<aux_ptr_t, tag>;

      /// Type of collection of auxiliary data for all main elements.
      using aux_coll_t = std::vector<aux_ptr_t>;

      /// Type of the source association.
      using assns_t = art::Assns<main_t, aux_t>;

      OneTo01Data(aux_coll_t&& data) : auxData(std::move(data)) {}

      /// Returns whether the element `i` is associated with auxiliary datum.
      bool has(std::size_t i) const { return get(i) == aux_ptr_t(); }

      /// Returns a copy of the pointer to data associated with element `i`.
      auxiliary_data_t get(std::size_t i) const { return auxiliary_data_t(auxData[i]); }

      /// Returns the range with the specified index (no check performed).
      auto operator[](std::size_t index) const -> decltype(auto)
      {
        static_assert(std::is_convertible<decltype(get(index)), auxiliary_data_t>(),
                      "Inconsistent data types.");
        return get(index);
      }

    private:
      aux_coll_t auxData; ///< Data associated to the main collection.

    }; // class OneTo01Data<>

  } // namespace details

  //@{
  /**
   * @brief Processes and returns an one-to-(zero/one) associated data object.
   * @tparam Tag the tag labelling this associated data
   *             (if omitted: second type of the association: `right_t`)
   * @tparam Assns type of association to be processed
   * @param assns association object to be processed
   * @param minSize minimum number of entries in the produced association data
   * @return a new `OneTo01Data` filled with associations from `tag`
   *
   * The content of the association object must fulfill the requirements of
   * @ref LArSoftProxyDefinitionOneToZeroOrOneSeqAssn "one-to-(zero or one) sequential association".
   * The `Assns` type is expected to be a `art::Assns` instance. At least,
   * the `Assns` type is required to have `left_t` and `right_t` definitions
   * representing respectively the main data type and the associated one, and
   * to respond to `begin()` and `end()` functions. The iterated object must
   * also respond to `std::get<0>()` with a `art::Ptr<left_t>` and to
   * `std::get<1>()` with a `art::Ptr<right_t>`.
   *
   * Elements in the main collection not associated with any object will present
   * an invalid _art_ pointer (`isNull()` true). If there is information for
   * less than `minSize` main objects, more records will be added to mark the
   * missing objects as not associated to anything.
   *
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * art::Assns<recob::Track, recob::Vertex> trackVertexAssns;
   * // ...
   * auto assData = proxy::makeOneTo01data(assns);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will have `assData` tagged as `recob::Vertex`.
   *
   * If `Tag` type is omitted, the class to the right side of the association
   * is used as tag.
   */
  template <typename Tag, typename Assns>
  auto makeOneTo01data(Assns const& assns, std::size_t minSize = 0);

  template <typename Assns>
  auto makeOneTo01data(Assns const& assns, std::size_t minSize = 0)
  {
    return makeOneTo01data<typename Assns::right_t>(assns, minSize);
  }
  //@}

  //@{
  /**
   * @brief Creates and returns an one-to-(zero/one) associated data object.
   * @tparam Tag the tag labelling this associated data
   *             (if omitted: second type of the association: `right_t`)
   * @tparam MainColl type of the main collection object
   * @tparam Assns type of the association object
   * @param mainColl the main collection object
   * @param assns association data object
   * @return a new `OneTo01Data` wrapping the information in `assns`
   * @see `makeOneTo01data(Assns const&, std::size_t)`
   *
   * This function operates like
   * `makeOneTo01data(Assns const&, std::size_t)`, where the size is
   * extracted from the main data collection.
   *
   * If `Tag` type is omitted, the class to the right side of the association
   * is used as tag.
   */
  template <typename Tag, typename MainColl, typename Assns>
  auto makeOneTo01data(MainColl const& mainColl, Assns const& assns)
  {
    return makeOneTo01data<Tag>(assns, mainColl.size());
  }

  template <typename MainColl, typename Assns>
  auto makeOneTo01data(MainColl const& mainColl, Assns const& assns)
  {
    return makeOneTo01data<typename Assns::right_t>(mainColl, assns);
  }
  //@}

} // namespace proxy

//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace proxy {

  namespace details {

    //--------------------------------------------------------------------------
    // Extends vector v with default-constructed data
    // and executes v[index]=value
    template <typename T>
    void extendAndAssign(std::vector<T>& v,
                         typename std::vector<T>::size_type index,
                         typename std::vector<T>::value_type const& value)
    {
      if (index >= v.size()) {
        v.reserve(index + 1);
        v.resize(index);
        v.push_back(value);
      }
      else
        v[index] = value;
    } // extendAndAssign()

    // Extends vector v with default-constructed data
    // and executes v[index]=move(value)
    template <typename T>
    void extendAndAssign(std::vector<T>& v,
                         typename std::vector<T>::size_type index,
                         typename std::vector<T>::value_type&& value)
    {
      if (index >= v.size()) {
        v.reserve(index + 1);
        v.resize(index);
        v.push_back(std::move(value));
      }
      else
        v[index] = std::move(value);
    } // extendAndAssign()

    //--------------------------------------------------------------------------
    template <std::size_t Key, std::size_t Data, typename Iter>
    auto associationOneToOneFullSequence(Iter begin, Iter end, std::size_t n)
    {
      //
      // Here we are actually not using the assumption that the keys are in
      // increasing order; which is just as good as long as we use a fast random
      // access container as STL vector.
      // We do assume the key side of the association to be valid, though.
      //
      using value_type = typename Iter::value_type;
      using data_t = std::tuple_element_t<Data, value_type>;
      std::vector<data_t> data(n); // all default-constructed
      for (auto it = begin; it != end; ++it) {
        auto const& keyPtr = std::get<Key>(*it);
        extendAndAssign(data, keyPtr.key(), std::get<Data>(*it));
      }
      return data;
    } // associationOneToOneFullSequence(Iter, Iter, std::size_t)

  } // namespace details

  //----------------------------------------------------------------------------
  //--- makeOneTo01data() implementation
  //----------------------------------------------------------------------------
  template <typename Tag, typename Assns>
  auto makeOneTo01data(Assns const& assns, std::size_t minSize /* = 0 */)
  {
    using Main_t = typename Assns::left_t;
    using Aux_t = typename Assns::right_t;
    using Metadata_t = lar::util::assns_metadata_t<Assns>;
    using AssociatedData_t = details::OneTo01Data<Main_t, Aux_t, Metadata_t, Tag>;

    using std::cbegin;
    using std::cend;
    return AssociatedData_t(
      details::associationOneToOneFullSequence<0U, 1U>(cbegin(assns), cend(assns), minSize));
  } // makeOneTo01data(assns)

  //----------------------------------------------------------------------------

} // namespace proxy

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_ONETO01DATA_H
