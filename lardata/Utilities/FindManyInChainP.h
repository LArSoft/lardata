/**
 * @file   lardata/Utilities/FindManyInChainP.h
 * @brief  Utility to navigate chains of associations.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 26, 2017
 *
 * This library is header-only.
 *
 */

#ifndef LARDATA_UTILITIES_FINDMANYINCHAINP_H
#define LARDATA_UTILITIES_FINDMANYINCHAINP_H

// framework
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

// C/C++ standard library
#include <cstdlib> // std::size_t
#include <initializer_list>
#include <utility> // std::forward()
#include <vector>

namespace lar {

  namespace details {

    template <typename T, typename R = void>
    struct enable_if_type_exists;

    template <typename T, typename R = void>
    using enable_if_type_exists_t = typename enable_if_type_exists<T, R>::type;

    template <typename H, typename R = void>
    using enable_if_is_handle_t = enable_if_type_exists_t<typename std::decay_t<H>::HandleTag, R>;

  } // namespace details

  /// Type for default tag in `FindManyInChainP` constructors.
  struct SameAsDataTag {};

  /// Value for default tag in `FindManyInChainP` constructors.
  constexpr SameAsDataTag SameAsData;

  /**
   * @brief  Query object collecting a list of associated objects.
   * @tparam Target type of objects to be fetched
   * @tparam Intermediate types of objects connecting to Target by association
   *
   * This query object collects information about all objects of type `Target`
   * associated to each specified object of type `Source`.
   * The `Source` type is implicitly specified in the constructor.
   * For example, each `recob::Shower` object is expected to be associated to
   * a number of `recob::Cluster` objects, and each one of these clusters must
   * be associated to `recob::Hit` objects. To retrieve all the hit objects
   * associated to a shower collection (_not recommended_):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto showers = event.getValidHandle<std::vector<recob::Shower>>(showerTag);
   * lar::FindManyInChainP<recob::Hit, recob::Cluster> showerToHits
   *   (showers, event, showerTag);
   *
   * for (std::size_t iShower = 0; iShower < showers.size(); ++iShower) {
   *
   *   recob::Shower const& shower = (*showers)[iShower];
   *
   *   decltype(auto) showerHits = showerToHits.at(iShower);
   *
   *   // ...
   *
   * } // for each shower
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this example, it is explicitly stated that the producer of the shower
   * associations to clusters is the same as the producer of the showers, as the
   * same input tag `showerTag` is used. It is also implicitly assumed that the
   * same producer which created the associated clusters is also responsible
   * for the creation of the associations between clusters and hits.
   *
   * The example shows the similarity with `art::FindManyP`. In fact, this
   * object has a similar interface, and a similar functionality too.
   * Nevertheless, prefer using the static member `find()` and then using
   * its result directly, as it has a more complete interface:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto showers = event.getValidHandle<std::vector<recob::Shower>>(showerTag);
   * std::size_t iShower = 0;
   * for (auto const& showerHits:
   *   lar::FindManyInChainP<recob::Hit, recob::Cluster>
   *     (showers, event, showerTag)
   *   )
   * {
   *   recob::Shower const& shower = (*showers)[iShower++];
   *   // all associated hits are already in showerHits
   *
   *   // ...
   *
   * } // for each shower
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * This is the recommended approach, as long as `lar::FindManyInChainP()` has
   * to be called only once.
   *
   * @note Due to the inability to retrieve provenance information in _gallery_,
   *       this class is not compatible with _gallery_.
   *
   */
  template <typename Target, typename... Intermediate>
  class FindManyInChainP {

  public:
    using Target_t = Target;                ///< Type of the associated objects.
    using TargetPtr_t = art::Ptr<Target_t>; ///< Pointer to associated objects.

    /// Type returned by `at()` method.
    using TargetPtrCollection_t = std::vector<TargetPtr_t>;

    /**
     * @brief Constructor: extracts target objects associated to all objects
     *        under the specified handle.
     * @tparam Source type of source: art Handle or collection of art pointers
     * @tparam Event type of event to be used (either _art_ or gallery `Event`)
     * @tparam InputTags a variable number of `art::InputTag` objects
     * @param source art Handle or collection of art pointers to source objects
     * @param event the event to read associations and objects from
     * @param tags input tags for each one of the required associations
     * @see find()
     *
     * This constructor finds the associated objects as in `find()`, and stored
     * the result. Access to it will be performed e.g. by the `at()` method.
     *
     */
    template <typename Source, typename Event, typename... InputTags>
    FindManyInChainP(Source&& source, Event const& event, InputTags... tags)
      : results(find(std::forward<Source>(source), event, tags...))
    {}

    /// Returns the number of `Source` objects we have information about.
    std::size_t size() const noexcept;

    /**
     * @brief Returns all the `Target` objects associated to specified object.
     * @param i index of the source object to query
     * @return a sequence of _art_ pointers to the associated `Target` objects
     * @throw std::out_of_range if the specified index is not valid
     *
     * The specified index matches the index of the element in the collection
     * this query object was constructed with. For example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto showers
     *   = event.getValidHandle<std::vector<recob::Shower>>(showerTag);
     * lar::FindManyInChainP<recob::Hit, recob::Cluster> showerToHits
     *   (shower, event, showerTag);
     *
     * for (std::size_t iShower = 0; iShower < showers.size(); ++iShower) {
     *
     *   recob::Shower const& shower = (*showers)[iShower];
     *
     *   decltype(auto) showerHits = showerToHits.at(iShower);
     *
     *   // ...
     *
     * } // for each shower
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * the `showerHits` will pertain the `iShower`-th `recob::Shower` (`shower`)
     * on each iteration of the loop.
     *
     * The returned sequence points to objects convertible to _art_ pointers to
     * `Target` (e.g., `art::Ptr<Target>`). The order of the objects in this
     * sequence is not defined. The same `Target` object _may_ appear more than
     * once if two intermediate objects are associated to the same target (or
     * intermediate object of tier closer to the target).
     * This sequence is of an unspecified type, guaranteed to support
     * `std::begin()`, `std::end()`, `size()` method and to be iterable in a
     * range-for loop.
     * For example, it may be a STL vector returned by constant reference.
     */
    TargetPtrCollection_t const& at(std::size_t i) const;

    /**
     * @brief Returns target objects associated to all objects contained in the
     *        specified source.
     * @tparam Source type of source: art Handle or collection of art pointers
     * @tparam Event type of event to be used (either _art_ or gallery `Event`)
     * @tparam InputTags a variable number of `art::InputTag` objects
     * @param source art Handle or collection of art pointers to source objects
     * @param event the event to read associations and objects from
     * @param tags input tags for each one of the required associations
     * @return collection of lists of pointers to associated objects, one for
     *         each source element, in the same order
     *
     * This methods returns a collection with an entry for each of the elements
     * pointed by the specified handle, in the same order as they are extracted
     * from the source collection.
     *
     * The input tag arguments must be convertible to `art::InputTag` type.
     * The first tag identifies the data product containing the associations
     * between the `Source` collection (the same pointed by `source`) and the
     * first `Intermediate` class. The second tag likewise points to the data
     * product containing the associations between the first and the second
     * `Intermediate` classes. The last tag, finally, points to the data product
     * containing the associations between the last `Intermediate` class and the
     * `Target` class.
     * If one tag is not specified, it is assumed that the same class that
     * produced the `Intermediate` (or `Source`) data product also produced the
     * association to the next `Intermediate` (or `Target`) data product.
     * Note that:
     * * if the tag is specified, _all_ associations are expected to come from
     *   that same tag; if the tag is _not_ specified, elements from different
     *   data products will be assumed to have different association data
     *   products too
     * * it is possible to "omit" the explicit specification of a tag by using
     *   the `lar::SameAsData` tag
     *
     */
    template <typename Source, typename Event, typename... InputTags>
    static std::vector<TargetPtrCollection_t> find(Source&& source,
                                                   Event const& event,
                                                   InputTags... tags);

  private:
    std::vector<TargetPtrCollection_t> results; ///< Stored results.

  }; // class FindManyInChainP<>

} // namespace lar

//------------------------------------------------------------------------------
//---  template implementation
//---
#include "FindManyInChainP.tcc" // expected in the same directory as this file

//------------------------------------------------------------------------------

#endif // LARDATA_UTILITIES_FINDMANYINCHAINP_H
