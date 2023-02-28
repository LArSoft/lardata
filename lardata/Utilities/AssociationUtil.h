/**
 * @file   AssociationUtil.h
 * @author brebel@fnal.gov
 * @brief  Utility object to perform functions of association
 *
 * @attention Please considering using the lightweight utility `art::PtrMaker`
 *            instead.
 *
 * This library provides a number of util::CreateAssn() functions;
 * for convenience, the ones supported as of January 2015 are listed here:
 *
 * -# `CreateAssn(art::Event&, std::vector<T> const&, art::Ptr<U> const&b, art::Assns<U,T> &, std::string, size_t)`
 *   one-to-one association, between a element of a vector (future data product)
 *   and a art-pointed object; allows for a instance name for the vector
 * -# `CreateAssn(art::Event&, std::vector<T> const&, art::Ptr<U> const&b, art::Assns<U,T> &, size_t)`
 *   one-to-one association, between a element of a vector (future data product)
 *   and a art-pointed object
 * -# `CreateAssn(art::Event&, art::Ptr<T> const&, art::Ptr<U> const&, art::Assns<U,T>&)`
 *   one-to-one association, between two art-pointed objects
 * -# `CreateAssn(art::Event&, std::vector<T> const&, art::PtrVector<U> const&, art::Assns<T,U>&, size_t)`
 *   one-to-many association, between a element of a vector (future data
 *   product) and a list of art-pointed objects in a art::PtrVector list
 * -# `CreateAssn(art::Event&, art::Ptr<T> const&, std::vector<art::Ptr<U>> const&, art::Assns<T,U>&)`
 *   one-to-many association, between an art-pointed object and all the elements
 *   in a std::vector of art-pointed objects
 * -# `CreateAssn(art::Event&, std::vector<T> const&, std::vector<art::Ptr<U>>&, art::Assns<T,U>&, size_t)`
 *   one-to-many association, between an element of a vector (future data
 *   product) and all the elements in a std::vector of art-pointed objects
 * -# `CreateAssn(art::Event&, std::vector<T> const&, std::vector<U> const&, art::Assns<T,U>&, size_t, size_t, size_t)`
 *   one-to-many association, between an element of a vector (future data
 *   product) and the elements in a subrange of a std::vector (also future data
 *   product)
 * -# `CreateAssn(art::Event&, art::Assns<T,U>&, size_t, Iter, Iter)`
 *   one-to-many association, between an element of a collection and the
 *   elements of another collection, whose indices are specified by the values
 *   in a subrange of a third collection (of indices)
 * -# @code CreateAssnD(art::Event&, art::Assns<T,U>&, size_t, size_t, typename art::Assns<T,U,D>::data_t const&) @endcode,
 *   @code CreateAssnD(art::Event&, art::Assns<T,U>&, size_t, size_t, typename art::Assns<T,U,D>::data_t&&) @endcode
 *   one-to-one association, between an element of a collection and the element
 *   of another collection, both specified by index, with additional data
 *
 * For all the associated objects, either side, that are not specified by
 * `art::Ptr`, the index of the object in its collection must be (or stay) the
 * same as the index in the final data product collection.
 *
 *
 * One-to-one associations
 * ------------------------
 *
 * the one (T)            | the other one (U) | special notes               | function
 * :--------------------: | :---------------: | --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------
 * element of std::vector | art pointer       | allows instance name for T  | CreateAssn(art::Event &, std::vector< T > const &, art::Ptr< U > const &, art::Assns< U, T > &, std::string, size_t)
 * element of std::vector | art pointer       |                             | CreateAssn(art::Event&, std::vector<T> const&, art::Ptr<U> const&, art::Assns<U,T> &, size_t)
 * art pointer            | art pointer       |                             | CreateAssn(art::Event&, art::Ptr<T> const&, art::Ptr<U> const&, art::Assns<U,T>&)
 * element by index       | element by index  | associates data too         | CreateAssnD(art::Event&, art::Ptr<T> const&, art::Ptr<U> const&, art::Assns<U,T,D>&, size_t, size_t, typename art::Assns<T,U,D>::data_t const&)
 * element by index       | element by index  | associates data too (moved) | CreateAssnD(art::Event&, art::Ptr<T> const&, art::Ptr<U> const&, art::Assns<U,T,D>&, size_t, size_t, typename art::Assns<T,U,D>::data_t&&)
 *
 *
 * One-to-many associations
 * -------------------------
 *
 * the one (T)            | the many (U)             | special notes              | function
 * :--------------------: | :----------------------: | -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------
 * element of std::vector | art::PtrVector           |                            | CreateAssn(art::Event&, std::vector<T> const&, art::PtrVector<U> const&, art::Assns<T,U>&, size_t)
 * art pointer            | std::vector<art::Ptr<U>> |                            | CreateAssn(art::Event&, art::Ptr<T> const&, std::vector<art::Ptr<U>> const&, art::Assns<T,U>&)
 * element of std::vector | std::vector<art::Ptr<U>> |                            | CreateAssn(art::Event&, std::vector<T> const&, std::vector<art::Ptr<U>>&, art::Assns<T,U>&, size_t)
 * element of std::vector | std::vector<U>           |                            | CreateAssn(art::Event&, std::vector<T> const&, std::vector<U> const&, art::Assns<T,U>&, size_t, size_t, size_t)
 * element by index       | range of indices         | does not need object lists | CreateAssn(art::Event&, art::Assns<T,U>&, size_t, Iter, Iter)
 *
 */

#ifndef ASSOCIATIONUTIL_H
#define ASSOCIATIONUTIL_H

// C/C++ standard libraries
#include <string>
#include <utility> // std::move()
#include <vector>

// framework libraries
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace util {

  // see https://cdcvs.fnal.gov/redmine/projects/art/wiki/Inter-Product_References
  // for information about using art::Assns

  /**
   * @brief Creates a single one-to-one association
   * @tparam T type of the new object to associate
   * @tparam U type of the object already in the data product or art::Ptr
   * @param evt reference to the current event
   * @param a vector of data products that are in, or will be put into, evt
   * @param b art::Ptr to the (new) object to be associated to the one in a
   * @param assn reference to association object where the new one will be put
   * @param a_instance name of the instance that will be used for a in evt
   * @param index index of the element in a to be associated with b (default: the last element)
   * @return whether the operation was successful (can it ever fail??)
   *
   * As example of usage: create a wire/raw digit association.
   * This code should live in the art::EDProduce::produce() method.
   * The raw::RawDigit product was created already by a DigitModuleLabel module.
   * The code is supposed to produce one recob::Wire for each existing
   * raw::RawDigit, and contextually associate the new wire to the source digit.
   * We are also assuming that there might be different RawDigit sets produced
   * by the same producer: we identify the one we care of by the string
   * spill_name and we create wires and associations with the same label
   * for convenience.
   *
   *     // this is the original list of digits, thawed from the event
   *     art::Handle< std::vector<raw::RawDigit>> digitVecHandle;
   *     evt.getByLabel(DigitModuleLabel, spill_name, digitVecHandle);
   *
   *     // the collection of wires that will be written as data product
   *     std::unique_ptr<std::vector<recob::Wire>> wirecol(new std::vector<recob::Wire>);
   *     // ... and an association set
   *     std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire>> WireDigitAssn
   *       (new art::Assns<raw::RawDigit,recob::Wire>);
   *
   *     for(size_t iDigit = 0; iDigit < digitVecHandle->size(); ++iDigit) {
   *       // turn the digit into a art::Ptr:
   *       art::Ptr<raw::RawDigit> digit_ptr(digitVecHandle, iDigit);
   *
   *       // store the wire in its final position in the data product;
   *       // the new wire is currently the last of the list
   *       wirecol->push_back(std::move(wire));
   *
   *       // add an association between the last object in wirecol
   *       // (that we just inserted) and digit_ptr
   *       if (!util::CreateAssn(*this, evt, *wirecol, digit_ptr, *WireDigitAssn, spill_name)) {
   *         throw art::Exception(art::errors::ProductRegistrationFailure)
   *           << "Can't associate wire #" << (wirecol->size() - 1)
   *           << " with raw digit #" << digit_ptr.key();
   *       } // if failed to add association
   *
   *     } // for digits
   *
   *     evt.put(std::move(wirecol), spill_name);
   *     evt.put(std::move(WireDigitAssn), spill_name);
   *
   */
  // MARK CreateAssn_01
  template <class T, class U>
  bool CreateAssn(art::Event& evt,
                  std::vector<T> const& a,
                  art::Ptr<U> const& b,
                  art::Assns<U, T>& assn,
                  std::string a_instance,
                  size_t index = UINT_MAX);

  /**
   * @brief Creates a single one-to-one association
   * @tparam T type of the new object to associate
   * @tparam U type of the object already in the data product or art::Ptr
   * @param evt reference to the current event
   * @param a vector of data products that are in, or will be put into, evt
   * @param b art::Ptr to the (new) object to be associated to the one in a
   * @param assn reference to association object where the new one will be put
   * @param index index of the element in a to be associated with b (default: the last element)
   * @return whether the operation was successful (can it ever fail??)
   *
   * The instance name of the product a will be in is assumed empty.
   * Example of usage:
   *
   *     // this is the original list of digits, thawed from the event
   *     art::Handle< std::vector<raw::RawDigit>> digitVecHandle;
   *     evt.getByLabel(DigitModuleLabel, digitVecHandle);
   *
   *     // the collection of wires that will be written as data product
   *     std::unique_ptr<std::vector<recob::Wire>> wirecol(new std::vector<recob::Wire>);
   *     // ... and an association set
   *     std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire>> WireDigitAssn
   *       (new art::Assns<raw::RawDigit,recob::Wire>);
   *
   *     for(size_t iDigit = 0; iDigit < digitVecHandle->size(); ++iDigit) {
   *       // turn the digit into a art::Ptr:
   *       art::Ptr<raw::RawDigit> digit_ptr(digitVecHandle, iDigit);
   *
   *       // store the wire in its final position in the data product;
   *       // the new wire is currently the last of the list
   *       wirecol->push_back(std::move(wire));
   *
   *       // add an association between the last object in wirecol
   *       // (that we just inserted) and digit_ptr
   *       if (!util::CreateAssn(*this, evt, *wirecol, digit_ptr, *WireDigitAssn)) {
   *         throw art::Exception(art::errors::ProductRegistrationFailure)
   *           << "Can't associate wire #" << (wirecol->size() - 1)
   *           << " with raw digit #" << digit_ptr.key();
   *       } // if failed to add association
   *
   *     } // for digits
   *
   *     evt.put(std::move(wirecol));
   *     evt.put(std::move(WireDigitAssn));
   *
   */
  // MARK CreateAssn_02
  template <class T, class U>
  inline bool CreateAssn(art::Event& evt,
                         std::vector<T> const& a,
                         art::Ptr<U> const& b,
                         art::Assns<U, T>& assn,
                         size_t index = UINT_MAX)
  {
    return CreateAssn(evt, a, b, assn, std::string(), index);
  }

  /**
   * @brief Creates a single one-to-one association
   * @tparam T type of one object to associate
   * @tparam U type of the other object to associate
   * @param evt reference to the current event
   * @param a art::Ptr to the first object in the association
   * @param b art::Ptr to the object to be associated to the one in a
   * @param assn reference to association object where the new one will be put
   * @return whether the operation was successful (can it ever fail??)
   *
   * This is the simplest way ever.
   * Neither the event not the producer references are used.
   */
  // MARK CreateAssn_03
  template <class T, class U>
  bool CreateAssn(art::Event& evt,
                  art::Ptr<T> const& a,
                  art::Ptr<U> const& b,
                  art::Assns<U, T>& assn);

  /**
   * @brief Creates a single one-to-many association
   * @tparam T type of the new object to associate
   * @tparam U type of the many objects already in the data product or art::Ptr
   * @param evt reference to the current event
   * @param a vector of data products that are in, or will be put into, evt
   * @param b art::PtrVector to the (new) objects to be associated to the one in a
   * @param assn reference to association object where the new one will be put
   * @param index index of the element in a to be associated with all the ones
   *             in b (default: the last element)
   * @return whether the operation was successful (can it ever fail??)
   *
   * A "one-to-many" association is actually a number of one-to-one
   * associations. If you want to keep the information of the order of the many,
   * you may have to use an association with a data member (the third template
   * parameter that we pretent not to exist).
   */
  // MARK CreateAssn_04
  template <class T, class U>
  bool CreateAssn(art::Event& evt,
                  std::vector<T> const& a,
                  art::PtrVector<U> const& b,
                  art::Assns<T, U>& assn,
                  size_t index = UINT_MAX);

  /**
   * @brief Creates a single one-to-many association
   * @tparam T type of the new object to associate
   * @tparam U type of the many objects already in the data product or art::Ptr
   * @param evt reference to the current event
   * @param a art::Ptr to the item to be associated with many
   * @param b vector to art::Ptr to the (new) objects to be associated to a
   * @param assn reference to association object where the new one will be put
   * @return whether the operation was successful (can it ever fail??)
   *
   * A "one-to-many" association is actually a number of one-to-one
   * associations. If you want to keep the information of the order of the many,
   * you may have to use an association with a data member (the third template
   * parameter that we pretent not to exist).
   */
  // method to create a 1 to many association, with the many being of type U
  // index is the location in the input std::vector<T> of the object you wish to
  // associate with the art::PtrVector<U>
  // MARK CreateAssn_05
  template <class T, class U>
  bool CreateAssn(art::Event& evt,
                  art::Ptr<T> const& a,
                  std::vector<art::Ptr<U>> const& b,
                  art::Assns<T, U>& assn);

  /**
   * @brief Creates a single one-to-many association
   * @tparam T type of the new object to associate
   * @tparam U type of the many objects already in the data product or art::Ptr
   * @param evt reference to the current event
   * @param a vector of data products that are in, or will be put into, evt
   * @param b vector to art::Ptr to the (new) objects to be associated to the one in a
   * @param assn reference to association object where the new one will be put
   * @param index index of the element in a to be associated with all the ones
   *             in b (default: the last element)
   * @return whether the operation was successful (can it ever fail??)
   *
   * A "one-to-many" association is actually a number of one-to-one
   * associations. If you want to keep the information of the order of the many,
   * you may have to use an association with a data member (the third template
   * parameter that we pretent not to exist).
   */
  // MARK CreateAssn_06
  template <class T, class U>
  bool CreateAssn(art::Event& evt,
                  std::vector<T> const& a,
                  std::vector<art::Ptr<U>> const& b,
                  art::Assns<T, U>& assn,
                  size_t index = UINT_MAX);

  /**
   * @brief Creates a single one-to-many association
   * @tparam T type of the new object to associate
   * @tparam U type of the many objects already in the data product or art::Ptr
   * @param evt reference to the current event
   * @param a vector of data products that are in, or will be put into, evt
   * @param b vector of the (new) objects to be associated to the one in a
   * @param assn reference to association object where the new one will be put
   * @param startU index in b of the first element to be associated to the one in a
   * @param endU index in b after the last element to be associated to the one in a
   * @param index index of the element in a to be associated with all the ones
   *             in b (default: the last element)
   * @return whether the operation was successful (can it ever fail??)
   *
   * Use this when the objects in b are not yet stored in the event and are in a
   * std::vector collection instead.
   *
   * The method gets the product id for those as well as for the element in a.
   * Also specify the range of entries to use from the std::vector collection of
   * U objects.
   *
   * A "one-to-many" association is actually a number of one-to-one
   * associations. If you want to keep the information of the order of the many,
   * you may have to use an association with a data member (the third template
   * parameter that we pretent not to exist).
   */
  // MARK CreateAssn_07
  template <class T, class U>
  bool CreateAssn(art::Event& evt,
                  std::vector<T> const& a,
                  std::vector<U> const& b,
                  art::Assns<T, U>& assn,
                  size_t startU,
                  size_t endU,
                  size_t index = UINT_MAX);

  /**
   * @brief Creates a single one-to-many association
   * @tparam T type of the new object to associate
   * @tparam U type of the many objects already in the data product or art::Ptr
   * @param evt reference to the current event
   * @param a vector of data products that are in, or will be put into, evt
   * @param b vector of the (new) objects to be associated to the one in a
   * @param assn reference to association object where the new one will be put
   * @param indices indices of the elements in b to be associated to the one in a
   * @param index index of the element in a to be associated with all the ones
   *             in b (default: the last element)
   * @return whether the operation was successful (can it ever fail??)
   *
   * Use this when the objects in b are not yet stored in the event and are in a
   * std::vector collection instead.
   *
   * The method gets the product id for those as well as for the element in a.
   * Also specify the entries to use from the std::vector collection of
   * U objects.
   *
   * A "one-to-many" association is actually a number of one-to-one
   * associations. If you want to keep the information of the order of the many,
   * you may have to use an association with a data member (the third template
   * parameter that we pretent not to exist).
   */
  // MARK CreateAssn_07a
  template <class T, class U>
  bool CreateAssn(art::Event& evt,
                  std::vector<T> const& a,
                  std::vector<U> const& b,
                  art::Assns<T, U>& assn,
                  std::vector<size_t> const& indices,
                  size_t index = UINT_MAX);

  /**
   * @brief Creates a single one-to-many association
   * @tparam T type of the new object to associate
   * @tparam U type of the many objects already in the data product or art::Ptr
   * @tparam Iter iterator to size_t-compatible elements
   * @param evt reference to the current event
   * @param assn reference to association object where the new one will be put
   * @param first_index index of the object of type T to be associated to all
   *        the others
   * @param from_second_index iterator pointing to the first of the indices
   *        of U type objects to be associated to the one of the first type
   * @param to_second_index iterator pointing after the last of the indices
   *        of U type objects to be associated to the one of the first type
   * @return whether the operation was successful (can it ever fail??)
   *
   * A "one-to-many" association is actually a number of one-to-one
   * associations. If you want to keep the information of the order of the many,
   * you may have to use an association with a data member (the third template
   * parameter that we pretent not to exist).
   *
   * Use this if the objects that have to be associated to the one of type T are
   * sparse, spread across a to-be-data-product, but you have a list of the
   * indices in the data product of the elements to associate to the one of type
   * T.
   * In other words, given that you have a data product "a" of type
   * `std::vector<T>` and a data product "b" of type `std::vector<U>`, this
   * `method creates an association between `a[first_index]` and
   * `b[*(from_second_index)], another between `a[first_index]` and
   * `b[*(from_second_index + 1)]`, etc.
   *
   * The surprising concept here is that you don't need to specify neither of
   * the collections of T or U elements. The data product is uniquely defined
   * by its type, producer, process and product label.
   * Here we assume that the type of the products are `std::vector<T>` and
   * `std::vector<U>`, and that the products have empty product labels,
   * and that the producer is prod for both of them.
   */
  // MARK CreateAssn_08
  template <typename T, typename U, typename Iter>
  bool CreateAssn(art::Event& evt,
                  art::Assns<T, U>& assn,
                  size_t first_index,
                  Iter from_second_index,
                  Iter to_second_index);

  //@{
  /**
   * @brief Creates a single one-to-one association with associated data
   * @tparam T type of the new object to associate
   * @tparam U type of the many objects already in the data product or art::Ptr
   * @tparam D type of the "metadata" coupled to this pair association
   * @param evt reference to the current event
   * @param assn reference to association object where the new one will be put
   * @param first_index index of the object of type T to be associated
   * @param second_index index of the object of type U to be associated
   * @param data "metadata" to be store in this association
   * @return whether the operation was successful (can it ever fail??)
   *
   * Use this if you want some metadata to travel together with the association.
   * An example may be the order of the second element within a list:
   *
   *     size_t a_index = 2;
   *     std::vector<size_t> b_indices{ 6, 9, 18, 12 };
   *     for (size_t i = 0; i < b_indices.size(); ++i)
   *       CreateAssn(prod, evt, assn, a_index, b_index[i], i);
   *
   * In this way, the association between the element #2 of "a" (a vector that
   * is not specified -- nor needed -- in this snippet of code) and the element
   * #18 will be remembered as being the third (metadata value of 2).
   * In this example metadata is of type `size_t` the association would be
   * declared as `art::Assn<A, B, size_t>`.
   * A FindMany query of that association might look like:
   *
   *     art::Handle<std::vector<A>> a_list; // read this from the event
   *
   *     art::FindMany<B, size_t> Query(a_list, event, ModuleLabel);
   *
   *     // search for the last of the objects associated to the third element:
   *     size_t a_index = 2; // this means third element
   *
   *     std::vector<size_t const*> const& metadata = Query.data(a_index);
   *     size_t largest_index = 0, last_item = 0;
   *     for (size_t iB = 0; iB < metadata.size(); ++iB) {
   *       if (largest_index >= *(metadata[iB])) continue;
   *       largest_index = *(metadata[iB]);
   *       last_item = iB;
   *     } // for iB
   *     B const& lastB = Query.at(last_item);
   *
   * In alternative, the elements and their metadata can be fetched
   * simultaneously with:
   *
   *     std::vector<art::Ptr<B>> const& Bs;
   *     std::vector<size_t const*> const& metadata;
   *
   *     size_t a_index = 2; // this means third element
   *     size_t nMatches = Query.get(a_index, Bs, metadata);
   *
   */
  // MARK CreateAssnD_01
  // MARK CreateAssnD_01a
  template <typename T, typename U, typename D>
  bool CreateAssnD(art::Event& evt,
                   art::Assns<T, U, D>& assn,
                   size_t first_index,
                   size_t second_index,
                   typename art::Assns<T, U, D>::data_t&& data);
  // MARK CreateAssnD_01b
  template <typename T, typename U, typename D>
  bool CreateAssnD(art::Event& evt,
                   art::Assns<T, U, D>& assn,
                   size_t first_index,
                   size_t second_index,
                   typename art::Assns<T, U, D>::data_t const& data);
  //@}

  // method to return all objects of type U that are not associated to
  // objects of type T. Label is the module label that would have produced
  // the associations and likely the objects of type T
  // this method assumes there is a one to many relationship between T and U
  // for example if you want to get all recob::Hits
  // that are not associated to recob::Clusters
  // std::vector<const recob::Hit*> hits = FindUNotAssociatedToU<recob::Cluster>(art::Handle<recob::Hit>, ...);
  template <class T, class U>
  std::vector<const U*> FindUNotAssociatedToT(art::Handle<U> b,
                                              art::Event const& evt,
                                              std::string const& label);

  // method to return all objects of type U that are not associated to
  // objects of type T. Label is the module label that would have produced
  // the associations and likely the objects of type T
  // this method assumes there is a one to many relationship between T and U
  // for example if you want to get all recob::Hits
  // that are not associated to recob::Clusters
  // std::vector< art::Ptr<recob::Hit> > hits = FindUNotAssociatedToTP<recob::Cluster>(art::Handle<recob::Hit>, ...);
  template <class T, class U>
  std::vector<art::Ptr<U>> FindUNotAssociatedToTP(art::Handle<U> b,
                                                  art::Event const& evt,
                                                  std::string const& label);

  // Methods make getting simple ART-independent association information.
  // --- GetAssociatedVectorOneI takes in a handle to an association, and a handle to a product on the event.
  //     The ouput is a vector of with the same number of entries as the handle to the product, containing an index
  //     to the location of one associated product in that product's collection.
  // --- GetAssociatedVectorOneP takes in a handle to an association, and a handle to a product on the event.
  //     The ouput is a vector of with the same number of entries as the handle to the product, containing a pointer
  //     to one associated product.
  // --- GetAssociatedVectorManyI takes in a handle to an association, and a handle to a product on the event.
  //     The ouput is a vector of with the same number of entries as the handle to the product, containing a vector
  //     of indices that give the locations of all associated products in those products' collection.
  // --- GetAssociatedVectorManyP takes in a handle to an association, and a handle to a product on the event.
  //     The ouput is a vector of with the same number of entries as the handle to the product, containing a vector
  //     of pointers to all associated products.

  template <class T, class U>
  std::vector<size_t> GetAssociatedVectorOneI(art::Handle<art::Assns<T, U>> h,
                                              art::Handle<std::vector<T>> index_p);
  template <class T, class U>
  std::vector<const U*> GetAssociatedVectorOneP(art::Handle<art::Assns<T, U>> h,
                                                art::Handle<std::vector<T>> index_p);

  template <class T, class U>
  std::vector<std::vector<size_t>> GetAssociatedVectorManyI(art::Handle<art::Assns<T, U>> h,
                                                            art::Handle<std::vector<T>> index_p);
  template <class T, class U>
  std::vector<std::vector<const U*>> GetAssociatedVectorManyP(art::Handle<art::Assns<T, U>> h,
                                                              art::Handle<std::vector<T>> index_p);

} // end namespace

//----------------------------------------------------------------------
// MARK CreateAssn_01
template <class T, class U>
bool util::CreateAssn(art::Event& evt,
                      std::vector<T> const& a,
                      art::Ptr<U> const& b,
                      art::Assns<U, T>& assn,
                      std::string a_instance,
                      size_t index /* = UINT_MAX */
)
{
  if (index == UINT_MAX) index = a.size() - 1;

  try {
    assn.addSingle(b, art::PtrMaker<T>{evt, a_instance}(index));
    return true;
  }
  catch (cet::exception& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

} // util::CreateAssn() [01]

//----------------------------------------------------------------------
// MARK CreateAssn_03
template <class T, class U>
bool util::CreateAssn(art::Event&,
                      art::Ptr<T> const& a,
                      art::Ptr<U> const& b,
                      art::Assns<U, T>& assn)
{

  try {
    assn.addSingle(b, a);
  }
  catch (cet::exception& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

  return true;
} // util::CreateAssn() [03]

//----------------------------------------------------------------------
// MARK CreateAssn_04
template <class T, class U>
bool util::CreateAssn(art::Event& evt,
                      std::vector<T> const& a,
                      art::PtrVector<U> const& b,
                      art::Assns<T, U>& assn,
                      size_t index /* = UINT_MAX */
)
{
  if (index == UINT_MAX) index = a.size() - 1;

  try {
    auto const aptr = art::PtrMaker<T>{evt}(index);
    for (art::Ptr<U> const& b_item : b)
      assn.addSingle(aptr, b_item);
  }
  catch (cet::exception& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

  return true;
} // util::CreateAssn() [04]

//----------------------------------------------------------------------
// MARK CreateAssn_05
template <class T, class U>
bool util::CreateAssn(art::Event&,
                      art::Ptr<T> const& a,
                      std::vector<art::Ptr<U>> const& b,
                      art::Assns<T, U>& assn)
{

  try {
    for (art::Ptr<U> const& b_item : b)
      assn.addSingle(a, b_item);
  }
  catch (cet::exception const& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

  return true;
} // util::CreateAssn() [05]

//----------------------------------------------------------------------
// MARK CreateAssn_06
template <class T, class U>
bool util::CreateAssn(art::Event& evt,
                      std::vector<T> const& a,
                      std::vector<art::Ptr<U>> const& b,
                      art::Assns<T, U>& assn,
                      size_t index /* = UINT_MAX */
)
{

  if (index == UINT_MAX) index = a.size() - 1;

  try {
    auto const aptr = art::PtrMaker<T>{evt}(index);
    for (art::Ptr<U> const& b_item : b)
      assn.addSingle(aptr, b_item);
  }
  catch (cet::exception& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

  return true;
} // util::CreateAssn() [06]

//----------------------------------------------------------------------
// MARK CreateAssn_07
template <class T, class U>
bool util::CreateAssn(art::Event& evt,
                      std::vector<T> const& a,
                      std::vector<U> const& /* b */,
                      art::Assns<T, U>& assn,
                      size_t startU,
                      size_t endU,
                      size_t index /* = UINT_MAX */
)
{

  if (index == UINT_MAX) index = a.size() - 1;

  try {
    auto const aptr = art::PtrMaker<T>{evt}(index);
    art::PtrMaker<U> const make_bptr{evt};
    for (size_t i = startU; i < endU; ++i) {
      assn.addSingle(aptr, make_bptr(i));
    }
  }
  catch (cet::exception& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

  return true;
} // util::CreateAssn() [07]

//----------------------------------------------------------------------
// MARK CreateAssn_07a
template <class T, class U>
bool util::CreateAssn(art::Event& evt,
                      std::vector<T> const& a,
                      std::vector<U> const& /* b */,
                      art::Assns<T, U>& assn,
                      std::vector<size_t> const& indices,
                      size_t index /* = UINT_MAX */
)
{

  if (index == UINT_MAX) index = a.size() - 1;

  try {
    auto const aptr = art::PtrMaker<T>{evt}(index);
    art::PtrMaker<U> const make_bptr{evt};
    for (size_t index : indices) {
      assn.addSingle(aptr, make_bptr(index));
    }
  }
  catch (cet::exception& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

  return true;
} // util::CreateAssn() [07a]

//----------------------------------------------------------------------
// MARK CreateAssn_08
template <typename T, typename U, typename Iter>
bool util::CreateAssn(art::Event& evt,
                      art::Assns<T, U>& assn,
                      size_t first_index,
                      Iter from_second_index,
                      Iter to_second_index)
{

  try {
    // we declare here that we want to associate the element first_index of the
    // (only) data product of type std::vector<T> with other objects.
    // This is the pointer to that element:
    auto const first_ptr = art::PtrMaker<T>{evt}(first_index);

    // we are going to associate that element in a with a number of elements
    // of the only data product of type std::vector<U>
    art::PtrMaker<U> const make_second_ptr{evt};
    for (; from_second_index != to_second_index; ++from_second_index) {
      assn.addSingle(first_ptr, make_second_ptr(*from_second_index));
    } // while
  }
  catch (cet::exception& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

  return true;
} // util::CreateAssn() [08]

//----------------------------------------------------------------------
// MARK CreateAssnD_01a
template <typename T, typename U, typename D>
bool util::CreateAssnD(art::Event& evt,
                       art::Assns<T, U, D>& assn,
                       size_t first_index,
                       size_t second_index,
                       typename art::Assns<T, U, D>::data_t&& data)
{

  try {
    // we declare here that we want to associate the element first_index of the
    // (only) data product of type std::vector<T> with the other object
    auto const first_ptr = art::PtrMaker<T>{evt}(first_index);

    // the same to associate the element second_index of the (only)
    // data product of type std::vector<U> with the first object.
    auto const second_ptr = art::PtrMaker<U>{evt}(second_index);

    assn.addSingle(first_ptr, second_ptr, std::move(data));
  }
  catch (cet::exception& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

  return true;
} // util::CreateAssnD() [01a]

template <typename T, typename U, typename D>
bool util::CreateAssnD(art::Event& evt,
                       art::Assns<T, U, D>& assn,
                       size_t first_index,
                       size_t second_index,
                       typename art::Assns<T, U, D>::data_t const& data)
{

  try {
    // we declare here that we want to associate the element first_index of the
    // (only) data product of type std::vector<T> with the other object
    auto const first_ptr = art::PtrMaker<T>{evt}(first_index);

    // the same to associate the element second_index of the (only)
    // data product of type std::vector<U> with the first object.
    auto const second_ptr = art::PtrMaker<U>{evt}(second_index);

    assn.addSingle(first_ptr, second_ptr, data);
  }
  catch (cet::exception& e) {
    mf::LogWarning("AssociationUtil")
      << "unable to create requested art:Assns, exception thrown: " << e;
    return false;
  }

  return true;
} // util::CreateAssnD() [01b]

//----------------------------------------------------------------------
template <class T, class U>
inline std::vector<const U*> util::FindUNotAssociatedToT(art::Handle<U> b,
                                                         art::Event const& evt,
                                                         std::string const& label)
{
  // Do a FindOne for type T for each object of type U
  // If the FindOne returns an invalid maybe ref, add the pointer
  // of object type U to the return vector

  std::vector<const U*> notAssociated;

  art::FindOne<T> const fa(b, evt, label);

  for (size_t u = 0; u < b->size(); ++u) {
    cet::maybe_ref<T const> t(fa.at(u));
    if (!t.isValid()) {
      art::Ptr<U> ptr(b, u);
      notAssociated.push_back(ptr.get());
    }
  }

  return notAssociated;
}

//----------------------------------------------------------------------
template <class T, class U>
inline std::vector<art::Ptr<U>> util::FindUNotAssociatedToTP(art::Handle<U> b,
                                                             art::Event const& evt,
                                                             std::string const& label)
{
  // Do a FindOneP for type T for each object of type U
  // If the FindOne returns an invalid maybe ref, add the pointer
  // of object type U to the return vector

  std::vector<art::Ptr<U>> notAssociated;

  art::FindOneP<T> const fa(b, evt, label);

  for (size_t u = 0; u < b->size(); ++u) {
    cet::maybe_ref<T const> t(fa.at(u));
    if (!t.isValid()) { notAssociated.emplace_back(b, u); }
  }

  return notAssociated;
}

template <class T, class U>
inline std::vector<size_t> util::GetAssociatedVectorOneI(art::Handle<art::Assns<T, U>> h,
                                                         art::Handle<std::vector<T>> index_p)
{
  std::vector<size_t> associated_index(index_p->size());
  for (auto const& pair : *h)
    associated_index.at(pair.first.key()) = pair.second.key();
  return associated_index;
}

template <class T, class U>
inline std::vector<const U*> util::GetAssociatedVectorOneP(art::Handle<art::Assns<T, U>> h,
                                                           art::Handle<std::vector<T>> index_p)
{
  std::vector<const U*> associated_pointer(index_p->size());
  for (auto const& pair : *h)
    associated_pointer.at(pair.first.key()) = &(*(pair.second));
  return associated_pointer;
}

template <class T, class U>
inline std::vector<std::vector<size_t>> util::GetAssociatedVectorManyI(
  art::Handle<art::Assns<T, U>> h,
  art::Handle<std::vector<T>> index_p)
{
  std::vector<std::vector<size_t>> associated_indices(index_p->size());
  for (auto const& pair : *h)
    associated_indices.at(pair.first.key()).push_back(pair.second.key());
  return associated_indices;
}

template <class T, class U>
inline std::vector<std::vector<const U*>> util::GetAssociatedVectorManyP(
  art::Handle<art::Assns<T, U>> h,
  art::Handle<std::vector<T>> index_p)
{
  std::vector<std::vector<const U*>> associated_pointers(index_p->size());
  for (auto const& pair : *h)
    associated_pointers.at(pair.first.key()).push_back(&(*(pair.second)));
  return associated_pointers;
}

//--------------------------------------------------------------------
// Functions to support unnecessary leading producer argument
//
//   These are legacy function signatures that accept a reference to a
//   producer module.  Although the signatures are still supported,
//   they should not be encouraged.  They, along with the other
//   CreateAssn(D) functions above, will be marked as [[deprecated]]
//   once art supports Assns::addMany (expected in art 3.04).
//--------------------------------------------------------------------

#include "art/Framework/Core/EDProducer.h"

namespace util {
  template <typename Producer, typename... Args>
  std::enable_if_t<std::is_base_of_v<art::EDProducer, Producer>, bool> CreateAssn(Producer const&,
                                                                                  Args&&... args)
  {
    return CreateAssn(std::forward<Args>(args)...);
  }

  template <typename Producer, typename... Args>
  std::enable_if_t<std::is_base_of_v<art::EDProducer, Producer>, bool> CreateAssnD(Producer const&,
                                                                                   Args&&... args)
  {
    return CreateAssnD(std::forward<Args>(args)...);
  }
}

#endif //ASSOCIATIONUTIL_H
