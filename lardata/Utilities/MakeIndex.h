/**
 * @file   MakeIndex.h
 * @brief  Procedures to create maps of object locations
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 13th, 2015
 */

#ifndef LARDATA_UTILITIES_MAKEINDEX_H
#define LARDATA_UTILITIES_MAKEINDEX_H

namespace util {

  /**
   * @brief Creates a map of indices from an existing collection
   * @tparam Coll type of the collection
   * @tparam KeyOf type of the extractor of the key
   * @param data the data collection
   * @param key_of instance of a functor extracting a key value from a datum
   * @return a vector with indices corresponding to the data keys
   *
   * This function maps the index of the items in data to an integral key
   * extracted from each item.
   * For example, if the items are wires and the key_of function extracts their
   * channel ID, the resulting vector will contain for each channel ID
   * the index in data of the wire with that channel ID.
   *
   * The key is converted into a unsigned integer (`size_t`).
   * If multiple items have the same key, the outcome for that key is undefined.
   * If no items has a specific key, the index of that key is assigned as
   * @code std::numeric_limits<size_t>::max() @endcode, i.e. an index larger
   * than the size of the original data collection.
   *
   * The returned vector is big enough to accommodate indices corresponding to
   * the keys of all the items in data. It may contain "holes" (that is, some
   * keys that have no corresponding items have a
   * @code std::numeric_limits<size_t>::max() @endcode value).
   * The memory allocated for the vector may be larger than necessary (if that
   * is a problem, `std::vector::shrink_to_fit()` can be used, but it may create
   * more problems than it solves).
   *
   */
  template <typename Coll, typename KeyOf>
  std::vector<size_t> MakeIndex(Coll const& data, KeyOf key_of = KeyOf())
  {

    // we start the index with the best guess that all the items will have
    // a unique key and they are contiguous:
    // the index would have the same size as the data
    std::vector<size_t> Index(data.size(), std::numeric_limits<size_t>::max());

    size_t min_size = 0; // minimum size needed to hold all keys

    size_t iDatum = 0;
    for (auto const& datum : data) {
      size_t key = size_t(key_of(datum));
      if (key >= min_size) min_size = key + 1;
      if (Index.size() <= key) {
        // make room for the entry: double the size
        Index.resize(std::max(key + 1, Index.size() * 2), std::numeric_limits<size_t>::max());
      } // if expand index
      Index[key] = iDatum;
      ++iDatum;
    } // for datum
    Index.resize(min_size);
    return Index;
  } // MakeIndex()

  /**
   * @brief Creates a map of objects from an existing collection
   * @tparam Coll type of the collection
   * @tparam KeyOf type of the extractor of the key
   * @param data the data collection
   * @param key_of instance of a functor extracting a key value from a datum
   * @return a vector with pointers to data corresponding to their keys
   *
   * This function maps the items in data to an integral key extracted from each
   * of them.
   * For example, if the items are wires and the key_of function extracts their
   * channel ID, the resulting vector will contain for each channel ID
   * the pointer to the wire with that channel ID.
   *
   * The key is converted into a unsigned integer (`size_t`).
   * If multiple items have the same key, the outcome for that key is undefined.
   * If no items has a specific key, the index of that key is assigned a
   * null pointer.
   *
   * The returned vector is big enough to accommodate pointers corresponding to
   * the keys of all the items in data. It may contain "holes" (that is, some
   * keys that have no corresponding items have a null pointer value).
   * The memory allocated for the vector may be larger than necessary (if that
   * is a problem, `std::vector::shrink_to_fit()` can be used, but it may create
   * more problems than it solves).
   *
   */
  template <typename Coll, typename KeyOf>
  auto MakeMap(Coll const& data, KeyOf key_of = KeyOf())
    -> std::vector<decltype(key_of(*(data.begin()))) const*>
  {
    using Mapped_t = decltype(key_of(*(data.begin())));
    using Ptr_t = Mapped_t const*;
    using Map_t = std::vector<Ptr_t>;

    // we start the index with the best guess that all the items will have
    // a unique key and they are contiguous:
    // the index would have the same size as the data
    Map_t Index(data.size(), nullptr);

    size_t min_size = 0; // minimum size needed to hold all keys

    for (auto const& datum : data) {
      size_t key = size_t(key_of(datum));
      if (key >= min_size) min_size = key + 1;
      if (Index.size() <= key) {
        // make room for the entry: double the size
        Index.resize(std::max(key + 1, Index.size() * 2), nullptr);
      } // if expand index
      Index[key] = &datum;
    } // for datum
    Index.resize(min_size);
    return Index;
  } // MakeMap()

} // namespace util

#endif // LARDATA_UTILITIES_MAKEINDEX_H
