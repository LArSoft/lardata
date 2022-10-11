/**
 * @file   CountersMap.h
 * @brief  Map of counters, stored compactly
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   August 22th, 2014
 */

#ifndef COUNTERSMAP_H
#define COUNTERSMAP_H

// interface include
#include <array>
#include <cstddef>    // std::ptrdiff_t
#include <functional> // std::less<>
#include <iterator>   // std::bidirectional_iterator_tag
#include <map>
#include <memory>      // std::allocator<>
#include <type_traits> // std::make_unsigned
#include <utility>     // std::pair<>

namespace lar {

  /// Returns true if the argument is a power of 2
  constexpr bool IsPowerOf2(unsigned long long int v)
  {
    return v & 1 ? v == 1 : IsPowerOf2(v >> 1);
  }

  /// Returns the position of the first set bit (0 for LSB)
  constexpr int LowestSetBit(unsigned long long int v);

  namespace details {

    /// Type of block of counters (just a STL array until SUBCOUNTERS are in)
    template <typename COUNTER, std::size_t NCounters>
    class CounterBlock : public std::array<COUNTER, NCounters> {
    public:
      using Counter_t = COUNTER;

      using Array_t = std::array<Counter_t, NCounters>; ///< type of base class
      using value_type = typename Array_t::value_type;

      /// Default constructor: initializes the array to 0
      CounterBlock() : Array_t() { Array_t::fill(0); }

      /// Convenience constructor: initializes all counters to 0, except one
      CounterBlock(size_t index, Counter_t value) : CounterBlock()
      {
        Array_t::operator[](index) = value;
      }

      void fill(const value_type& value) { Array_t::fill(value); }
      void fill(size_t begin, size_t n, const value_type& value)
      {
        std::fill(Array_t::data() + begin, Array_t::data() + begin + n, value);
      }

    }; // CounterBlock

    template <typename KEY, typename COUNTER, size_t SIZE>
    struct CountersMapTraits {

      using Key_t = KEY;         ///< Type of counter key in the map.
      using Counter_t = COUNTER; ///< Type of the single counter.

      /// Number of counters in one counter block.
      static constexpr size_t NCounters = SIZE;

      /// Type of counter block actually stored.
      using CounterBlock_t = CounterBlock<Counter_t, NCounters>;

      /// General type of map (no special allocator specified).
      using PlainBaseMap_t = std::map<Key_t, CounterBlock_t, std::less<Key_t>>;

      /// Type of allocator for the plain map.
      using DefaultAllocator_t = typename PlainBaseMap_t::allocator_type;

      /// Base type of map, allowing a custom allocator.
      template <typename Alloc>
      using BaseMap_t = std::map<Key_t, CounterBlock_t, std::less<Key_t>, Alloc>;

      /// Type of value in the map.
      using MapValue_t = typename PlainBaseMap_t::value_type;

    }; // struct CountersMapTraits

  } // namespace details

  /**
   * @brief Map storing counters in a compact way
   * @param KEY the type of the key of the counters map
   * @param COUNTER the type of a basic counter (can be signed or unsigned)
   * @param BLOCKSIZE the number of counters in a cluster
   * @param ALLOC allocator for the underlying STL map
   * @param SUBCOUNTERS split each counter in subcounters (not implemented yet)
   *
   * This class is designed for the need of a vast number of counters with
   * a integral numerical key, when the counter keys are usually clustered.
   * This can be more or less efficient than a straight counters STL map
   * according to how often the counters are clustered (if that does not happen
   * often, CountersMap can have considerable memory overhead).
   *
   * Counters are allocated in contiguous blocks of BLOCKSIZE counters.
   * The selling point here is that a map node has some overhead (typically
   * at least 3 pointers) and dynamically allocating it costs a lot (40 bytes
   * have been observed). If you need a counter with a range of 100 (1 byte),
   * that's far from optimal in many aspects (memory allocation also takes
   * time).
   *
   * The requirement of the key to be numerical is so that the concept of
   * "next counter" is well defined and we can store contiguous counters
   * in a fixed structure.
   *
   * <h3>Subcounters</h3>
   * The idea behind subcounters is that you migt want to split a counter into
   * subcounters to save memory if the maximum counter value is smaller than
   * the range of the counter type.
   * <strong>The implementation is delayed</strong> since in the end the same
   * effect can be achieved by using a small counter type (e.g. signed char),
   * unless the range is smaller that 16 (or 4, or 2), in which case the
   * character can be split into bits. That will cause some overhead for
   * increment and decrement instructions.
   */
  template <typename KEY,
            typename COUNTER,
            size_t SIZE,
            typename ALLOC =
              typename details::CountersMapTraits<KEY, COUNTER, SIZE>::DefaultAllocator_t,
            unsigned int SUBCOUNTERS = 1>
  class CountersMap {
    static_assert(SUBCOUNTERS == 1, "subcounters not implemented yet");
    static_assert(IsPowerOf2(SIZE), "the size of the cluster of counters must be a power of 2");

    /// Set of data types pertaining this counter.
    using Traits_t = details::CountersMapTraits<KEY, COUNTER, SIZE>;

  public:
    using Key_t = KEY;         ///< type of counter key in the map
    using Counter_t = COUNTER; ///< type of the single counter
    using Allocator_t = ALLOC; ///< type of the single counter

    /// This class
    using CounterMap_t = CountersMap<KEY, COUNTER, SIZE, ALLOC, SUBCOUNTERS>;

    /// Number of counters in one counter block
    static constexpr size_t NCounters = Traits_t::NCounters;

    /// Number of subcounters in one counter block
    static constexpr size_t NSubcounters = NCounters * SUBCOUNTERS;

    /// Type of the subcounter (that is, the actual counter)
    using SubCounter_t = Counter_t;

    using CounterBlock_t = typename Traits_t::CounterBlock_t;

    /// Type of the map used in the implementation
    using BaseMap_t = typename Traits_t::template BaseMap_t<typename std::allocator_traits<
      Allocator_t>::template rebind_alloc<typename Traits_t::MapValue_t>>;

    /*
    /// Iterator through the allocated elements
    using const_iterator = double_fwd_const_iterator<
      typename BaseMap_t::const_iterator,
      PairSecond<typename BaseMap_t::value_type>
      >;
    */

    ///@{
    /// @name STL type definitions
    using mapped_type = SubCounter_t;
    using allocator_type = Allocator_t; ///< type of the single counter
    // TODO others are missing
    ///@}

    using value_type = std::pair<const Key_t, SubCounter_t>;

    /// Iterator through the counters (shown as value_type)
    class const_iterator;

    /// Default constructor (empty map)
    CountersMap() {}

    /// Constructor, specifies an allocator
    CountersMap(Allocator_t alloc) : BaseMap_t(alloc) {}

    /// Read-only access to an element; returns 0 if no counter is present
    SubCounter_t operator[](Key_t key) const { return GetSubCounter(key); }

    /**
     * @brief Sets the specified counter to a count
     * @param key key of the counter to be increased
     * @param value count value to be set
     * @return new value of the counter
     *
     * Given the complex underlying structure (especially when subcounters are
     * involved), a special method is used to set the value of a counter.
     * This is equivalent to map's operator[] in non-constant version, but the
     * creation of a new counter is explicit.
     */
    SubCounter_t set(Key_t key, SubCounter_t value);

    /**
     * @brief Increments by 1 the specified counter
     * @param key key of the counter to be increased
     * @return new value of the counter
     */
    SubCounter_t increment(Key_t key);

    /**
     * @brief Decrements by 1 the specified counter
     * @param key key of the counter to be decreased
     * @return new value of the counter
     */
    SubCounter_t decrement(Key_t key);

    ///@{
    /// @name Iterators (experimental)
    // this stuff needs serious work to be completed

    /// Returns an iterator to the begin of the counters
    const_iterator begin() const;

    /// Returns an iterator past-the-end of the counters
    const_iterator end() const;
    ///@}

    /// Returns whether the map has no counters
    bool empty() const { return counter_map.empty(); }

    /// Returns the number of allocated counters
    typename std::make_unsigned<Key_t>::type n_counters() const
    {
      return counter_map.size() * NSubcounters;
    }

    /**
     * @brief Returns whether the counters in this map are equivalent to another
     * @param to a STL map
     * @param first_difference key of the first difference
     * @return whether the counters are equivalent
     *
     * The maps are considered equal if all keys in each are present in the
     * other, and their counters have the same value, with one exception:
     * counters in CountersMap can not to exist in the other map, but only if
     * their count is 0 (these should be plenty of them since counters are
     * created in blocks).
     *
     * If there is a key in one map but not in the other, first_difference
     * stores that key; if a counter is present in both maps but with a
     * different count, first_difference reports the key of that counter.
     * If no difference is found, first_difference is left untouched.
     */
    template <typename OALLOC>
    bool is_equal(const std::map<Key_t, SubCounter_t, std::less<Key_t>, OALLOC>& to,
                  Key_t& first_difference) const;

    /**
     * @brief Returns whether the counters in this map are equivalent to another
     * @param to a STL map
     * @return whether the counters are equivalent
     */
    template <typename OALLOC>
    bool is_equal(const std::map<Key_t, SubCounter_t, std::less<Key_t>, OALLOC>& to) const
    {
      Key_t dummy;
      return is_equal(to, dummy);
    }

  protected:
    using BlockKey_t = Key_t; ///< type of block key
    using CounterIndex_t = typename std::make_unsigned<Key_t>::type;
    ///< type of index in the block

    /// Structure with the index of the counter, split as needed
    struct CounterKey_t {
      BlockKey_t block;       ///< key of the counter block
      CounterIndex_t counter; ///< index of the counter in the block
      // there should be a subcluster index here too

      CounterKey_t() : block(0), counter(0) {}

      /// Constructor from a pair
      CounterKey_t(BlockKey_t major, CounterIndex_t minor) : block(major), counter(minor) {}

      /// Initialize from a mangled key
      CounterKey_t(Key_t key)
        : //    CounterKey_t(key >> MinorKeyBits, key & MinorKeyMask) {}
        CounterKey_t(key & ~MinorKeyMask, key & MinorKeyMask)
      {}

      /// Returns the full key
      //  Key_t Key() const { return (block << MinorKeyBits) + counter; }
      Key_t Key() const { return block + counter; }

      /// Conversion of this key into a number (@see Key())
      operator Key_t() const { return Key(); }

      CounterKey_t& operator++()
      {
        if (++counter == MinorKeyRange) next_block();
        return *this;
      }
      CounterKey_t& operator--()
      {
        if (counter-- == 0) {
          prev_block();
          counter = MinorKeyRange - 1;
        }
        return *this;
      } // operator--()
      CounterKey_t& operator++(int)
      {
        CounterKey_t old(*this);
        this->operator++();
        return old;
      }
      CounterKey_t& operator--(int)
      {
        CounterKey_t old(*this);
        this->operator--();
        return old;
      }

      /// Skips to the beginning of this block
      CounterKey_t& start_block()
      {
        counter = 0;
        return *this;
      }

      /// Skips to the beginning of the previous block
      CounterKey_t& prev_block()
      {
        block -= MinorKeyRange;
        return start_block();
      }

      /// Skips to the beginning of the next block
      CounterKey_t& next_block()
      {
        block += MinorKeyRange;
        return start_block();
      }

      /// Number of values of the minor key
      static constexpr Key_t MinorKeyRange = NSubcounters;

      /// Number of bits for the minor key
      static constexpr Key_t MinorKeyBits = LowestSetBit(MinorKeyRange);

      /// Bit mask for the minor key
      static constexpr Key_t MinorKeyMask = MinorKeyRange - 1;

      static_assert(MinorKeyBits > 0, "Wrong COUNTER value for lar::CountersMap");

    }; // struct CounterKey_t

    BaseMap_t counter_map; ///< the actual data structure for counters

    /// Returns the value of the counter at the specified split key
    Counter_t GetCounter(CounterKey_t key) const;

    /// Returns the value of the subcounter at the specified split key
    // Since subcounters are not implemented yet, this is a no-op
    SubCounter_t GetSubCounter(CounterKey_t key) const { return GetCounter(key); }

    /// Returns the counter at the specified split key
    Counter_t& GetOrCreateCounter(CounterKey_t key);

    /// Returns a split key corresponding to the specified key
    static CounterKey_t SplitKey(Key_t key) { return key; }

    /// Creates a const_iterator (useful in derived classes)
    static const_iterator make_const_iterator(typename BaseMap_t::const_iterator it, size_t ix)
    {
      return {it, ix};
    }

  private:
    /// Sets the specified counter to a value (no check on value range)
    SubCounter_t unchecked_set(CounterKey_t key, SubCounter_t delta);

    /// Adds a delta to the specified counter (no check on underflow/overflow)
    SubCounter_t unchecked_add(CounterKey_t key, SubCounter_t delta);

  }; // class CountersMap

} // namespace lar

//------------------------------------------------------------------------------

namespace lar {

  namespace details {
    /// Internally used by LowestSetBit
    inline constexpr int LowestSetBitScaler(unsigned long long int v, int b)
    {
      return (v & 1) ? b : LowestSetBitScaler(v >> 1, b + 1);
    }

    namespace counters_map {
      static constexpr bool bDebug = true;

    } // namespace counters_map
  }   // namespace details

  inline constexpr int LowestSetBit(unsigned long long int v)
  {
    return (v == 0) ? -1 : details::LowestSetBitScaler(v, 0);
  }

  // providing "definitions" of the static constants;
  // these are required if the code is going to take the address of these
  // constants, which is most often due to the use of "const size_t&" in some
  // function, where "size_t" is a template argument type
  // (or else the reference would not be needed).
  // I am not providing the same for the protected and private members
  // (just because of laziness).
  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  constexpr size_t CountersMap<K, C, S, A, SUB>::NCounters;

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  constexpr size_t CountersMap<K, C, S, A, SUB>::NSubcounters;

  // CountersMap<>::const_iterator does not fully implement the STL iterator
  // interface, since it does not implement operator-> () (for technical reason:
  // the value does not actually exist and it does not have an address),
  // in addition to std::swap().
  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  class CountersMap<K, C, S, A, SUB>::const_iterator : public std::bidirectional_iterator_tag {
    friend CountersMap<K, C, S, A, SUB>;

  public:
    using value_type = typename CounterMap_t::value_type; ///< value type: pair
    using difference_type = std::ptrdiff_t;
    using pointer = const value_type*;
    using reference = const value_type&;
    using iterator_category = std::bidirectional_iterator_tag;

    using iterator_type = CounterMap_t::const_iterator; ///< this type

    /// Default constructor
    const_iterator() = default;

    /// Access to the pointed pair
    value_type operator*() const { return {key(), iter->second[index]}; }

    iterator_type& operator++()
    {
      if (++index >= NSubcounters) {
        ++iter;
        index = 0;
      }
      return *this;
    } // operator++
    iterator_type operator++(int)
    {
      iterator_type old(*this);
      this->operator++();
      return old;
    }

    iterator_type& operator--()
    {
      if (index == 0) {
        --iter;
        index = NSubcounters - 1;
      }
      else
        --index;
      return *this;
    } // operator--()
    iterator_type operator--(int)
    {
      iterator_type old(*this);
      this->operator--();
      return old;
    }

    bool operator==(const iterator_type& as) const
    {
      return (iter == as.iter) && (index == as.index);
    }
    bool operator!=(const iterator_type& as) const
    {
      return (iter != as.iter) || (index != as.index);
    }

    /// Returns the key of the pointed item as a CounterKey_t
    CounterKey_t key() const { return {iter->first, index}; }

  protected:
    typename BaseMap_t::const_iterator iter;
    ///< iterator to the block of counters
    CounterIndex_t index; ///< index of the counted in the subblock

    /// Private constructor (from a map iterator and a block index)
    const_iterator(typename BaseMap_t::const_iterator it, size_t ix) : iter(it), index(ix) {}

  }; // CountersMap<>::const_iterator

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  inline typename CountersMap<K, C, S, A, SUB>::SubCounter_t CountersMap<K, C, S, A, SUB>::set(
    Key_t key,
    SubCounter_t value)
  {
    return unchecked_set(CounterKey_t(key), value);
  }

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  inline typename CountersMap<K, C, S, A, SUB>::SubCounter_t
  CountersMap<K, C, S, A, SUB>::increment(Key_t key)
  {
    return unchecked_add(CounterKey_t(key), +1);
  }

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  inline typename CountersMap<K, C, S, A, SUB>::SubCounter_t
  CountersMap<K, C, S, A, SUB>::decrement(Key_t key)
  {
    return unchecked_add(CounterKey_t(key), -1);
  }

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  inline typename CountersMap<K, C, S, A, SUB>::const_iterator CountersMap<K, C, S, A, SUB>::begin()
    const
  {
    return const_iterator{counter_map.begin(), 0};
  }

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  inline typename CountersMap<K, C, S, A, SUB>::const_iterator CountersMap<K, C, S, A, SUB>::end()
    const
  {
    return const_iterator{counter_map.end(), 0};
  }

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  template <typename OALLOC>
  bool CountersMap<K, C, S, A, SUB>::is_equal(
    const std::map<Key_t, SubCounter_t, std::less<Key_t>, OALLOC>& to,
    Key_t& first_difference) const
  {
    // this function is not optimised
    using CompMap_t = const std::map<Key_t, SubCounter_t, std::less<Key_t>, OALLOC>;
    typename CompMap_t::const_iterator to_iter = to.begin(), to_end = to.end();
    for (const typename const_iterator::value_type& p : *this) {

      if (to_iter != to_end) { // we have still counters to find
        // if counter is already beyond the next non-empty one froml STL map,
        // then we are missing some
        if (p.first > to_iter->first) {
          first_difference = to_iter->first;
          return false;
        }
        //    while (p.first > to_iter->first) {
        //      ++nMissingKeys;
        //      std::cout << "ERROR missing key " << to_iter->first << std::endl;
        //      if (++to_iter == to_end) break;
        //    } // while
        //  } // if
        //
        //  if (to_iter != to_end) { // we have still counters to find
        if (p.first == to_iter->first) {
          // if the counter is in SLT map, the two counts must match
          //  std::cout << "  " << p.first << " " << p.second << std::endl;
          if (to_iter->second != p.second) {
            //  std::cout << "ERROR wrong counter value " << p.second
            //    << ", expected " << to_iter->second << std::endl;
            //  ++nMismatchValue;
            first_difference = to_iter->first;
            return false;
          }
          ++to_iter; // done with it
        }
        else if (p.first < to_iter->first) {
          // if the counter is not in STL map, then it must be 0
          if (p.second != 0) {
            //  ++nExtraKeys;
            //  std::cout << "ERROR extra key " << p.first << " (" << p.second << ")"
            //    << std::endl;
            first_difference = p.first;
            return false;
          }
          //  else {
          //    std::cout << "  " << p.first << " " << p.second << " (not in STL)"
          //      << std::endl;
          //  }
        }
      }
      else { // if we are at the end of compared map
        // no more keys in STL map
        if (p.second != 0) {
          //  ++nExtraKeys;
          //  std::cout << "ERROR extra key " << p.first << " (" << p.second << ")"
          //    << std::endl;
          first_difference = p.first;
          return false;
        }
      }
    } // for element in map

    return true;
  } // CountersMap<>::is_equal()

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  inline typename CountersMap<K, C, S, A, SUB>::Counter_t CountersMap<K, C, S, A, SUB>::GetCounter(
    CounterKey_t key) const
  {
    typename BaseMap_t::const_iterator iBlock = counter_map.find(key.block);
    return (iBlock == counter_map.end()) ? 0 : iBlock->second[key.counter];
  } // CountersMap<>::GetCounter() const

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  inline typename CountersMap<K, C, S, A, SUB>::Counter_t&
  CountersMap<K, C, S, A, SUB>::GetOrCreateCounter(CounterKey_t key)
  {
    return counter_map[key.block][key.counter];
  }

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  typename CountersMap<K, C, S, A, SUB>::SubCounter_t CountersMap<K, C, S, A, SUB>::unchecked_set(
    CounterKey_t key,
    SubCounter_t value)
  {
    // take it low level here
    // first get the iterator to the block
    typename BaseMap_t::iterator iBlock = counter_map.lower_bound(key.block);
    if (iBlock != counter_map.end()) {
      if (iBlock->first == key.block) { // sweet, we have that counter already
        return iBlock->second[key.counter] = value;
      }
    }
    // no counter there yet: create one;
    // hint to insert before the block in the position we jave already found
    // (this is optimal in STL map for C++11);
    // create a block using the special constructor;
    // the temporary value created here is moved to the map
    counter_map.insert(iBlock, {key.block, {key.counter, value}});
    return value;
  } // CountersMap<>::unchecked_add() const

  template <typename K, typename C, size_t S, typename A, unsigned int SUB>
  typename CountersMap<K, C, S, A, SUB>::SubCounter_t CountersMap<K, C, S, A, SUB>::unchecked_add(
    CounterKey_t key,
    SubCounter_t delta)
  {
    // take it low level here
    // first get the iterator to the block
    typename BaseMap_t::iterator iBlock = counter_map.lower_bound(key.block);
    if (iBlock != counter_map.end()) {
      if (iBlock->first == key.block) { // sweet, we have that counter already
        return iBlock->second[key.counter] += delta;
      }
    }
    // no counter there yet: create one;
    // hint to insert before the block in the position we jave already found
    // (this is optimal in STL map for C++11);
    // create a block using the special constructor;
    // the temporary value created here is moved to the map
    counter_map.insert(iBlock, {key.block, {key.counter, delta}});
    return delta;
  } // CountersMap<>::unchecked_add() const

} // namespace lar

#endif // BULKALLOCATOR_H
