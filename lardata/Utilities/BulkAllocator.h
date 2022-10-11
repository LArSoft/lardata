/**
 * @file   BulkAllocator.h
 * @brief  Memory allocator for large amount of (small) objects
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   August 17th, 2014
 */

//--- BEGIN issue #19494 -------------------------------------------------------
/// @bug BulkAllocator.h is currently broken; see issue #19494.
// We are leaving it here because, being a header, it will not bother unless
// explicitly invoked. Note that there is a unit test for it too.
#error("BulkAllocator.h is currently broken; see issue #19494.")
//--- END issue #19494 ---------------------------------------------------------

#ifndef BULKALLOCATOR_H
#define BULKALLOCATOR_H

// interface include
#include <cstdlib>   // std::free
#include <memory>    // std::allocator<>, std::unique_ptr<>
#include <stdexcept> // std::logic_error

namespace lar {
  /// Namespace hiding implementation details
  namespace details {
    /// Namespace specific to bulk allocator
    namespace bulk_allocator {
      template <typename T>
      class BulkAllocatorBase;
    } // namespace bulk_allocator
  }   // namespace details

  /// Exception thrown when BulkAllocator-specific allocation errors happen
  class memory_error : public std::bad_alloc {
  public:
    memory_error() : std::bad_alloc() {}
    memory_error(const char* message) : std::bad_alloc(), msg(message) {}

    virtual const char* what() const throw() override { return msg; }

  private:
    const char* msg = nullptr;
  }; // class memory_error

  /**
   * @brief Aggressive allocator reserving a lot of memory in advance
   * @param T type being allocated
   *
   * This allocator appropriates memory in large chunks of GetChunkSize()
   * elements of type T. The memory will never be deleted! (but read further)
   *
   * @note With C++17, an allocator called `std::pmr::monotonic_buffer_resource`
   *       is available that seems to have pretty much the same functionality as
   *       this one (but STL quality).
   *       When C++17 is adopted by LArSoft (and the supported compilers have a
   *       complete enough support for it), the users of `BulkAllocator` should
   *       migrate to that one. Note that the interface is different, and
   *       probably the way to use it is also different.
   *
   * <h3>Deletion policy</h3>
   *
   * This allocator does not release not reuse deallocated memory. This design
   * choice is meant to reflect a specific use case where a large amount of
   * elements is created and then used, and the created object is fairly static.
   * Tracking freed memory fragments takes time and more memory, and reusing
   * them too.
   * Nevertheless, the allocator has a user count; when no user is present,
   * all the memory is deallocated. This can be convenient, or disastrous:
   * remember that the elements of a container can (or just might) not survive
   * after the container is destroyed. Clearing the container will not trigger
   * this self-distruction; if you are completely sure that no other container
   * is currently using the same allocator, you can explicitly Free() its
   * memory.
   *
   * <h3>One allocator for them all</h3>
   *
   * Since STL containers do not necessarely store their allocator but they
   * may create it with a default constructor, allocators should be formally
   * stateless, and every newly-created allocator should be equivalent
   * (or else a copy of an allocator will not know what the original has
   * allocated already).
   *
   * This is implemented hiding a singleton in the allocator (as a static
   * member). Each allocator type has its own singleton, i.e., a
   * BulkAllocator<int> does not share memory with a BulkAllocator<double>,
   * but all BulkAllocator<int> share.
   */
  template <typename T>
  class BulkAllocator : public std::allocator<T> {
  public:
    using BaseAllocator_t = std::allocator<T>;

    // import types from the STL allocator
    typedef typename BaseAllocator_t::size_type size_type;
    typedef typename BaseAllocator_t::difference_type difference_type;

    typedef typename BaseAllocator_t::value_type value_type;

    typedef typename BaseAllocator_t::pointer pointer;
    typedef typename BaseAllocator_t::const_pointer const_pointer;

    typedef typename BaseAllocator_t::reference reference;
    typedef typename BaseAllocator_t::const_reference const_reference;

    template <typename U>
    struct rebind {
      typedef BulkAllocator<U> other;
    };

    /// Default constructor: uses the default chunk size
    BulkAllocator() noexcept : BulkAllocator(GetChunkSize(), false) {}

    /// Constructor: sets chunk size and optionally allocates the first chunk
    BulkAllocator(size_type ChunkSize, bool bPreallocate = false) noexcept
    {
      CreateGlobalAllocator(ChunkSize, bPreallocate);
    }

    /// Copy constructor: default
    BulkAllocator(const BulkAllocator& a) noexcept = default;

    /// Move  constructor: default
    BulkAllocator(BulkAllocator&& a) noexcept = default;

    /// General copy constructor; currently, it does not preallocate
    template <class U>
    BulkAllocator(const BulkAllocator<U>& a) noexcept : BaseAllocator_t(a)
    {
      SetChunkSize(a.GetChunkSize());
    }

    /// Copy assignment: default
    BulkAllocator& operator=(const BulkAllocator& a) = default;

    /// Move assignment: default
    BulkAllocator& operator=(BulkAllocator&& a) = default;

    /// Destructor; memory is freed only if no allocators are left around
    ~BulkAllocator() { GlobalAllocator.RemoveUser(); }

    /// Allocates memory for n elements
    pointer allocate(size_type n, const void* /* hint */ = 0);

    /// Frees n elements at p
    void deallocate(pointer p, size_type n);

    /// Releases all the allocated memory: dangerous!
    static void Free() { GlobalAllocator.Free(); }

    /// Returns the chunk size of the underlying global allocator
    static size_type GetChunkSize() { return GlobalAllocator.GetChunkSize(); }

    /// Sets chunk size of global allocator; only affects future allocations!
    static void SetChunkSize(size_type ChunkSize) { GlobalAllocator.SetChunkSize(ChunkSize); }

  private:
    typedef details::bulk_allocator::BulkAllocatorBase<T>
      SharedAllocator_t; ///< shared allocator type

    /// Makes sure we have a valid "global allocator"
    void CreateGlobalAllocator(size_type ChunkSize, bool bPreallocate = false);

    /// The allocator shared by all instances of this object
    static SharedAllocator_t GlobalAllocator;

  }; // class BulkAllocator<>

} // namespace lar

//------------------------------------------------------------------------------
#include <algorithm> // std::max()
#include <array>
#include <iostream>
#include <string>
#include <typeinfo>
#include <vector>
#ifdef __GNUG__
#include <cxxabi.h>
#endif // __GNUG__

namespace lar {
  namespace details {
    ///@{
    /**
     * @brief Demangles the name of a type
     * @tparam T type to be demangled
     * @param [anonymous] parameter to determine the type
     * @return a string with the demangled name
     *
     * This function relies on GCC ABI; if there is no GCC, no demangling
     * happens.
     * One version of this function takes no parameters, and the type must be
     * specified explicitly in the call. The other takes one parameter, that
     * is not actually used but allows the compiler to understand which type
     * to use. The following usese are equivalent:
     * @code
     * std::vector<int> v;
     * std::cout << demangle<std::vector<int>>() << std::endl;
     * std::cout << demangle(v) << std::endl;
     * @endcode
     */
    template <typename T>
    std::string demangle()
    {
      std::string name = typeid(T).name();
#ifdef __GNUG__
      int status; // some arbitrary value to eliminate the compiler warning
      std::unique_ptr<char, void (*)(void*)> res{
        abi::__cxa_demangle(name.c_str(), NULL, NULL, &status), std::free};
      return (status == 0) ? res.get() : name;
#else  // !__GNUG__
      return name;
#endif // ?__GNUG__
    }  // demangle()

    template <typename T>
    inline std::string demangle(const T&)
    {
      return demangle<T>();
    }
    ///@}

    namespace bulk_allocator {
      constexpr bool bDebug = false;

      /// A simple reference counter, keep track of a number of users.
      class ReferenceCounter {
      public:
        typedef unsigned int Counter_t; ///< type of user counter

        /// Returns whether there are currently users
        bool hasUsers() const { return counter > 0; }

        /// Returns the number of registered users
        Counter_t Count() const { return counter; }

        /// Adds a user to the users count
        void AddUser() { ++counter; }

        /// Removed a user to the users count; returns false if no user yet
        bool RemoveUser()
        {
          if (!counter) return false;
          --counter;
          return true;
        }

      private:
        Counter_t counter = 0;
      }; // class ReferenceCounter

      /**
       * @brief A class managing a memory pool
       *
       * The management policy is to allocate *big* chunks of memory.
       * Memory is never freed, until the last user is removed (which is
       * responsibility of the caller), this object is destroyed of Free() is
       * explicitly called.
       * The amount of memory on
       *
       * This class has a users counter. The count must be explicitly handled by
       * the caller.
       */
      template <typename T>
      class BulkAllocatorBase : public ReferenceCounter {
      public:
        typedef size_t size_type;
        typedef T value_type;
        typedef T* pointer;

        /// Constructor; preallocates memory if explicitly requested
        BulkAllocatorBase(size_type NewChunkSize = DefaultChunkSize, bool bPreallocate = false);

        /// Destructor: releases all the memory pool (@see Free())
        ~BulkAllocatorBase() { Free(); }

        /// Releases the pool of memory; all pointer to it become invalid
        void Free();

        /// Returns a pointer to memory for n new values of type T
        pointer Get(size_type n);

        /// Releases memory pointed by the specified pointer (but it does not).
        void Release(pointer) {}

        /// Add a new pool user with the current parameters
        void AddUser() { ReferenceCounter::AddUser(); }

        /// Add a new pool user with new parameters
        void AddUser(size_type NewChunkSize, bool bPreallocate = false);

        /// Removed a user to the users count; if no user is left, free the pool
        bool RemoveUser();

        /// Returns the total number of entries in the pool
        size_type AllocatedCount() const;

        /// Returns the total number of used entries in the pool
        size_type UsedCount() const;

        /// Returns the total number of unused entries in the pool
        size_type FreeCount() const;

        /// Returns the number of memory pool chunks allocated
        size_type NChunks() const { return MemoryPool.size(); }

        /// Returns an array equivalent to { UsedCount(), FreeCount() }
        std::array<size_type, 2> GetCounts() const;

        /// Sets the chunk size for the future allocations
        void SetChunkSize(size_type NewChunkSize, bool force = false);

        /// Returns the current chunk size
        size_type GetChunkSize() const { return ChunkSize; }

        /// Preallocates a chunk of the current ChunkSize;
        /// @see Preallocate(size_type)
        void Preallocate() { Preallocate(ChunkSize); }

      private:
        typedef std::allocator<T> Allocator_t;
        typedef typename Allocator_t::difference_type difference_type;

        Allocator_t allocator; ///< the actual allocator we use

        /// Internal memory chunk; like a std::vector, but does not construct
        class MemoryChunk_t {
        public:
          Allocator_t* allocator; ///< reference to the allocator to be used

          pointer begin = nullptr; ///< start of the pool
          pointer end = nullptr;   ///< end of the pool
          pointer free = nullptr;  ///< first unused element of the pool

          ///< Constructor: allocates memory
          MemoryChunk_t(Allocator_t& alloc, size_type n) : allocator(&alloc)
          {
            begin = n ? allocator->allocate(n) : nullptr;
            end = begin + n;
            free = begin;
          }                                             // MemoryChunk_t()
          MemoryChunk_t(const MemoryChunk_t&) = delete; ///< Can't copy
          MemoryChunk_t(MemoryChunk_t&&);               ///< Move constructor

          ~MemoryChunk_t() { allocator->deallocate(begin, size()); }

          MemoryChunk_t& operator=(const MemoryChunk_t&) = delete;
          ///< Can't assign
          MemoryChunk_t& operator=(MemoryChunk_t&&); ///< Move assignment

          /// Returns the number of elements in this pool
          size_type size() const { return end - begin; }

          /// Returns the number of free elements in this pool
          size_type available() const { return end - free; }

          /// Returns the number of used elements in this pool
          size_type used() const { return free - begin; }

          /// Returns whether the chunk is full
          bool full() const { return !available(); }

          /// Returns a pointer to a free item, or nullptr if none is available
          pointer get() { return (free != end) ? free++ : nullptr; }

          /// Returns a pointer to n free items, or nullptr if not available
          pointer get(size_t n)
          {
            pointer ptr = free;
            if ((free += n) <= end) return ptr;
            free = ptr;
            return nullptr;
          }

        private:
          ///< Default constructor (does nothing)
          MemoryChunk_t() : allocator(nullptr) {}

        }; // class MemoryChunk_t

        typedef std::vector<MemoryChunk_t> MemoryPool_t;

        size_type ChunkSize;     ///< size of the chunks to add
        MemoryPool_t MemoryPool; ///< list of all memory chunks; first is free

        /// Default chunk size (default: 10000)
        static size_type DefaultChunkSize;

        /// Preallocates a chunk of the given size; allocates if free space < n
        void Preallocate(size_type n);

      }; // class BulkAllocatorBase<>

      template <typename T>
      BulkAllocatorBase<T>::MemoryChunk_t::MemoryChunk_t(MemoryChunk_t&& from)
      {
        std::swap(allocator, from.allocator);
        std::swap(begin, from.begin);
        std::swap(end, from.end);
        std::swap(free, from.free);
      } // BulkAllocatorBase<T>::MemoryChunk_t::MemoryChunk_t(MemoryChunk_t&&)

      template <typename T>
      typename BulkAllocatorBase<T>::MemoryChunk_t& BulkAllocatorBase<T>::MemoryChunk_t::operator=(
        MemoryChunk_t&& from)
      {
        std::swap(allocator, from.allocator);
        std::swap(begin, from.begin);
        std::swap(end, from.end);
        std::swap(free, from.free);
        return *this;
      } // BulkAllocatorBase<T>::MemoryChunk_t::operator= (MemoryChunk_t&&)

      template <typename T>
      typename BulkAllocatorBase<T>::size_type BulkAllocatorBase<T>::DefaultChunkSize = 10000;

      template <typename T>
      BulkAllocatorBase<T>::BulkAllocatorBase(size_type NewChunkSize,
                                              bool bPreallocate /* = false */)
        : ChunkSize(NewChunkSize), MemoryPool()
      {
        Preallocate(bPreallocate ? ChunkSize : 0);
        if (bDebug) {
          std::cout << "BulkAllocatorBase[" << ((void*)this) << "] created for type "
                    << demangle<value_type>() << " with chunk size " << GetChunkSize() << " x"
                    << sizeof(value_type) << " byte => " << (GetChunkSize() * sizeof(value_type))
                    << " bytes/chunk" << std::endl;
        } // if debug
      }   // BulkAllocatorBase<T>::BulkAllocatorBase()

      template <typename T>
      void BulkAllocatorBase<T>::Free()
      {
        if (bDebug) {
          std::cout << "BulkAllocatorBase[" << ((void*)this) << "] freeing " << NChunks()
                    << " memory chunks with " << AllocatedCount() << " elements" << std::endl;
        } // if debug
        MemoryPool.clear();
      } // BulkAllocatorBase<T>::Free()

      template <typename T>
      bool BulkAllocatorBase<T>::RemoveUser()
      {
        ReferenceCounter::RemoveUser();
        if (hasUsers()) return true;
        Free();
        return false;
      } // BulkAllocatorBase<T>::RemoveUser()

      template <typename T>
      void BulkAllocatorBase<T>::AddUser(size_type NewChunkSize, bool bPreallocate /* = false */)
      {
        AddUser();
        SetChunkSize(NewChunkSize);
        Preallocate(bPreallocate ? ChunkSize : 0);
      } // BulkAllocatorBase<T>::AddUser(size_type, bool )

      template <typename T>
      void BulkAllocatorBase<T>::Preallocate(size_type n)
      {
        if (MemoryPool.empty() || (MemoryPool.front().available() < n))
          MemoryPool.emplace_back(allocator, n);
      } // BulkAllocatorBase<T>::RemoveUser()

      template <typename T>
      typename BulkAllocatorBase<T>::size_type BulkAllocatorBase<T>::AllocatedCount() const
      {
        size_type n = 0;
        for (const auto& chunk : MemoryPool)
          n += chunk.size();
        return n;
      } // AllocatedCount()

      template <typename T>
      typename BulkAllocatorBase<T>::size_type BulkAllocatorBase<T>::UsedCount() const
      {
        size_type n = 0;
        for (const auto& chunk : MemoryPool)
          n += chunk.used();
        return n;
      } // BulkAllocatorBase<T>::UsedCount()

      template <typename T>
      std::array<typename BulkAllocatorBase<T>::size_type, 2> BulkAllocatorBase<T>::GetCounts()
        const
      {
        // BUG the double brace syntax is required to work around clang bug 21629
        // (https://bugs.llvm.org/show_bug.cgi?id=21629)
        std::array<BulkAllocatorBase<T>::size_type, 2> stats = {{0U, 0U}};
        for (const auto& chunk : MemoryPool) {
          stats[0] += chunk.used();
          stats[1] += chunk.available();
        } // for
        return stats;
      } // BulkAllocatorBase<T>::GetCounts()

      template <typename T>
      void BulkAllocatorBase<T>::SetChunkSize(size_type NewChunkSize, bool force /* = false */)
      {
        if ((GetChunkSize() == NewChunkSize) && !force) return;
        if (bDebug) {
          std::cout << "BulkAllocatorBase[" << ((void*)this) << "]"
                    << " changing chunk size: " << GetChunkSize() << " => " << NewChunkSize << ": x"
                    << sizeof(value_type) << " byte => " << (NewChunkSize * sizeof(value_type))
                    << " bytes/chunk" << std::endl;
        }
        ChunkSize = NewChunkSize;
      } // BulkAllocatorBase<T>::SetChunkSize()

      template <typename T>
      inline typename BulkAllocatorBase<T>::pointer BulkAllocatorBase<T>::Get(size_type n)
      {
        if (n == 0) return nullptr;
        // get the free pointer from the latest chunk (#0)
        pointer ptr = MemoryPool.front().get(n);
        if (ptr) return ptr;
        // no free element left in that chunk:
        // - create a new one in the first position of the pool (move the rest)
        // - return the pointer from the new pool
        //   (if it fails, it fails... but how could that happen?)
        if (bDebug) {
          std::array<size_type, 2> stats = GetCounts();
          std::cout << "BulkAllocatorBase[" << ((void*)this) << "] allocating "
                    << std::max(ChunkSize, n) << " more elements (on top of the current "
                    << (stats[0] + stats[1]) << " elements, " << stats[1] << " unused)"
                    << std::endl;
        } // if debug
        return MemoryPool.emplace(MemoryPool.begin(), allocator, std::max(ChunkSize, n))->get(n);
      } // BulkAllocatorBase<T>::Get()

    } // namespace bulk_allocator
  }   // namespace details

  template <typename T>
  typename BulkAllocator<T>::SharedAllocator_t BulkAllocator<T>::GlobalAllocator;

  template <typename T>
  void BulkAllocator<T>::CreateGlobalAllocator(size_type ChunkSize, bool bPreallocate /* = false */)
  {
    GlobalAllocator.AddUser(ChunkSize, bPreallocate);
  } // BulkAllocator<T>::CreateGlobalAllocator()

  template <typename T>
  inline typename BulkAllocator<T>::pointer BulkAllocator<T>::allocate(size_type n,
                                                                       const void* /* hint = 0 */)
  {
    return GlobalAllocator.Get(n);
  }

  template <typename T>
  inline void BulkAllocator<T>::deallocate(pointer p, size_type n)
  {
    return GlobalAllocator.Release(p);
  } // BulkAllocator<T>::deallocate()
} // namespace lar

#endif // BULKALLOCATOR_H
