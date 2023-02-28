/**
 * @file   FindAllP.h
 * @brief  More association queries
 * @date   December 18th, 2014
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 * The code here is implementation details, living in lar::util::details
 * namespace.
 */

#ifndef FINDALLP_H
#define FINDALLP_H 1

// C/C++ standard libraries
#include <climits>    // CHAR_BIT
#include <functional> // std::hash<>
#include <stdexcept>  // std::out_of_range
#include <unordered_map>
#include <vector>

// framework libraries
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

/// LArSoft-specific namespace
namespace lar {

  /// LArSoft utility namespace
  namespace util {

    /// LArSoft utility implementation details
    namespace details {

      /// Hash functions for art and larsoft objects
      template <typename T>
      struct hash {
        using result_type = size_t;
        using argument_type = T;

        result_type operator()(argument_type const& v) const;

      }; // class hash<>

      /** **********************************************************************
       * @brief A class holding many associations between objects
       * @tparam Source type of the object we use as query key (indexing object)
       * @tparam Dest type of the object we query for
       *
       * This class is conceptually related to the query object FindOneP.
       * This object is a cache of possible query results of the type:
       * which Dest object is associated to this specific Src object?
       * The cache is structured so that only one Dest object is known for each
       * Src.
       */
      template <typename Source, typename Dest>
      class UniqueAssociationCache {
      public:
        using Source_t = Source;
        using Dest_t = Dest;

        using SourcePtr_t = art::Ptr<Source_t>;
        using DestPtr_t = art::Ptr<Dest_t>;

        /// type for a cache of dest products for a given source product ID
        using InProductCache_t = std::vector<DestPtr_t>;

        /// type of hash functor for the cache
        using Hash_t = details::hash<art::ProductID>;

        /// type for the complete cache, keyed by source product ID
        using Cache_t = std::unordered_map<art::ProductID, InProductCache_t, Hash_t>;

        Cache_t AssnCache; ///< association cache, keyed by product ID and index

        /// Constructor: an empty cache
        UniqueAssociationCache() = default;

        /**
         * @brief Returns the specified element of the cache
         * @param src art pointer to the object we want the association of
         * @return the requested element, or a null pointer if not found
         */
        DestPtr_t operator[](SourcePtr_t const& src) const
        {
          return AssnCache[src.id()][src.key()];
        }

        /// Empties the cache
        void clear() { AssnCache.clear(); }

        size_t NProductIDs() const { return AssnCache.size(); }

      }; // class UniqueAssociationCache<>

      /** **********************************************************************
       * @brief Query object reading *all* the associations between two classes
       * @tparam Source type of the object we use as query key (indexing object)
       * @tparam Dest type of the object we query for
       *
       * When assigned an event, this object reads all the associations from
       * Source type classes to Dest type classes in the event, and stores
       * their information in a map to track a Dest object from its Source one.
       * In fact, it assumes that only one Dest object is associated,
       * event-wise, to each single Source object.
       */
      template <typename Source, typename Dest>
      class FindAllP {
      public:
        using Source_t = Source;
        using Dest_t = Dest;

        /// Default constructor: empty query, read information with Read()
        FindAllP() = default;

        /// Constructor: reads all associations from the specified event
        FindAllP(art::Event& event) : FindAllP() { Read(event); }

        /// Constructor: reads one association from the specified event
        FindAllP(art::Event& event, art::InputTag assnTag) : FindAllP() { Read(event, assnTag); }

        /**
         * @brief Returns the object associated to the specified one
         * @param src a art pointer to the source object
         * @return a pointer to the associated object, or a null pointer if none
         */
        art::Ptr<Dest_t> const& operator[](art::Ptr<Dest_t> const& src) const;

        /// Returns whether there are associations from objects in product id
        bool hasProduct(art::ProductID const& id) const;

        /// Returns if there are associations from objects with product as ptr
        bool hasProduct(art::Ptr<Source_t> const& ptr) const;

        /**
         * @brief Reads all the associations from the event
         * @throw art::Exception if multiple dest objects are found for one
         *   source object
         */
        unsigned int Read(art::Event& event);

        /**
         * @brief Reads the specified association from the event
         * @param event the event to read the associations from
         * @param assnTag the input tag for the association
         * @throw art::Exception if multiple dest objects are found for one
         *   source object
         *
         * The input tag for the association is usually simply a string with the
         * name of the module that produced the association, and often the same
         * module has also produced the source objects as well.
         */
        unsigned int Read(art::Event& event, art::InputTag const& assnTag);

        /**
         * @brief Reads the specified association from the event
         * @param event the event to read the associations from
         * @param assnTag the input tag for the association
         * @throw art::Exception if multiple dest objects are found for one
         *   source object
         *
         * The existing associations already in cache are not removed.
         *
         * The input tag for the association is usually simply a string with the
         * name of the module that produced the association, and often the same
         * module has also produced the source objects as well.
         */
        unsigned int Add(art::Event& event, art::InputTag const& assnTag);

      protected:
        using Assns_t = art::Assns<Source_t, Dest_t>;

        /// Type of the cache
        using Cache_t = UniqueAssociationCache<Source_t, Dest_t>;

        Cache_t cache; ///< set of associations, keyed by product ID and key

        /// Adds all associations in the specified handle; returns their number
        unsigned int Merge(art::Handle<Assns_t>& handle);
      }; // class FindAllP<>

      /** **********************************************************************
       * @brief Resizes a vector to a size power of 2, with a minimum size
       * @param v the vector to be resized
       * @param min_size the minimal size of the vector
       *
       * The vector v is resized, any new element is default-initialized.
       * The target size is the smallest power of 2 not smaller than min_size,
       * or 0 if min_size is 0.
       */
      template <typename T>
      void ResizeToPower2(std::vector<T>& v, size_t min_size);

    } // namespace details
  }   // namespace util
} // namespace lar

//******************************************************************************
//***  Template implementation
//******************************************************************************

namespace lar {
  namespace util {
    namespace details {
      //------------------------------------------------------------------------
      //---  FindAllP

      template <typename Source, typename Dest>
      auto FindAllP<Source, Dest>::operator[](art::Ptr<Dest_t> const& src) const
        -> art::Ptr<Dest_t> const&
      {
        // we expect a missing match to be exceptional
        try {
          return cache.AssnCache.at(src.id()).at(src.key());
        }
        catch (std::out_of_range) {
          return {};
        }
      } // FindAllP<>::operator[]

      template <typename Source, typename Dest>
      inline bool FindAllP<Source, Dest>::hasProduct(art::ProductID const& id) const
      {
        return cache.AssnCache.count(id) > 0;
      }

      template <typename Source, typename Dest>
      inline bool FindAllP<Source, Dest>::hasProduct(art::Ptr<Source_t> const& ptr) const
      {
        return hasProduct(ptr.id());
      }

      template <typename Source, typename Dest>
      unsigned int FindAllP<Source, Dest>::Read(art::Event& event)
      {

        // read all the associations between source and destination class types
        //std::vector<art::Handle<Assns_t>> assns_list;
        //event.getManyByType(assns_list);
        auto assns_list = event.getMany<Assns_t>();

        MF_LOG_DEBUG("FindAllP") << "Read(): read " << assns_list.size() << " association sets";

        unsigned int count = 0;
        // parse all the associations, and translate them into a local cache
        for (art::Handle<Assns_t> handle : assns_list)
          count += Merge(handle);

        MF_LOG_DEBUG("FindAllP") << "Read " << count << " associations for " << cache.NProductIDs()
                                 << " product IDs";

        return count;
      } // FindAllP::Read(Event)

      template <typename Source, typename Dest>
      unsigned int FindAllP<Source, Dest>::Read(art::Event& event, art::InputTag const& assnTag)
      {
        cache.clear();
        return Add(event, assnTag);
      } // FindAllP::Read(Event, InputTag)

      template <typename Source, typename Dest>
      unsigned int FindAllP<Source, Dest>::Add(art::Event& event, art::InputTag const& assnTag)
      {

        // read the association between source and destination class types
        art::Handle<Assns_t> handle;
        if (!event.getByLabel(assnTag, handle)) {
          throw art::Exception(art::errors::ProductNotFound)
            << "no association found with input tag '" << assnTag << "'";
        }

        return Merge(handle);
      } // FindAllP::Add(Event, InputTag)

      template <typename Source, typename Dest>
      unsigned int FindAllP<Source, Dest>::Merge(art::Handle<Assns_t>& handle)
      {
        // product ID of the last source object; initialized invalid
        art::ProductID LastProductID = art::Ptr<Source_t>().id();
        typename Cache_t::InProductCache_t const* AssnsList = nullptr;

        unsigned int count = 0;

        MF_LOG_DEBUG("FindAllP") << "Merge(): importing " << handle->size() << " associations from "
                                 << handle.provenance();

        for (auto const& assn : *handle) {
          // assn is a std::pair<art::Ptr<Source_t>, art::Ptr<Dest_t>>
          art::Ptr<Source_t> const& src = assn.first;

          if (src.isNull()) {
            MF_LOG_ERROR("FindAllP")
              << "Empty pointer found in association " << handle.provenance();
            continue; // this should not happen
          }

          art::Ptr<Dest_t> const& dest = assn.second;

          // if we have changed product (that should be fairly rare),
          // update the running pointers
          if (src.id() != LastProductID) {
            LastProductID = src.id();
            AssnsList = &(cache.AssnCache[LastProductID]);

            // if the list is empty, it means we have just created it!
            if (AssnsList->empty()) {
              // allocate enough space to accomodate all the associations,
              // (provided that source IDs are sequencial);
              // in fact typically all the associations in the same handle
              // have the same product ID
              ResizeToPower2(*AssnsList, handle->size());
            }
          } // if different product ID

          // make sure there is enough room in the vector
          typename art::Ptr<Source_t>::key_type key = src.key();
          if (key >= AssnsList->size()) ResizeToPower2(*AssnsList, key + 1);

          // store the association to dest
          art::Ptr<Dest_t>& dest_cell = (*AssnsList)[key];
          if (dest_cell.isNonnull() && (dest_cell != dest)) {
            throw art::Exception(art::errors::InvalidNumber)
              << "Object Ptr" << src << " is associated with at least two objects: " << dest
              << " and " << dest_cell;
          }
          dest_cell = dest;
          ++count;
        } // for all associations in a list

        MF_LOG_DEBUG("FindAllP") << "Merged " << count << " associations from "
                                 << handle.provenance();
        return count;
      } // FindAllP::Merge()

      //------------------------------------------------------------------------
      template <>
      auto hash<art::ProductID>::operator()(argument_type const& id) const -> result_type
      {
        // make sure we have enough bits in result_type;
        // if not, we need a more clever algorithm
        //static_assert(
        //  sizeof(id.processIndex()) + sizeof(id.productIndex())
        //    <= sizeof(result_type),
        //  "hash return type not large enough for hashing art::ProductID"
        //);
        // stack the process and product IDs in one integer
        //return result_type(
        //(id.processIndex() << sizeof(id.productIndex() * CHAR_BIT))
        //+ id.productIndex()
        //);
        return result_type(id.value());
      } // hash<art::ProductID>::operator()

      //------------------------------------------------------------------------
      template <typename T>
      void ResizeToPower2(std::vector<T>& v, size_t min_size)
      {
        if (min_size == 0) {
          v.clear();
          return;
        }
        size_t new_size = 1;
        while (new_size < min_size)
          new_size *= 2;

        v.resize(new_size);
      } // ResizeToPower2()

    } // namespace details
  }   // namespace util
} // namespace lar

//------------------------------------------------------------------------------
#endif // FINDALLP_H 1
