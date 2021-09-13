/**
 * @file   lardata/Utilities/FindManyInChainP.tcc
 * @brief  Template implementation for `FindManyInChainP.h`.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 26, 2017
 *
 */

#ifndef LARDATA_UTILITIES_FINDMANYINCHAINP_TCC
#define LARDATA_UTILITIES_FINDMANYINCHAINP_TCC

#ifndef LARDATA_UTILITIES_FINDMANYINCHAINP_H
#error "FindManyInChainP.tcc must not be included directly. Include FindManyInChainP.h instead."
#endif // LARDATA_UTILITIES_FINDMANYINCHAINP_H

// framework
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Utilities/Exception.h"

// C/C++ standard library
#include <map>
#include <iterator> // std::begin(), std::cbegin(), std::distance()...
#include <tuple> // std::tuple_cat(), ...
#include <algorithm> // std::lower_bound(), std::sort(), std::transform()...
#include <utility> // std::pair<>, std::move(), std::declval()...
#include <type_traits> // std::decay_t<>, std::enable_if_t<>, ...


//------------------------------------------------------------------------------
//---  template implementation
//---
namespace lar {

  namespace details {

    //--------------------------------------------------------------------------
    // this is a copy of things in cetlib; should be replaced by the original
    // if they get exposed
    template <class T, class R /* = void */>
    struct enable_if_type_exists { using type = R; };

    //--------------------------------------------------------------------------
    /// Trait evaluating to true if H is a art Handle type.
    template <typename H, typename = void>
    struct is_handle: public std::false_type {};

    template <typename H>
    struct is_handle<H, enable_if_is_handle_t<H>>: public std::true_type {};

    template <typename H>
    constexpr bool is_handle_v = is_handle<H>::value;

    //--------------------------------------------------------------------------
    template <typename T>
    struct is_art_ptr_impl: public std::false_type {};

    template <typename T>
    struct is_art_ptr_impl<art::Ptr<T>>: public std::true_type {};

    /// Whether the specified type is an art::Ptr.
    template <typename T>
    using is_art_ptr = is_art_ptr_impl<std::decay_t<T>>;


    //--------------------------------------------------------------------------

    /// Hosts as `type` the type in position `N` (first is `0`).
    // this might be implemented with std::tuple_element;
    // here we get control on error messages
    template <unsigned int N, typename... Args>
    struct get_type;

    template <unsigned int N, typename... Args>
    using get_type_t = typename get_type<N, Args...>::type;


    template <unsigned int N, typename First, typename... Others>
    struct get_type<N, First, Others...> {
      static_assert
        (N <= sizeof...(Others), "No type in the specified position.");
      using type = get_type_t<N-1, Others...>;
    }; // get_type<>

    template <typename First, typename... Others>
    struct get_type<0U, First, Others...> { using type = First; };


    //--------------------------------------------------------------------------

    template <unsigned int NCopies>
    struct TupleAppender {
      template <typename Tuple, typename T>
      static auto append(Tuple&& t, T const& value)
        {
          return TupleAppender<(NCopies - 1)>::append(
            std::tuple_cat(std::forward<Tuple>(t), std::make_tuple(value)),
            value
            );
        }
    }; // TupleAppender<>

    template <>
    struct TupleAppender<0> {
      template <typename Tuple, typename T>
      static decltype(auto) append(Tuple&& t, T const&) { return std::forward<Tuple>(t); }
    }; // TupleAppender<0>


    template <unsigned int NCopies, typename T, typename Tuple>
    auto appendToTuple(Tuple&& t, T value)
      { return TupleAppender<NCopies>::append(std::forward<Tuple>(t), value); }


    //--------------------------------------------------------------------------
    /// Returns the input tag of the product identified by the handle.
    template <typename Handle>
    art::InputTag tagFromHandle(Handle&& handle)
      { return handle.provenance()->inputTag(); }

    /// Returns the input tag of the product identified by `id`.
    template <typename Data, typename Event>
    art::InputTag tagFromProductID
      (art::ProductID const& id, Event const& event)
    {
      // This is not efficient for repeated queries--perhaps a map can
      // be created that maps ProductID to input tag.
      auto pd = event.getProductDescription(id);
      if (pd != nullptr) {
        return pd->inputTag();
      }
      throw art::Exception(art::errors::ProductNotFound)
        << "Couldn't find data product with product ID " << id << "\n";
    } // tagFromProductID()


    //--------------------------------------------------------------------------
    template <typename Value>
    class SimpleDataIndex {
        public:
      using Value_t = Value;

      /// Constructor: indexes pointers to data in the specified collection.
      template <typename BeginIter, typename EndIter>
      SimpleDataIndex(BeginIter begin, EndIter end) { init(begin, end); }

      /// Constructor: indexes pointers to data in the specified collection.
      template <typename Coll>
      SimpleDataIndex(Coll& data)
        { using std::begin; using std::end; init(begin(data), end(data)); }

      /// Returns a pointer to the matched data, nullptr if none.
      std::size_t const* operator() (Value_t const& key) const
        {
          auto const e = index.cend();
          auto const it = std::lower_bound(index.cbegin(), e, key, comparer());
          return ((it == e) || (key < it->first))? nullptr: &(it->second);
        }


        private:
      using IndexElement_t = std::pair<Value_t const*, std::size_t>;

      std::vector<IndexElement_t> index; ///< The index.

      template <typename BeginIter, typename EndIter>
      void init(BeginIter begin, EndIter end)
        {
          index.reserve(std::distance(begin, end)); // hope it's fast
          std::size_t i = 0;
          for (auto it = begin; it != end; ++it, ++i)
            index.emplace_back(&*it, i);
          std::sort(index.begin(), index.end(), comparer());
        }

      struct KeyComparer {
        bool operator() (Value_t const& a, Value_t const& b) const
          { return a < b; }
        bool operator() (IndexElement_t  const& a, Value_t const& b) const
          { return *(a.first) < b; }
        bool operator() (Value_t const& a, IndexElement_t const& b) const
          { return a < *(b.first); }
        bool operator() (IndexElement_t const& a, IndexElement_t const& b) const
          { return *(a.first) < *(b.first); }
      }; // struct KeyComparer

      static auto comparer() { return KeyComparer(); }

    }; // class SimpleDataIndex


    template <typename Key, typename Value>
    class DataIndex {
        public:
      using Key_t = Key;
      using Value_t = Value;
      using ValuePtr_t = std::add_pointer_t<Value_t>;

      /// Constructor: indexes pointers to data in the specified collection.
      template <typename BeginIter, typename EndIter, typename KeyExtractor>
      DataIndex(BeginIter begin, EndIter end, KeyExtractor&& getKey)
        { init(begin, end, std::forward<KeyExtractor>(getKey)); }

      /// Constructor: indexes pointers to data in the specified collection.
      template <typename Coll, typename KeyExtractor>
      DataIndex(Coll& data, KeyExtractor&& getKey)
        {
          using std::begin;
          using std::end;
          init(begin(data), end(data), std::forward<KeyExtractor>(getKey));
        }

      /// Returns a pointer to the matched data, nullptr if none.
      ValuePtr_t operator() (Key_t const& key) const
        {
          auto const e = index.cend();
          auto const it = std::lower_bound(index.cbegin(), e, key, comparer());
          return ((it == e) || (key < it->first))? nullptr: it->second;
        }


        private:
      using IndexElement_t = std::pair<Key_t, Value_t*>;

      std::vector<IndexElement_t> index; ///< The index.

      template <typename BeginIter, typename EndIter, typename KeyExtractor>
      void init(BeginIter begin, EndIter end, KeyExtractor&& getKey)
        {
          index.reserve(std::distance(begin, end)); // hope it's fast
          for (auto it = begin; it != end; ++it)
            index.emplace_back(getKey(*it), &*it);
          std::sort(index.begin(), index.end(), comparer());
        }

      struct KeyComparer {
        bool operator() (Key_t const& a, Key_t const& b) const
          { return a < b; }
        bool operator() (IndexElement_t  const& a, Key_t const& b) const
          { return a.first < b; }
        bool operator() (Key_t const& a, IndexElement_t const& b) const
          { return a < b.first; }
        bool operator() (IndexElement_t const& a, IndexElement_t const& b) const
          { return a.first < b.first; }
      }; // struct KeyComparer

      static auto comparer() { return KeyComparer(); }

    }; // class DataIndex


    /// Returns a new index for data using the keys extracted by getKey.
    template <typename Coll, typename KeyExtractor>
    auto makeIndex(Coll& data, KeyExtractor&& getKey)
      {
        using Value_t = typename Coll::value_type;
        using Key_t = std::decay_t<decltype(getKey(std::declval<Value_t>()))>;
        return
          DataIndex<Key_t, Value_t>(data, std::forward<KeyExtractor>(getKey));
      }

    /// Returns a new index for data in range using the data itself as key.
    template <typename BeginIter, typename EndIter>
    auto makeSimpleIndex(BeginIter begin, EndIter end)
      {
        using Value_t = std::decay_t<decltype(*begin)>;
        return SimpleDataIndex<Value_t>(begin, end);
      }

    /// Returns a new index for data using the data itself as key.
    template <typename Coll>
    auto makeSimpleIndex(Coll& data)
      {
        using Value_t = typename Coll::value_type;
        return SimpleDataIndex<Value_t>(data);
      }


    //--------------------------------------------------------------------------
    //---  finders implementation
    //---

    namespace AssociationFinderBase {

      /// Information on a single source: its index in the source list, and
      /// "key" (it can be any sortable object).
      template <typename Key = std::size_t>
      struct SourceIDinfo_t: public std::tuple<Key, std::size_t> {
        using Base_t = std::tuple<Key, std::size_t>;
        using Key_t = Key;

        /// Returns the key of the art pointer to this source element.
        Key_t const& key() const { return std::get<0>(*this); }
        /// Returns the index of this source element within the source list.
        std::size_t position() const { return std::get<1>(*this); }

        using Base_t::Base_t; // inherit constructor

        /// Comparison functor: compares art pointer keys.
        struct Comparer_t {
          bool operator()
            (SourceIDinfo_t const& a, SourceIDinfo_t const& b) const
            { return a.key() < b.key(); }
          bool operator() (SourceIDinfo_t const& a, Key_t const& b) const
            { return a.key() < b; }
          bool operator() (Key_t const& a, SourceIDinfo_t const& b) const
            { return a < b.key(); }
        }; // struct SourceIDinfo_t::Comparer_t

        static Comparer_t comparer() { return {}; }
      }; // struct SourceIDinfo_t

      template <typename Key = std::size_t>
      struct SourceVector_t: private std::vector<SourceIDinfo_t<Key>> {
        using Base_t = std::vector<SourceIDinfo_t<Key>>;
        using Info_t = typename Base_t::value_type;
        using Key_t = typename Info_t::Key_t;

        using const_iterator = typename Base_t::const_iterator;

        using Base_t::begin;
        using Base_t::end;
        using Base_t::size;
        using Base_t::empty;
        using Base_t::reserve;

        struct MatchConstIterator_t: public std::pair<const_iterator, bool> {
          using Base_t = std::pair<const_iterator, bool>;
          using Base_t::Base_t; // inherit all constructors
          auto const& operator->() const { return std::get<0>(*this); }
          operator bool() const { return std::get<1>(*this); }
          bool operator!() const { return !this->operator bool(); }
        }; // struct MatchConstIterator_t

        /// Default constructor: starts empty.
        SourceVector_t() = default;

        /// Constructor: from the elements between the specified iterators.
        template <typename BeginIter, typename EndIter>
        SourceVector_t(BeginIter b, EndIter e) { assign(b, e); }

        /// Sorts the vector according to keys.
        void sort()
          { std::sort(Base_t::begin(), Base_t::end(), Info_t::comparer()); }

        /// Constructs a new element at the end of the vector.
        void emplace_back(Key_t const& key, std::size_t pos)
          { Base_t::emplace_back(key, pos); }

        /// Returns iterator to element with key, if it exists.
        /// Requires the vector to be sort()ed.
        /// If key is not present it may return end() or the wrong element.
        auto find(Key_t const& key) const
          {
            auto const e = end();
            auto it = std::lower_bound(begin(), e, key, Info_t::comparer());
            return MatchConstIterator_t(it, (it != e) && (it->key() == key));
          }

        /// Resets the content to map a container delimited by the iterators.
        template <typename BeginIter, typename EndIter>
        void assign(BeginIter b, EndIter e)
          {
            Base_t::clear();
            std::size_t i = 0;
            for (auto it = b; it != e; ++it, ++i) emplace_back(it->key(), i);
            sort();
          }

      }; // struct SourceVector_t


      /// Returns a tuple of NTags elements, including all specified tags and
      /// then copies of the default one to make it N.
      template <unsigned int NTags, typename DefaultTag, typename... InputTags>
      auto makeTagsTuple(DefaultTag defTag, InputTags... tags) {
        return appendToTuple<(NTags - sizeof...(InputTags))>
          (std::make_tuple(tags...), defTag);
      } // makeTagsTuple()


      /// A named pair of a key and all its connections.
      /// Comparisons are based solely on key.
      template <typename A, typename B>
      struct Connections: private std::pair<A, std::vector<B>> {
        using Base_t = std::pair<A, std::vector<B>>;
        using Key_t = A;
        using ConnectedItem_t = B;
        using Connections_t = Connections<A, B>;
        using Connected_t = std::tuple_element_t<1, Base_t>;
        using Base_t::Base_t; // import pair constructors

        /// Additional constructor: key with no connected items.
        Connections(Key_t const& key): Base_t(key, {}) {}

        Key_t const& key() const { return std::get<0>(*this); }
        Connected_t const& connectedTo() const { return std::get<1>(*this); }
        Connected_t& connectedTo() { return std::get<1>(*this); }

        struct Comparer_t {
          bool operator() (Key_t const& a, Key_t const& b) const
            { return a < b; }
          bool operator() (Key_t const& a, Connections_t const& b) const
            { return a < b.key(); }
          bool operator() (Connections_t const& a, Key_t const& b) const
            { return a.key() < b; }
          bool operator() (Connections_t const& a, Connections_t const& b) const
            { return a.key() < b.key(); }
        }; // Comparer_t

        static Comparer_t comparer() { return {}; }

        /// Compares the connection lists by key only.
        bool operator< (Connections_t const& b) const
          { return comparer()(*this, b); }

      }; // struct Connections<>


      struct KeysFromTag {};
      struct KeysFromIndexTag {};
      struct KeysFromOtherListTag {};
      constexpr KeysFromTag KeysFrom;
      constexpr KeysFromIndexTag KeysFromIndex;
      constexpr KeysFromOtherListTag KeysFromOtherList;

      /// A connections list.
      template <typename Key, typename Connected>
      class ConnectionList: private std::vector<Connections<Key, Connected>> {
        using Connections_t = Connections<Key, Connected>;
        using Data_t = std::vector<Connections_t>;

        Data_t& data() { return *this; }
        Data_t const& data() const { return *this; }

          public:
        /// Type of an item connected to the key.
        using ConnectedItem_t = typename Connections_t::ConnectedItem_t;
        /// Type of list of connected items.
        using Connected_t = typename Connections_t::Connected_t;

        // import types from the base class
        using value_type = typename Data_t::value_type;


        /// Initializes an empty connection list.
        ConnectionList() = default;

        /// Initializes with no connection for each of the keys from b to e.
        template <typename BeginIter, typename EndIter>
        ConnectionList(KeysFromTag, BeginIter b, EndIter e)
          {
            data().reserve(std::distance(b, e)); // hope this is fast
            std::transform(b, e, std::back_inserter(data()),
              [](auto const& key) { return Connections_t(key, {}); }
              );
          }

        /// Initializes with no connection for each of the keys from a list.
        template <typename Other>
        ConnectionList
          (KeysFromOtherListTag, ConnectionList<Key, Other> const& source)
          {
            data().reserve(source.size()); // hope this is fast
            using std::cbegin;
            using std::cend;
            std::transform(
              cbegin(source), cend(source), std::back_inserter(data()),
              [](auto const& conn) { return Connections_t(conn.key(), {}); }
              );
          }

        /// Initializes with no connection for each of the generated keys.
        template <typename KeyMaker>
        ConnectionList(KeysFromIndexTag, std::size_t n, KeyMaker keyMaker)
          {
            data().reserve(n);
            for (std::size_t i = 0; i < n; ++i)
              data().emplace_back(keyMaker(i));
          }

        // import methods from the base class (not the constructors)
        using Data_t::size;
        using Data_t::operator[];
        using Data_t::clear;
        using Data_t::cbegin;
        using Data_t::begin;
        using Data_t::cend;
        using Data_t::end;

        /// @{
        /// @name Special methods

        void addConnectionAt(std::size_t index, ConnectedItem_t const& item)
          { data()[index].connectedTo().push_back(item); }
        void addConnectionAt(std::size_t index, ConnectedItem_t&& item)
          { data()[index].connectedTo().push_back(std::move(item)); }

        /// Returns all connected objects on a vector, one element per key.
        /// As a R-value reference method, the current content is lost
        std::vector<Connected_t> values() &&
          {
            std::vector<Connected_t> result;
            result.reserve(data().size());
            std::transform(
              data().begin(), data().end(), std::back_inserter(result),
              [](auto& conn){ return std::move(conn.connectedTo()); }
              );
            return result;
          } // values()

        /// @}

      }; // struct ConnectionList<>


      /// Class handling a connection list.
      template <typename Source, typename Target>
      class PtrConnectionManager {

          public:
        using Manager_t = PtrConnectionManager<Source, Target>;
        using Source_t = Source;
        using Target_t = Target;
        using SourcePtr_t = art::Ptr<Source_t>;
        using TargetPtr_t = art::Ptr<Target_t>;
        using ConnectionList_t = ConnectionList<SourcePtr_t, TargetPtr_t>;

        /// Constructor: steals the connections from the argument.
        PtrConnectionManager(ConnectionList_t&& data): data(std::move(data)) {}


        /// Returns all connected objects as a single connection of pointers.
        std::vector<TargetPtr_t> allConnected() const
          {

            // this might as well be implemented with ranges

            std::vector<TargetPtr_t> result;

            using std::cbegin;
            using std::cend;

            for (std::size_t i = 0; i < data.size(); ++i) {
              auto const& connected = data[i].connectedTo();
              result.insert(result.end(), cbegin(connected), cend(connected));
            } // for
            return result;

          } // allConnected()

        /// Returns a new connection list connecting each of our sources to
        /// FurtherPtr_t objects, using shared TargetPtr_t to map connections.
        template <typename FurtherPtr_t>
        ConnectionList<SourcePtr_t, FurtherPtr_t> join
          (ConnectionList<TargetPtr_t, FurtherPtr_t>&& other) const
        {
          // internally, we'll call the current data A -> B, other data B -> C,
          // and we are after A -> C
          using APtr_t = SourcePtr_t;
        //  using BPtr_t = TargetPtr_t;
          using CPtr_t = FurtherPtr_t;

          auto const& AB = data;
          auto&& BC = other;
          ConnectionList<APtr_t, CPtr_t> AC(KeysFromOtherList, data);

          // the match object returns the iterator to the list connected to the
          // specified key; it does not own the data structure;
          auto const match
            = makeIndex(BC, [](auto const& elem){ return elem.key(); });

          for (std::size_t iA = 0; iA < AB.size(); iA++) {
            auto const& connAB = AB[iA]; // connections for this A
            auto const& BsForA = connAB.connectedTo();

            auto& CsForA = AC[iA].connectedTo(); // connected list to fill (empty)

            // for each B connected to this A:
            for (auto const& key: BsForA) {

              auto iConnBC = match(key);
              if (!iConnBC) continue; // no C's associated to the matching B

              // matched! move or copy the while list of Cs
              auto& CsForB = iConnBC->connectedTo(); // (this list)
              if (CsForA.empty()) CsForA = std::move(CsForB); // move
              else CsForA.insert(CsForA.end(), cbegin(CsForB), cend(CsForB));

            } // for connected Bs
          } // for iA

          return AC;
        } // join()

          private:
        ConnectionList_t data; ///< All connection information.

      }; // PtrConnectionManager<>

      /// Helper function to created the right type of manager from a temporary.
      template <typename A, typename B>
      auto makeConnectionManager(ConnectionList<A, B>&& data)
        {
          return PtrConnectionManager
              <typename A::value_type, typename B::value_type>(std::move(data));
        }

    } // namespace AssociationFinderBase


    template <typename Target>
    struct AssociationFinder {

      using Target_t = Target;
      using TargetPtr_t = art::Ptr<Target_t>;

      /// Type of result of finders: one list of associated targets per source.
      template <typename Source>
      using Result_t
        = AssociationFinderBase::ConnectionList<art::Ptr<Source>, TargetPtr_t>;

      template <typename Source>
      using Assns_t = art::Assns<Source, Target>;

      template <typename Assns>
      using AssnKey_t = typename Assns::left_t;

      template <typename Assns>
      using AssnKeyPtr_t = art::Ptr<AssnKey_t<Assns>>;

      template <typename Assns>
      using ResultFromAssns_t = Result_t<AssnKey_t<Assns>>;

      template <typename Handle>
      using ValueFromHandle_t = typename Handle::element_type::value_type;

      template <typename Handle>
      using AssnsFromHandle_t = Assns_t<ValueFromHandle_t<Handle>>;

      template <typename Handle>
      using ResultFromHandle_t = Result_t<art::Ptr<ValueFromHandle_t<Handle>>>;

      /// Helper finding a associations from a single complete data product.
      template <typename Assns>
      static ResultFromAssns_t<Assns> findFromDataProduct(
        art::ProductID const& sourceID, std::size_t nSources,
        Assns const& assns
        )
      {
        using SourcePtr_t = AssnKeyPtr_t<Assns>;

        using namespace AssociationFinderBase;

        // as many lists in result as sources, keys created from source ID
        ResultFromAssns_t<Assns> result(
          KeysFromIndex, nSources,
          [id=sourceID](std::size_t i){ return SourcePtr_t(id, i, nullptr); }
          );

        // get the associations, following the content of the assns data product
        for (decltype(auto) assn: assns) {
          auto const& sourcePtr = assn.first;

          // does this association contain a pointer with an ID different than
          // the one we are looking for? (it may mean our assumption that the
          // same module produced source and its association is wrong)
          if (sourcePtr.id() != sourceID) continue;

          // we follow the assumption that the data product is complete with
          // nSources elements, therefore no pointer can exist with a larger
          // key:
          assert(sourcePtr.key() < nSources);

          // push target pointer into the result of the matched source
          result.addConnectionAt(sourcePtr.key(), assn.second);

        } // for

        return result;
      } // findFromDataProduct()


      /// Helper finding a associations from a single list.
      template <
        typename Source, typename IndexBegin, typename IndexEnd, typename Event
        >
      static Result_t<Source> findSingle(
        art::ProductID const& sourceID,
        IndexBegin sbegin, IndexEnd send,
        Assns_t<Source> const& assns
        )
      {
        // type of the source object (hidden in PtrColl)
        using Source_t = Source;
        using SourcePtr_t = art::Ptr<Source_t>;
        using namespace AssociationFinderBase;

        std::size_t const nSources = std::distance(sbegin, send);
        Result_t<Source_t> result(nSources);

        // These are the source "pointers" we still have to find; they are all
        // on the same product ID. We keep track of the original position too
        // (unordered map may or may not be better...);
        // given that these are all from the same product ID, the product ID is
        // not saved. SourceVector_t provides log(N) lookup.
        SourceVector_t<std::size_t> sourceInfos(sbegin, send);

        // get the association, following the content of the assns data product
        for (decltype(auto) assn: assns) {
          SourcePtr_t const& sourcePtr = assn.first;

          // does this association contain a pointer with an ID different than
          // the one we are looking for? (it may mean our assumption that the
          // same module produced source and its association is wrong)
          if (sourcePtr.id() != sourceID) continue;

          // is this pointer interesting?
          auto iSrcInfo = sourceInfos.find(sourcePtr.key());
          if (!iSrcInfo) continue; // it's not it

          // match! push target pointer into the result of the matched source
          result[iSrcInfo->position()].connectedTo().push_back(assn.second);

        } // for

        return result;
      } // findSingle()


      /// Helper finding a single degree of associations from the tag of each
      /// source art pointer.
      template <typename PtrCollBegin, typename PtrCollEnd, typename Event>
      static auto findWithRange(
        PtrCollBegin sbegin, PtrCollEnd send, Event const& event,
        art::InputTag const& tag
      ) {
        /*
         * The strategy of this implementation is:
         *
         * 1. collect all source art pointers, sorted by key for faster lookup
         * 2. parse all the associated pairs
         *     1. if the source pointer of a pair is in the list of
         *        interesting source pointers, push the target pointer of the
         *        pair into the results for this source
         *
         */

        using SourcePtr_t = std::decay_t<decltype(*sbegin)>;
        static_assert(is_art_ptr<SourcePtr_t>(),
          "Collection for AssociationFinder (FindManyInChainP) is not art pointers!"
          );

        // type of the source object (hidden in PtrColl)
        using Source_t = typename SourcePtr_t::value_type;
        using namespace AssociationFinderBase;

        Result_t<Source_t> result(KeysFrom, sbegin, send);
      //  std::size_t const nSources = result.size();

        // use this index for fast lookup of the sources
        auto match = makeSimpleIndex(sbegin, send);

        // fetch the association data product
        auto const& assns
          = *(event.template getValidHandle<Assns_t<Source_t>>(tag));

        // get the association, following the content of the assns data product
        for (decltype(auto) assn: assns) {
          SourcePtr_t const& sourcePtr = assn.first;

          // is this pointer interesting?
          auto iPos = match(sourcePtr);
          if (!iPos) continue; // it's not it

          // match! push target pointer into the result of the matched source
          result[*iPos].connectedTo().push_back(assn.second);

        } // for

        return result;
      } // findWithRange()


      /// Helper finding a single degree of associations from the tag of each
      /// source art pointer.
      template <typename PtrCollBegin, typename PtrCollEnd, typename Event>
      static auto findWithRange
        (PtrCollBegin sbegin, PtrCollEnd send, Event const& event)
      {
        /*
         * The strategy of this implementation is:
         *
         * 1. collect all the source art pointers, grouped by product ID
         *    (and sorted by key for faster lookup)
         * 2. for each interesting product ID:
         *    1. fetch the association collection; this is assumed to have been
         *       created with the same input tag as the source product
         *    2. parse all the associated pairs
         *       1. if the source pointer of a pair is in the list of
         *          interesting source pointers, push the target pointer of the
         *          pair into the results for this source
         *
         * The maximum complexity of this algorithm is N log(M), where M is no
         * larger than the maximum number of source pointers with a single
         * product ID and N is the number of associations in each association
         * data product.
         *
         */

        using SourcePtr_t = std::decay_t<decltype(*sbegin)>;
        static_assert(is_art_ptr<SourcePtr_t>(),
          "Collection for AssociationFinder (FindManyInChainP) is not art pointers!"
          );

        // type of the source object (hidden in PtrColl)
        using Source_t = typename SourcePtr_t::value_type;
        using namespace AssociationFinderBase;

        Result_t<Source_t> result(KeysFrom, sbegin, send);

        // These are the source pointers we still have to find,
        // organised by product ID; we keep track of the original position too
        // (unordered map may or may not be better...);
        // given that these are all from the same product ID, the product ID is
        // not saved.
        // Also, for fast lookup the lists are sorted by pointer key.
        std::map<art::ProductID, SourceVector_t<std::size_t>>
          sourcesLeft;
        std::size_t iSource = 0;
        for (auto it = sbegin; it != send; ++it, ++iSource)
          sourcesLeft[it->id()].emplace_back(it->key(), iSource);
        for (auto& sourcesWithinID: sourcesLeft)
          sourcesWithinID.second.sort();

        // look for all sources in each product ID
        for (auto const& sourcesWithinID: sourcesLeft) {

          art::ProductID sourceID = sourcesWithinID.first;
          auto const& sourceInfos = sourcesWithinID.second;

          // need to get the association between source and target,
          // as produced by the same producer that produced the source itself
          art::InputTag tag
            = tagFromProductID<std::vector<Source_t>>(sourceID, event);

          // fetch the association data product
          auto const& assns
            = *(event.template getValidHandle<Assns_t<Source_t>>(tag));

          // get the association, following the content of the assns data product
          for (decltype(auto) assn: assns) {
            SourcePtr_t const& sourcePtr = assn.first;

            // does this association contain a pointer with an ID different than
            // the one we are looking for? (it may mean our assumption that the
            // same module produced source and its association is wrong)
            if (sourcePtr.id() != sourceID) continue;

            // is this pointer interesting?
            auto iSrcInfo = sourceInfos.find(sourcePtr.key());
            if (!iSrcInfo) continue; // it's not it

            // match! push target pointer into the result of the matched source
            result[iSrcInfo->position()].connectedTo().push_back(assn.second);

          } // for

        } // for sourcesWithinID

        return result;
      } // findWithRange()


    }; // class AssociationFinder



    /// Helper finding a single degree of associations from a specified tag.
    template <
      typename Target, typename Handle, typename Event,
      typename = std::enable_if_t<is_handle_v<Handle>>
      >
    auto findAssociations
      (Handle&& handle, Event const& event, art::InputTag const& assnsTag)
    {
      using Finder = AssociationFinder<Target>;
      using Assns_t
        = typename Finder::template AssnsFromHandle_t<std::decay_t<Handle>>;

      auto const& assns = *(event.template getValidHandle<Assns_t>(assnsTag));

      // FIXME: this implementation is NOT canvas-compatible (Handle::id())
      return Finder::findFromDataProduct(handle.id(), handle->size(), assns);
    } // findAssociations(Handle, tag)


    /// Helper finding a single degree of associations from any tag.
    template<
      typename Target, typename Handle, typename Event,
      typename = std::enable_if_t<is_handle_v<Handle>>
      >
    auto findAssociations(Handle handle, Event const& event, lar::SameAsDataTag)
    {
      /*
       * By definition, here all the interesting sources share the same product
       * ID, which simplifies the implementation.
       */

      // need to get the association between source and target,
      // as produced by the same producer that produced the source itself
      // FIXME: the implementation of tagFromHandle is NOT canvas-compatible
      art::InputTag tag = tagFromHandle(handle);

      return findAssociations<Target>
        (std::forward<Handle>(handle), event, tag);
    } // findAssociations(Handle, same tag)


    template<
      typename Target, typename PtrColl, typename Event,
      typename = std::enable_if_t<!is_handle_v<PtrColl>>
      >
    auto findAssociations
      (PtrColl const& coll, Event const& event, lar::SameAsDataTag)
    {
      using std::cbegin;
      using std::cend;
      return AssociationFinder<Target>::findWithRange
        (cbegin(coll), cend(coll), event);
    }


    /// Helper finding a single degree of associations from a specified tag.
    template <
      typename Target, typename PtrColl, typename Event,
      typename = std::enable_if_t<!is_handle_v<PtrColl>>
      >
    auto findAssociations
      (PtrColl const& coll, Event const& event, art::InputTag const& tag)
    {
      using std::cbegin;
      using std::cend;
      return AssociationFinder<Target>::findWithRange
        (cbegin(coll), cend(coll), event, tag);
    } // findAssociations(Handle, tag)


    //--------------------------------------------------------------------------

    /**
     * @tparam Target type to be associated to the source data
     * @tparam Tier tier: 0 -> get target, 1 -> get leftmost intermediate, ...
     * @tparam Intermediate intermediate types, leftmost is closest to `Target`
     */
    template <typename Target, unsigned int Tier, typename... Intermediate>
    struct FindManyInChainPimpl {

      /// Total number of tiers (original source + all intermediates).
      static constexpr unsigned int Tiers = sizeof...(Intermediate) + 1;

      static_assert
        (Tier < Tiers, "Invalid tier: there is no type in position #N");

      /// Type used for the association to source.
      using Intermediate_t = get_type_t<(Tier - 1), Intermediate...>;

      template <typename Source, typename Event, typename InputTags>
      static auto find
        (Source&& source, Event const& event, InputTags const& tags)
      {
        static_assert(std::tuple_size<InputTags>() == Tiers,
          "Wrong number of input tags specified");

        // number of the intermediate being addressed
        constexpr auto nTag = Tiers - 1 - Tier;

        // find the associations between source and the next tier:
        // Source <==> Intermediate;
        // result is a ConnectionList object, that we store in a "manager"
        auto iq = makeConnectionManager(
          findAssociations<Intermediate_t>
            (source, event, std::get<nTag>(tags))
          );

        // collapse the result for input into the next tier;
        // this will result in a long sequence of art pointers
        auto const intermediateData = iq.allConnected();

        // process the next tier: Intermediate <==> Target;
        // this is also a ConnectionList
        auto oq
          = FindManyInChainPimpl<Target, (Tier - 1U), Intermediate...>::find
            (intermediateData, event, tags);

        // combine the result into jet another connection list:
        // Source <==> Intermediate (+) Intermediate <==> Target
        return iq.join(std::move(oq)); // AssociationFinderBase::combine(iq, oq);

      } // find()


    }; // FindManyInChainPimpl<>


    // Specialization for tier 0
    // (target association to first intermediate in list)
    template <typename Target, typename... Intermediate>
    struct FindManyInChainPimpl<Target, 0U, Intermediate...> {

      /// Total number of tiers (original source + all intermediates).
      static constexpr unsigned int Tiers = sizeof...(Intermediate) + 1;

      template <typename Source, typename Event, typename InputTags>
      static auto find
        (Source&& source, Event const& event, InputTags const& tags)
      {
        static_assert(std::tuple_size<InputTags>() == Tiers,
          "Wrong number of input tags specified");
        // find the associations between the last intermediate and the target
        // (using the last tag)
        return findAssociations<Target>
          (source, event, std::get<Tiers - 1>(tags));
      } // find()

    }; // FindManyInChainPimpl



    //--------------------------------------------------------------------------

  } // namespace details

} // namespace lar

//------------------------------------------------------------------------------
//
// Implementation note:
// - FindManyInChainP is the front-end utility with minimal template arguments;
//     it sets up and runs FindManyInChainPimpl...
// - FindManyInChainPimpl is the implementation of recursion where associations
//     are extracted for the results of the previous level ("tier") of
//     associated objects; it juggles with tags and parameters, and relies on
//     findAssociations()...
// - findAssociations() implements the extraction of a single association tier,
//     from a (flattened) collection of input objects
//
//------------------------------------------------------------------------------
template <typename Target, typename... Intermediate>
template <typename Source, typename Event, typename... InputTags>
auto lar::FindManyInChainP<Target, Intermediate...>::find
  (Source&& source, Event const& event, InputTags... tags)
  -> std::vector<TargetPtrCollection_t>
{
  constexpr auto Tiers = sizeof...(Intermediate) + 1U;

  // create a parameter pack with one tag per association
  auto const allTags
    = details::AssociationFinderBase::makeTagsTuple<Tiers>
    (SameAsData, std::forward<InputTags>(tags)...);

  return details::FindManyInChainPimpl<Target, (Tiers - 1), Intermediate...>
    ::find(source, event, allTags).values();

} // lar::FindManyInChainP<Target, Intermediate...>::FindManyInChainP()


//------------------------------------------------------------------------------
template <typename Target, typename... Intermediate>
std::size_t
lar::FindManyInChainP<Target, Intermediate...>::size() const noexcept
  { return results.size(); }


//------------------------------------------------------------------------------
template <typename Target, typename... Intermediate>
typename lar::FindManyInChainP<Target, Intermediate...>::TargetPtrCollection_t const&
lar::FindManyInChainP<Target, Intermediate...>::at(std::size_t i) const
  { return results.at(i); }


//------------------------------------------------------------------------------


#endif // LARDATA_UTILITIES_FINDMANYINCHAINP_TCC

// Local variables:
// mode: c++
// End:
