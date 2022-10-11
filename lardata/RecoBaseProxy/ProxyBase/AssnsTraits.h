/**
 * @file   lardata/RecoBaseProxy/ProxyBase/AssnsTraits.h
 * @brief  Traits for _art_ associations.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_ASSNSTRAITS_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_ASSNSTRAITS_H

// LArSoft libraries
#include "larcorealg/CoreUtils/MetaUtils.h" // util::is_not_same<>

// framework libraries
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"

// C/C++ standard
#include <utility> // std::pair

namespace lar {

  namespace util {

    //--- BEGIN Traits for art associations ------------------------------------
    /**
     * @defgroup ArtAssociationsMetaprogramming Traits for art associations
     * @ingroup Metaprogramming
     *
     * A small set of metaprogramming classes is meant to hide the small
     * differences in interface currently in `art::Assns` (_art_ 2.10).
     *
     * These classes are all contained in the `lar::util` namespace.
     */
    /// @{

    //--------------------------------------------------------------------------
    // was proxy::details::hasMetadata
    /// Trait: `value` true if `Assns` (association or its node) has metadata.
    template <typename Assns>
    struct assns_has_metadata;

    // was: proxy::details::AssnWithMetadata_v
    /// Trait: true if `Assns` (association or its node) has metadata.
    template <typename Assns>
    constexpr bool assns_has_metadata_v = assns_has_metadata<Assns>::value;

    // was: proxy::details::AssnsMetadataTypeStruct, proxy::details::AssnNodeMetadataType
    /// Trait: `type` is metadata in `Assns` (association or its node).
    template <typename Assns>
    struct assns_metadata_type;

    // was: proxy::details::AssnsMetadata_t, proxy::details::AssnNodeMetadata_t
    /// Trait: type of metadata in `Assns` (association or its node).
    template <typename Assns>
    using assns_metadata_t = typename assns_metadata_type<Assns>::type;

    //--------------------------------------------------------------------------
    // FIXME simplify this code if issue #18769 is accepted
    // was: proxy::details::AssnsIteratorTypeStruct
    /// Trait: `type` is iterator of `Assns`.
    template <typename Assns>
    struct assns_iterator_type;

    /// Trait: type of iterator of `Assns`.
    // was: proxy::details::AssnsIterator_t
    template <typename Assns>
    using assns_iterator_t = typename assns_iterator_type<Assns>::type;

    //--------------------------------------------------------------------------
    // was: proxy::details::AssnsNodeTraits
    /**
     * @brief Data types for the specified association type (or its node).
     * @tparam Assns _art_ association type, or its node
     *
     * The trait class is expected to provide the following types:
     * * `left_t`: type at the left side of the association
     * * `right_t`: type at the right side of the association
     * * `data_t`: type of data bound to the association
     * * `leftptr_t`: _art_ pointer to the left side
     * * `rightptr_t`: _art_ pointer to the right side
     * * `dataptr_t`: pointer to the bound metadata
     * * `assns_t`: type of the _art_ association object
     * * `art_assns_node_t`: node in the associations list, representing a
     *      single connection between two objects and its metadata
     *
     * and the following constant value:
     *
     * * `hasMetadata`: shortcut to know whether this node supports any metadata
     *
     * The trait class is not defined for types other than `art::Assns` and its
     * node type.
     */
    template <typename Assns>
    struct assns_traits;

    //--------------------------------------------------------------------------

    ///@}
    //--- END Traits for art associations --------------------------------------

  } // namespace util

} // namespace lar

//------------------------------------------------------------------------------
//---  template implementation
//------------------------------------------------------------------------------
namespace lar {

  namespace util {

    namespace details {

      // was: proxy::details::isAssnMetadata
      template <typename T>
      using isAssnMetadata = ::util::is_not_same<T, void>;

      template <typename Assns>
      struct node_of;

      template <typename Assns>
      using node_of_t = typename node_of<Assns>::type;

      template <typename L, typename R, typename D>
      struct node_of<art::Assns<L, R, D>> {
        using type = art::AssnsNode<L, R, D>;
      };

      template <typename L, typename R>
      struct node_of<art::Assns<L, R, void>> {
        using type = std::pair<art::Ptr<L>, art::Ptr<R>>;
      };

    } // namespace details

    //--------------------------------------------------------------------------
    template <typename L, typename R, typename D>
    struct assns_metadata_type<art::Assns<L, R, D>> {
      using type = D;
    };

    template <typename L, typename R, typename D>
    struct assns_metadata_type<art::AssnsNode<L, R, D>> {
      using type = D;
    };

    template <typename L, typename R>
    struct assns_metadata_type<std::pair<art::Ptr<L>, art::Ptr<R>>> {
      using type = void;
    };

    //--------------------------------------------------------------------------
    template <typename Assns>
    struct assns_has_metadata : details::isAssnMetadata<assns_metadata_t<Assns>> {};

    //--------------------------------------------------------------------------
    template <typename L, typename R, typename D>
    struct assns_iterator_type<art::Assns<L, R, D>> {
      using type = typename art::Assns<L, R, D>::const_iterator;
    }; // struct assns_iterator_type

    template <typename L, typename R>
    struct assns_iterator_type<art::Assns<L, R, void>> {
      using type = typename art::Assns<L, R, void>::const_iterator;
    };

    //--------------------------------------------------------------------------
    // was proxy::details::AssnsNodeTraitsBase
    template <typename L, typename R, typename D>
    struct assns_traits<art::Assns<L, R, D>> {

      using left_t = L;                     ///< Type at the left side of the association.
      using right_t = R;                    ///< Type at the right side of the association.
      using data_t = D;                     ///< Type of data bound to the association.
      using leftptr_t = art::Ptr<left_t>;   ///< Art pointer to the left side.
      using rightptr_t = art::Ptr<right_t>; ///< Art pointer to the right side.
      using dataptr_t = data_t const*;      ///< Pointer to the bound metadata.

      using assns_t = art::Assns<L, R, D>; ///< Type of the association.

      /// Type of the association iterator.
      using assns_iterator_t = typename assns_iterator_type<assns_t>::type;

      /// Type of a node (element) in the association.
      using art_assns_node_t = details::node_of_t<assns_t>;

      /// Shortcut to know whether this node supports any metadata,
      static constexpr bool hasMetadata = details::isAssnMetadata<data_t>();

    }; // struct assns_traits<Assns<L,R,D>>

    template <typename L, typename R, typename D>
    struct assns_traits<art::AssnsNode<L, R, D>> : public assns_traits<art::Assns<L, R, D>> {};

    template <typename L, typename R>
    struct assns_traits<std::pair<art::Ptr<L>, art::Ptr<R>>>
      : public assns_traits<art::Assns<L, R>> {};

    //--------------------------------------------------------------------------

  } // namespace util

} // namespace lar

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_ASSNSTRAITS_H
