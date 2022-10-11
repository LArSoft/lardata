/**
 * @file   lardata/RecoBaseProxy/ProxyBase/AssnsNodeAsTuple.h
 * @brief  Specializations of STL tuple utilities for `art::AssnsNode`.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @see    lardata/RecoBaseProxy/ProxyBase.h
 *
 * This library is header-only.
 */

#ifndef LARDATA_RECOBASEPROXY_PROXYBASE_ASSNSNODEASTUPLE_H
#define LARDATA_RECOBASEPROXY_PROXYBASE_ASSNSNODEASTUPLE_H

// framework libraries
#include "canvas/Persistency/Common/Assns.h" // art::AssnsNode
#include "canvas/Persistency/Common/Ptr.h"

// C/C++ standard
#include <cstdlib> // std::size_t
#include <tuple>   // std::tuple_element_t<>, std::get()

// FIXME simplify this code if issue #18769 is accepted
namespace std {

  //----------------------------------------------------------------------------
  //--- specializations of std::tuple interface for art::AssnsNode
  //----------------------------------------------------------------------------

  // specialize for indices 0, 1, and 2; for all others, it's an incomplete type
  template <typename L, typename R, typename D>
  class tuple_element<0U, art::AssnsNode<L, R, D>> {
  public:
    using type = art::Ptr<L>;
  };

  template <typename L, typename R, typename D>
  class tuple_element<1U, art::AssnsNode<L, R, D>> {
  public:
    using type = art::Ptr<R>;
  };

  template <typename L, typename R, typename D>
  class tuple_element<2U, art::AssnsNode<L, R, D>> {
  public:
    using type = D const*;
  };

  //----------------------------------------------------------------------------
  template <std::size_t I, typename L, typename R, typename D>
  constexpr std::tuple_element_t<I, art::AssnsNode<L, R, D>>& get(
    art::AssnsNode<L, R, D>& t) noexcept;

  template <std::size_t I, typename L, typename R, typename D>
  constexpr std::tuple_element_t<I, art::AssnsNode<L, R, D>>&& get(
    art::AssnsNode<L, R, D>&& t) noexcept;

  template <std::size_t I, typename L, typename R, typename D>
  constexpr std::tuple_element_t<I, art::AssnsNode<L, R, D>> const& get(
    const art::AssnsNode<L, R, D>& t) noexcept;

  template <std::size_t I, typename L, typename R, typename D>
  constexpr std::tuple_element_t<I, art::AssnsNode<L, R, D>> const&& get(
    const art::AssnsNode<L, R, D>&& t) noexcept;

  template <class T, typename L, typename R, typename D>
  constexpr T& get(art::AssnsNode<L, R, D>& t) noexcept;

  template <class T, typename L, typename R, typename D>
  constexpr T&& get(art::AssnsNode<L, R, D>&& t) noexcept;

  template <class T, typename L, typename R, typename D>
  constexpr const T& get(const art::AssnsNode<L, R, D>& t) noexcept;

  template <class T, typename L, typename R, typename D>
  constexpr const T&& get(const art::AssnsNode<L, R, D>&& t) noexcept;

  //----------------------------------------------------------------------------

} // namespace std

//------------------------------------------------------------------------------
//--- implementation
//---
namespace util {
  namespace details {

    //
    // support utilities for implementation (they live outside std namespace)
    //

    template <std::size_t I, typename L, typename R, typename D>
    struct AssnsNodeGetter; // incomplete type, except for specializations...

    template <typename L, typename R, typename D>
    struct AssnsNodeGetter<0U, L, R, D> {

      using AssnsNode_t = art::AssnsNode<L, R, D>;
      using Element_t = std::tuple_element_t<0U, AssnsNode_t>;

      static constexpr Element_t& get(AssnsNode_t& node) noexcept { return node.first; }

      static constexpr Element_t const& get(AssnsNode_t const& node) noexcept { return node.first; }

      static constexpr Element_t&& get(AssnsNode_t&& node) noexcept
      {
        return std::move(node.first);
      }

      static constexpr Element_t const&& get(AssnsNode_t const&& node) noexcept
      {
        return std::move(node.first);
      }

    }; // struct AssnsNodeGetter<0U>

    template <typename L, typename R, typename D>
    struct AssnsNodeGetter<1U, L, R, D> {

      using AssnsNode_t = art::AssnsNode<L, R, D>;
      using Element_t = std::tuple_element_t<1U, AssnsNode_t>;

      static constexpr Element_t& get(AssnsNode_t& node) noexcept { return node.second; }

      static constexpr Element_t const& get(AssnsNode_t const& node) noexcept
      {
        return node.second;
      }

      static constexpr Element_t&& get(AssnsNode_t&& node) noexcept
      {
        return std::move(node.second);
      }

      static constexpr Element_t const&& get(AssnsNode_t const&& node) noexcept
      {
        return std::move(node.second);
      }

    }; // struct AssnsNodeGetter<1U>

    template <typename L, typename R, typename D>
    struct AssnsNodeGetter<2U, L, R, D> {

      using AssnsNode_t = art::AssnsNode<L, R, D>;
      using Element_t = std::tuple_element_t<2U, AssnsNode_t>;

      static constexpr Element_t& get(AssnsNode_t& node) noexcept { return node.data; }

      static constexpr Element_t const& get(AssnsNode_t const& node) noexcept { return node.data; }

      static constexpr Element_t&& get(AssnsNode_t&& node) noexcept { return std::move(node.data); }

      static constexpr Element_t const&& get(AssnsNode_t const&& node) noexcept
      {
        return std::move(node.data);
      }

    }; // struct AssnsNodeGetter<2U>

  } // namespace details
} // namespace util

namespace std {

  //----------------------------------------------------------------------------
  //--- implementation of specializations of std::get() for art::AssnsNode
  //----------------------------------------------------------------------------

  template <std::size_t I, typename L, typename R, typename D>
  constexpr std::tuple_element_t<I, art::AssnsNode<L, R, D>>& get(
    art::AssnsNode<L, R, D>& node) noexcept
  {
    return util::details::AssnsNodeGetter<I, L, R, D>::get(node);
  }

  template <std::size_t I, typename L, typename R, typename D>
  constexpr std::tuple_element_t<I, art::AssnsNode<L, R, D>>&& get(
    art::AssnsNode<L, R, D>&& node) noexcept
  {
    return util::details::AssnsNodeGetter<I, L, R, D>::get(std::move(node));
  }

  template <std::size_t I, typename L, typename R, typename D>
  constexpr std::tuple_element_t<I, art::AssnsNode<L, R, D>> const& get(
    art::AssnsNode<L, R, D> const& node) noexcept
  {
    return util::details::AssnsNodeGetter<I, L, R, D>::get(node);
  }

  template <std::size_t I, typename L, typename R, typename D>
  constexpr std::tuple_element_t<I, art::AssnsNode<L, R, D>> const&& get(
    art::AssnsNode<L, R, D> const&& node) noexcept
  {
    return util::details::AssnsNodeGetter<I, L, R, D>::get(std::move(node));
  }

  // not implemented yet:
  template <class T, typename L, typename R, typename D>
  constexpr T& get(art::AssnsNode<L, R, D>& t) noexcept;

  template <class T, typename L, typename R, typename D>
  constexpr T&& get(art::AssnsNode<L, R, D>&& t) noexcept;

  template <class T, typename L, typename R, typename D>
  constexpr const T& get(const art::AssnsNode<L, R, D>& t) noexcept;

  template <class T, typename L, typename R, typename D>
  constexpr const T&& get(const art::AssnsNode<L, R, D>&& t) noexcept;

  //----------------------------------------------------------------------------

} // namespace std

// -----------------------------------------------------------------------------

#endif // LARDATA_RECOBASEPROXY_PROXYBASE_ASSNSNODEASTUPLE_H
