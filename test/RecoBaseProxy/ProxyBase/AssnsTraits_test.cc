/**
 * @file   AssnsTraits_test.cc
 * @brief  Unit tests on _art_ association trait utilities.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   April 17, 2018
 *
 * The relevant components of this test are static.
 */

// LArSoft libraries
#include "lardata/RecoBaseProxy/ProxyBase/AssnsTraits.h"

// framework libraries
#include "canvas/Persistency/Common/Assns.h"

// C/C++ standard libraries
#include <type_traits> // std::is_same<>, std::decay_t
#include <utility>     // std::declval()

// -----------------------------------------------------------------------------
// types used for the test
using PlainAssns = art::Assns<int, double>;
using MetaAssns = art::Assns<long int, float, char>;

// iterator types
using PlainAssnsIter = std::decay_t<decltype(std::declval<PlainAssns>().begin())>;
using MetaAssnsIter = std::decay_t<decltype(std::declval<MetaAssns>().begin())>;

// node types
using PlainAssnsNode = std::decay_t<decltype(*(std::declval<PlainAssnsIter>()))>;
using MetaAssnsNode = std::decay_t<decltype(*(std::declval<MetaAssnsIter>()))>;

// -----------------------------------------------------------------------------
// lar::util::assns_metadata_type
//
static_assert(std::is_same<lar::util::assns_metadata_type<PlainAssns>::type, void>(),
              "assns_metadata_type<PlainAssns> should have been void! ");
static_assert(std::is_same<lar::util::assns_metadata_type<PlainAssnsNode>::type, void>(),
              "assns_metadata_type<PlainAssns> node should have been void! ");
static_assert(std::is_same<lar::util::assns_metadata_type<MetaAssns>::type, char>(),
              "assns_metadata_type<MetaAssns> should have been char! ");
static_assert(std::is_same<lar::util::assns_metadata_type<MetaAssnsNode>::type, char>(),
              "assns_metadata_type<MetaAssns> node should have been char! ");

//
// lar::util::assns_metadata_t
//
static_assert(std::is_same<lar::util::assns_metadata_t<PlainAssns>, void>(),
              "assns_metadata_type_t<PlainAssns> should have been void! ");
static_assert(std::is_same<lar::util::assns_metadata_t<PlainAssnsNode>, void>(),
              "assns_metadata_type_t<PlainAssns> node should have been void! ");
static_assert(std::is_same<lar::util::assns_metadata_t<MetaAssns>, char>(),
              "assns_metadata_type_t<MetaAssns> should have been char! ");
static_assert(std::is_same<lar::util::assns_metadata_t<MetaAssnsNode>, char>(),
              "assns_metadata_type_t<MetaAssns> node should have been char! ");

// -----------------------------------------------------------------------------
// lar::util::assns_has_metadata
//
static_assert(!lar::util::assns_has_metadata<PlainAssns>(),
              "assns_has_metadata<PlainAssns> should have been false! ");
static_assert(!lar::util::assns_has_metadata<PlainAssnsNode>(),
              "assns_has_metadata<PlainAssns> node should have been false! ");
static_assert(lar::util::assns_has_metadata<MetaAssns>(),
              "assns_has_metadata<MetaAssns> should have been true! ");
static_assert(lar::util::assns_has_metadata<MetaAssnsNode>(),
              "assns_has_metadata<MetaAssns> node should have been true! ");

//
// lar::util::assns_has_metadata_v
//
static_assert(!lar::util::assns_has_metadata_v<PlainAssns>,
              "assns_has_metadata<PlainAssns> should have been false! ");
static_assert(!lar::util::assns_has_metadata_v<PlainAssnsNode>,
              "assns_has_metadata<PlainAssns> node should have been false! ");
static_assert(lar::util::assns_has_metadata_v<MetaAssns>,
              "assns_has_metadata<MetaAssns> should have been true! ");
static_assert(lar::util::assns_has_metadata_v<MetaAssnsNode>,
              "assns_has_metadata<MetaAssns> node should have been true! ");

//------------------------------------------------------------------------------
// lar::util::assns_traits
//
using PlainAssnsTraits = lar::util::assns_traits<PlainAssns>;
static_assert(std::is_same<PlainAssnsTraits::left_t, int>(),
              "left element of PlainAssns should have been `int`!");
static_assert(std::is_same<PlainAssnsTraits::right_t, double>(),
              "right element of PlainAssns should have been `double`!");
static_assert(std::is_same<PlainAssnsTraits::data_t, void>(),
              "metadata of PlainAssns should have been `void`!");
static_assert(std::is_same<PlainAssnsTraits::leftptr_t, art::Ptr<int>>(),
              "left element pointer of PlainAssns should have been `art::Ptr<int>`!");
static_assert(std::is_same<PlainAssnsTraits::rightptr_t, art::Ptr<double>>(),
              "right element pointer of PlainAssns should have been `art::Ptr<double>`!");
static_assert(std::is_same<PlainAssnsTraits::dataptr_t, void const*>(),
              "metadata pointer of PlainAssns should have been `void const*`!");
static_assert(std::is_same<PlainAssnsTraits::assns_t, PlainAssns>(),
              "PlainAssns traits should report `PlainAssns` as association class!");
static_assert(std::is_same<PlainAssnsTraits::assns_iterator_t, PlainAssnsIter>(),
              "PlainAssns traits should report `PlainAssnsIter` as association iterator!");
static_assert(std::is_same<PlainAssnsTraits::art_assns_node_t, PlainAssnsNode>(),
              "PlainAssns traits should report `PlainAssnsNode` as association node!");
static_assert(!PlainAssnsTraits::hasMetadata,
              "PlainAssns should have reported not to have metadata!");

using MetaAssnsTraits = lar::util::assns_traits<MetaAssns>;
static_assert(std::is_same<MetaAssnsTraits::left_t, long int>(),
              "left element of MetaAssns should have been `long int`!");
static_assert(std::is_same<MetaAssnsTraits::right_t, float>(),
              "right element of MetaAssns should have been `float`!");
static_assert(std::is_same<MetaAssnsTraits::data_t, char>(),
              "metadata of MetaAssns should have been `char`!");
static_assert(std::is_same<MetaAssnsTraits::leftptr_t, art::Ptr<long int>>(),
              "left element pointer of MetaAssns should have been `art::Ptr<long int>`!");
static_assert(std::is_same<MetaAssnsTraits::rightptr_t, art::Ptr<float>>(),
              "right element pointer of MetaAssns should have been `art::Ptr<float>`!");
static_assert(std::is_same<MetaAssnsTraits::dataptr_t, char const*>(),
              "metadata pointer of MetaAssns should have been `char const*`!");
static_assert(std::is_same<MetaAssnsTraits::assns_t, MetaAssns>(),
              "MetaAssns traits should report `MetaAssns` as association class!");
static_assert(std::is_same<MetaAssnsTraits::assns_iterator_t, MetaAssnsIter>(),
              "MetaAssns traits should report `MetaAssnsIter` as association iterator!");
static_assert(std::is_same<MetaAssnsTraits::art_assns_node_t, MetaAssnsNode>(),
              "MetaAssns traits should report `MetaAssnsNode` as association node!");
static_assert(MetaAssnsTraits::hasMetadata, "MetaAssns should have reported to have metadata!");

//------------------------------------------------------------------------------
// lar::util::assns_iterator_type
//
static_assert(std::is_same<lar::util::assns_iterator_type<PlainAssns>::type, PlainAssnsIter>(),
              "assns_metadata_type<PlainAssns> should have been void! ");
static_assert(std::is_same<lar::util::assns_iterator_type<MetaAssns>::type, MetaAssnsIter>(),
              "assns_metadata_type<MetaAssns> should have been char! ");

//------------------------------------------------------------------------------
// lar::util::assns_iterator_t
//
static_assert(std::is_same<lar::util::assns_iterator_t<PlainAssns>, PlainAssnsIter>(),
              "assns_metadata_type<PlainAssns> should have been void! ");
static_assert(std::is_same<lar::util::assns_iterator_t<MetaAssns>, MetaAssnsIter>(),
              "assns_metadata_type<MetaAssns> should have been char! ");

//==============================================================================
// we need a main program because every test that compiles must run
int main() {}

//------------------------------------------------------------------------------
