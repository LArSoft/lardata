/**
 * @file    getManyByRegexTag_test.cc
 * @brief   Unit test for `lar::util::getManyByRegexTag()`.
 * @see     `lardata/ArtDataHelper/GetManyByRegexTag.h`
 *
 * See http://www.boost.org/libs/test for the Boost test library home page.
 *
 * `getManyByRegexTag()` is a function template on the event type `Evt`,
 * duck-typed to a very small interface:
 *   - `Evt::template HandleT<T>`      a product handle type,
 *   - `e.template getInputTags<T>(selector)`  the input tags of type `T`
 *                                     accepted by an `art::SelectorBase`,
 *   - `e.template getHandle<T>(tag)`  a handle, exposing `isValid()`.
 *
 * That lets us validate the real filtering logic (including the real
 * `lar::util::detail::RegexTagSelector`) against a lightweight mock event,
 * with no art job required.
 */

// Boost libraries
#define BOOST_TEST_MODULE (getManyByRegexTag_test)
#include "boost/test/unit_test.hpp"

// LArSoft libraries
#include "lardata/ArtDataHelper/GetManyByRegexTag.h"

// framework libraries
#include "art/Framework/Principal/SelectorBase.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Persistency/Provenance/BranchType.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/TypeLabel.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/TypeID.h"
#include "fhiclcpp/ParameterSetID.h"

// C/C++ standard libraries
#include <regex> // std::regex_error
#include <set>
#include <string>
#include <utility> // std::move
#include <vector>

namespace {

  //----------------------------------------------------------------------------
  /// A stand-in data product; its contents are irrelevant to the algorithm.
  struct MyProd {};

  /**
   * @brief Minimal handle fulfilling the interface `getManyByRegexTag()` uses.
   *
   * The algorithm only ever calls `isValid()`; the `tag` is carried along so
   * the test can inspect *which* products were selected.
   */
  template <typename T>
  struct MockHandle {
    art::InputTag tag{};
    bool valid{true};
    bool isValid() const { return valid; }
  };

  /**
   * @brief Minimal event duck-typed to the `Evt` concept the algorithm needs.
   *
   * Configure it with `add(tag, valid)`: `valid == false` models a product
   * whose tag is advertised by `getInputTags<T>()` but whose `getHandle<T>()`
   * comes back invalid (e.g. dropped on input) and must be filtered out.
   */
  class MockEvent {
  public:
    template <typename T>
    using HandleT = MockHandle<T>;

    void add(art::InputTag tag, bool valid = true) { products_.push_back({std::move(tag), valid}); }

    /// Like art, hands the selector a `BranchDescription` for each product.
    template <typename T>
    std::vector<art::InputTag> getInputTags(art::SelectorBase const& selector) const
    {
      std::vector<art::InputTag> tags;
      for (auto const& p : products_) {
        art::BranchDescription const bd{
          art::InEvent,
          art::TypeLabel{art::TypeID{typeid(T)}, p.tag.instance(), /* supportsView = */ false},
          p.tag.label(),
          fhicl::ParameterSetID{},
          art::ProcessConfiguration{p.tag.process(), fhicl::ParameterSetID{}, ""}};
        if (selector.match(bd)) tags.push_back(p.tag);
      }
      return tags;
    }

    template <typename T>
    HandleT<T> getHandle(art::InputTag const& t) const
    {
      for (auto const& p : products_) {
        if (p.tag == t) return MockHandle<T>{p.tag, p.valid};
      }
      return MockHandle<T>{t, false};
    }

  private:
    struct Entry {
      art::InputTag tag;
      bool valid;
    };
    std::vector<Entry> products_;
  };

  //----------------------------------------------------------------------------
  /// Encoded ("label:instance:process") tags of the selected handles, as a set.
  template <typename Handles>
  std::set<std::string> selected(Handles const& handles)
  {
    std::set<std::string> s;
    for (auto const& h : handles)
      s.insert(h.tag.encode());
    return s;
  }

  /// Event shared by the matching test cases.
  MockEvent makeEvent()
  {
    MockEvent e;
    e.add({"gaushit", "", "reco1"});                        // 0
    e.add({"gaushit", "", "reco2"});                        // 1
    e.add({"gaushitStage2", "", "reco1"});                  // 2
    e.add({"fasthit", "", "reco1"});                        // 3
    e.add({"gaushit", "uncalib", "reco1"});                 // 4
    e.add({"brokenhit", "", "reco1"}, /* valid = */ false); // 5: filtered out
    return e;
  }

} // local namespace

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(EmptyTagMatchesEveryValidProduct)
{
  auto const e = makeEvent();

  // Empty fields => ".*" => match anything. The only product dropped is the
  // one whose handle is invalid ("brokenhit").
  auto const handles = lar::util::getManyByRegexTag<MyProd>(e, art::InputTag{});

  BOOST_TEST(handles.size() == 5u);
  std::set<std::string> const expected{"gaushit::reco1",
                                       "gaushit::reco2",
                                       "gaushitStage2::reco1",
                                       "fasthit::reco1",
                                       "gaushit:uncalib:reco1"};
  BOOST_TEST((selected(handles) == expected));
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(LabelExactMatchIsAnchored)
{
  auto const e = makeEvent();

  // std::regex_match is a *full* match: "gaushit" must NOT select
  // "gaushitStage2".
  auto const handles = lar::util::getManyByRegexTag<MyProd>(e, art::InputTag{"gaushit", "", ""});

  std::set<std::string> const expected{"gaushit::reco1", "gaushit::reco2", "gaushit:uncalib:reco1"};
  BOOST_TEST((selected(handles) == expected));
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(LabelWildcardMatch)
{
  auto const e = makeEvent();

  auto const handles = lar::util::getManyByRegexTag<MyProd>(e, art::InputTag{"gaushit.*", "", ""});

  // Now "gaushitStage2" is included as well; "fasthit" still is not.
  std::set<std::string> const expected{
    "gaushit::reco1", "gaushit::reco2", "gaushitStage2::reco1", "gaushit:uncalib:reco1"};
  BOOST_TEST((selected(handles) == expected));
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(InstanceFilter)
{
  auto const e = makeEvent();

  // Only the product carrying a (non-empty) matching instance is selected.
  auto const handles = lar::util::getManyByRegexTag<MyProd>(e, art::InputTag{"", "uncalib", ""});

  std::set<std::string> const expected{"gaushit:uncalib:reco1"};
  BOOST_TEST((selected(handles) == expected));
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(ProcessFilter)
{
  auto const e = makeEvent();

  auto const handles = lar::util::getManyByRegexTag<MyProd>(e, art::InputTag{"", "", "reco2"});

  std::set<std::string> const expected{"gaushit::reco2"};
  BOOST_TEST((selected(handles) == expected));
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(AllThreeFieldsMustMatch)
{
  auto const e = makeEvent();

  // label "gaushit" (anchored) AND process "reco1".
  auto const handles =
    lar::util::getManyByRegexTag<MyProd>(e, art::InputTag{"gaushit", "", "reco1"});

  std::set<std::string> const expected{"gaushit::reco1", "gaushit:uncalib:reco1"};
  BOOST_TEST((selected(handles) == expected));
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(InvalidHandlesAreFilteredOut)
{
  auto const e = makeEvent();

  // "brokenhit" matches the pattern but its handle is invalid => dropped.
  auto const handles = lar::util::getManyByRegexTag<MyProd>(e, art::InputTag{"brokenhit", "", ""});

  BOOST_TEST(handles.empty());
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(NoMatchReturnsEmpty)
{
  auto const e = makeEvent();

  // NOTE: GetManyByRegexTag.h documents "@throws std::runtime_error if no
  // product of type T matches", but the implementation returns an empty
  // vector instead. This test pins the *actual* behavior. If the intended
  // contract is to throw, change this to BOOST_CHECK_THROW and fix the
  // implementation accordingly.
  auto const handles =
    lar::util::getManyByRegexTag<MyProd>(e, art::InputTag{"does_not_exist", "", ""});

  BOOST_TEST(handles.empty());
}

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(InvalidRegexThrows)
{
  auto const e = makeEvent();

  // An unbalanced bracket is not a valid std::regex.
  BOOST_CHECK_THROW(lar::util::getManyByRegexTag<MyProd>(e, art::InputTag{"[", "", ""}),
                    std::regex_error);
}

//------------------------------------------------------------------------------
