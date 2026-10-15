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
 *                                     grouped by process, newest first,
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
   * Constructed with the process history, oldest process first.
   * Configure it with `add(tag, valid)`: `valid == false` models a product
   * whose tag is advertised by `getInputTags<T>()` but whose `getHandle<T>()`
   * comes back invalid (e.g. dropped on input) and must be filtered out.
   */
  class MockEvent {
  public:
    template <typename T>
    using HandleT = MockHandle<T>;

    explicit MockEvent(std::vector<std::string> processHistory)
      : processHistory_{std::move(processHistory)}
    {}

    void add(art::InputTag tag, bool valid = true) { products_.push_back({std::move(tag), valid}); }

    /**
     * Like art (`Principal::findGroups()`), hands the selector a
     * `BranchDescription` for each product, walking the process history in
     * reverse: the tags come out grouped by process, newest process first,
     * independently of the order products were `add()`-ed.
     */
    template <typename T>
    std::vector<art::InputTag> getInputTags(art::SelectorBase const& selector) const
    {
      std::vector<art::InputTag> tags;
      for (auto proc = processHistory_.crbegin(); proc != processHistory_.crend(); ++proc) {
        for (auto const& p : products_) {
          if (p.tag.process() != *proc) continue;
          art::BranchDescription const bd{
            art::InEvent,
            art::TypeLabel{art::TypeID{typeid(T)}, p.tag.instance(), /* supportsView = */ false},
            p.tag.label(),
            fhicl::ParameterSetID{},
            art::ProcessConfiguration{p.tag.process(), fhicl::ParameterSetID{}, ""}};
          if (selector.match(bd)) tags.push_back(p.tag);
        }
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
    std::vector<std::string> processHistory_;
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

  /// Shorthand: encoded tags selected from `e` by the pattern `tag`.
  std::set<std::string> select(MockEvent const& e, art::InputTag const& tag)
  {
    return selected(lar::util::getManyByRegexTag<MyProd>(e, tag));
  }

  /// The event of the examples in the `getManyByRegexTag()` documentation.
  MockEvent makeEvent()
  {
    MockEvent e{{"reco1", "reco2"}}; // reco2 is the newest
    // older process first, to show insertion order is irrelevant
    e.add({"gaushit", "", "reco1"});
    e.add({"gaushit", "uncalib", "reco1"});
    e.add({"gaushitStage2", "", "reco1"});
    e.add({"gaushit", "", "reco2"});
    e.add({"fasthit", "", "reco2"});
    return e;
  }

} // local namespace

//------------------------------------------------------------------------------
// Pattern matching: one test per row of the documentation table.
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(EmptyProcessPicksNewestProcess)
{
  // gaushit::reco1 also matches, but reco2 is newer.
  BOOST_TEST((select(makeEvent(), {"gaushit", "", ""}) == std::set<std::string>{"gaushit::reco2"}));
}

BOOST_AUTO_TEST_CASE(LabelExactMatchIsAnchored)
{
  // "gaushit" must NOT select "gaushitStage2" (std::regex_match is a full match).
  BOOST_TEST(
    (select(makeEvent(), {"gaushit", "", "reco1"}) == std::set<std::string>{"gaushit::reco1"}));
}

BOOST_AUTO_TEST_CASE(LabelWildcardMatch)
{
  std::set<std::string> const expected{"gaushit::reco1", "gaushitStage2::reco1"};
  BOOST_TEST((select(makeEvent(), {"gaushit.*", "", "reco1"}) == expected));
}

BOOST_AUTO_TEST_CASE(OlderProcessProductsAreIgnored)
{
  // gaushitStage2 only exists in reco1: it is dropped, because reco2 has a match.
  BOOST_TEST(
    (select(makeEvent(), {"gaushit.*", "", ""}) == std::set<std::string>{"gaushit::reco2"}));
}

BOOST_AUTO_TEST_CASE(AnyInstanceNewestProcess)
{
  // gaushit:uncalib only exists in reco1: it is dropped, because reco2 has a match.
  BOOST_TEST(
    (select(makeEvent(), {"gaushit", ".*", ""}) == std::set<std::string>{"gaushit::reco2"}));
}

BOOST_AUTO_TEST_CASE(EmptyInstanceMatchesOnlyEmptyInstance)
{
  // "" selects gaushit::reco1 only; ".*" adds the "uncalib" instance.
  BOOST_TEST(
    (select(makeEvent(), {"gaushit", "", "reco1"}) == std::set<std::string>{"gaushit::reco1"}));
  std::set<std::string> const expected{"gaushit::reco1", "gaushit:uncalib:reco1"};
  BOOST_TEST((select(makeEvent(), {"gaushit", ".*", "reco1"}) == expected));
}

BOOST_AUTO_TEST_CASE(NewestProcessWithAMatch)
{
  // reco2 is newer but has no match: reco1 is the newest process *with a match*.
  BOOST_TEST((select(makeEvent(), {"gaushit", "uncalib", ""}) ==
              std::set<std::string>{"gaushit:uncalib:reco1"}));
}

BOOST_AUTO_TEST_CASE(AnyLabelNewestProcess)
{
  std::set<std::string> const expected{"gaushit::reco2", "fasthit::reco2"};
  BOOST_TEST((select(makeEvent(), {".*", "", ""}) == expected));
}

BOOST_AUTO_TEST_CASE(ProcessRegex)
{
  // "reco.*" matches both processes: same as an empty process pattern.
  std::set<std::string> const expected{"gaushit::reco2", "fasthit::reco2"};
  BOOST_TEST((select(makeEvent(), {".*", "", "reco.*"}) == expected));
}

BOOST_AUTO_TEST_CASE(AllEmptyPatternsSelectNothing)
{
  // An empty label pattern matches only an empty label: no product has one.
  BOOST_TEST(select(makeEvent(), art::InputTag{}).empty());
}

//------------------------------------------------------------------------------
// Validity of the handles.
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(InvalidHandlesAreFilteredOut)
{
  MockEvent e{{"reco1"}};
  e.add({"gaushit", "", "reco1"});
  e.add({"brokenhit", "", "reco1"}, /* valid = */ false);

  BOOST_TEST((select(e, {".*", "", ""}) == std::set<std::string>{"gaushit::reco1"}));
}

BOOST_AUTO_TEST_CASE(InvalidNewestProductDoesNotFallBack)
{
  // The newest matching process is chosen before checking validity:
  // the valid reco1 product is not used in place of the invalid reco2 one.
  MockEvent e{{"reco1", "reco2"}};
  e.add({"gaushit", "", "reco1"});
  e.add({"gaushit", "", "reco2"}, /* valid = */ false);

  BOOST_TEST(select(e, {"gaushit", "", ""}).empty());
}

//------------------------------------------------------------------------------
// Edge cases.
//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(NoMatchReturnsEmpty)
{
  BOOST_TEST(select(makeEvent(), {"does_not_exist", "", ""}).empty());
}

BOOST_AUTO_TEST_CASE(InvalidRegexThrows)
{
  // An unbalanced bracket is not a valid std::regex.
  BOOST_CHECK_THROW(select(makeEvent(), {"[", "", ""}), std::regex_error);
}

//------------------------------------------------------------------------------
