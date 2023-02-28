////////////////////////////////////////////////////////////////////////
// Class:       AssnAnalyzer
// Plugin Type: analyzer (art v2_05_00)
// File:        AssnAnalyzer_module.cc
//
// Generated at Fri Dec  9 00:12:59 2016 by Saba Sehrish using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/AssnsAlgorithms.h" // art::for_each_group()
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardata/Utilities/ForEachAssociatedGroup.h"

#include <algorithm> // std::copy()
#include <iterator>  // std::inserter()
#include <set>

class AssnAnalyzer;

class AssnAnalyzer : public art::EDAnalyzer {
public:
  typedef std::vector<int> intvec_t;
  typedef std::vector<std::string> strvec_t;
  typedef art::Assns<std::string, int> strintAssns_t;

  explicit AssnAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AssnAnalyzer(AssnAnalyzer const&) = delete;
  AssnAnalyzer(AssnAnalyzer&&) = delete;
  AssnAnalyzer& operator=(AssnAnalyzer const&) = delete;
  AssnAnalyzer& operator=(AssnAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  art::InputTag fInputLabel;

  std::set<std::string> fEnabledTests;

  void for_each_associated_group_test(art::Event const& e) const;
  void associated_groups_test(art::Event const& e) const;
  void associated_groups_with_left_test(art::Event const& e) const;
};

namespace {

  bool starts_with(std::string const& s, std::string const& key)
  {
    return s.substr(0, key.length()) == key;
  }

} // local namespace

AssnAnalyzer::AssnAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer(p), fInputLabel(p.get<art::InputTag>("input_label"))
{
  auto enableTests = p.get<std::vector<std::string>>("enableTests");
  if (enableTests.empty()) {
    fEnabledTests = {"forEachAssociatedGroup", "associatedGroups", "associatedGroupsWithLeft"};
  }
  else {
    std::copy(
      enableTests.begin(), enableTests.end(), std::inserter(fEnabledTests, fEnabledTests.begin()));
  }
}

void AssnAnalyzer::analyze(art::Event const& e)
{
  if (fEnabledTests.count("forEachAssociatedGroup"))
    AssnAnalyzer::for_each_associated_group_test(e);
  if (fEnabledTests.count("associatedGroups")) AssnAnalyzer::associated_groups_test(e);
  if (fEnabledTests.count("associatedGroupsWithLeft"))
    AssnAnalyzer::associated_groups_with_left_test(e);
}

void AssnAnalyzer::for_each_associated_group_test(art::Event const& e) const
{
  typedef typename art::Assns<int, std::string> istr_assns;
  auto const& int_to_str_assns = *e.getValidHandle<istr_assns>(fInputLabel);
  auto vs = strvec_t{"one", "one-a", "two", "two-a", "three", "three-a"};

  strvec_t strvec;
  auto strings = [&strvec](auto strs) {
    for (auto s = begin(strs); s != end(strs); ++s) {
      std::cout << *s << std::flush << " \"" << **s << "\"" << std::endl;
      strvec.push_back(**s);
    }
  };

  art::for_each_group(int_to_str_assns, strings);

  //strings should be same as vs
  for (auto k = 0; k < 6; ++k) {
    if (strvec[k] != vs[k]) {
      throw art::Exception(art::errors::LogicError)
        << "String #" << k << " expected to be '" << vs[k] << "', got '" << strvec[k]
        << "' instead!\n";
    }
  }

} // for_each_associated_group_test()

void AssnAnalyzer::associated_groups_test(art::Event const& e) const
{
  // this is the exact same test as for_each_associated_group_test(),
  // but written with an explicit loop

  typedef typename art::Assns<int, std::string> istr_assns;
  auto const& int_to_str_assns = *e.getValidHandle<istr_assns>(fInputLabel);
  auto vs = strvec_t{"one", "one-a", "two", "two-a", "three", "three-a"};

  strvec_t strvec;
  for (auto strs : util::associated_groups(int_to_str_assns)) {
    for (art::Ptr<std::string> const& s : strs) {
      std::cout << s << std::flush << " \"" << *s << "\"" << std::endl;
      strvec.push_back(*s);
    }
  } // for associated groups

  //strings should be same as vs
  for (auto k = 0; k < 6; ++k) {
    if (strvec[k] != vs[k]) {
      throw art::Exception(art::errors::LogicError)
        << "String #" << k << " expected to be '" << vs[k] << "', got '" << strvec[k]
        << "' instead!\n";
    }
  }

} // associated_groups_test()

void AssnAnalyzer::associated_groups_with_left_test(art::Event const& e) const
{
  // this is the exact same test as associated_groups_test(),
  // but passing around also the key

  typedef typename art::Assns<int, std::string> istr_assns;
  auto const& int_to_str_assns = *e.getValidHandle<istr_assns>(fInputLabel);
  std::vector<std::pair<int, std::string>> vs = {
    {1, "one"}, {1, "one-a"}, {2, "two"}, {2, "two-a"}, {3, "three"}, {3, "three-a"}};

  std::vector<std::pair<int, std::string>> strvec;
  for (auto const& group : util::associated_groups_with_left(int_to_str_assns)) {
    // user code here:
    auto const& key = std::get<0>(group);  // group.first also ok
    auto const& strs = std::get<1>(group); // group.second also ok

    std::cout << "#" << (*key) << " (" << key << ")" << std::endl;
    for (art::Ptr<std::string> const& s : strs) {
      std::cout << " - " << s << " \"" << *s << "\"" << std::endl;
      strvec.emplace_back(*key, *s);
    }
  } // for associated groups

  //strings should be same as vs
  for (auto k = 0; k < 6; ++k) {
    std::string const& s = strvec[k].second;
    int key = 0; // (unknown)
    if (starts_with(s, "one"))
      key = 1;
    else if (starts_with(s, "two"))
      key = 2;
    else if (starts_with(s, "three"))
      key = 3;

    if (key != vs[k].first) {
      throw art::Exception(art::errors::LogicError)
        << "String #" << k << " expected to have key '" << vs[k].first << "', got '" << key
        << "' instead!\n";
    }

    if (s != vs[k].second) {
      throw art::Exception(art::errors::LogicError)
        << "String #" << k << " expected to be '" << vs[k].second << "', got '" << s
        << "' instead!\n";
    }
  }

} // associated_groups_test()

DEFINE_ART_MODULE(AssnAnalyzer)
