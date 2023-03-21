/**
 * @file   associated_groups_test.cc
 * @brief  Unit test for `util::associated_groups()`.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   February 7th, 2017
 *
 */

// LArSoft libraries
#include "lardata/Utilities/ForEachAssociatedGroup.h"

// framework libraries
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/ProductID.h"

// Boost libraries
#define BOOST_TEST_MODULE (PointIsolationAlg_test)
#include "boost/test/unit_test.hpp"

// C/C++ standard libraries
#include <array>
#include <iostream>

//------------------------------------------------------------------------------
// ROOT libraries
#include "TClassEdit.h"
#include "TInterpreter.h"
#include "TROOT.h" // gROOT

// C/C++ standard libraries
#include <string>
#include <typeinfo>

template <typename T>
TClass* QuickGenerateTClass()
{

  // magic! this interpreter call is needed before GetNormalizedName() is called
  TInterpreter* interpreter = gROOT->GetInterpreter();

  // demangle name of type T
  int err; // we'll ignore errors
  char* classNameC = TClassEdit::DemangleTypeIdName(typeid(T), err);

  // "normalise" it
  std::string normalizedClassName;
  TClassEdit::GetNormalizedName(normalizedClassName, classNameC);

  // clean up
  std::free(classNameC);

  // generate and register the TClass; load it and be silent.
  return interpreter->GenerateTClass(normalizedClassName.c_str(), kTRUE, kTRUE);

} // QuickGenerateTClass()

//------------------------------------------------------------------------------
void AssociatedGroupsTest()
{

  // types used in the association (they actually do not matter)
  struct TypeA {};
  struct TypeB {};

  using MyAssns_t = art::Assns<TypeA, TypeB>;

  // art::Assns constructor tries to have ROOT initialise its streamer, which
  // requires a TClass instance which is not present at this time.
  // This trick creates that TClass.
  QuickGenerateTClass<MyAssns_t>();

  using Index_t = art::Ptr<TypeA>::key_type;

  // association description: B's for each A
  std::array<std::pair<Index_t, std::vector<Index_t>>, 3U> expected;
  expected[0] = {0, {0, 3, 6}};
  expected[1] = {1, {2, 4, 6}};
  expected[2] = {3, {8, 10, 12, 13}};
  art::ProductID aPID{5}, bPID{12};

  // fill the association
  MyAssns_t assns;
  for (auto const& pair : expected) {
    auto const& aIndex = pair.first;
    for (auto const& bIndex : pair.second) {
      assns.addSingle({aPID, aIndex, nullptr}, {bPID, bIndex, nullptr});
    } // for bIndex
  }   // for pair

  std::vector<std::vector<Index_t>> results;
  for (auto Bs : util::associated_groups(assns)) {
    std::cout << "Association group #" << results.size() << ":" << std::endl;
    results.emplace_back();
    std::vector<Index_t>& myBs = results.back();
    for (art::Ptr<TypeB> const& B : Bs) {
      std::cout << "  " << B << std::endl;
      myBs.push_back(B.key());
    }
  } // for associated groups

  //strings should be same as vs
  std::cout << "Starting check..." << std::endl;
  BOOST_TEST(results.size() == expected.size());
  for (size_t i = 0; i < results.size(); ++i) {
    auto const& Bs = results[i];
    auto const& expectedBs = expected[i].second;
    BOOST_TEST_MESSAGE("  element #" << i << ", A=" << expected[i].first);
    BOOST_TEST(Bs.size() == expectedBs.size());
    for (size_t j = 0; j < expectedBs.size(); ++j) {
      BOOST_TEST_MESSAGE("    assn #" << j);
      BOOST_TEST(Bs[j] == expectedBs[j]);
    } // for j
  }   // for results

} // associated_groups_test()

//------------------------------------------------------------------------------
//--- tests
//
BOOST_AUTO_TEST_CASE(AssociatedGroupsTestCase)
{
  AssociatedGroupsTest();
} // AssociatedGroupsTestCase
