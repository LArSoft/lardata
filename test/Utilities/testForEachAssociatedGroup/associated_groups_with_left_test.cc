/**
 * @file   associated_groups_with_left_test.cc
 * @brief  Unit test for `util::associated_groups_with_left()`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 14, 2021
 *
 */

// LArSoft libraries
#include "lardata/Utilities/ForEachAssociatedGroup.h"
#include "larcorealg/CoreUtils/enumerate.h"

// framework libraries
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"

// Boost libraries
#define BOOST_TEST_MODULE ( PointIsolationAlg_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// C/C++ standard libraries
#include <array>
#include <iostream>


//------------------------------------------------------------------------------
// ROOT libraries
#include "TROOT.h" // gROOT
#include "TInterpreter.h"
#include "TClassEdit.h"

// C/C++ standard libraries
#include <string>
#include <typeinfo>

template <typename T>
TClass* QuickGenerateTClass() {

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
void associated_groups_with_left_test() {

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
  expected[0] = { 0, { 0, 3, 6 } };
  expected[1] = { 1, { 2, 4, 6 } };
  expected[2] = { 3, { 8, 10, 12, 13 } };
  art::ProductID aPID{ 5 }, bPID{ 12 };

  // fill the association
  MyAssns_t assns;
  for (auto const& pair: expected) {
    auto const& aIndex = pair.first;
    for (auto const& bIndex: pair.second) {
      assns.addSingle({ aPID, aIndex, nullptr }, { bPID, bIndex, nullptr });
    } // for bIndex
  } // for pair

  std::vector<std::pair<Index_t, std::vector<Index_t>>> results;
  for (auto const& [ A, Bs ]: util::associated_groups_with_left(assns)) {
    std::cout << "Association group #" << results.size() << ":" << std::endl;
    std::vector<Index_t> myBs;
    myBs.reserve(Bs.size());
    for(art::Ptr<TypeB> const& B: Bs) {
      std::cout << "  " << B << std::endl;
      myBs.push_back(B.key());
    }
    results.emplace_back(A.key(), std::move(myBs));
  } // for associated groups

  //strings should be same as vs
  std::cout << "Starting check..." << std::endl;
  BOOST_CHECK_EQUAL(results.size(), expected.size());
  for (auto const& [ i, result, expected ]: util::enumerate(results, expected))
  {
    auto const& [ A, Bs ] = result;
    auto const& [ expectedA, expectedBs ] = expected;
    BOOST_TEST_MESSAGE("  element #" << i << ", A=" << expectedA);
    BOOST_CHECK_EQUAL(A, expectedA);
    BOOST_CHECK_EQUAL_COLLECTIONS
      (Bs.cbegin(), Bs.cend(), expectedBs.cbegin(), expectedBs.cend());
  } // for results

} // associated_groups_with_left_test()


//------------------------------------------------------------------------------
//--- tests
//
BOOST_AUTO_TEST_CASE(AssociatedGroupsWithLeftTestCase) {
  associated_groups_with_left_test();
} // AssociatedGroupsWithLeftTestCase
