/**
 * @file   lardata/test/Utilities/CollectionView_test.cc
 * @brief  Unit test for `CollectionView` class.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   August 3rd, 2017
 *
 * This is a Boost unit test with no specific configuration.
 */

// LArSoft libraries
#include "lardata/Utilities/CollectionView.h"

// Boost libraries
#define BOOST_TEST_MODULE ( CollectionView_test )
#include <boost/test/unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK()

// C/C++ standard libraries
#include <deque>
#include <list>
#include <forward_list>
#include <numeric> // std::iota()
#include <sstream>
#include <iostream>


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(VectorTestCase) {
  //
  // test on contiguous access collection
  //
  std::vector<int> c{ 3, 4, 5 };
  auto const cbegin = c.cbegin();
  auto const cend = c.cend();

  auto cv = lar::makeCollectionView(cbegin, cend);

  BOOST_CHECK_EQUAL(cv.empty(), c.empty());
  BOOST_CHECK_EQUAL(cv.size(), c.size());

  //
  // iterators
  //
  BOOST_CHECK((cv.cbegin() == c.cbegin()));
  BOOST_CHECK((cv.cend() == c.cend()));
  BOOST_CHECK((cv.crbegin() == c.crbegin()));
  BOOST_CHECK((cv.crend() == c.crend()));

  //
  // elements
  //
  BOOST_CHECK_EQUAL(&(cv.front()), &(c.front()));
  BOOST_CHECK_EQUAL(cv.front(), c.front());
  BOOST_CHECK_EQUAL(&(cv.back()), &(c.back()));
  BOOST_CHECK_EQUAL(cv.back(), c.back());

  //
  // data
  //
  BOOST_CHECK_EQUAL(cv.data(), c.data());

  //
  // range-for iteration
  //
  auto ic = cbegin;
  size_t i = 0;
  for (auto const& d: cv) {
    BOOST_CHECK_EQUAL(d, *ic);
    BOOST_CHECK_EQUAL(cv[i], d);
    BOOST_CHECK_EQUAL(&(cv[i]), &(c[i]));
    BOOST_CHECK_EQUAL(cv.at(i), d);
    BOOST_CHECK_EQUAL(&(cv.at(i)), &(c.at(i)));

    ++i;
    ++ic;
  } // for
  BOOST_CHECK((ic == cend));

} // BOOST_AUTO_TEST_CASE(VectorTestCase)


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(DequeTestCase) {
  //
  // test on random access collection
  //
  std::deque<int> c{ 3, 4, 5 };
  auto const cbegin = c.cbegin();
  auto const cend = c.cend();

  auto cv = lar::makeCollectionView(cbegin, cend);

  BOOST_CHECK_EQUAL(cv.empty(), c.empty());
  BOOST_CHECK_EQUAL(cv.size(), c.size());

  //
  // iterators
  //
  BOOST_CHECK((cv.cbegin() == c.cbegin()));
  BOOST_CHECK((cv.cend() == c.cend()));
  BOOST_CHECK((cv.crbegin() == c.crbegin()));
  BOOST_CHECK((cv.crend() == c.crend()));

  //
  // elements
  //
  BOOST_CHECK_EQUAL(&(cv.front()), &(c.front()));
  BOOST_CHECK_EQUAL(cv.front(), c.front());
  BOOST_CHECK_EQUAL(&(cv.back()), &(c.back()));
  BOOST_CHECK_EQUAL(cv.back(), c.back());

  //
  // range-for iteration
  //
  auto ic = cbegin;
  size_t i = 0;
  for (auto const& d: cv) {
    BOOST_CHECK_EQUAL(d, *ic);
    BOOST_CHECK_EQUAL(cv[i], d);
    BOOST_CHECK_EQUAL(&(cv[i]), &(c[i]));
    BOOST_CHECK_EQUAL(cv.at(i), d);
    BOOST_CHECK_EQUAL(&(cv.at(i)), &(c.at(i)));

    ++i;
    ++ic;
  } // for
  BOOST_CHECK((ic == cend));

} // BOOST_AUTO_TEST_CASE(DequeTestCase)


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(ListTestCase) {
  //
  // test on bidirectional access collection
  //
  std::list<int> c{ 3, 4, 5 };
  auto const cbegin = c.cbegin();
  auto const cend = c.cend();

  auto cv = lar::makeCollectionView(cbegin, cend);

  BOOST_CHECK_EQUAL(cv.empty(), c.empty());
  BOOST_CHECK_EQUAL(cv.size(), c.size());

  //
  // iterators
  //
  BOOST_CHECK((cv.cbegin() == c.cbegin()));
  BOOST_CHECK((cv.cend() == c.cend()));
  BOOST_CHECK((cv.crbegin() == c.crbegin()));
  BOOST_CHECK((cv.crend() == c.crend()));

  //
  // elements
  //
  BOOST_CHECK_EQUAL(&(cv.front()), &(c.front()));
  BOOST_CHECK_EQUAL(cv.front(), c.front());
  BOOST_CHECK_EQUAL(&(cv.back()), &(c.back()));
  BOOST_CHECK_EQUAL(cv.back(), c.back());

  //
  // range-for iteration
  //
  auto ic = cbegin;
  for (auto const& d: cv) {
    BOOST_CHECK_EQUAL(d, *ic);

    ++ic;
  } // for
  BOOST_CHECK((ic == cend));

} // BOOST_AUTO_TEST_CASE(ListTestCase)


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(ForwardListTestCase) {
  //
  // test on forward access collection
  //
  std::forward_list<int> c{ 3, 4, 5 };
  auto const cbegin = c.cbegin();
  auto const cend = c.cend();

  auto cv = lar::makeCollectionView(cbegin, cend);

  BOOST_CHECK_EQUAL(cv.empty(), c.empty());

  //
  // iterators
  //
  BOOST_CHECK((cv.cbegin() == c.cbegin()));
  BOOST_CHECK((cv.cend() == c.cend()));

  //
  // elements
  //
  BOOST_CHECK_EQUAL(&(cv.front()), &(c.front()));
  BOOST_CHECK_EQUAL(cv.front(), c.front());

  //
  // range-for iteration
  //
  auto ic = cbegin;
  for (auto const& d: cv) {
    BOOST_CHECK_EQUAL(d, *ic);

    ++ic;
  } // for
  BOOST_CHECK((ic == cend));

} // BOOST_AUTO_TEST_CASE(ForwardListTestCase)


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(DocumentationTestCase) {

  std::ostringstream out;

  /* The promises:
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * std::vector<int> range(5);
   * std::iota(range.begin(), range.end(), 1); // { 1, 2, 3, 4, 5 }
   *
   * for (int d: lar::wrapCollectionIntoView(range)) {
   *   std::cout << d << " ";
   * }
   * std::cout << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * which will print "1 2 3 4 5 ".
   */
  out.str("");

  std::vector<int> range(5);
  std::iota(range.begin(), range.end(), 1); // { 1, 2, 3, 4, 5 }

  for (int d: lar::wrapCollectionIntoView(range)) {
    out << d << " ";
  }
  std::cout << out.str() << std::endl;

  BOOST_CHECK_EQUAL(out.str(), "1 2 3 4 5 ");

  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * decltype(auto) view = lar::wrapCollectionIntoView(range);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  {
    out.str("");
    decltype(auto) view = lar::wrapCollectionIntoView(range);

    for (int d: view) {
      out << d << " ";
    }

    BOOST_CHECK_EQUAL(out.str(), "1 2 3 4 5 ");
  }

  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const& view = lar::wrapCollectionIntoView(range);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  {
    out.str("");
    auto const& view = lar::wrapCollectionIntoView(range);

    for (int d: view) {
      out << d << " ";
    }

    BOOST_CHECK_EQUAL(out.str(), "1 2 3 4 5 ");
  }

  /* The promise:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> v(10);
   * std::iota(v.begin(), v.end(), 0); // { 0, 1, ..., 9 }
   *
   * for (int d: lar::makeCollectionView(v.begin() + 4, v.begin() + 7)) {
   *   std::cout << d << " ";
   * }
   * std::cout << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will print "4 5 6 ".
   */
  out.str("");
  std::vector<int> v(10);
  std::iota(v.begin(), v.end(), 0); // { 0, 1, ..., 9 }

  for (int d: lar::makeCollectionView(v.begin() + 4, v.begin() + 7)) {
    out << d << " ";
  }
  std::cout << out.str() << std::endl;

  BOOST_CHECK_EQUAL(out.str(), "4 5 6 ");

  /* The promise:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * class IntVector {
   *    using vector_t = std::vector<int>;
   *
   *    vector_t data;
   *
   *      public:
   *    IntVector(vector_t&& data): data(std::move(data)) {}
   *
   *    auto begin() const -> decltype(auto) { return data.cbegin(); }
   *    auto end() const -> decltype(auto) { return data.cend(); }
   *
   * }; // struct IntVector
   *
   * using IntViewBase_t = lar::CollectionView<IntVector>;
   *
   * struct MyCollection: public IntViewBase_t {
   *   MyCollection(std::vector<int>&& data) : IntViewBase_t(std::move(data)) {}
   * }; // class MyCollection
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  out.str("");
  {
    std::vector<int> v_data(10);
    std::iota(v_data.begin(), v_data.end(), 0); // { 0, 1, ..., 9 }

    class IntVector {
       using vector_t = std::vector<int>;

       vector_t data;

         public:
       IntVector(vector_t&& data): data(std::move(data)) {}

       auto begin() const -> decltype(auto) { return data.cbegin(); }
       auto end() const -> decltype(auto) { return data.cend(); }

    }; // struct IntVector

    using IntViewBase_t = lar::CollectionView<IntVector>;

    struct MyCollection: public IntViewBase_t {
      MyCollection(std::vector<int>&& data) : IntViewBase_t(std::move(data)) {}
    }; // class MyCollection

    MyCollection v(std::move(v_data));

    for (int d: v) {
      out << d << " ";
    }
    std::cout << out.str() << std::endl;

    BOOST_CHECK_EQUAL(out.str(), "0 1 2 3 4 5 6 7 8 9 ");
  }


} // BOOST_AUTO_TEST_CASE(DocumentationTestCase)
