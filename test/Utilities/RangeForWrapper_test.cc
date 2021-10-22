/**
 * @file   RangeForWrapper_test.cc
 * @brief  Test for RangeForWrapper utilities
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 29, 2016
 * @see    RangeForWrapper.h
 *
 * The test is run with no arguments.
 *
 */

// LArSoft libraries
#include "lardata/Utilities/RangeForWrapper.h"

// Boost libraries
#define BOOST_TEST_MODULE ( RangeForWrapper_test )
#include <boost/test/unit_test.hpp>

// C/C++ standard libraries
#include <iterator>
#include <numeric> // std::accumulate()
#include <limits> // std::numeric_limits<>
#include <type_traits> // std::is_same<>, std::is_lvalue_reference<>, ...
#include <initializer_list>


//------------------------------------------------------------------------------
//--- make up a class with begin and end iterators with different types
//---

template <typename Value>
class base_iterator {
  using traits_t = std::iterator_traits<Value*>;
    public:
  using difference_type   = typename traits_t::difference_type;
  using value_type        = typename traits_t::value_type;
  using pointer           = typename traits_t::pointer;
  using reference         = typename traits_t::reference;
  using iterator_category = typename traits_t::iterator_category;

  pointer ptr;

  base_iterator(pointer ptr = nullptr): ptr(ptr) {}

  reference operator* () const { return *ptr; }
  pointer operator-> () const { return ptr; }

  base_iterator& operator++ () { ++ptr; return *this; }
  base_iterator& operator-- () { --ptr; return *this; }

  reference operator[] (difference_type offset) const { return ptr[offset]; }

  difference_type operator- (base_iterator const& other) const
    { return ptr - other.ptr; }

};

template <typename ValueL, typename ValueR>
bool operator!=(base_iterator<ValueL> const& a, base_iterator<ValueR> const& b)
  { return a.ptr != b.ptr; }

struct begin_iterator_t : base_iterator<int> {};
struct end_iterator_t : base_iterator<int> {};

struct begin_const_iterator_t : base_iterator<int const> {};
struct end_const_iterator_t : base_iterator<int const> {};

using wrapper_type = util::details::RangeForWrapperIterator<begin_iterator_t, end_iterator_t>;
using const_wrapper_type = util::details::RangeForWrapperIterator<begin_const_iterator_t, end_const_iterator_t>;
BOOST_TEST_DONT_PRINT_LOG_VALUE(wrapper_type)
BOOST_TEST_DONT_PRINT_LOG_VALUE(const_wrapper_type)

struct Data {
  std::vector<int> data;

  using begin_iterator = begin_iterator_t;
  using end_iterator = end_iterator_t;
  using begin_const_iterator = begin_const_iterator_t;
  using end_const_iterator = end_const_iterator_t;

  Data(std::initializer_list<int> data): data(data) {}

  begin_const_iterator begin() const { return { &*data.cbegin() }; }
  begin_iterator begin() { return { &*data.begin() }; }
  end_const_iterator end() const { return { &*data.cend() }; }
  end_iterator end() { return { &*data.end() }; }

  bool empty() const { return data.empty(); }
  auto size() const { return data.size(); }
  auto operator[](std::size_t index) -> decltype(auto)
    { return data[index]; }
  auto operator[](std::size_t index) const -> decltype(auto)
    { return data[index]; }

}; // class Data

template <typename T>
T copy(T const& v) { return v; }


//------------------------------------------------------------------------------

template <typename DataColl>
void const_test(DataColl const& data, int expected_total)
{
  static_assert(
    !std::is_lvalue_reference
      <decltype(copy(data) | util::range_for)>::value,
    "util::range_for on a rvalue should return a rvalue"
    );

  int total = 0;
  for (int d: data | util::range_for) total += d;

  BOOST_TEST(total == expected_total);

  // from a temporary
  total = 0;
  for (int d: copy(data) | util::range_for) total += d;

  BOOST_TEST(total == expected_total);

} // const_test()


template <typename DataColl>
void test(DataColl& data, int expected_total) {
  static_assert(
    !std::is_lvalue_reference
      <decltype(copy(data) | util::range_for)>::value,
    "util::range_for on a rvalue should return a rvalue"
    );

  //
  // from a lvalue
  //
  int total = 0;
  for (int d: data | util::range_for) total += d;

  BOOST_TEST(total == expected_total);

  //
  // from a rvalue (temporary)
  //
  total = 0;
  for (int d: copy(data) | util::range_for) total += d;

  BOOST_TEST(total == expected_total);

  //
  // from a temporary, which is changed
  //
  total = 0;
  for (int& d: copy(data) | util::range_for) total += (d *= 3);

  BOOST_TEST(total == 3 * expected_total);

  // original value is still unchanged
  total = 0;
  for (int d: data | util::range_for) total += d;

  BOOST_TEST(total == expected_total);

  //
  // from a lvalue, which is changed
  //
  for (int& d: data | util::range_for) d *= 3;

  total = 0;
  for (int d: data | util::range_for) total += d;

  BOOST_TEST(total == 3 * expected_total);

} // test()


//-----------------------------------------------------------------------------
// iterator tests

template <typename Iter, typename RefIter>
void iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  // these asserts do not completely test the concepts;
  // for example, the copy-constructable concept requires that the value of
  // the copy-constructed object and the one it was constructed from are
  // "equivalent", while std::is_copy_constructable only guarantees that
  // copy construction can happen, not testing the equivalence.
  static_assert
    (std::is_move_constructible<Iter>(), "Iterator concept violation");
  static_assert
    (std::is_copy_constructible<Iter>(), "Iterator concept violation");
  static_assert
    (std::is_copy_assignable<Iter>(), "Iterator concept violation");
  static_assert
    (std::is_destructible<Iter>(), "Iterator concept violation");
  // C++17
//  static_assert
//    (std::is_swappable<Iter>(), "Iterator concept violation");

  const bool isEnd = refIter == refEnd;
  const bool isSingular = (iter == Iter());
  const bool isDereferenciable = !isEnd && !isSingular;
  const bool isLast = (std::next(refIter) == refEnd);

  Iter ia = iter;
  BOOST_TEST(bool(ia == iter));

  if (isDereferenciable) {
    //BOOST_TEST(*iter == *refIter);
    BOOST_TEST((*iter == *refIter));
  }
  if (!isLast && !isEnd && !isSingular) {
    //BOOST_TEST(*++iter == *++refIter);
    BOOST_TEST((*++iter == *++refIter));
  }

} // iterator_tests()


template <typename Iter, typename RefIter>
void const_input_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  const bool isEnd = refIter == refEnd;
  const bool isSingular = (iter == Iter());
  const bool isNormal = !isEnd && !isSingular;
//  const bool isDereferenciable = isNormal;
  const bool isLast = (std::next(refIter) == refEnd);

  static_assert(std::is_same<
    typename std::iterator_traits<Iter>::reference,
    decltype(Iter(iter).operator*())
    >(), "Inconsistent return type for dereference operator");

  if (!isEnd) {
    auto ia = iter;
    if (!isSingular) BOOST_TEST(++ia != iter);
  }

  if (isNormal) {
    auto ia = iter, ib = iter;
    ia++; ++ib;
    BOOST_TEST(ia == ib);
    if (!isLast) {
      auto ia = iter;
      BOOST_TEST(*(ia++) == *iter);
    }
  }

  {
    auto ia = iter;
    if (!isEnd)
      (void)ia++;
  }


} // const_input_iterator_tests()

template <typename Iter, typename RefIter>
void input_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  const_input_iterator_tests(iter, refIter, refEnd);
} // input_iterator_tests()


template <typename Iter, typename RefIter>
void const_output_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  const bool isEnd = refIter == refEnd;
  const bool isSingular = (iter == Iter());
  const bool isNormal = !isEnd && !isSingular;
//  const bool isDereferenciable = isNormal;
//  const bool isLast = (std::next(refIter) == refEnd);

  static_assert(std::is_same<
    typename std::iterator_traits<Iter>::reference,
    decltype(Iter(iter).operator*())
    >(), "Inconsistent return type for dereference operator");

  if (!isEnd) {
    auto ia = iter;
    auto addr = &ia;
    if (!isSingular) {
      BOOST_TEST((++ia != iter));
      BOOST_TEST((&++ia == addr));
    }
  }

  if (isNormal) {
    auto ia = iter, ib = iter;
    ia++; ++ib;
    BOOST_TEST(ia == ib);
  }

  {
    auto ia = iter;
    if (!isEnd) (void)ia++;
  }


} // const_output_iterator_tests()

template <typename Iter, typename RefIter>
void output_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  const bool isEnd = refIter == refEnd;

  const_output_iterator_tests(iter, refIter, refEnd);

  if (!isEnd) {
    auto value = *iter;

    auto newValue = std::numeric_limits<typename Iter::value_type>::max();

    *iter = newValue;

    // we couldn't check the result since an output operator is not expected
    // to return anything useful with *iter; but since we have already used it
    // to store the old value...
    BOOST_TEST((*iter == newValue));

    *iter = value;
    BOOST_TEST((*iter == value));

    auto ia = iter, ib = iter;

    *ia++ = newValue;
    ++ib;
    BOOST_TEST((*iter == newValue));
    BOOST_TEST(ia == ib);

    *iter = value;
    BOOST_TEST((*iter == value));
  }

} // output_iterator_tests()


template <typename Iter, typename RefIter>
void const_forward_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  const bool isEnd = refIter == refEnd;
  const bool isSingular = (iter == Iter());
  const bool isNormal = !isEnd && !isSingular;
  const bool isDereferenciable = isNormal;

  static_assert(std::is_default_constructible<Iter>(),
    "Forward iterator concept violation");

  using traits_t = std::iterator_traits<Iter>;
  using dereference_t = std::remove_reference_t<decltype(*iter)>;
  constexpr bool isConst = std::is_const<dereference_t>();

  static_assert(
    std::is_same<
      typename traits_t::reference,
      std::add_lvalue_reference_t<
        std::conditional_t<isConst,
          std::add_const_t<typename traits_t::value_type>,
          typename traits_t::value_type
        >
      >
    >(),
    "Forward iterator concept violation"
    );

  if (isDereferenciable) {
    auto ia = iter, ib = iter;

    BOOST_TEST((ia == ib));

    BOOST_TEST((&*ia == &*ib));

    BOOST_TEST((++ia == ++ib));

  }

} // const_forward_iterator_tests()

template <typename Iter, typename RefIter>
void forward_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  const_forward_iterator_tests(iter, refIter, refEnd);
} // forward_iterator_tests()


template <typename Iter, typename RefIter>
void const_bidirectional_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  const bool isEnd = refIter == refEnd;
//  const bool isLast = std::next(refIter) == refEnd;
  const bool isSingular = (iter == Iter());
  const bool isNormal = !isEnd && !isSingular;
  const bool isDereferenciable = isNormal;

  auto ia = isEnd? iter: std::next(iter);
  auto iaRef = isEnd? refIter: std::next(refIter);

  auto ib = ia;
  --ib;
  BOOST_TEST((ib != ia));
  BOOST_TEST((*ib == *std::prev(iaRef)));
  BOOST_TEST((++ib == ia));
  if (isDereferenciable) BOOST_TEST((*ib == *iaRef));

  ib--;
  BOOST_TEST((ib != ia));
  BOOST_TEST((*ib == *std::prev(iaRef)));
  BOOST_TEST((++ib == ia));
  if (isDereferenciable) BOOST_TEST((*ib == *iaRef));

} // const_bidirectional_iterator_tests()

template <typename Iter, typename RefIter>
void bidirectional_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  const_bidirectional_iterator_tests(iter, refIter, refEnd);
} // bidirectional_iterator_tests()


template <typename Iter, typename RefIter>
void const_random_access_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
} // const_random_access_iterator_tests()

template <typename Iter, typename RefIter>
void random_access_iterator_tests
  (Iter iter, RefIter refIter, RefIter refEnd)
{
  const_random_access_iterator_tests(iter, refIter, refEnd);
} // random_access_iterator_tests()


//
// dispatching
//
template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::input_iterator_tag)
{
  iterator_tests(iter, refIter, refEnd);
  const_input_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<const>(input iterator)


template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<!IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::input_iterator_tag tag)
{
  iterator_tests(iter, refIter, refEnd);
  input_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<non-const>(input iterator)


template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::output_iterator_tag)
{
  iterator_tests(iter, refIter, refEnd);
  const_output_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<const>(output iterator)


template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<!IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::output_iterator_tag tag)
{
  iterator_tests(iter, refIter, refEnd);
  output_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<non-const>(output iterator)


template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::forward_iterator_tag)
{
  iterator_test_impl<IsConst>(iter, refIter, refEnd, std::input_iterator_tag{});
  iterator_test_impl<IsConst>(iter, refIter, refEnd, std::output_iterator_tag{});
  const_forward_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<const>(random access iterator)


template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<!IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::forward_iterator_tag tag)
{
  iterator_test_impl<IsConst>(iter, refIter, refEnd, std::input_iterator_tag{});
  iterator_test_impl<IsConst>(iter, refIter, refEnd, std::output_iterator_tag{});
  forward_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<non-const>(forward iterator)


template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::bidirectional_iterator_tag)
{
  iterator_test_impl<IsConst>(iter, refIter, refEnd, std::forward_iterator_tag{});
  const_random_access_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<const>(forward iterator)


template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<!IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::bidirectional_iterator_tag tag)
{
  iterator_test_impl<IsConst>(iter, refIter, refEnd, std::forward_iterator_tag{});
  bidirectional_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<non-const>(bidirectional iterator)


template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::random_access_iterator_tag)
{
  iterator_test_impl<IsConst>(iter, refIter, refEnd, std::bidirectional_iterator_tag{});
  const_bidirectional_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<const>(bidirectional iterator)


template <bool IsConst, typename Iter, typename RefIter>
std::enable_if_t<!IsConst> iterator_test_impl
  (Iter iter, RefIter refIter, RefIter refEnd, std::random_access_iterator_tag tag)
{
  iterator_test_impl<IsConst>(iter, refIter, refEnd, std::bidirectional_iterator_tag{});
  random_access_iterator_tests(iter, refIter, refEnd);
} // iterator_test_impl<non-const>(random access iterator)


template <typename Iter, typename RefIter>
void iterator_test(Iter iter, RefIter refIter, RefIter refEnd) {
  //
  // dispatcher
  //
  using traits_t = std::iterator_traits<Iter>;
  constexpr bool IsConst
    = std::is_const<std::remove_reference_t<decltype(*iter)>>();

  iterator_test_impl<IsConst>
    (iter, refIter, refEnd, typename traits_t::iterator_category{});

} // iterator_test()


template <bool IsConst>
void RangeForWrapperIteratorStandardsTest() {

  using base_reference_container_t = std::vector<int>;
  using reference_container_t = std::conditional_t
    <IsConst, base_reference_container_t const, base_reference_container_t>;
  reference_container_t vdata = { 2, 3, 4 };
  using basic_container_t = Data;

  static_assert(std::is_same<basic_container_t::begin_iterator::reference, int&>());
  static_assert(std::is_same<typename std::iterator_traits<basic_container_t::begin_iterator>::reference, int&>());

  static_assert(std::is_same<basic_container_t::begin_const_iterator::reference, std::add_const_t<int>&>());
  static_assert(std::is_same<typename std::iterator_traits<basic_container_t::begin_const_iterator>::reference, std::add_const_t<int>&>());

  using container_t
    = std::conditional_t<IsConst, basic_container_t const, basic_container_t>;
  //
  // non-const iterator interface (iterators may still be constant)
  //
  container_t data = { 2, 3, 4 };
  auto range = data | util::range_for;
  using std::begin;
  using std::end;
  auto rbegin = begin(range);
  auto rend = end(range);

  auto vbegin = vdata.begin();
  auto vend = vdata.end();

  BOOST_TEST((unsigned)std::distance(rbegin, rend) == vdata.size());

  iterator_test(rbegin, vbegin, vend);
  iterator_test(rend, vend, vend);


  //
  // const iterator interface
  //
  using std::cbegin;
  using std::cend;
  auto rcbegin = cbegin(range);
  auto rcend = cend(range);

  auto vcbegin = vdata.cbegin();
  auto vcend = vdata.cend();

  BOOST_TEST((unsigned)std::distance(rcbegin, rcend) == vdata.size());

  iterator_test(rcbegin, vcbegin, vcend);
  iterator_test(rcend, vcend, vcend);

  //
  // extra access (partial support for random access)
  //
  BOOST_TEST(range.size() == data.size());
  BOOST_TEST(range.empty() == data.empty());
  for (std::size_t i = 0; i < data.size(); ++i) {
    BOOST_TEST(range[i] == data[i]);
  }

} // RangeForWrapperIteratorStandardsTest()


//-----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(RangeForWrapperSameIterator_test) {

  std::vector vdata { 2, 3, 4 };

  // static checks
  static_assert(std::is_same<
      std::decay_t<decltype(vdata)>,
      std::decay_t<decltype(vdata | util::range_for)>
    >::value,
    "util::range_for should be pass-through!"
    );
  static_assert(
    std::is_lvalue_reference<decltype(vdata | util::range_for)>::value,
    "Pass-through on a lvalue should return a lvalue reference"
    );

  static_assert(
    !std::is_lvalue_reference
      <decltype(copy(vdata) | util::range_for)>::value,
    "Pass-through on a rvalue should return a rvalue (or its reference)"
    );

  BOOST_TEST(&vdata == &(vdata | util::range_for));

  auto const expected_total = std::accumulate(vdata.begin(), vdata.end(), 0);

  const_test(vdata, expected_total);

  test(vdata, expected_total);

} // BOOST_AUTO_TEST_CASE(RangeForWrapperSameIterator_test_test)


BOOST_AUTO_TEST_CASE(RangeForWrapperDifferentIterator_test) {

  std::vector vdata { 2, 3, 4 };

  Data data { 2, 3, 4 };

  auto expected_total = std::accumulate(vdata.begin(), vdata.end(), 0);

  // static check
  static_assert(!std::is_same<
      std::decay_t<decltype(data)>,
      std::decay_t<decltype(data | util::range_for)>
    >::value,
    "util::range_for should generate a wrapper!"
    );

//  for (double d: data); // this should fail compilation

  const_test(data, expected_total);

  test(data, expected_total);

} // BOOST_AUTO_TEST_CASE(RangeForWrapperDifferentIterator_test_test)


BOOST_AUTO_TEST_CASE(RangeForWrapperIteratorStandardsTestCase) {

  RangeForWrapperIteratorStandardsTest<false>(); // mutable range test
  RangeForWrapperIteratorStandardsTest<true>(); // constant range test

} // BOOST_AUTO_TEST_CASE(RangeForWrapperIteratorStandardsTestCase)
