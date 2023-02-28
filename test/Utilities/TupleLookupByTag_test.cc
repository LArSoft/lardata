/**
 * @file   TupleLookupByTag_test.cc
 * @brief  Unit tests on `TupleLookupByTag.h` utilities.
 * @author Gianluca Petrillo
 * @date   August 17, 2018
 *
 * Most of the tests are static and their failure will trigger compilation
 * errors.
 */

#include "lardata/Utilities/TupleLookupByTag.h"

// LArSoft utilities
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

// C/C++ standard libraries
#include <cassert>
#include <cstddef> // std::size_t
#include <memory>  // std::addressof()
#include <tuple>
#include <type_traits>
#include <utility> // std::move(), std::forward()

//
// * in the following tests, static tests that are expected to fail
//   (with compilation failure) are commented out
// * error messages are not very smart; anyway the compiler will report the
//   line number where the assertion failed
//

//------------------------------------------------------------------------------
//---  Implementation detail tests
//------------------------------------------------------------------------------
//
// test objects preparation
//
template <typename Tag>
struct TestTagged {
  using myTag = Tag;
  int value;
  TestTagged(int value) : value(value) {}
};
struct TestTagA {};
struct TestTagB {};
struct TestTagC {};
using TestTaggedA = TestTagged<TestTagA>;
using TestTaggedB = TestTagged<TestTagB>;
using TestTaggedC = TestTagged<TestTagC>;

using TestTuple_t = std::tuple<int, char, int>;
using TestTaggedTuple_t = std::tuple<TestTaggedA, TestTaggedB, TestTaggedA>;

template <typename Tagged>
struct TestExtractTag {
  using type = typename Tagged::myTag;
};

// static_assert(std::is_same<TestExtractTag<int>::type, TestTagA>(), "Bug :-O");
static_assert(std::is_same<TestExtractTag<TestTaggedA>::type, TestTagA>(), "Bug :-O");
static_assert(std::is_same<TestExtractTag<TestTaggedB>::type, TestTagB>(), "Bug :-O");
static_assert(std::is_same<TestExtractTag<TestTaggedC>::type, TestTagC>(), "Bug :-O");

//
// util::self_type
//
static_assert(std::is_same<util::self_type<int>::type, int>(), "Bug :'(");
static_assert(std::is_same<util::self_type<TestTaggedA>::type, TestTaggedA>(), "Bug :'(");

static_assert(std::is_same<util::self_t<int>, int>(), "Bug :'(");
static_assert(std::is_same<util::self_t<TestTaggedA>, TestTaggedA>(), "Bug :'(");
static_assert(!std::is_same<util::self_t<TestTaggedA>, TestTaggedB>(), "Bug :'(");

//
// util::typelist_element_type
//
static_assert(std::is_same<util::typelist_element_type<0U, int, int&, int>::type, int>(), "Bug :O");
static_assert(std::is_same<util::typelist_element_type<1U, int, int&, int>::type, int&>(),
              "Bug :O");
static_assert(std::is_same<util::typelist_element_type<2U, int, int&, int>::type, int>(), "Bug :O");
// static_assert( std::is_same<util::typelist_element_type<3U, int, int&, int>::type, void>(), "Bug :O");

static_assert(std::is_same<util::typelist_element_t<0U, int, int&, int>, int>(), "Bug :O");
static_assert(std::is_same<util::typelist_element_t<1U, int, int&, int>, int&>(), "Bug :O");
static_assert(std::is_same<util::typelist_element_t<2U, int, int&, int>, int>(), "Bug :O");
// static_assert( std::is_same<util::typelist_element_t<3U, int, int&, int>, void>(), "Bug :O");

//
// util::details::index_of_extracted_type_checked,
// util::details::index_of_extracted_type_impl
//
static_assert(
  util::details::index_of_extracted_type_checked<util::self_type, int, 3, 3, TestTuple_t>() == 3,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked<util::self_type, TestTaggedA, 3, 3, TestTaggedTuple_t>() == 3,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked<TestExtractTag, TestTaggedA, 3, 3, TestTaggedTuple_t>() == 3,
  "Bug :'(");

static_assert(
  util::details::index_of_extracted_type_checked<util::self_type, int, 3, 2, TestTuple_t>() == 2,
  "Bug :'(");
static_assert(
  util::details::index_of_extracted_type_impl<util::self_type, int, 2, int, TestTuple_t>() == 2,
  "Bug :'(");
static_assert(
  util::details::index_of_extracted_type_checked<util::self_type, int, 3, 1, TestTuple_t>() == 2,
  "Bug :'(");
static_assert(
  util::details::index_of_extracted_type_impl<util::self_type, int, 1, char, TestTuple_t>() == 2,
  "Bug :'(");
static_assert(
  util::details::index_of_extracted_type_checked<util::self_type, int, 3, 0, TestTuple_t>() == 0,
  "Bug :'(");
static_assert(
  util::details::index_of_extracted_type_impl<util::self_type, int, 0, int, TestTuple_t>() == 0,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked<util::self_type, TestTaggedA, 3, 2, TestTaggedTuple_t>() == 2,
  "Bug :'(");
static_assert(util::details::index_of_extracted_type_impl<util::self_type,
                                                          TestTaggedA,
                                                          2,
                                                          TestTaggedA,
                                                          TestTaggedTuple_t>() == 2,
              "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked<util::self_type, TestTaggedA, 3, 1, TestTaggedTuple_t>() == 2,
  "Bug :'(");
static_assert(util::details::index_of_extracted_type_impl<util::self_type,
                                                          TestTaggedA,
                                                          1,
                                                          TestTaggedB,
                                                          TestTaggedTuple_t>() == 2,
              "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked<util::self_type, TestTaggedA, 3, 0, TestTaggedTuple_t>() == 0,
  "Bug :'(");
static_assert(util::details::index_of_extracted_type_impl<util::self_type,
                                                          TestTaggedA,
                                                          0,
                                                          TestTaggedA,
                                                          TestTaggedTuple_t>() == 0,
              "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked<TestExtractTag, TestTagA, 3, 2, TestTaggedTuple_t>() == 2,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_impl<TestExtractTag, TestTagA, 2, TestTagA, TestTaggedTuple_t>() == 2,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked<TestExtractTag, TestTagA, 3, 1, TestTaggedTuple_t>() == 2,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_impl<TestExtractTag, TestTagA, 1, TestTagB, TestTaggedTuple_t>() == 2,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked<TestExtractTag, TestTagA, 3, 0, TestTaggedTuple_t>() == 0,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_impl<TestExtractTag, TestTagA, 0, TestTagA, TestTaggedTuple_t>() == 0,
  "Bug :'(");

//
// util::details::index_of_extracted_type_checked_after
//
static_assert(
  util::details::index_of_extracted_type_checked_after<util::self_type, int, 3, 3, TestTuple_t>() ==
    3,
  "Bug :'(");
static_assert(
  util::details::index_of_extracted_type_checked_after<util::self_type, int, 3, 2, TestTuple_t>() ==
    3,
  "Bug :'(");
static_assert(
  util::details::index_of_extracted_type_checked_after<util::self_type, int, 3, 1, TestTuple_t>() ==
    2,
  "Bug :'(");
static_assert(
  util::details::index_of_extracted_type_checked_after<util::self_type, int, 3, 0, TestTuple_t>() ==
    2,
  "Bug :'(");
static_assert(util::details::index_of_extracted_type_checked_after<util::self_type,
                                                                   TestTaggedA,
                                                                   3,
                                                                   3,
                                                                   TestTaggedTuple_t>() == 3,
              "Bug :'(");
static_assert(util::details::index_of_extracted_type_checked_after<util::self_type,
                                                                   TestTaggedA,
                                                                   3,
                                                                   2,
                                                                   TestTaggedTuple_t>() == 3,
              "Bug :'(");
static_assert(util::details::index_of_extracted_type_checked_after<util::self_type,
                                                                   TestTaggedA,
                                                                   3,
                                                                   1,
                                                                   TestTaggedTuple_t>() == 2,
              "Bug :'(");
static_assert(util::details::index_of_extracted_type_checked_after<util::self_type,
                                                                   TestTaggedA,
                                                                   3,
                                                                   0,
                                                                   TestTaggedTuple_t>() == 2,
              "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked_after<TestExtractTag, TestTagA, 3, 3, TestTaggedTuple_t>() ==
    3,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked_after<TestExtractTag, TestTagA, 3, 2, TestTaggedTuple_t>() ==
    3,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked_after<TestExtractTag, TestTagA, 3, 1, TestTaggedTuple_t>() ==
    2,
  "Bug :'(");
static_assert(
  util::details::
      index_of_extracted_type_checked_after<TestExtractTag, TestTagA, 3, 0, TestTaggedTuple_t>() ==
    2,
  "Bug :'(");

//
// util::details::index_of_type_base
//
static_assert(util::details::index_of_type_base<util::self_type, int, TestTuple_t>() == 0,
              "Bug :'(");
static_assert(util::details::index_of_type_base<util::self_type, char, TestTuple_t>() == 1,
              "Bug :'(");
static_assert(util::details::index_of_type_base<util::self_type, void, TestTuple_t>() == 3,
              "Bug :'(");
static_assert(
  util::details::index_of_type_base<util::self_type, TestTaggedA, TestTaggedTuple_t>() == 0,
  "Bug :'(");
static_assert(
  util::details::index_of_type_base<util::self_type, TestTaggedB, TestTaggedTuple_t>() == 1,
  "Bug :'(");
static_assert(
  util::details::index_of_type_base<util::self_type, TestTaggedC, TestTaggedTuple_t>() == 3,
  "Bug :'(");
static_assert(util::details::index_of_type_base<TestExtractTag, TestTagA, TestTaggedTuple_t>() == 0,
              "Bug :'(");
static_assert(util::details::index_of_type_base<TestExtractTag, TestTagB, TestTaggedTuple_t>() == 1,
              "Bug :'(");
static_assert(util::details::index_of_type_base<TestExtractTag, TestTagC, TestTaggedTuple_t>() == 3,
              "Bug :'(");

//
// util::details::index_of_type_helper
//
// static_assert(util::details::index_of_type_helper<util::self_type, int,         TestTuple_t      >::value == 0, "Bug :'(");
static_assert(util::details::index_of_type_helper<util::self_type, char, TestTuple_t>::value == 1,
              "Bug :'(");
// static_assert(util::details::index_of_type_helper<util::self_type, void,        TestTuple_t      >::value == 3, "Bug :'(");
// static_assert(util::details::index_of_type_helper<util::self_type, TestTaggedA, TestTaggedTuple_t>::value == 0, "Bug :'(");
static_assert(
  util::details::index_of_type_helper<util::self_type, TestTaggedB, TestTaggedTuple_t>::value == 1,
  "Bug :'(");
// static_assert(util::details::index_of_type_helper<util::self_type, TestTaggedC, TestTaggedTuple_t>::value == 3, "Bug :'(");
// static_assert(util::details::index_of_type_helper<TestExtractTag,  TestTagA,    TestTaggedTuple_t>::value == 0, "Bug :'(");
static_assert(
  util::details::index_of_type_helper<TestExtractTag, TestTagB, TestTaggedTuple_t>::value == 1,
  "Bug :'(");
// static_assert(util::details::index_of_type_helper<TestExtractTag,  TestTagC,    TestTaggedTuple_t>::value == 3, "Bug :'(");

//
// util::index_of_extracted_type
//
// static_assert(util::index_of_extracted_type<util::self_type, int,         TestTuple_t      >() == 0, "Bug :'(");
static_assert(util::index_of_extracted_type<util::self_type, char, TestTuple_t>() == 1, "Bug :'(");
// static_assert(util::index_of_extracted_type<util::self_type, void,        TestTuple_t      >() == 3, "Bug :'(");
// static_assert(util::index_of_extracted_type<util::self_type, TestTaggedA, TestTaggedTuple_t>() == 0, "Bug :'(");
static_assert(util::index_of_extracted_type<util::self_type, TestTaggedB, TestTaggedTuple_t>() == 1,
              "Bug :'(");
// static_assert(util::index_of_extracted_type<util::self_type, TestTaggedC, TestTaggedTuple_t>() == 3, "Bug :'(");
// static_assert(util::index_of_extracted_type<TestExtractTag,  TestTagA,    TestTaggedTuple_t>() == 0, "Bug :'(");
static_assert(util::index_of_extracted_type<TestExtractTag, TestTagB, TestTaggedTuple_t>() == 1,
              "Bug :'(");
// static_assert(util::index_of_extracted_type<TestExtractTag,  TestTagC,    TestTaggedTuple_t>() == 3, "Bug :'(");

//
// util::index_of_type
//
// static_assert(util::index_of_type<int,         TestTuple_t      >() == 0, "Bug :'(");
static_assert(util::index_of_type<char, TestTuple_t>() == 1, "Bug :'(");
// static_assert(util::index_of_type<void,        TestTuple_t      >() == 3, "Bug :'(");
// static_assert(util::index_of_type<TestTaggedA, TestTaggedTuple_t>() == 0, "Bug :'(");
static_assert(util::index_of_type<TestTaggedB, TestTaggedTuple_t>() == 1, "Bug :'(");
// static_assert(util::index_of_type<TestTaggedC, TestTaggedTuple_t>() == 3, "Bug :'(");

//
// util::has_extracted_type
//
static_assert(util::has_extracted_type<util::self_type, int, TestTuple_t>(), "Bug :'(");
static_assert(util::has_extracted_type<util::self_type, char, TestTuple_t>(), "Bug :'(");
static_assert(!util::has_extracted_type<util::self_type, void, TestTuple_t>(), "Bug :'(");
static_assert(util::has_extracted_type<util::self_type, TestTaggedA, TestTaggedTuple_t>(),
              "Bug :'(");
static_assert(util::has_extracted_type<util::self_type, TestTaggedB, TestTaggedTuple_t>(),
              "Bug :'(");
static_assert(!util::has_extracted_type<util::self_type, TestTaggedC, TestTaggedTuple_t>(),
              "Bug :'(");
static_assert(util::has_extracted_type<TestExtractTag, TestTagA, TestTaggedTuple_t>(), "Bug :'(");
static_assert(util::has_extracted_type<TestExtractTag, TestTagB, TestTaggedTuple_t>(), "Bug :'(");
static_assert(!util::has_extracted_type<TestExtractTag, TestTagC, TestTaggedTuple_t>(), "Bug :'(");

//
// util::has_type
//
static_assert(util::has_type<int, TestTuple_t>(), "Bug :'(");
static_assert(util::has_type<char, TestTuple_t>(), "Bug :'(");
static_assert(!util::has_type<void, TestTuple_t>(), "Bug :'(");
static_assert(util::has_type<TestTaggedA, TestTaggedTuple_t>(), "Bug :'(");
static_assert(util::has_type<TestTaggedB, TestTaggedTuple_t>(), "Bug :'(");
static_assert(!util::has_type<TestTaggedC, TestTaggedTuple_t>(), "Bug :'(");

//
// util::type_is_in
//
static_assert(util::type_is_in<int, int, char, int>(), "Buuug.");
static_assert(util::type_is_in<char, int, char, int>(), "Buuug.");
static_assert(!util::type_is_in<void, int, char, int>(), "Buuug.");

//
// util::details::count_type_in_list
//
static_assert(util::count_type_in_list<int, int, char, int>() == 2, "Buuug.");
static_assert(util::count_type_in_list<char, int, char, int>() == 1, "Buuug.");
static_assert(util::count_type_in_list<void, int, char, int>() == 0, "Buuug.");

//
// util::details::count_type_in_tuple
//
static_assert(util::count_type_in_tuple<int, std::tuple<int, char, int>>() == 2, "Buuug.");
static_assert(util::count_type_in_tuple<char, std::tuple<int, char, int>>() == 1, "Buuug.");
static_assert(util::count_type_in_tuple<void, std::tuple<int, char, int>>() == 0, "Buuug.");

//
// util::has_duplicate_types_impl
//
static_assert(
  util::details::has_duplicate_types_impl<std::tuple<int, char, int>, int, char, void>(),
  "Buuug.");
static_assert(util::details::has_duplicate_types_impl<std::tuple<int, char, int>, int, char>(),
              "Buuug.");
static_assert(util::details::has_duplicate_types_impl<std::tuple<int, char, int>, int, void>(),
              "Buuug.");
static_assert(
  util::details::has_duplicate_types_impl<std::tuple<int, char, int>, int, char, void>(),
  "Buuug.");
static_assert(!util::details::has_duplicate_types_impl<std::tuple<int, char, int>, char, void>(),
              "Buuug.");
static_assert(util::details::has_duplicate_types_impl<std::tuple<int, char, int>, int>(), "Buuug.");
static_assert(!util::details::has_duplicate_types_impl<std::tuple<int, char, int>, char>(),
              "Buuug.");
static_assert(!util::details::has_duplicate_types_impl<std::tuple<int, char, int>, void>(),
              "Buuug.");

//
// util::has_duplicate_types
//
static_assert(util::has_duplicate_types<std::tuple<int, char, int>>(), "Buuug.");
static_assert(!util::has_duplicate_types<std::tuple<int, char, long>>(), "Buuug.");

//
// util::details::count_tags_in_list
//
static_assert(
  util::details::count_tags_in_list<TestExtractTag, TestTaggedB, TestTagA, TestTagB, TestTagA>() ==
    1,
  "Buuug.");
static_assert(
  util::details::count_tags_in_list<TestExtractTag, TestTaggedA, TestTagA, TestTagB, TestTagA>() ==
    2,
  "Buuug.");
static_assert(
  util::details::count_tags_in_list<TestExtractTag, TestTaggedC, TestTagA, TestTagB, TestTagA>() ==
    0,
  "Buuug.");

//
// util::count_extracted_types
//
static_assert(util::count_extracted_types<util::self_type, int, TestTuple_t>() == 2, "Bug :'(");
static_assert(util::count_extracted_types<util::self_type, char, TestTuple_t>() == 1, "Bug :'(");
static_assert(util::count_extracted_types<util::self_type, void, TestTuple_t>() == 0, "Bug :'(");
static_assert(util::count_extracted_types<util::self_type, TestTaggedA, TestTaggedTuple_t>() == 2,
              "Bug :'(");
static_assert(util::count_extracted_types<util::self_type, TestTaggedB, TestTaggedTuple_t>() == 1,
              "Bug :'(");
static_assert(util::count_extracted_types<util::self_type, TestTaggedC, TestTaggedTuple_t>() == 0,
              "Bug :'(");
static_assert(util::count_extracted_types<TestExtractTag, TestTagA, TestTaggedTuple_t>() == 2,
              "Bug :'(");
static_assert(util::count_extracted_types<TestExtractTag, TestTagB, TestTaggedTuple_t>() == 1,
              "Bug :'(");
static_assert(util::count_extracted_types<TestExtractTag, TestTagC, TestTaggedTuple_t>() == 0,
              "Bug :'(");

//
// util::count_types
//
static_assert(util::count_types<int, TestTuple_t>() == 2, "Bug :'(");
static_assert(util::count_types<char, TestTuple_t>() == 1, "Bug :'(");
static_assert(util::count_types<void, TestTuple_t>() == 0, "Bug :'(");
static_assert(util::count_types<TestTaggedA, TestTaggedTuple_t>() == 2, "Bug :'(");
static_assert(util::count_types<TestTaggedB, TestTaggedTuple_t>() == 1, "Bug :'(");
static_assert(util::count_types<TestTaggedC, TestTaggedTuple_t>() == 0, "Bug :'(");

namespace my {

  // A class supporting tuple-type operations.
  template <typename... Data>
  class MyTuple {
    using tuple_t = std::tuple<Data...>;
    tuple_t data;

  public:
    template <std::size_t I>
    using element_type = std::tuple_element_t<I, tuple_t>;

    MyTuple(Data&&... data) : data(std::forward<Data>(data)...) {}
    MyTuple(tuple_t&& data) : data(std::move(data)) {}

    template <std::size_t I>
    auto get() -> decltype(auto)
    {
      return std::get<I>(data);
    }

    template <std::size_t I>
    auto get() const -> decltype(auto)
    {
      return std::get<I>(data);
    }

    template <typename T>
    auto get() -> decltype(auto)
    {
      return std::get<T>(data);
    }

    template <typename T>
    auto get() const -> decltype(auto)
    {
      return std::get<T>(data);
    }

    static constexpr std::size_t tuple_size() { return std::tuple_size<tuple_t>(); }

  }; // MyTuple

  template <typename... Data>
  auto make_my_tuple(Data&&... data)
  {
    return MyTuple<Data...>{std::forward<Data>(data)...};
  }

  template <std::size_t I, typename... T>
  auto get(MyTuple<T...> const& t) -> decltype(auto)
  {
    return t.template get<I>();
  }

  template <typename Target, typename... T>
  auto get(MyTuple<T...> const& t) -> decltype(auto)
  {
    return t.template get<Target>();
  }

} // namespace my

namespace std {
  // specialization of tuple-like operations for the new type

  template <std::size_t I, typename... T>
  class tuple_element<I, my::MyTuple<T...>> {
  public:
    using type = typename my::MyTuple<T...>::template element_type<I>;
  }; // class tuple_element<MyTuple>

  template <typename... T>
  class tuple_size<my::MyTuple<T...>>
    : public std::integral_constant<std::size_t, my::MyTuple<T...>::tuple_size()> {};

} // namespace std

//
// test types using the custom tuple
//
template <typename Tag, typename Payload = void>
struct TaggedType : public TaggedType<Tag> {
  Payload data;

  TaggedType(Payload data) : data(data) {}
}; // TaggedType

template <typename Tag>
struct TaggedType<Tag, void> {
  using tag = Tag;
};

using TagA = util::TagN<0>;
using TagB = util::TagN<1>;
using TagC = util::TagN<2>;

void testMakeTagged()
{

  struct MyData {
    int content = 5;
  };

  struct MyStonedData : public MyData, public lar::UncopiableAndUnmovableClass {};

  MyData lightData;        // moveable, copiable
  MyStonedData heavyStone; // unmoveable, uncopiable

  decltype(auto) lightDataTagged = util::makeTagged<TagA>(lightData);
  decltype(auto) heavyStoneTagged = util::makeTagged<TagA>(heavyStone);
  decltype(auto) lightDataCopyTagged = util::makeTagged<TagA>(MyData(lightData));

  static_assert(std::is_lvalue_reference<decltype(util::makeTagged<TagA>(lightData))>(),
                "makeTagged(moveable lvalue)   does not produce a reference");
  static_assert(std::is_lvalue_reference<decltype(util::makeTagged<TagA>(heavyStone))>(),
                "makeTagged(unmoveable lvalue) does not produce a reference");
  static_assert(!std::is_lvalue_reference<decltype(util::makeTagged<TagA>(MyData(lightData)))>(),
                "makeTagged(rvalue) produces a (lvalue) reference");
  static_assert(!std::is_rvalue_reference<decltype(util::makeTagged<TagA>(MyData(lightData)))>(),
                "makeTagged(rvalue) produces a (rvalue) reference");

  assert(std::addressof(lightDataTagged) == std::addressof(lightData));
  assert(std::addressof(heavyStoneTagged) == std::addressof(heavyStone));
  assert(std::addressof(lightDataCopyTagged) != std::addressof(lightData));

} // testMakeTagged()

int main()
{

  //
  // test data
  //
  auto const testTuple = std::make_tuple<TestTaggedA, TestTaggedB, TestTaggedA>(1, 2, 3);
  assert((std::get<0>(testTuple).value == 1));
  assert((std::get<1>(testTuple).value == 2));
  assert((std::get<2>(testTuple).value == 3));
  using std::get;
  //  static_assert((index_of_extracted_type<TestExtractTag, TestTagA, TestTaggedTuple_t>(testTuple)>::value == 0), "Bug!");
  static_assert(util::index_of_extracted_type<TestExtractTag, TestTagB, TestTaggedTuple_t>() == 1,
                "Bug!");
  assert(get<1U>(testTuple).value == 2);
  assert((get<util::index_of_extracted_type<TestExtractTag, TestTagB, TestTaggedTuple_t>::value>(
            testTuple)
            .value) == 2);
  //  assert((util::getByExtractedType<TestExtractTag, TestTagA>(testTuple).value == 1));
  assert((util::getByExtractedType<TestExtractTag, TestTagB>(testTuple).value == 2));
  //  assert((util::getByExtractedType<TestExtractTag, TestTagC>(testTuple).value == 3));

  using DataA = TaggedType<TagA, int>;
  using DataB = TaggedType<TagB, int>;
  using DataC = TaggedType<TagC, char>;

  auto data = my::make_my_tuple<DataA, DataC, DataB>(64, 'b', 66);
  auto dataWithDupl = my::make_my_tuple<DataA, DataC, DataA>(64, 'b', 66);

  //
  // traditional std::get()
  //
  using std::get;
  static_assert(std::is_same<std::decay_t<decltype(get<0U>(data))>, DataA>(), "Unexpected type 1");
  static_assert(std::is_same<std::decay_t<decltype(get<1U>(data))>, DataC>(), "Unexpected type 2");
  static_assert(std::is_same<std::decay_t<decltype(get<2U>(data))>, DataB>(), "Unexpected type 3");
  static_assert(std::is_same<std::decay_t<decltype(get<0U>(dataWithDupl))>, DataA>(),
                "Unexpected type 1 (dupl)");
  static_assert(std::is_same<std::decay_t<decltype(get<1U>(dataWithDupl))>, DataC>(),
                "Unexpected type 2 (dupl)");
  static_assert(std::is_same<std::decay_t<decltype(get<2U>(dataWithDupl))>, DataA>(),
                "Unexpected type 3 (dupl)");

  //
  // traditional std::get() (by type)
  //
  using std::get;
  static_assert(std::is_same<std::decay_t<decltype(get<DataA>(data))>, DataA>(),
                "Unexpected type 1");
  static_assert(std::is_same<std::decay_t<decltype(get<DataC>(data))>, DataC>(),
                "Unexpected type 2");
  static_assert(std::is_same<std::decay_t<decltype(get<DataB>(data))>, DataB>(),
                "Unexpected type 3");
  //  static_assert(std::is_same<std::decay_t<decltype(get<DataA>(dataWithDupl))>, DataA>(), "Unexpected type 1 (dupl)"); // does not compile: duplicate types!
  static_assert(std::is_same<std::decay_t<decltype(get<DataC>(dataWithDupl))>, DataC>(),
                "Unexpected type 2 (dupl)");
  //  static_assert(std::is_same<std::decay_t<decltype(get<DataA>(dataWithDupl))>, DataA>(), "Unexpected type 3 (dupl)"); // does not compile: duplicate types!

  //
  // traditional std::tuple_element()
  //
  using std::get;
  static_assert(std::is_same<std::tuple_element_t<0U, decltype(data)>, DataA>(),
                "Unexpected type 1");
  static_assert(std::is_same<std::tuple_element_t<1U, decltype(data)>, DataC>(),
                "Unexpected type 2");
  static_assert(std::is_same<std::tuple_element_t<2U, decltype(data)>, DataB>(),
                "Unexpected type 3");
  static_assert(std::is_same<std::tuple_element_t<0U, decltype(dataWithDupl)>, DataA>(),
                "Unexpected type 1 (dupl)");
  static_assert(std::is_same<std::tuple_element_t<1U, decltype(dataWithDupl)>, DataC>(),
                "Unexpected type 2 (dupl)");
  static_assert(std::is_same<std::tuple_element_t<2U, decltype(dataWithDupl)>, DataA>(),
                "Unexpected type 3 (dupl)");

  //
  // traditional std::tuple_size
  //
  static_assert(std::tuple_size<decltype(data)>() == 3U, "Unexpected tuple size");
  static_assert(std::tuple_size<decltype(dataWithDupl)>() == 3U, "Unexpected tuple size (dupl)");

  //
  // util::index_of_type()
  //
  static_assert(util::index_of_type<DataA, decltype(data)>() == 0U, "Unexpected type 1");
  static_assert(util::index_of_type<DataC, decltype(data)>() == 1U, "Unexpected type 2");
  static_assert(util::index_of_type<DataB, decltype(data)>() == 2U, "Unexpected type 3");
  //  static_assert(util::index_of_type<DataA, decltype(dataWithDupl)>() == 0U, "Unexpected type 1 (dupl)");
  static_assert(util::index_of_type<DataC, decltype(dataWithDupl)>() == 1U,
                "Unexpected type 2 (dupl)");
  //  static_assert(util::index_of_type<DataA, decltype(dataWithDupl)>() == 2U, "Unexpected type 3 (dupl)");

  //
  // util::index_of_tag()
  //
  static_assert(util::index_of_tag<TagA, decltype(data)>() == 0U, "Unexpected tagged type 1");
  static_assert(util::index_of_tag<TagC, decltype(data)>() == 1U, "Unexpected tagged type 2");
  static_assert(util::index_of_tag<TagB, decltype(data)>() == 2U, "Unexpected tagged type 3");
  //  static_assert(util::index_of_tag<TagA, decltype(dataWithDupl)>() == 0U, "Unexpected tagged type 1 (dupl)");
  static_assert(util::index_of_tag<TagC, decltype(dataWithDupl)>() == 1U,
                "Unexpected tagged type 2 (dupl)");
  //  static_assert(util::index_of_tag<TagA, decltype(dataWithDupl)>() == 2U, "Unexpected tagged type 3 (dupl)");

  //
  // util::type_with_tag()
  //
  static_assert(std::is_same<util::type_with_tag_t<TagA, decltype(data)>, DataA>(),
                "Unexpected tagged type 1");
  static_assert(std::is_same<util::type_with_tag_t<TagC, decltype(data)>, DataC>(),
                "Unexpected tagged type 2");
  static_assert(std::is_same<util::type_with_tag_t<TagB, decltype(data)>, DataB>(),
                "Unexpected tagged type 3");
  //  static_assert(std::is_same<util::type_with_tag_t<TagA, decltype(dataWithDupl)>, DataA>(), "Unexpected tagged type 1 (dupl)");
  static_assert(std::is_same<util::type_with_tag_t<TagC, decltype(dataWithDupl)>, DataC>(),
                "Unexpected tagged type 2 (dupl)");
  //  static_assert(std::is_same<util::type_with_tag_t<TagA, decltype(dataWithDupl)>, DataA>(), "Unexpected tagged type 3 (dupl)");

  //
  // util::has_type()
  //
  static_assert(util::has_type<DataA, decltype(data)>(), "Unexpected type 1");
  static_assert(util::has_type<DataC, decltype(data)>(), "Unexpected type 2");
  static_assert(util::has_type<DataB, decltype(data)>(), "Unexpected type 3");
  static_assert(util::has_type<DataA, decltype(dataWithDupl)>(), "Unexpected type 1 (dupl)");
  static_assert(util::has_type<DataC, decltype(dataWithDupl)>(), "Unexpected type 2 (dupl)");
  static_assert(!util::has_type<DataB, decltype(dataWithDupl)>(), "Unexpected type 3 (dupl)");

  //
  // util::has_tag()
  //
  static_assert(util::has_tag<TagA, decltype(data)>(), "Unexpected tagged type 1");
  static_assert(util::has_tag<TagC, decltype(data)>(), "Unexpected tagged type 2");
  static_assert(util::has_tag<TagB, decltype(data)>(), "Unexpected tagged type 3");
  static_assert(util::has_tag<TagA, decltype(dataWithDupl)>(), "Unexpected tagged type 1 (dupl)");
  static_assert(util::has_tag<TagC, decltype(dataWithDupl)>(), "Unexpected tagged type 2 (dupl)");
  static_assert(!util::has_tag<TagB, decltype(dataWithDupl)>(), "Unexpected tagged type 3 (dupl)");

  //
  // util::count_types()
  //
  static_assert(util::count_types<DataA, decltype(data)>() == 1, "Unexpected type 1");
  static_assert(util::count_types<DataC, decltype(data)>() == 1, "Unexpected type 2");
  static_assert(util::count_types<DataB, decltype(data)>() == 1, "Unexpected type 3");
  static_assert(util::count_types<DataA, decltype(dataWithDupl)>() == 2,
                "Unexpected type 1 (dupl)");
  static_assert(util::count_types<DataC, decltype(dataWithDupl)>() == 1,
                "Unexpected type 2 (dupl)");
  static_assert(util::count_types<DataB, decltype(dataWithDupl)>() == 0,
                "Unexpected type 3 (dupl)");

  //
  // util::count_tags()
  //
  static_assert(util::count_tags<TagA, decltype(data)>() == 1, "Unexpected type 1");
  static_assert(util::count_tags<TagC, decltype(data)>() == 1, "Unexpected type 2");
  static_assert(util::count_tags<TagB, decltype(data)>() == 1, "Unexpected type 3");
  static_assert(util::count_tags<TagA, decltype(dataWithDupl)>() == 2, "Unexpected type 1 (dupl)");
  static_assert(util::count_tags<TagC, decltype(dataWithDupl)>() == 1, "Unexpected type 2 (dupl)");
  static_assert(util::count_tags<TagB, decltype(dataWithDupl)>() == 0, "Unexpected type 3 (dupl)");

  //
  // util::has_duplicate_types()
  //
  static_assert(!util::has_duplicate_types<decltype(data)>(), "Type has duplicate tags!");
  static_assert(util::has_duplicate_types<decltype(dataWithDupl)>(), "Type has no duplicate tags");

  //
  // util::has_duplicate_tags()
  //
  static_assert(!util::has_duplicate_tags<decltype(data)>(), "Type has duplicate tags!");
  static_assert(util::has_duplicate_tags<decltype(dataWithDupl)>(), "Type has no duplicate tags");

  //
  // util::getByTag()
  //
  static_assert(std::is_same<std::decay_t<decltype(util::getByTag<TagA>(data).data)>, int>(),
                "Unexpected type 1");
  static_assert(std::is_same<std::decay_t<decltype(util::getByTag<TagC>(data).data)>, char>(),
                "Unexpected type 2");
  static_assert(std::is_same<std::decay_t<decltype(util::getByTag<TagB>(data).data)>, int>(),
                "Unexpected type 3");
  //  static_assert(std::is_same<std::decay_t<decltype(util::getByTag<TagA>(dataWithDupl).data)>, int >(), "Unexpected type 1 (dupl)"); // does not compile: duplicate types!
  static_assert(
    std::is_same<std::decay_t<decltype(util::getByTag<TagC>(dataWithDupl).data)>, char>(),
    "Unexpected type 2 (dupl)");
  //  static_assert(std::is_same<std::decay_t<decltype(util::getByTag<TagA>(dataWithDupl).data)>, int >(), "Unexpected type 3 (dupl)"); // does not compile: duplicate types!

  assert((util::getByTag<TagA>(data)).data == 64);
  assert((util::getByTag<TagC>(data)).data == 'b');
  assert((util::getByTag<TagB>(data)).data == 66);
  //  assert((util::getByTag<TagA>(dataWithDupl)).data ==  64); // does not compile: duplicate types!
  assert((util::getByTag<TagC>(dataWithDupl)).data == 'b');
  //  assert((util::getByTag<TagA>(dataWithDupl)).data ==  66); // does not compile: duplicate types!

  //
  // makeTagged()
  //
  testMakeTagged();

  return 0;
} // main()
