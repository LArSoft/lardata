/**
 * @file   Dereference_test.cc
 * @brief  Test for Dereference.h utilities.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 23rd, 2015
 *
 * A static test is performed to verify that the compilation succeeds with the
 * types int, int*, unique_ptr<int> and a custom smart pointer to int.
 */

// static tests
#include <type_traits>
#include <memory>

// Boost libraries
/*
 * Boost Magic: define the name of the module;
 * and do that before the inclusion of Boost unit test headers
 * because it will change what they provide.
 * Among the those, there is a main() function and some wrapping catching
 * unhandled exceptions and considering them test failures, and probably more.
 * This also makes fairly complicate to receive parameters from the command line
 * (for example, a random seed).
 */
#define BOOST_TEST_MODULE ( StatCollector_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()
#include <boost/test/tools/floating_point_comparison.hpp> // BOOST_CHECK_CLOSE()

// library to be tested:
#include "lardata/Utilities/Dereference.h"


//------------------------------------------------------------------------------
// custom "smart pointer" type
template <typename T>
struct MyPtr {
	using element_type = T;
	MyPtr(T* p = nullptr): ptr(p) {}

	T& operator* () { return *ptr; }

	T* ptr;
}; // struct MyPtr<>


//******************************************************************************
//*** static (compilation-time) tests
//***
// test details::is_type<>
static_assert(lar::util::details::is_type<decltype(*(std::unique_ptr<int>()))>::value,
  "is_type<*(unique_ptr<int>())> not working");


// test details::has_dereference_class
static_assert(!lar::util::details::has_dereference_class<int>::value,
  "details::has_dereference_class<int> is not working");
static_assert(lar::util::details::has_dereference_class<int*>::value,
  "details::has_dereference_class<int*> is not working");
static_assert(lar::util::details::has_dereference_class<MyPtr<int>>::value,
  "details::has_dereference_class<MyPtr<int>> is not working");
static_assert(lar::util::details::has_dereference_class<std::unique_ptr<int>>::value,
  "details::has_dereference_class<unique_ptr<int>> is not working");


// test details::dereferenced_type; should always have type = int&
static_assert(std::is_same<
    typename lar::util::details::dereferenced_type<int, false>::type,
    int
  >::value,
  "details::dereferenced_type<int> is not working"
  );
static_assert(std::is_same<
    typename lar::util::details::dereferenced_type<int*, true>::type,
    int&
  >::value,
  "details::dereferenced_type<int*> is not working"
  );
static_assert(std::is_same<
    typename lar::util::details::dereferenced_type<MyPtr<int>, true>::type,
    int&
  >::value,
  "details::dereferenced_type<MyPtr<int>> is not working"
  );
static_assert(std::is_same<
    typename lar::util::details::dereferenced_type<std::unique_ptr<int>, true>::type,
    int&
  >::value,
  "details::dereferenced_type<unique_ptr<int>> is not working"
  );


// test dereference_class pointer type; should always be int*
static_assert(std::is_same<
    lar::util::details::dereference_class<int, false>::reference_type,
    int&
  >::value,
  "details::dereference_class<int> not working"
  );
static_assert(std::is_same<
    lar::util::details::dereference_class<int*, true>::reference_type,
    int&
  >::value,
  "details::dereference_class<int*> not working"
  );
static_assert(std::is_same<
    lar::util::details::dereference_class<MyPtr<int>, true>::reference_type,
    int&
  >::value,
  "details::dereference_class<MyPtr<int>> not working"
  );
static_assert(std::is_same<
    lar::util::details::dereference_class<std::unique_ptr<int>, true>::reference_type,
    int&
  >::value,
  "details::dereference_class<unique_ptr<int>> not working"
  );


// test dereference_class pointer type; should always be int*
static_assert(std::is_same<
    lar::util::details::make_pointer_class<int, false>::pointer_type,
    int*
  >::value,
  "details::make_pointer_class<int> not working"
  );
static_assert(std::is_same<
    lar::util::details::make_pointer_class<int*, true>::pointer_type,
    int*
  >::value,
  "details::make_pointer_class<int*> not working"
  );
static_assert(std::is_same<
    lar::util::details::make_pointer_class<MyPtr<int>, true>::pointer_type,
    int*
  >::value,
  "details::make_pointer_class<MyPtr<int>> not working"
  );
static_assert(std::is_same<
    lar::util::details::make_pointer_class<std::unique_ptr<int>, true>::pointer_type,
    int*
  >::value,
  "details::make_pointer_class<unique_ptr<int>> not working"
  );


// test dereferenced_type type; should always be int
static_assert(std::is_same<
    typename lar::util::dereferenced_type<int>::type,
    int
  >::value,
  "dereferenced_type<int> not working"
  );
static_assert(std::is_same<
    typename lar::util::dereferenced_type<int*>::type,
    int&
  >::value,
  "dereferenced_type<int*> not working"
  );
static_assert(std::is_same<
    typename lar::util::dereferenced_type<MyPtr<int>>::type,
    int&
  >::value,
  "dereferenced_type<MyPtr<int>> not working"
  );
static_assert(std::is_same<
    typename lar::util::dereferenced_type<std::unique_ptr<int>>::type,
    int&
  >::value,
  "dereferenced_type<unique_ptr<int>> not working"
  );


//******************************************************************************
//***  testing starts here;
//***  still mostly a compilation test
template <typename T>
void test() {

  T value = T(17);
  T* cptr = &value;
  MyPtr<T> my_ptr(&value);
  auto uptr = std::make_unique<T>(value);

  T* ptr;
  ptr = lar::util::make_pointer(uptr);
  BOOST_CHECK_EQUAL(*ptr, value);
  ptr = lar::util::make_pointer(my_ptr);
  BOOST_CHECK_EQUAL(*ptr, value);
  ptr = lar::util::make_pointer(cptr);
  BOOST_CHECK_EQUAL(*ptr, value);
  ptr = lar::util::make_pointer(value);
  BOOST_CHECK_EQUAL(*ptr, value);

  BOOST_CHECK_EQUAL(lar::util::dereference(uptr), value);
  BOOST_CHECK_EQUAL(lar::util::dereference(my_ptr), value);
  BOOST_CHECK_EQUAL(lar::util::dereference(cptr), value);
  BOOST_CHECK_EQUAL(lar::util::dereference(value), value);

} // test<>()


//******************************************************************************
BOOST_AUTO_TEST_CASE(TestInt) {
  test<int>();
}

BOOST_AUTO_TEST_CASE(TestConstInt) {
  test<const int>();
}
