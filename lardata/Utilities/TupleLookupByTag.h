/**
 * @file   TupleLookupByTag.h
 * @brief  Utilities to address elements of a tuple-like class by tag.
 * @author Gianluca Petrillo
 * @date   August 17, 2018
 *
 * A "tagged" class is a class containing a type named "tag". The class is
 * tagged by it. A container of objects, tuple-like, can be accessed by
 * specifying the type of the tag. For example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::tuple<TaggedItem1, TaggedItem2, TaggedItem3> data;
 *
 * auto tagBdata = util::getByTag<TagB>(data);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will assign tagBdata with one of the three components of the `data` objects,
 * the one whose type has `tag` equal to `TagB`.
 *
 * The type being processed must be "tuple-like" (in the example above, the type
 * of `data` is in fact a tuple). This means that the type must respond to
 * `std::tuple_element` and `std::tuple_size`.
 *
 * Beside these utilities, equivalent utilities are exposed that allow a
 * different definition of the tag (via a class "returning" the tag of its
 * only template argument), and that operate on the types directly rather than
 * on the tags (equivalent to define the tag as the object itself).
 *
 * All the type trait utilities follow the standard C++ convention: the ones
 * "returning" a type have that type in a `type` member, and the ones returning
 * a value have that value in a `value` member, and their instances can be
 * implicitly converted to that value. Aliases ending with `_t` and `_v` are
 * also provided.
 *
 * This is a header-only library providing mostly type trait information.
 *
 */

#ifndef LARDATA_UTILITIES_TUPLELOOKUPBYTAG_H
#define LARDATA_UTILITIES_TUPLELOOKUPBYTAG_H

// LArSoft libraries
#include "larcorealg/CoreUtils/MetaUtils.h"

// C/C++ standard libraries
#include <cstddef> // std::size_t
#include <tuple>
#include <type_traits>

/**
 * @addtogroup Utilities General utilities
 * @brief General programming utilities.
 *
 * @{
 */
/**
 * @namespace util
 * @brief Namespace for general, not LArSoft-specific utilities.
 */
/**
 * @addtogroup Metaprogramming General utilities for metaprogramming
 * @brief General utilities for use with templates and metaprogramming.
 */
/**
 * @}
 */

namespace util {

  namespace details {

    //--------------------------------------------------------------------------
    // forward declaration of implementation details
    //--------------------------------------------------------------------------
    //---  General utilities
    //--------------------------------------------------------------------------
    template <typename Target, typename... T>
    struct count_type_in_list_impl;

    template <typename Target, typename... T>
    struct type_is_in_impl;

    //--------------------------------------------------------------------------
    //---  Courtesy traits
    //--------------------------------------------------------------------------
    template <typename Tuple>
    struct has_duplicate_types_unwrapper;

    //--------------------------------------------------------------------------
    //---  Tagging
    //--------------------------------------------------------------------------
    template <typename Tagged, typename = void>
    struct TagExtractorImpl;

    //--------------------------------------------------------------------------

  } // namespace details

  //--- BEGIN Metaprogramming --------------------------------------------------
  /// @ingroup MetaprogrammingBase
  /// @{

  /// Returns how many of the types in `T` exactly match `Target`.
  template <typename Target, typename... T>
  using count_type_in_list = details::count_type_in_list_impl<Target, T...>;

  /// Direct access to the value in `count_type_in_list`.
  template <typename Target, typename... T>
  constexpr unsigned int count_type_in_list_v = count_type_in_list<Target, T...>();

  /// Returns the `N` type of the type list.
  template <std::size_t N, typename... T>
  using typelist_element_type = std::tuple_element<N, std::tuple<T...>>;

  /// Direct access to the value in `typelist_element_type`.
  template <std::size_t N, typename... T>
  using typelist_element_t = typename typelist_element_type<N, T...>::type;

  /// Holds whether the `Target` type is among the ones in the `T` pack.
  template <typename Target, typename... T>
  using type_is_in = details::type_is_in_impl<Target, T...>;

  /// Direct access to the value in `type_is_in`.
  template <typename Target, typename... T>
  constexpr bool type_is_in_v = type_is_in<Target, T...>();

  /// Holds whether the `Target` type is element of the specified `std::tuple`.
  template <typename Target, typename Tuple>
  struct count_type_in_tuple;

  /// @}
  //--- END Metaprogramming ----------------------------------------------------

  //----------------------------------------------------------------------------
  /**
   * @defgroup MetaprogrammingGeneral General utility traits
   * @brief Traits of general utility.
   * @ingroup Metaprogramming
   * @{
   *
   * @details These traits are used in the implementation of the tag-related
   * traits, but they are general enough that are deemed worth being exposed.
   */

  //----------------------------------------------------------------------------
  /**
   * @brief Returns type `TargetClass<U...>` from a `SrcTuple<T...>`.
   * @tparam SrcTuple a tuple-like type enclosing a list of types
   * @tparam Extractor class extracting the `U` type from each `T`
   * @tparam TargetClass the template type of the class to be returned
   * @see `extract_to_tuple_type_t`, `to_tuple`
   *
   * The input template argument `SrcTuple` is an object supporting
   * `std::tuple_size` and `std::tuple_element` to return all its contained
   * types.
   *
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename... T>
   * struct MyFirstClass {
   *   // ...
   * };
   *
   * template <typename... T>
   * struct MySecondClass {
   *   // ...
   * };
   *
   * using DataClass = MyFirstDataClass<int, char, std::string>;
   * using PointerDataClass = util::extract_to_tuple_type
   *   <DataClass, std::add_pointer, MySecondClass>;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will yield a `PointerDataClass` as
   * `MySecondClass<int*, char*, std::string*>` (note that in this partial
   * example `MyFirstClass` is not yet a tuple-like type).
   *
   * The extractor is a trait whose first element is the type to be transformed,
   * and whose member `type` holds the transformed type.
   *
   * @note The extractor must be able to return a type for each and every type
   *       in `SrcTuple`.
   *
   */
  template <typename SrcTuple,
            template <typename T, typename...>
            class Extractor,
            template <typename...> class TargetClass = std::tuple>
  struct extract_to_tuple_type;

  /// Direct access to the type in `extract_to_tuple_type`.
  template <typename SrcTuple,
            template <typename T, typename...>
            class Extractor,
            template <typename...> class TargetClass = std::tuple>
  using extract_to_tuple_type_t =
    typename extract_to_tuple_type<SrcTuple, Extractor, TargetClass>::type;

  /// `extract_to_tuple_type` with no `T`-to-`U` type transformation.
  /// @see `extract_to_tuple_type`, `to_tuple_t`
  template <typename Tuple, template <typename...> class TargetClass = std::tuple>
  using to_tuple = extract_to_tuple_type<Tuple, self_type, TargetClass>;

  /// Direct access to the type in `to_tuple`.
  template <typename Tuple, template <typename...> class TargetClass = std::tuple>
  using to_tuple_t = typename to_tuple<Tuple, TargetClass>::type;

  //----------------------------------------------------------------------------
  /**
   * @brief Returns the index of the element in `Tuple` with the specified type.
   * @tparam Extractor extract the candidate type from a `Tuple` element
   * @tparam Target the type being sought
   * @tparam Tuple tuple-like data structure containing types.
   * @see `index_of_extracted_type_v`, `has_extracted_type`,
   *      `count_extracted_types`, `has_type`,
   *
   * Given a tuple-like type `Tuple`, this traits returns the index of the one
   * element type which "contains" the target type `Target`.
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct ExtractValueType { using type = typename T::value_type; };
   *
   * using tuple_t = std::tuple<std::vector<std::string>, std::vector<int>>;
   *
   * constexpr std::size_t IntegersIndex
   *   = util::index_of_extracted_type<ExtractValueType, int, tuple_t>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `IntegersIndex` will hold value `1`, pointing to the container in `tuple_t`
   * with `int` as `value_type`.
   *
   * If the target type is not present, or if it is present more than once, a
   * compilation error will ensue.
   *
   * @note Currently there is no equivalent trait to ask for the index of the
   *       second or following type; this can be implemented on request.
   */
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  struct index_of_extracted_type;

  /// Direct access to the value in `index_of_extracted_type`.
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  constexpr std::size_t index_of_extracted_type_v =
    index_of_extracted_type<Extractor, Target, Tuple>();

  /// Like `index_of_extracted_type`, but querying the element types directly.
  /// @see `index_of_extracted_type`, `index_of_type_v`, `has_type`,
  ///      `count_types`
  template <typename Target, typename Tuple>
  using index_of_type = index_of_extracted_type<self_type, Target, Tuple>;

  /// Direct access to the value in `index_of_type`.
  template <typename Target, typename Tuple>
  constexpr std::size_t index_of_type_v = index_of_type<Target, Tuple>();

  //----------------------------------------------------------------------------
  /**
   * @brief Returns the element type in `Tuple` with the specified type.
   * @tparam Extractor extract the candidate type from a `Tuple` element
   * @tparam Target the type being sought
   * @tparam Tuple tuple-like data structure containing types.
   * @see `type_with_extracted_type_t`, `index_of_extracted_type`,
   *      `has_extracted_type`, `count_extracted_types`, `has_type`,
   *
   * Given a tuple-like type `Tuple`, this traits returns the one element type
   * which "contains" the target type `Target`.
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct ExtractValueType { using type = typename T::value_type; };
   *
   * using tuple_t = std::tuple<std::vector<std::string>, std::vector<int>>;
   *
   * using int_container_t
   *   = util::type_with_extracted_type<ExtractValueType, int, tuple_t>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `int_container_t` will be `std::vector<int>`, as the type in `tuple_t`
   * with `int` as `value_type`.
   *
   * If the target type is not present, or if it is present more than once, a
   * compilation error will ensue.
   *
   * @note Currently there is no equivalent trait to ask for the second or
   *       following type; this can be implemented on request.
   */
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  using type_with_extracted_type =
    std::tuple_element<index_of_extracted_type_v<Extractor, Target, Tuple>, Tuple>;

  /// Direct access to the value in `type_with_extracted_type`.
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  using type_with_extracted_type_t =
    typename type_with_extracted_type<Extractor, Target, Tuple>::type;

  //----------------------------------------------------------------------------
  /**
   * @brief Trait holding whether an element in `Tuple` type contains `Target`.
   * @tparam Extractor trait exposing the target type in an element
   * @tparam Target the target type to be found
   * @tparam Tuple the tuple-like type to be investigated
   * @see `has_extracted_type_v`, `index_of_extracted_type`,
   *      `count_extracted_types`, `has_type`
   *
   * Given a tuple-like type `Tuple`, this trait returns whether there is at
   * least an element type which "contains" the target type `Target`.
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct ExtractValueType { using type = typename T::value_type; };
   *
   * using tuple_t
   *   = std::tuple<std::vector<float>, std::vector<int>, std::list<int>>;
   *
   * constexpr bool hasIntegerContainers
   *   = util::has_extracted_type<ExtractValueType, int, tuple_t>();
   * constexpr bool hasDoubleContainers
   *   = util::has_extracted_type<ExtractValueType, double, tuple_t>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `hasIntegerContainers` will hold value `true`, meaning there are two types
   * with `int` as `value_type`. Likewise, `hasDoubleContainers` will be
   * `false`.
   *
   * @note Before C++17, the following construct is invalid:
   *       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   *       using tuple_t
   *         = std::tuple<std::vector<float>, std::vector<int>, std::list<int>>;
   *       tuple_t data {
   *         { 0.5, 1.5 },
   *         { 0, 2, 4 },
   *         { 1, 3, 5 }
   *         };
   *       if (util::has_extracted_type<ExtractValueType, int, tuple_t>())
   *         return util::getByExtractedType<ExtractValueType, int>(data);
   *       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *       because `util::getByExtractedType<ExtractValueType, double>(data)`
   *       needs to be compilable, but it is not since it has more than one
   *       `int` extracted types. This is solved in C++17 by the use of a
   *       "constexpr if".
   */
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  struct has_extracted_type;

  /// Direct access to the value in `has_extracted_type`.
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  constexpr bool has_extracted_type_v = has_extracted_type<Extractor, Target, Tuple>();

  /// Like `has_extracted_type`, but querying the element types directly.
  /// @see `has_extracted_type`, `index_of_type`, `count_types`, `has_type_v`
  template <typename Target, typename Tuple>
  using has_type = has_extracted_type<self_type, Target, Tuple>;

  /// Direct access to the value in `has_type`.
  template <typename Target, typename Tuple>
  constexpr bool has_type_v = has_type<Target, Tuple>();

  //----------------------------------------------------------------------------
  /**
   * @brief Returns the value of the element containing the specified type.
   * @tparam Extractor trait extracting the target type from an element
   * @tparam Target the type sought
   * @tparam Tuple a tuple-like type to get values from
   * @param data the tuple to get data from
   * @return the value of the element containing the specified `Target` type
   * @see `count_extracted_types`, `index_of_extracted_type`,
   *      `has_extracted_type`
   *
   * The value of the selected element in `data` is returned, as in `std::get`.
   * The value is taken from the element of data which contains the target type
   * `Target`. The target type is extracted from each element in turn, using the`
   * `Extractor` trait.
   *
   * This function has the same limitations as `util::index_of_extracted_type`.
   *
   * In the example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct ExtractValueType { using type = typename T::value_type; };
   *
   * auto myData = std::make_tuple(
   *   std::vector<float>{ 1.5, 2.5 },
   *   std::vector<int>{ 0, 2, 4 },
   *   std::list<int>{ 1, 3, 5 }
   *   );
   *
   * decltype(auto) realData
   *   = util::getByExtractedType<ExtractValueType, float>(data);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `realData` will be a reference `std::vector<float>&` to the element
   * containing `{ 1.5, 2.5 }`. Also note that in this example requesting the
   * data of value type `int` will cause a compilation error since there is more
   * than one element with that value type.
   */
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  auto getByExtractedType(Tuple const& data) -> decltype(auto);

  /**
   * @brief Traits holding whether elements of `Tuple` have duplicate types.
   * @tparam Extractor trait to extract the type to check from each element
   * @tparam Tuple tuple-type object to check
   * @see `count_extracted_types_v`, `index_of_extracted_type`,
   *      `has_extracted_type`,
   *      `has_duplicate_types`, `has_duplicate_extracted_types_v`
   *
   * This trait holds whether there are multiple occurrences of elements in
   * `Tuple` which have the same extracted type.
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct ExtractValueType { using type = typename T::value_type; };
   *
   * using tuple_t
   *   = std::tuple<std::vector<float>, std::vector<int>, std::list<int>>;
   *
   * constexpr bool hasDuplicate
   *   = util::has_duplicate_extracted_types<ExtractValueType, tuple_t>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `hasDuplicate` is `true` because there are two elements in `tuple_t` with
   * `int` as `value_type` (that is, the type delivered by `ExtractValueType`).
   */
  template <template <typename T, typename...> class Extractor, typename Tuple>
  struct has_duplicate_extracted_types
    : public details::has_duplicate_types_unwrapper<extract_to_tuple_type_t<Tuple, Extractor>> {};

  /// Direct access to the value in `has_duplicate_extracted_types`.
  template <template <typename T, typename...> class Extractor, typename Tuple>
  constexpr bool has_duplicate_extracted_types_v =
    has_duplicate_extracted_types<Extractor, Tuple>();

  /// Like `has_duplicate_extracted_types`, but on the element types directly.
  /// @see `has_duplicate_extracted_types`, `index_of_type`, `count_types`,
  ///      `has_types`, `has_duplicate_types_v`
  template <typename Tuple>
  using has_duplicate_types = has_duplicate_extracted_types<self_type, Tuple>;

  /// Direct access to the value in `has_duplicate_types`.
  template <typename Tuple>
  constexpr bool has_duplicate_types_v = has_duplicate_types<Tuple>();

  /**
   * @brief Counts the elements of a tuple-like type containing a `Target` type.
   * @tparam Extractor trait exposing the target type in an element
   * @tparam Target the target type to be counted
   * @tparam Tuple the tuple-like type to be investigated
   * @see `count_extracted_types_v`, `index_of_extracted_type`,
   *      `has_extracted_type`, `count_types`
   *
   * Given a tuple-like type `Tuple`, this trait returns how many of its
   * elements "contain" the target type `Target`.
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct ExtractValueType { using type = typename T::value_type; };
   *
   * using tuple_t
   *   = std::tuple<std::vector<float>, std::vector<int>, std::list<int>>;
   *
   * constexpr unsigned int nIntegerContainers
   *   = util::count_extracted_types<ExtractValueType, int, tuple_t>();
   * constexpr unsigned int nDoubleContainers
   *   = util::count_extracted_types<ExtractValueType, double, tuple_t>();
   * constexpr unsigned int nFloatContainers
   *   = util::count_extracted_types<ExtractValueType, float, tuple_t>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `nIntegerContainers` will hold value `2`, meaning there are two types
   * with `int` as `value_type`. Likewise, `nDoubleContainers` will be `1`
   * and `nFloatContainers` will be `0`.
   */
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  struct count_extracted_types;

  /// Direct access to the value in `count_extracted_types`.
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  constexpr unsigned int count_extracted_types_v =
    count_extracted_types<Extractor, Target, Tuple>();

  /// Counts the number of `Target` elements in the specified `Tuple`.
  /// @see `count_extracted_types`, `index_of_extracted_type_v`,
  ///      `has_extracted_type_v`, `count_types_v`
  template <typename Target, typename Tuple>
  using count_types = count_extracted_types<self_type, Target, Tuple>;

  /// Direct access to the value in `count_extracted_types`.
  template <typename Target, typename Tuple>
  constexpr unsigned int count_types_v = count_types<Target, Tuple>();

  /// @}

  /**
   * @defgroup MetaprogrammingTagged Tag-related traits
   * @brief  Traits for types with a `tag`.
   * @ingroup Metaprogramming
   *
   * @details Tag-related traits operate on "tagged" types. A tagged type is a
   * type which contains a `tag` type definition, and that type is the tag type.
   *
   * In the examples, the types used are defined as:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct TagA {};
   * struct TagB {};
   * struct TagC {};
   * struct TagD {};
   *
   * using IntTaggedA     = util::add_tag_t<std::vector<int>,    TagA>;
   * using DoubleTaggedB  = util::add_tag_t<std::vector<double>, TagB>;
   * using StringTaggedC  = util::add_tag_t<std::string,         TagC>;
   * using ComplexTaggedA = util::add_tag_t<std::complex<float>, TagA>;
   *
   * using Tuple_t     = std::tuple<IntTaggedA, DoubleTaggedB, StringTaggedC>;
   * using DuplTuple_t = std::tuple<IntTaggedA, DoubleTaggedB, ComplexTaggedA>;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Note that `DuplTuple_t` has two elements with the same tag (`TagA`).
   *
   * @note All the traits and function here will generate a compilation error
   *       if any of the elements is not tagged.
   *
   * @{
   */

  /**
   * @brief A type with a specified tag.
   * @tparam T the base type being tagged
   * @tparam Tag the tag to be assigned to the new type
   * @see `add_tag`, `TagExtractor`
   *
   * The new type `TaggedType<T, Tag>` inherits from `T` (including its
   * constructors) and adds (and hides existing) `tag` type with value `Tag`.
   * This type is suitable to be used in the tag-related traits.
   *
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct VectorTag {};
   * using TaggedVector = TaggedType<std::vector<int>, VectorTag>;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename T, typename Tag>
  struct TaggedType : public T {

    // Forward all constructors
    template <typename... Args>
    TaggedType(Args&&... args) : T(std::forward<Args>(args)...)
    {}

    using tag = Tag;       ///< Tag of this object.
    using tagged_type = T; ///< Type of the object which was tagged.

  }; // struct TaggedType

  /// Trait holding a type derived from `T`, tagged with `Tag`.
  /// @see `TaggedType`
  template <typename T, typename Tag>
  struct add_tag {
    using type = TaggedType<std::remove_reference_t<T>, Tag>;
  };

  /// Direct access to the type contained in `add_tag`.
  /// @see `TaggedType`
  template <typename T, typename Tag>
  using add_tag_t = typename add_tag<T, Tag>::type;

  /// Trait holding the type contained in a `TaggedType` (or the type itself).
  template <typename Tagged>
  struct remove_tag {
    using type = Tagged;
  };

  template <typename T, typename Tag>
  struct remove_tag<TaggedType<T, Tag>> {
    using type = typename TaggedType<T, Tag>::tagged_type;
  };

  /// Direct access to the type contained in `remove_tag`.
  template <typename Tagged>
  using remove_tag_t = typename remove_tag<Tagged>::type;

  /**
   * @brief "Converts" `obj` to an object with tag `Tag`.
   * @tparam Tag tag to be added to the object
   * @tparam T type of the object to be tagged (implicitly deduced)
   * @param obj (l-value) reference to the object to be tagged
   * @return a reference to `obj`, reinterpreted as tagged
   * @see `TaggedType`, `add_tag`
   *
   * The returned object is the same as `obj`, reinterpreted as a different type
   * derived from `T` and tagged with `Tag`.
   */
  template <typename Tag, typename T>
  auto makeTagged(T& obj) -> decltype(auto)
  {
    return static_cast<add_tag_t<T, Tag>&>(obj);
  }

  /**
   * @brief "Converts" `obj` to an object with tag `Tag`.
   * @tparam Tag tag to be added to the object
   * @tparam T type of the object to be tagged (implicitly deduced)
   * @param obj (l-value) constant reference to the object to be tagged
   * @return a reference to `obj`, reinterpreted as tagged
   * @see `TaggedType`, `add_tag`
   *
   * The returned object is the same as `obj`, reinterpreted as a different type
   * derived from `T` and tagged with `Tag`.
   */
  template <typename Tag, typename T>
  auto makeTagged(T const& obj) -> decltype(auto)
  {
    return static_cast<add_tag_t<T const, Tag> const&>(obj);
  }

  /**
   * @brief "Converts" `obj` to an object with tag `Tag`.
   * @tparam Tag tag to be added to the object
   * @tparam T type of the object to be tagged (implicitly deduced)
   * @param obj (r-value) reference to the object to be tagged
   * @return a reference to `obj`, reinterpreted as tagged
   * @see `TaggedType`, `add_tag`
   *
   * The returned object is an object of a new type, derived from `T` , with
   * data copied from `obj`, and tagged with `Tag`.
   * The argument object may be a temporary.
   */
  template <typename Tag, typename T>
  auto makeTagged(T const&& obj) -> decltype(auto)
  {
    return add_tag_t<T, Tag>(obj); /* copy, since it's constant */
  }

  /**
   * @brief "Converts" `obj` to an object with tag `Tag`.
   * @tparam Tag tag to be added to the object
   * @tparam T type of the object to be tagged (implicitly deduced)
   * @param obj the object to be tagged
   * @return a reference to `obj`, reinterpreted as tagged
   * @see `TaggedType`, `add_tag`
   *
   * The returned object is the same as `obj`, reinterpreted as a different type
   * derived from `T` and tagged with `Tag`.
   * The returned object is a temporary of the new type, whose content is moved
   * (`std::move()`) from the argument object `obj`.
   */
  template <typename Tag, typename T>
  auto makeTagged(T&& obj) -> decltype(auto)
  {
    return add_tag_t<T, Tag>(std::move(obj));
  }

  /// "Converts" a tagged type back to its original type.
  template <typename Tagged>
  auto removeTag(Tagged& tagged) -> decltype(auto)
  {
    return static_cast<remove_tag_t<Tagged>&>(tagged);
  }

  /// "Converts" a tagged type back to its original type.
  template <typename Tagged>
  auto removeTag(Tagged const& tagged) -> decltype(auto)
  {
    return static_cast<remove_tag_t<Tagged> const&>(tagged);
  }

  /// "Converts" a tagged type back to its original type.
  template <typename Tagged>
  auto removeTag(Tagged const&& tagged) -> decltype(auto)
  {
    return static_cast<remove_tag_t<Tagged> const&&>(tagged);
  }

  /// "Converts" a tagged type back to its original type.
  template <typename Tagged>
  auto removeTag(Tagged&& tagged) -> decltype(auto)
  {
    return static_cast<remove_tag_t<Tagged>&&>(tagged);
  }

  /// Tag class parametrized by a sequence of numbers.
  template <std::size_t...>
  struct TagN {};

  /**
   * @brief Extracts the tag from a type.
   * @tparam Tagged type of tagged type
   *
   * This trait holds the tag type of `Tagged` (that is, `Tagged::tag`) as
   * `type` member.
   * If `Tagged` type has no `tag` type member, a compilation error will ensue.
   */
  template <typename Tagged>
  struct TagExtractor : public details::TagExtractorImpl<Tagged> {};

  /// Trait holding the tag of `Tagged` as `type`.
  template <typename Tagged>
  using tag_of = TagExtractor<Tagged>;

  /// Direct access to the type in `tag_of`.
  template <typename Tagged>
  using tag_of_t = typename tag_of<Tagged>::type;

  /// Returns a tuple with all the tags from `SrcTuple`.
  template <typename SrcTuple>
  using extract_tags = extract_to_tuple_type<SrcTuple, TagExtractor>;

  /// Direct access to the type in `extract_tags`.
  template <typename SrcTuple>
  using extract_tags_t = typename extract_tags<SrcTuple>::type;

  /**
   * @brief Trait holding the index of the element of `Tuple` with tag `Tag`.
   * @tparam Tag the sought tag
   * @tparam Tuple the tuple-like type holding the elements to check
   * @see `index_of_tag_v`, `count_tags`, `has_duplicate_tags`, `has_tag`
   *
   * Given a tuple-like type `Tuple`, this traits returns the index of the one
   * element type tagged with `Tag`.
   * If the target type is not present, or if it is present more than once, a
   * compilation error will ensue.
   *
   * For example (see above for the definitions):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * constexpr std::size_t TagAindex = util::index_of_tag<TagA, Tuple_t>();
   * constexpr std::size_t TagBindex = util::index_of_tag<TagB, Tuple_t>();
   * constexpr std::size_t TagBduplIndex
   *   = util::index_of_tag<TagB, DuplTuple_t>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `TagAindex` will hold value `0`, pointing to the container in `Tuple_t`
   * of type `IntTaggedA`, while `TagBindex` will be `1` and `TagBduplIndex`
   * will be also `1`. Instead, the expression
   * `util::index_of_tag<TagD, Tuple_t>()` would not compile because no element
   * in `Tuple_t` is tagged with `TagD`, and the expression
   * `util::index_of_tag<TagA, DuplTuple_t>()` would not compile because two
   * elements of `DuplTuple_t` are tagged `TagA`.
   *
   * @note Currently there is no equivalent trait to ask for the index of the
   *       second or following type, allowing for duplicate tags; this can be
   *       implemented on request.
   */
  template <typename Tag, typename Tuple>
  using index_of_tag = index_of_extracted_type<TagExtractor, Tag, Tuple>;

  /// Direct access to the value in `index_of_tag`.
  template <typename Tag, typename Tuple>
  constexpr std::size_t index_of_tag_v = index_of_tag<Tag, Tuple>();

  /**
   * @brief Trait holding the type of the element of `Tuple` with tag `Tag`.
   * @tparam Tag the sought tag
   * @tparam Tuple the tuple-like type holding the elements to check
   * @see `type_with_tag_t`, `index_of_tag`, `count_tags`, `has_duplicate_tags`,
   *      `has_tag`
   *
   * Given a tuple-like type `Tuple`, this traits returns the one element type
   * tagged with `Tag`.
   * If the target type is not present, or if it is present more than once, a
   * compilation error will ensue.
   *
   * For example (see above for the definitions):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using TagA_t = typename util::type_with_tag<TagA, Tuple_t>::type;
   * using TagB_t = typename util::type_with_tag<TagB, Tuple_t>::type;
   * using TagBdupl_t = typename util::type_with_tag<TagB, DuplTuple_t>::type;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `TagA_t` will be `IntTaggedA`, the type in `Tuple_t` with `TagA`, while
   * `TagB_t` will be `DoubleTaggedB` and `TagBdupl_t` will be also
   * `DoubleTaggedB`. Instead, the type `util::type_with_tag<TagD, Tuple_t>`
   * would not compile because no element in `Tuple_t` is tagged with `TagD`,
   * and the type `util::type_with_tag<TagA, DuplTuple_t>()` would not compile
   * because two elements of `DuplTuple_t` are tagged `TagA`.
   *
   * @note Currently there is no equivalent trait to ask for the second or
   *       following type, allowing for duplicate tags; this can be implemented
   *       on request.
   */
  template <typename Tag, typename Tuple>
  using type_with_tag = type_with_extracted_type<TagExtractor, Tag, Tuple>;

  /// Direct access to the value in `type_with_tag`.
  template <typename Tag, typename Tuple>
  using type_with_tag_t = typename type_with_tag<Tag, Tuple>::type;

  /**
   * @brief Trait informing if there are elements in `Tuple` with tag `Tag`.
   * @tparam Tag the sought tag
   * @tparam Tuple the tuple-like type holding the elements to check
   * @see `index_of_tag`, `count_tags`, `has_duplicate_tags`, `has_tag_v`
   *
   * Given a tuple-like type `Tuple`, this traits returns whether there is at
   * least one element type tagged with `Tag`.
   *
   * For example (see above for the definitions):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * constexpr bool hasTagA     = util::index_of_tag<TagA, Tuple_t    >();
   * constexpr bool hasTagB     = util::index_of_tag<TagB, Tuple_t    >();
   * constexpr bool hasTagC     = util::index_of_tag<TagC, Tuple_t    >();
   * constexpr bool hasTagD     = util::index_of_tag<TagD, Tuple_t    >();
   * constexpr bool hasTagAdupl = util::index_of_tag<TagA, DuplTuple_t>();
   * constexpr bool hasTagBdupl = util::index_of_tag<TagB, DuplTuple_t>();
   * constexpr bool hasTagCdupl = util::index_of_tag<TagC, DuplTuple_t>();
   * constexpr bool hasTagDdupl = util::index_of_tag<TagD, DuplTuple_t>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `hasTagA`, `hasTagB` and `hasTagC`, will be `true` and `hasTagD` will be
   * `false`. Likewise, `hasTagAdupl` and `hasTagBdupl` will be `true`, while
   * `hasTagCdupl` and `hasTagDdupl` will be `false`.
   */
  template <typename Tag, typename Tuple>
  using has_tag = has_extracted_type<TagExtractor, Tag, Tuple>;

  /// Direct access to the value in `has_tag`.
  template <typename Tag, typename Tuple>
  constexpr bool has_tag_v = has_tag<Tag, Tuple>();

  /**
   * @brief Trait counting the elements in `Tuple` with tag `Tag`.
   * @tparam Tag the sought tag
   * @tparam Tuple the tuple-like type holding the elements to check
   * @see `index_of_tag`, `count_tags_v`, `has_duplicate_tags`, `has_tag`
   *
   * Given a tuple-like type `Tuple`, this traits returns the number of element
   * types tagged with `Tag`.
   *
   * For example (see above for the definitions):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * constexpr unsigned int nTagA     = util::count_tags<TagA, Tuple_t    >();
   * constexpr unsigned int nTagB     = util::count_tags<TagB, Tuple_t    >();
   * constexpr unsigned int nTagC     = util::count_tags<TagC, Tuple_t    >();
   * constexpr unsigned int nTagD     = util::count_tags<TagD, Tuple_t    >();
   * constexpr unsigned int nTagAdupl = util::count_tags<TagA, DuplTuple_t>();
   * constexpr unsigned int nTagBdupl = util::count_tags<TagB, DuplTuple_t>();
   * constexpr unsigned int nTagCdupl = util::count_tags<TagC, DuplTuple_t>();
   * constexpr unsigned int nTagDdupl = util::count_tags<TagD, DuplTuple_t>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `nTagA`, `nTagB` and `nTagC`, will be `1` and `nTagD` will be `0`.
   * Likewise, `nTagAdupl` will be `1`, `nTagBdupl` will be `2`, while
   * `nTagCdupl` and `nTagDdupl` will be `0`.
   */
  template <typename Tag, typename Tuple>
  using count_tags = count_extracted_types<TagExtractor, Tag, Tuple>;

  /// Direct access to the value in `count_tags`.
  template <typename Tag, typename Tuple>
  constexpr unsigned int count_tags_v = count_tags<Tag, Tuple>();

  /**
   * @brief Trait reporting if multiple elements in `Tuple` have the same tag.
   * @tparam Tuple the tuple-like type holding the elements to check
   * @see `index_of_tag`, `count_tags`, `has_duplicate_tags_v`, `has_tag`
   *
   * Given a tuple-like type `Tuple`, this traits returns whether any of the
   * tags in the elements appears more than once.
   *
   * For example (see above for the definitions):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * constexpr bool hasDuplTags     = util::has_duplicate_tags<Tuple_t    >();
   * constexpr bool hasDuplTagsDupl = util::has_duplicate_tags<DuplTuple_t>();
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `hasDuplTags`, will be `false` and `hasDuplTagsDupl` will be `true`.
   */
  template <typename Tuple>
  using has_duplicate_tags = has_duplicate_extracted_types<TagExtractor, Tuple>;

  /// Direct access to the value in `has_duplicate_tags`.
  template <typename Tuple>
  constexpr bool has_duplicate_tags_v = has_duplicate_tags<Tuple>();

  /**
   * @brief Returns the object with the specified tag.
   * @tparam Tag the sought tag
   * @tparam Tuple the tuple-like type holding the data
   * @see `index_of_tag`, `count_tags`, `has_duplicate_tags`, `has_tag`
   *
   * This function operates in a fashion similar to `std::get()`, where instead
   * of specifying the index or type of the object to retrieve, the tag is
   * specified.
   *
   * For example (see above for the definitions):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * Tuple_t data(
   *   { 0, 1, 2 },   // std::vector<int> assigned to TagA
   *   { 0.5, 1.5 },  // std::vector<double> assigned to TagB
   *   "middle point" // std::string assigned to TagC
   *   );
   * decltype(auto) TagBdata = util::getByTag<TagB>(data);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * `TagBdata` will be a reference to a `std::vector<double>` in `data`,
   * currently with values `{ 0.5, 1.5 }".
   * The attempt to use `util::getByTag()` on an argument with duplicate tags
   * (like `DuplTuple_t`) will ensue a compilation error.
   */
  template <typename Tag, typename Tuple>
  auto getByTag(Tuple const& data) -> decltype(auto)
  {
    return getByExtractedType<TagExtractor, Tag>(data);
  }

  /// @}

} // namespace util

//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
namespace util {

  namespace details {

    //--------------------------------------------------------------------------
    //
    // CAUTION: prolonged exposition to this code may result into loss of sight.
    //
    //--------------------------------------------------------------------------
    //--- implementation details (most are not documented)
    //--------------------------------------------------------------------------
    //--- General utilities
    //--------------------------------------------------------------------------
    template <typename Target, typename... T>
    struct count_type_in_list_impl : public std::integral_constant<unsigned int, 0U> {};

    template <typename Target, typename First, typename... Others>
    struct count_type_in_list_impl<Target, First, Others...>
      : public std::integral_constant<unsigned int,
                                      count_type_in_list_impl<Target, Others...>::value> {};

    template <typename Target, typename... Others>
    struct count_type_in_list_impl<Target, Target, Others...>
      : public std::integral_constant<unsigned int,
                                      count_type_in_list_impl<Target, Others...>::value + 1> {};

    //--------------------------------------------------------------------------
    //--- Courtesy traits
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //---  Tagging
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    template <typename Target, typename... T>
    struct type_is_in_impl : public std::false_type {};

    template <typename Target, typename First, typename... Others>
    struct type_is_in_impl<Target, First, Others...>
      : public std::integral_constant<bool, type_is_in_impl<Target, Others...>::value> {};

    template <typename Target, typename... Others>
    struct type_is_in_impl<Target, Target, Others...> : public std::true_type {};

    //--------------------------------------------------------------------------
    /**
     * @brief Implementation for `extract_to_tuple_type`.
     * @tparam TargetClass as in `extract_to_tuple_type`
     * @tparam Extractor as in `extract_to_tuple_type`
     * @tparam Tuple as in `extract_to_tuple_type`
     * @tparam I index of the type being extracted at this iteration
     * @tparam N total number of types to be extracted
     * @tparam T types already extracted
     *
     * The recursive implementation adds the type extracted from the `I`-th
     * `Tuple` element to the list of T.
     *
     * The total number of iterations needs to be passed as argument since the
     * end-of-iteration specialization needs to define `I` as
     * `std::tuple_size_v<Tuple>`, but C++ (at least up to C++14) does not allow
     * for a typed template argument to be so complicate ().
     * This may change in the future,
     */
    template <template <typename...> class TargetClass,
              template <typename T, typename...>
              class Extractor,
              typename Tuple,
              std::size_t I,
              std::size_t N,
              typename... T>
    struct extract_to_tuple_type_impl {
      using type = typename extract_to_tuple_type_impl<
        TargetClass,
        Extractor,
        Tuple,
        (I + 1),
        N,
        T...,
        typename Extractor<std::tuple_element_t<I, Tuple>>::type>::type;
    }; // extract_to_tuple_type_impl

    /// End-of-recursion specialization: packs all types already extracted into the target type.
    template <template <typename...> class TargetClass,
              template <typename T, typename...>
              class Extractor,
              typename Tuple,
              std::size_t N,
              typename... T>
    struct extract_to_tuple_type_impl<TargetClass, Extractor, Tuple, N, N, T...> {
      using type = TargetClass<T...>;
    }; // extract_to_tuple_type_impl<N>

    //--------------------------------------------------------------------------

    // TODO try alternative simple implementation:
    template <typename SrcTuple,
              template <typename T, typename...>
              class Extractor,
              template <typename...>
              class TargetClass,
              std::size_t... I>
    struct extract_to_tuple_type_impl_simple {
      using type = TargetClass<typename Extractor<std::tuple_element_t<I, SrcTuple>>::type...>;
    }; // extract_to_tuple_type_impl<N>

    //--------------------------------------------------------------------------

    // Part of implementation for `index_of_extracted_type`.
    template <template <typename T, typename...> class Extractor,
              typename Target,
              std::size_t I,
              typename Elem,
              typename Tuple>
    struct index_of_extracted_type_impl;

    // Part of implementation for `index_of_extracted_type`.
    template <template <typename T, typename...> class Extractor,
              typename Target,
              std::size_t N,
              std::size_t I,
              typename Tuple>
    struct index_of_extracted_type_checked
      : public index_of_extracted_type_impl<
          Extractor,
          Target,
          I,
          typename Extractor<std::tuple_element_t<I, Tuple>>::type,
          Tuple> {};

    // Part of implementation for `index_of_extracted_type`.
    template <template <typename T, typename...> class Extractor,
              typename Target,
              std::size_t N,
              typename Tuple>
    struct index_of_extracted_type_checked<Extractor, Target, N, N, Tuple>
      : public std::integral_constant<std::size_t, N> {};
    template <template <typename T, typename...> class Extractor,
              typename Target,
              std::size_t I,
              typename Elem,
              typename Tuple>
    struct index_of_extracted_type_impl
      : public index_of_extracted_type_checked<Extractor,
                                               Target,
                                               std::tuple_size<Tuple>::value,
                                               (I + 1),
                                               Tuple> {};

    // Part of implementation for `index_of_extracted_type`.
    template <template <typename T, typename...> class Extractor,
              typename Target,
              std::size_t I,
              typename Tuple>
    struct index_of_extracted_type_impl<Extractor, Target, I, Target, Tuple>
      : public std::integral_constant<std::size_t, I> {
      static constexpr std::size_t N = std::tuple_size<Tuple>();
      static_assert(I < N, "Internal logic error.");
    };

    // Part of implementation for `index_of_extracted_type`.
    template <template <typename T, typename...> class Extractor,
              typename Target,
              std::size_t N,
              std::size_t After,
              typename Tuple>
    struct index_of_extracted_type_checked_after
      : public index_of_extracted_type_checked<Extractor, Target, N, (After + 1), Tuple> {};

    // Part of implementation for `index_of_extracted_type`.
    template <template <typename T, typename...> class Extractor,
              typename Target,
              std::size_t N,
              typename Tuple>
    struct index_of_extracted_type_checked_after<Extractor, Target, N, N, Tuple>
      : public std::integral_constant<std::size_t, N> {};

    // Part of implementation for `index_of_extracted_type`.
    // This implementation relies on std::tuple_size and std::tuple_element;
    // an implementation assuming Tuple to be std::tuple would be more efficient...
    template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
    struct index_of_type_base
      : public std::integral_constant<std::size_t,
                                      index_of_extracted_type_checked<Extractor,
                                                                      Target,
                                                                      std::tuple_size<Tuple>::value,
                                                                      0U,
                                                                      Tuple>::value> {};

    // Part of implementation for `index_of_extracted_type`.
    template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
    struct index_of_type_helper {
      static constexpr std::size_t N = std::tuple_size<Tuple>::value;
      static constexpr std::size_t value = index_of_type_base<Extractor, Target, Tuple>();

      static_assert(value < N, "The specified tuple does not have the sought type.");
      static_assert(index_of_extracted_type_checked_after<Extractor, Target, N, value, Tuple>() >=
                      N,
                    "The specified tuple has more than one element with the sought type.");
    }; // struct index_of_type_helper

    //--------------------------------------------------------------------------
    // Part of implementation of `has_duplicate_types`.
    template <typename Tuple, typename... T>
    struct has_duplicate_types_impl : public std::false_type {};

    template <typename Tuple, typename First, typename... Others>
    struct has_duplicate_types_impl<Tuple, First, Others...>
      : public std::integral_constant<bool,
                                      (count_type_in_tuple<First, Tuple>::value > 1U) ||
                                        has_duplicate_types_impl<Tuple, Others...>::value> {};

    // Part of implementation of `has_duplicate_types`.
    template <typename... T>
    struct has_duplicate_types_unwrapper<std::tuple<T...>>
      : public has_duplicate_types_impl<std::tuple<T...>, T...> {};

    //--------------------------------------------------------------------------
    // Part of implementation of `count_extracted_types` (might be exposed...).
    template <template <typename T, typename...> class Extractor, typename Target, typename... Tags>
    struct count_tags_in_list
      : public count_type_in_list<typename Extractor<Target>::type, Tags...> {};

    template <template <typename T, typename...> class Extractor,
              typename Target,
              typename... TagTuple>
    struct count_tags_in_tuple;

    template <template <typename T, typename...> class Extractor, typename Target, typename... Tags>
    struct count_tags_in_tuple<Extractor, Target, std::tuple<Tags...>>
      : public count_type_in_tuple<typename Extractor<Target>::type, Tags...> {};

    //--------------------------------------------------------------------------
    //--- Tag traits
    //--------------------------------------------------------------------------
    template <typename Tagged, typename /* = void */>
    struct TagExtractorImpl {
      static_assert(always_false_type<Tagged>(), "This type is not tagged.");
    };

    template <typename Tagged>
    struct TagExtractorImpl<
      Tagged,
      std::enable_if_t<always_true_type<typename std::remove_reference_t<Tagged>::tag>::value>> {
      using type = typename std::remove_reference_t<Tagged>::tag;
    };

    //--------------------------------------------------------------------------

  } // namespace details

  //----------------------------------------------------------------------------
  //--- Implementation of the exposed traits
  //----------------------------------------------------------------------------
  //--- General utilities
  //----------------------------------------------------------------------------
  template <typename Target, typename Tuple>
  struct count_type_in_tuple {
    static_assert(always_false_type<Tuple>(),
                  "count_type_in_tuple requires an instance of std::tuple");
  };

  // specialization: only works for std::tuples
  template <typename Target, typename... T>
  struct count_type_in_tuple<Target, std::tuple<T...>> : public count_type_in_list<Target, T...> {};

  //----------------------------------------------------------------------------
  //--- Courtesy traits
  //----------------------------------------------------------------------------
  template <typename SrcTuple,
            template <typename T, typename...>
            class Extractor,
            template <typename...>
            class TargetClass /* = std::tuple */
            >
  struct extract_to_tuple_type {
  private:
    static constexpr std::size_t N = std::tuple_size<SrcTuple>();

  public:
    using type =
      typename details::extract_to_tuple_type_impl<TargetClass, Extractor, SrcTuple, 0U, N>::type;
  }; // extract_to_tuple_type

  //----------------------------------------------------------------------------
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  struct index_of_extracted_type
    : public std::integral_constant<
        std::size_t,
        details::index_of_type_helper<Extractor, Target, Tuple>::value> {};

  //----------------------------------------------------------------------------
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  struct has_extracted_type
    : public std::integral_constant<bool,
                                    (details::index_of_type_base<Extractor, Target, Tuple>::value <
                                     std::tuple_size<Tuple>::value)> {};

  //----------------------------------------------------------------------------
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  struct count_extracted_types
    : public count_type_in_tuple<Target, extract_to_tuple_type_t<Tuple, Extractor>> {};

  //----------------------------------------------------------------------------
  template <template <typename T, typename...> class Extractor, typename Target, typename Tuple>
  auto getByExtractedType(Tuple const& data) -> decltype(auto)
  {
    using std::get;
    return get<index_of_extracted_type_v<Extractor, Target, Tuple>>(data);
  }

  //----------------------------------------------------------------------------
  //--- Tagging utilities
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
} // namespace util

//------------------------------------------------------------------------------

#endif // LARDATA_UTILITIES_TUPLELOOKUPBYTAG_H
