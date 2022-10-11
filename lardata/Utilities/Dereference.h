/**
 * @brief  Functor to dereference an object if the object is a pointer
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 23rd, 2015
 *
 */

#ifndef LARDATA_UTILITIES_DEREFERENCE_H
#define LARDATA_UTILITIES_DEREFERENCE_H 1

// C/C++ standard libraries
#include <type_traits>

/// LArSoft namespace
namespace lar {
  /// LArSoft utilities namespace
  namespace util {
    /// Implementation of utility details
    namespace details {

      /// Class compiling only if type T exists (then, it's std::true_type)
      template <typename T>
      struct is_type : public std::true_type {};

      //@{
      /**
       * @brief Class defining whether the specified type can be dereferenced
       * @tparam T the type to be tested for dereferenciation
       * @tparam Enable triggers the conditional compilation
       *
       * Usage:
       *
       *     static_assert(!has_dereference_class<int>::value, "error");
       *     static_assert(has_dereference_class<int*>::value, "error");
       *
       *
       * Implementation details
       * -----------------------
       *
       * The class is implemented by two definitions:
       * 1. a general one, always equivalent to std::false_type
       * 2. a special one, always equivalent to std::true_type, defined only
       *    if the type T can be dereferenced (that is, if *(T()) is valid)
       *
       * The selection mechanism will always select the second form if
       * available, since it is more specialized.
       */
      template <typename T, typename Enable = void>
      struct has_dereference_class : public std::false_type {};

      template <typename T>
      struct has_dereference_class<
        T,
        typename std::enable_if<is_type<decltype(*(T()))>::value, void>::type>
        : public std::true_type {};
      //@}

      //@{
      /**
       * @brief Class holding the type dereferenced from an object of type T
       * @tparam T the type to be tested
       * @tparam CanDereference whether T can be dereferenced or not
       * @see lar::util::dereferenced_type
       *
       * This class can be used as:
       *
       *      static_assert(std::is_same<
       *        typename dereferenced_type<int, false>::type,
       *        int
       *        >::value,
       *        "error");
       *      static_assert(std::is_same<
       *        typename dereferenced_type<int*, true>::type,
       *        int&
       *        >::value,
       *        "error");
       *
       * The type is contained in the `type` member of the class.
       * The `type` is precisely what is obtained by dereferencing T (that is
       * often a reference to a type).
       * Note that the second parameter must express correctly whether the first
       * type can be dereferenced or not. Clearly, this class is not very useful
       * by itself, since we have to tell it how to do the trick
       * (lar::util::dereferenced_type puts everything together).
       * It is used in association with has_dereference_class.
       *
       * The implementation is the usual two-definition one:
       * 1. generic class, always with the same type T
       *    (that we associate with CanDereference false)
       * 2. specialized class, with CanDereference explicitly set to true,
       *    that has the type *T
       *
       */
      template <typename T, bool CanDereference>
      struct dereferenced_type {
        //  using type = typename std::add_lvalue_reference<T>::type;
        using type = T;
      }; // dereferenced_type <>

      template <typename T>
      struct dereferenced_type<T, true> {
        using type = decltype(*T());
      }; // dereferenced_type<T, true>
      //@}

      //@{
      /**
       * @brief Functor returning the dereferenced value of the argument
       * @tparam T type of the argument
       * @tparam CanDereference whether T can be dereferenced or not
       *
       * The functor defines a call operator returning:
       * 1. a reference to the value pointed by the argument (CanDereference)
       * 2. or, a reference to the argument itself (!CanDereference).
       *
       * The behaviour is therefore determined by the CanDereference parameter.
       * Note that the second parameter must express correctly whether the first
       * type can be dereferenced or not. Clearly, this class is not very useful
       * by itself, since we have to tell it how to do the trick.
       * It is used in association with has_dereference_class.
       *
       * This class is state-less.
       */
      template <typename T, bool CanDereference>
      struct dereference_class {
        using argument_type = T;
        using reference_type = typename std::add_lvalue_reference<
          typename dereferenced_type<T, CanDereference>::type>::type;

        reference_type operator()(argument_type& ref) const { return ref; }
      }; // dereference_class<T, bool>

      template <typename T>
      struct dereference_class<T, true> {
        using argument_type = T;
        using reference_type =
          typename std::add_lvalue_reference<typename dereferenced_type<T, true>::type>::type;

        reference_type operator()(argument_type& ref) const { return *ref; }
      }; // dereference_class<T, true>
      //@}

      //@{
      /**
       * @brief Functor returning the pointer to a value in the argument
       * @tparam T type of the argument
       * @tparam CanDereference whether T can be dereferenced or not
       *
       * The functor defines a call operator returning:
       * 1. a pointer to the value pointed by the argument (CanDereference), or
       * 2. a pointer to the argument itself (!CanDereference).
       *
       * The behaviour is therefore determined by the CanDereference parameter.
       * Note that the second parameter must express correctly whether the first
       * type can be dereferenced or not. Clearly, this class is not very useful
       * by itself, since we have to tell it how to do the trick.
       * It is used in association with has_dereference_class.
       *
       * This class is state-less.
       */
      template <typename T, bool CanDereference>
      struct make_pointer_class {
        using argument_type = T;
        using pointer_type =
          typename std::add_pointer<typename dereferenced_type<T, CanDereference>::type>::type;

        pointer_type operator()(argument_type& ref) const { return &ref; }
      }; // make_pointer_class<T, bool>

      template <typename T>
      struct make_pointer_class<T, true> {
        using argument_type = T;
        using pointer_type =
          typename std::add_pointer<typename dereferenced_type<T, true>::type>::type;

        pointer_type operator()(argument_type& ref) const { return &*ref; }
      }; // make_pointer_class<T, true>
      //@}

    } // namespace details

    /** ************************************************************************
     * @brief Class defining the dereferenced type of the specified type
     * @tparam T the type to be tested for dereferenciation
     *
     * Usage:
     *
     *     static_assert(
     *       std::is_same<typename lar::util::dereferenced_type<int>::type, int>,
     *       "error"
     *       );
     *     static_assert(
     *       std::is_same<typename lar::util::dereferenced_type<int*>::type, int&>,
     *       "error"
     *       );
     *
     * The type is contained in the `type` member of the class.
     * The `type` is precisely what is obtained by dereferencing T (that is
     * often a reference to a type).
     *
     * The interaction with const `T` types and with constant dereference
     * operators (`T::operator* () const`) has not been investigated.
     */
    template <typename T>
    struct dereferenced_type
      : public details::dereferenced_type<T, details::has_dereference_class<T>::value> {};

    /** ************************************************************************
     * @brief Returns the value pointed by the argument, or the argument itself
     * @tparam T the type of the argument
     * @param v the value to be dereferenced
     * @return the value v points to, if any, or v itself
     *
     * This function allows the use of the same template code to process both
     * pointers and pointed values.
     * For example:
     *
     *     template <typename T>
     *     std::vector<int> extract_int(std::vector<T> const& from) {
     *       std::vector<int> v;
     *       v.reserve(from.size());
     *       std::transform(from.begin(), from.end(), std::back_inserter(v),
     *         util::dereference);
     *       return v;
     *     }
     *
     *     int value = 10;
     *     std::vector<int> v(10, value);
     *     std::vector<int*> v_ptr(10, &value);
     *
     *     extract_int(v);
     *     extract_int(v_ptr);
     *
     * shows that the same function can be used on vectors of integers and of
     * their pointers.
     */
    template <typename T>
    inline
      typename details::dereference_class<T,
                                          details::has_dereference_class<T>::value>::reference_type
      dereference(T& v)
    {
      return details::dereference_class<T, details::has_dereference_class<T>::value>()(v);
    }

    /** ************************************************************************
     * @brief Returns a pointer to the value of argument, or the argument itself
     * @tparam T the type of the argument
     * @param v the value to be taken the pointer of
     * @return a pointer to v, or v itself if it is already a C pointer
     *
     * This function allows the use of the same template code to process both
     * pointers and pointed values.
     * For example:
     *
     *     template <typename T>
     *     std::vector<int*> extract_int(std::vector<T> const& from) {
     *       std::vector<int*> v;
     *       v.reserve(from.size());
     *       std::transform(from.begin(), from.end(), std::back_inserter(v),
     *         util::make_pointer);
     *       return v;
     *     }
     *
     *     int value = 10;
     *     std::vector<int> v(10, value);
     *     std::vector<int*> v_ptr(10, &value);
     *
     *     extract_int(v);
     *     extract_int(v_ptr);
     *
     * shows that the same function can be used on vectors of integers and of
     * their pointers.
     */
    template <typename T>
    inline
      typename details::make_pointer_class<T,
                                           details::has_dereference_class<T>::value>::pointer_type
      make_pointer(T& v)
    {
      return details::make_pointer_class<T, details::has_dereference_class<T>::value>()(v);
    }

  } // namespace util

} // namespace lar

#endif // LARDATA_UTILITIES_DEREFERENCE_H
