/**
 * @file    hexfloat.h
 * @brief   Helper to support output of real numbers in base 16
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    April 29, 2016
 *
 * @note This is going to be obsolete and deprecated as soon as GCC supports
 * `std::hexfloat` from C++ standard (gcc 5.1).
 *
 * @note The manipulator `lar::hexfloat` is *not* equivalent to `std::hexfloat`
 * in that the former takes one argument and formats only that argument,
 * while the second takes no argument and all the following floats are formatted
 * in base 16.
 */

#ifndef LARDATA_RECOBASEART_DUMPERS_HEXFLOAT_H
#define LARDATA_RECOBASEART_DUMPERS_HEXFLOAT_H 1

// C/C++ standard libraries
#include <cstdio>      // std::snprintf()
#include <iosfwd>      // std::ostream
#include <type_traits> // std::is_same<>
#include <utility>     // std::forward()

namespace lar {

  namespace details {

    template <typename T>
    struct OptionalHexFloatFormatter {
    public:
      using real_t = T;

      // sprintf() does not support floats;
      // it will probably do some bad thing with them
      static_assert(!std::is_same<std::decay_t<T>, float>::value,
                    "OptionalHexFloatFormatter does not support float values");

      OptionalHexFloatFormatter(real_t v, bool start_active = true) : active(start_active), value(v)
      {}

      /// Prints the value set at construction
      template <typename Stream>
      Stream& operator()(Stream&& os) const
      {
        write(std::forward<Stream>(os), value);
        return os;
      }

      /// Prints the specified value into the specified stream
      template <typename Stream>
      void write(Stream&& os, real_t v) const
      {
        if (active)
          write_hexfloat(std::forward<Stream>(os), v);
        else
          write_standard(std::forward<Stream>(os), v);
      }

      /// Prints the specified value into the specified stream
      template <typename Stream>
      static void write_hexfloat(Stream&& os, real_t v)
      {
        constexpr auto buf_size = 8 * sizeof(real_t) + 1;
        char buf[buf_size];
        std::snprintf(buf, buf_size, "%+24.14a", v);
        os << buf;
      }

      /// Prints the specified value into the specified stream
      template <typename Stream>
      static void write_standard(Stream&& os, real_t v)
      {
        os << v;
      }

    private:
      bool active;  ///< whether we are writing in base 16
      real_t value; ///< the value to be printed

    }; // OptionalHexFloatFormatter

    template <typename T>
    std::ostream& operator<<(std::ostream& os, details::OptionalHexFloatFormatter<T> fmt)
    {
      return fmt(os);
    }

  } // namespace details

  /**
    * @brief Helper for formatting floats in base 16
    *
    * Example of use:
    *
    *     lar::OptionalHexFloat hexfloat;
    *     constexpr double value = 0.375;
    *
    *     std::cout << "Hex: " << hexfloat(value) << std::endl;
    *
    *     hexfloat.disable();
    *
    *     std::cout << "Dec: " << hexfloat(value) << std::endl;
    *
    * The first printout is expected to be in base 16, the second one in base
    * 10.
    *
    */
  class OptionalHexFloat {

  public:
    /// Constructor: if start_active is true, it will print floats in base 16
    OptionalHexFloat(bool start_active = true) : active(start_active) {}

    /// Returns whether base 16 printing is enabled
    bool enabled() const { return active; }

    /// Enables base 16 printing (or disables it if enable is false)
    void enable(bool enable = true) { active = enable; }

    /// Disables base 16 printing
    void disable() { active = false; }

    /// Returns an object that knows what to do with an output stream
    template <typename T>
    auto operator()(T value) const
    {
      return formatter_t<T>(value, active);
    }

    /// Returns an object that knows what to do with an output stream
    template <typename T>
    auto operator()(bool this_active, T value) const
    {
      return formatter_t<T>(this_active, value);
    }

  private:
    template <typename T>
    using formatter_t = details::OptionalHexFloatFormatter<T>;

    bool active; ///< whether we are writing in base 16

  }; // OptionalHexFloat

} // namespace lar

#endif // LARDATA_RECOBASEART_DUMPERS_HEXFLOAT_H
