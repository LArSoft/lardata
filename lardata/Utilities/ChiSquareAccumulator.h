/**
 * @file   lardata/Utilities/ChiSquareAccumulator.h
 * @brief  Computes a simple &chi;&sup2; sum from data and a expectation function.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 26, 2018
 *
 * This is a header-only library.
 */

#ifndef LARDATA_UTILITIES_CHISQUAREACCUMULATOR_H
#define LARDATA_UTILITIES_CHISQUAREACCUMULATOR_H

// C/C++ standard libraries
#include <utility> // std::move(), std::forward()

namespace lar {
  namespace util {

    /**
     * @brief Computes a &chi;&sup2; from expectation function and data points.
     * @tparam F type of the function
     * @tparam T type of data
     *
     * The formula used is the simple
     * @f$ \chi^{2} = \sum_{i} \left(\frac{y_{i} - e(x_{i})}{\sigma_{i}}\right)^{2} @f$
     * with each observed point being @f$ ( x_{i}, y_{i} \pm \sigma_{i} ) @f$
     * and with @f$ e() @f$ the function describing the expectation (e.g. a fit
     * result).
     *
     * The parameter `F` must be usable as a unary functor, that is it must
     * accept a single argument convertible from type `Data_t` (that is `T`),
     * and return a value convertible back to type `Data_t`.
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * double const a =  2.0;
     * double const b = -1.0;
     * auto f = [a,b](double x){ return a + b * x; };
     * lar::util::ChiSquareAccumulator<decltype(f)> chiSquare;
     *
     * chiSquare.add(0.0, 1.0, 0.5); // add ( 0 ; 1.0 +/- 0.5 )
     * chiSquare.add(1.0, 1.0, 0.5); // add ( 1 ; 1.0 +/- 0.5 )
     * chiSquare.add(2.0, 1.0, 0.5); // add ( 2 ; 1.0 +/- 0.5 )
     *
     * double const chi2value = chiSquare();
     * int degreesOfFreedom = int(chiSquare.N()) - 3;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will check three observations against the prediction of `2 - x`,
     * returning a `chi2value` of `8.0` and a `degreesOfFreedom` of `0`
     * (note that the `3` degrees are manually subtracted).
     */
    template <typename F, typename T = double>
    class ChiSquareAccumulator {
    public:
      using Function_t = F; ///< Type of function for the expectation.
      using Data_t = T;     ///< Type of parameter and observed values.

      //@{
      /**
       * @brief Constructor: uses the specified expectation function.
       * @param expected expectation function
       *
       * The expectation function domain must be a single dimension of type
       * `Data_t`.
       */
      ChiSquareAccumulator(Function_t const& expected) : fExpected(expected) {}
      // @}

      // --- BEGIN -- Access to results ----------------------------------------
      /// @name Access to results
      /// @{

      // @{
      /// Returns the value of &chi;&sup2; currently accumulated.
      Data_t chiSquare() const { return fChiSq; }
      Data_t operator()() const { return chiSquare(); }
      operator Data_t() const { return chiSquare(); }
      //@}

      /// Returns the number of added points (it's not degrees of freedom yet!).
      unsigned int N() const { return fN; }

      /// Returns the expected value for the specified parameter.
      Data_t expected(Data_t x) const { return fExpected(x); }

      /// @}
      // --- END -- Access to results ------------------------------------------

      // --- BEGIN -- Data manipulation ----------------------------------------
      /// @name Data manipulation
      /// @{

      /**
       * @brief Adds a data point to the &chi;&sup2;.
       * @param x parameter
       * @param y observed data with the `x` parameter
       *
       * The &chi;&sup2; is increased by @f$ \left(y - e(x)\right)^{2} @f$
       * where _e_ is the expectation function (`expected()`).
       * The observed values are considered to have nominal uncertainty `1`.
       */
      void add(Data_t x, Data_t y)
      {
        fChiSq += sqr(y - expected(x));
        ++fN;
      }

      /**
       * @brief Adds a data point to the &chi;&sup2;.
       * @param x parameter
       * @param y observed data with the `x` parameter
       * @param s uncertainty on the observed data (default: `1.0`)
       *
       * The &chi;&sup2; is increased by @f$ \left(\frac{y - e(x)}{s}\right)^{2} @f$
       * where _e_ is the expectation function (`expected()`).
       */
      void add(Data_t x, Data_t y, Data_t s)
      {
        fChiSq += sqr(z(y, expected(x), s));
        ++fN;
      }

      /// Resets all the counts, starting from no data.
      void clear()
      {
        fChiSq = Data_t{0};
        fN = 0U;
      }

      /// @}
      // --- END -- Data manipulation ------------------------------------------

    private:
      unsigned int fN = 0U; ///< Number of data entries.
      Data_t fChiSq = 0.0;  ///< Accumulated &chi;&sup2; value.

      Function_t fExpected; ///< Function for the expectation.

      /// Normal variable.
      static Data_t z(Data_t x, Data_t mu, Data_t sigma) { return (x - mu) / sigma; }

      /// The usual square function.
      static Data_t sqr(Data_t v) { return v * v; }

    }; // ChiSquareAccumulator<>

    //--------------------------------------------------------------------------
    /**
     * @brief Creates a `ChiSquareAccumulator` object with the specified function.
     * @tparam F type of function (deduced from `e`)
     * @param e expectation function
     * @return a `ChiSquareAccumulator<F>` instance with specified expectation
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto zero = [](double){ return 0.0; }; // expectation function
     * auto chiSquare = lar::util::makeChiSquareAccumulator(zero);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * declare `chiSquare` in a way equivalent to:
     * `lar::util::ChiSquareAccumulator<decltype(zero)> chiSquare(zero)`.
     */
    template <typename F>
    auto makeChiSquareAccumulator(F&& e)
    {
      return lar::util::ChiSquareAccumulator<F>(std::forward<F>(e));
    }

    /**
     * @brief Creates a `ChiSquareAccumulator` object with the specified function.
     * @tparam T type of data (default: `double`)
     * @tparam F type of function (deduced from `e`)
     * @param e expectation function
     * @return a `ChiSquareAccumulator<F,T>` instance with specified expectation
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto zero = [](float){ return 0.0F; }; // expectation function
     * auto chiSquare = lar::util::makeChiSquareAccumulator<float>(zero);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * declare `chiSquare` in a way equivalent to:
     * `lar::util::ChiSquareAccumulator<decltype(zero), float> chiSquare(zero)`.
     */
    template <typename T, typename F>
    auto makeChiSquareAccumulator(F&& e)
    {
      return lar::util::ChiSquareAccumulator<F, T>(std::forward<F>(e));
    }

    //--------------------------------------------------------------------------

  } // namespace util
} // namespace lar

#endif // LARDATA_UTILITIES_CHISQUAREACCUMULATOR_H
