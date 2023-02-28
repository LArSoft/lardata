/**
 * @file    SimpleFits.h
 * @brief   Classes performing simple fits
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 31st, 2015
 *
 * Currently includes:
 *  - LinearFit
 *  -
 *
 */

#ifndef SIMPLEFITS_H
#define SIMPLEFITS_H 1

// C/C++ standard libraries
#include <algorithm> // std::for_each()
#include <array>
#include <cmath>     // std::sqrt()
#include <iterator>  // std::begin(), std::end()
#include <ostream>   // std::endl
#include <stdexcept> // std::range_error
#include <tuple>
#include <type_traits> // std::enable_if<>, std::is_const<>

#include "lardata/Utilities/FastMatrixMathHelper.h" // lar::util::details::FastMatrixOperations
#include "lardataalg/Utilities/StatCollector.h"     // lar::util::identity

namespace lar {
  namespace util {

    namespace details {

      /// Class providing data collection for the simple polynomial fitters
      template <typename T, unsigned int D>
      class FitDataCollector {

        /**
         * @details
         * The "extractor" is just a C++ trick to replace a statement like:
         *
         *     (N == 0)? weight.Weights(): sums.template SumN<N>()
         *
         * where N is a constexpr number. The problem is that SumN<0> is not
         * valid (in our case, it raises a static assertion). The analysis of
         * the full statement shows that it would not be a problem, since the
         * compiler could remove the ternary operator since it always knows
         * whether N is 0 or not. Matter of fact, the compiler can also choose
         * to be dumb and compile the whole thing, including a SumN<0>.
         *
         * The trick is in the form of a class that gets specialized for N = 0.
         * A long list of limitations in C++ (up to 2014 at least) forces the
         * implementation in a struct.
         */
        template <unsigned int Power, unsigned int N>
        struct SumExtractor {
          static T Sum(WeightTracker<T> const&, DataTracker<Power, T> const& sums)
          {
            return sums.template SumN<N>();
          }
        }; // SumExtractor

        // specialization
        template <unsigned int Power>
        struct SumExtractor<Power, 0U> {
          static T Sum(WeightTracker<T> const& weight, DataTracker<Power, T> const&)
          {
            return weight.Weights();
          }
        }; // SumExtractor<>

      public:
        /// Degree of the fit
        static constexpr unsigned int Degree = D;
        static constexpr unsigned int NParams = Degree + 1;

        using Data_t = T; ///< type of the data

        /// type of measurement without uncertainty
        using Measurement_t = std::tuple<Data_t, Data_t>;

        /// type of measurement with uncertainty
        using MeasurementAndUncertainty_t = std::tuple<Data_t, Data_t, Data_t>;

        // default constructor, destructor and all the rest

        /// @{
        /// @name Add elements

        /**
         * @brief Adds one entry with specified x, y and uncertainty
         * @param x value of x
         * @param y value of y
         * @param sy value of uncertainty on y (1 by default)
         * @return whether the point was added
         *
         * If the uncertainty is exactly 0, the entry is ignored and not added.
         */
        bool add(Data_t x, Data_t y, Data_t sy = Data_t(1.0));

        /**
         * @brief Adds one entry with specified x, y and uncertainty
         * @param value the ( x ; y ) pair
         * @param sy value of uncertainty on y (1 by default)
         * @return whether the point was added
         *
         * If the uncertainty is exactly 0, the entry is ignored and not added.
         */
        bool add(Measurement_t value, Data_t sy = Data_t(1.0))
        {
          return add(std::get<0>(value), std::get<1>(value), sy);
        }

        /**
         * @brief Adds one entry with specified x, y and uncertainty
         * @param value ( x ; y ; sy ), sy being the uncertainty on y
         * @return whether the point was added
         *
         * If the uncertainty is exactly 0, the entry is ignored and not added.
         */
        bool add(MeasurementAndUncertainty_t value)
        {
          return add(std::get<0>(value), std::get<1>(value), std::get<2>(value));
        }

        /**
         * @brief Adds measurements from a sequence, with no uncertainty
         * @tparam Iter forward iterator to the pairs to be added
         * @param begin iterator pointing to the first element to be added
         * @param end iterator pointing after the last element to be added
         *
         * The value pointed by the iterator must be a tuple with types compatible
         * with Measurement_t.
         */
        template <typename Iter>
        void add_without_uncertainty(Iter begin, Iter end)
        {
          add_without_uncertainty(begin, end, identity());
        }

        /**
         * @brief Adds measurements from a sequence with no uncertainty
         * @tparam Iter forward iterator to the elements to be added
         * @tparam Pred a predicate to extract the element from iterator value
         * @param begin iterator pointing to the first element to be added
         * @param end iterator pointing after the last element to be added
         * @param extractor the predicate extracting the value to be inserted
         *
         * The predicate is required to react to a call like with:
         *
         *     Measurement_t Pred::operator() (typename Iter::value_type);
         *
         */
        template <typename Iter, typename Pred>
        void add_without_uncertainty(Iter begin, Iter end, Pred extractor);

        /**
         * @brief Adds all measurements from a container, with no uncertainty
         * @tparam Cont type of container of the elements to be added
         * @tparam Pred a predicate to extract the element from iterator value
         * @param cont container of the elements to be added
         * @param extractor the predicate extracting the value to be inserted
         *
         * The predicate is required to react to a call like with:
         *
         *     Measurement_t Pred::operator() (typename Cont::value_type);
         *
         * The container must support the range-based for loop syntax, that is
         * is must have std::begin<Cont>() and std::end<Cont>() defined.
         */
        template <typename Cont, typename Pred>
        void add_without_uncertainty(Cont cont, Pred extractor)
        {
          add_without_uncertainty(std::begin(cont), std::end(cont), extractor);
        }

        /**
         * @brief Adds all measurements from a container, with no uncertainty
         * @tparam Cont type of container of the elements to be added
         * @param cont container of the elements to be added
         *
         * The container must support the range-based for loop syntax, that is
         * is must have std::begin<Cont>() and std::end<Cont>() defined.
         * The value in the container must be convertible to the Measurement_t
         * type.
         */
        template <typename Cont>
        void add_without_uncertainty(Cont cont)
        {
          add_without_uncertainty(std::begin(cont), std::end(cont));
        }

        /**
         * @brief Adds measurements with uncertainties from a sequence
         * @tparam VIter forward iterator to the elements to be added
         * @tparam UIter forward iterator to the uncertainties
         * @tparam VPred predicate to extract the measurement from iterator value
         * @tparam UPred predicate to extract the uncertainty from iterator value
         * @param begin_value iterator to the first measurement to be added
         * @param end_value iterator after the last measurement to be added
         * @param begin_uncertainty iterator to the uncertainty of first measurement
         * @param value_extractor predicate extracting the measurement to be inserted
         * @param uncertainty_extractor predicate extracting the uncertainty
         * @return number of points added
         *
         * Each element is added with the uncertainty pointed by the matching
         * element in the list pointed by begin_weight: the @f$ y @f$ measurement
         * of `*(begin_value)` will have uncertainty `*(begin_uncertainty)`, the
         * next measurement `*(begin_value + 1)` will have uncertainty
         * `*(begin_uncertainty + 1)`, etc.
         *
         * The predicates are required to react to a call like with:
         *
         *     Measurement_t VPred::operator() (typename VIter::value_type);
         *     Data_t UPred::operator() (typename UIter::value_type);
         *
         * Points with zero or infinite uncertainty are ignored and not added.
         */
        template <typename VIter, typename UIter, typename VPred, typename UPred = identity>
        unsigned int add_with_uncertainty(VIter begin_value,
                                          VIter end_value,
                                          UIter begin_uncertainty,
                                          VPred value_extractor,
                                          UPred uncertainty_extractor = UPred());

        /**
         * @brief Adds measurements with uncertainties from a sequence
         * @tparam Iter forward iterator to MeasurementAndUncertainty_t to be added
         * @param begin iterator pointing to the first measurement to be added
         * @param end iterator pointing after the last measurement to be added
         * @return number of points added
         *
         * The value pointed by the iterator must be a MeasurementAndUncertainty_t
         * or something with elements convertible to its own (Data_t triplet).
         * For more complicate structures, use the version with two predicates
         * (using the uncertainty iterator the same as the measurement iterator).
         *
         * Points with zero or infinite uncertainty are ignored and not added.
         */
        template <typename Iter>
        unsigned int add_with_uncertainty(Iter begin, Iter end);

        /**
         * @brief Adds measurements with uncertainties from a container
         * @tparam Cont type of container of MeasurementAndUncertainty_t elements
         * @param cont container of MeasurementAndUncertainty_t to be added
         * @return number of points added
         *
         * The values in the container must be tuples with each element
         * convertible to Data_t.
         *
         * Points with zero or infinite uncertainty are ignored and not added.
         */
        template <typename Cont>
        unsigned int add_with_uncertainty(Cont cont)
        {
          return add_with_uncertainty(std::begin(cont), std::end(cont));
        }

        ///@}

        /// Clears all the statistics
        void clear();

        /// @{
        /// @name Statistic retrieval

        /// Returns the number of entries added
        int N() const { return s2.N(); }

        /**
         * @brief Returns an average of the uncertainties
         * @return the uncertainty average
         * @throws std::range_error if no entry was added
         *
         * The average is the square root of the harmonic average of the variances
         * (that is, the errors squared):
         * @f$ \bar{s}^{-2} = \frac{1}{N} \sum_{i=1}^{N} s_{y}^{-2} @f$
         */
        Data_t AverageUncertainty() const { return WeightToUncertainty(s2.AverageWeight()); }

        /// Returns the square of the specified value
        template <typename V>
        static constexpr V sqr(V const& v)
        {
          return v * v;
        }

        /// Returns the weighted sum of x^n
        Data_t XN(unsigned int n) const { return (n == 0) ? s2.Weights() : x.Sum(n); }

        /// Returns the weighted sum of x^n y
        Data_t XNY(unsigned int n) const { return (n == 0) ? y.Weights() : xy.Sum(n); }

        /// Returns the weighted sum of x^N
        template <unsigned int N>
        Data_t XN() const
        {
          return SumExtractor<decltype(x)::Power, N>::Sum(s2, x);
        }

        /// Returns the weighted sum of x^N y
        template <unsigned int N>
        Data_t XNY() const
        {
          return SumExtractor<decltype(xy)::Power, N>::Sum(y, xy);
        }

        /// Returns the weighted sum of y^2
        Data_t Y2() const { return y2.template SumN<1>(); }

        /// @}

        /// Prints the statistics into a stream
        template <typename Stream>
        void Print(Stream& out) const;

        /// Transforms an uncertainty into a weight (@f$ s^{-2} @f$)
        static Data_t UncertaintyToWeight(Data_t s) { return Data_t(1.) / sqr(s); }

        /// Transforms a weight back to an uncertainty (@f$ s^{-1/2} @f$)
        static Data_t WeightToUncertainty(Data_t w) { return Data_t(1.) / std::sqrt(w); }

      protected:
        WeightTracker<Data_t> s2;          ///< accumulator for uncertainty
        DataTracker<Degree * 2, Data_t> x; ///< accumulator for variable x^k
        WeightTracker<Data_t> y;           ///< accumulator for y
        DataTracker<1, Data_t> y2;         ///< accumulator for y2
        DataTracker<Degree, Data_t> xy;    ///< accumulator for variable xy

      }; // class FitDataCollector<>

      template <typename T, unsigned int D>
      inline std::ostream& operator<<(std::ostream& out, FitDataCollector<T, D> const& stats)
      {
        stats.Print(out);
        return out;
      }

      /// Base class providing data collection for the simple polynomial fitters
      template <typename T, unsigned int D>
      class SimplePolyFitterDataBase {
        using Collector_t = FitDataCollector<T, D>; ///< class storing input

      public:
        /// Degree of the fit
        static constexpr unsigned int Degree = Collector_t::Degree;

        using Data_t = typename Collector_t::Data_t; ///< type of the data

        /// type of measurement without uncertainty
        using Measurement_t = typename Collector_t::Measurement_t;

        /// type of measurement with uncertainty
        using MeasurementAndUncertainty_t = typename Collector_t::MeasurementAndUncertainty_t;

        // default constructor, destructor and all the rest

        /// @{
        /// @name Add elements
        /// @see FitDataCollector

        bool add(Data_t x, Data_t y, Data_t sy = Data_t(1.0)) { return stats.add(x, y, sy); }

        bool add(Measurement_t value, Data_t sy = Data_t(1.0)) { return stats.add(value, sy); }

        bool add(MeasurementAndUncertainty_t value) { return stats.add(value); }

        template <typename Iter>
        void add_without_uncertainty(Iter begin, Iter end)
        {
          stats.add_without_uncertainty(begin, end);
        }

        template <typename Iter, typename Pred>
        void add_without_uncertainty(Iter begin, Iter end, Pred extractor)
        {
          stats.add_without_uncertainty(begin, end, extractor);
        }

        template <typename Cont, typename Pred>
        void add_without_uncertainty(Cont cont, Pred extractor)
        {
          stats.add_without_uncertainty(cont, extractor);
        }

        template <typename Cont>
        void add_without_uncertainty(Cont cont)
        {
          stats.add_without_uncertainty(cont);
        }

        template <typename VIter, typename UIter, typename VPred, typename UPred = identity>
        unsigned int add_with_uncertainty(VIter begin_value,
                                          VIter end_value,
                                          UIter begin_uncertainty,
                                          VPred value_extractor,
                                          UPred uncertainty_extractor = UPred())
        {
          return stats.add_with_uncertainty(
            begin_value, end_value, begin_uncertainty, value_extractor, uncertainty_extractor);
        }

        template <typename Iter>
        unsigned int add_with_uncertainty(Iter begin, Iter end)
        {
          return stats.add_with_uncertainty(begin, end);
        }

        template <typename Cont>
        unsigned int add_with_uncertainty(Cont cont)
        {
          return stats.add_with_uncertainty(cont);
        }

        ///@}

        /// Clears all the statistics
        void clear() { stats.clear(); }

        /// @{
        /// @name Statistic retrieval
        /// @see FitDataCollector

        int N() const { return stats.N(); }

        Data_t AverageUncertainty() const { return stats.AverageUncertainty(); }

        /// @}

        /// Prints the collected statistics into a stream
        template <typename Stream>
        void PrintStats(Stream& out) const
        {
          stats.Print(out);
        }

        /// Returns the square of the specified value
        template <typename V>
        static constexpr V sqr(V const& v)
        {
          return Collector_t::sqr(v);
        }

      protected:
        Collector_t stats; ///< statistics collected from fit data input

        /// Returns the weighted sum of x^n
        Data_t XN(unsigned int n) const { return stats.XN(n); }

        /// Returns the weighted sum of x^n y
        Data_t XNY(unsigned int n) const { return stats.XNY(n); }

      }; // class SimplePolyFitterDataBase<>

      /// Simple fitter abstract interface
      template <typename T, unsigned int N>
      class SimpleFitterInterface {
      public:
        /// Number of parameters in the fit
        static constexpr unsigned int NParams = N;

        using Data_t = T; ///< type of the data

        using MatrixOps = FastMatrixOperations<Data_t, NParams>;

        /// type of set of fit parameters
        using FitParameters_t = std::array<Data_t, NParams>;

        /// type of matrix for covariance (a std::array)
        using FitMatrix_t = typename MatrixOps::Matrix_t;

        /// Virtual destructor: compiler's default
        virtual ~SimpleFitterInterface() = default;

        /**
         * @brief Returns if the fit has valid results
         * @return if the fit has valid results
         *
         * The fit has no valid results if:
         * 1. insufficient data has been add()ed (no more than the fit Degree)
         * 2. if input points are vertically aligned
         *
         * Note that checking point 2 is expensive in terms of time.
         */
        virtual bool isValid() const = 0;

        /// @{
        /// @name Fitting

        /**
         * @brief Computes and returns all the parameters of the fit result
         * @return the full set of parameters of the fit
         * @throws std::runtime_error if there is no unique solution
         */
        virtual FitParameters_t FitParameters() const = 0;

        /**
         * @brief Computes and returns all the parameter errors of the fit result
         * @return the full set of parameter errors of the fit
         * @throws std::runtime_error if there is no unique solution
         */
        virtual FitParameters_t FitParameterErrors() const = 0;

        /**
         * @brief Computes and returns all the covariance matrix of the fit result
         * @return the the covariance matrix of the fit
         * @throws std::runtime_error if there is no unique solution
         *
         * The matrix is symmetric, and stored in a linear way (say, first row,
         * then second row).
         */
        virtual FitMatrix_t FitParameterCovariance() const = 0;

        /**
         * @brief Returns the parameter n of the fit result
         * @param n degree of the parameter; must be no larger than Degree
         * @return the parameter of the fit, in y/x^n units
         * @throws std::runtime_error if there is no unique solution
         */
        virtual Data_t FitParameter(unsigned int n) const { return FitParameters()[n]; }

        /**
         * @brief Returns the error on parameter n of the fit result
         * @param n degree of the parameter; must be no larger than Degree
         * @return the error on the parameter of the fit, in y/x^n units
         * @throws std::runtime_error if there is no unique solution
         */
        virtual Data_t FitParameterError(unsigned int n) const
        {
          return std::sqrt(FitParameterCovariance()[n * (NParams + 1)]);
        }

        /**
         * @brief Returns the @f$ \chi^{2} @f$ of the fit
         * @return the @f$ \chi^{2} @f$ of the fit (not divided by NDF())
         */
        virtual Data_t ChiSquare() const = 0;

        /**
         * @brief Returns the degrees of freedom in the determination of the fit
         * @return the degrees of freedom in the determination of the fit
         *
         * The return value may be 0 or negative if insufficient points have been
         * added.
         */
        virtual int NDF() const = 0;

        /// @}

        /**
         * @brief Fills the specified parameters
         * @param params the fitted values of the parameters
         * @param Xmat the matrix of the x^n/s^2 sums
         * @param Smat the covariance matrix
         * @param det the determinant of Xmat
         * @return true if the fit is valid (i.e. if a unique solution exists)
         */
        virtual bool FillResults(FitParameters_t& params,
                                 FitMatrix_t& Xmat,
                                 Data_t& det,
                                 FitMatrix_t& Smat) const = 0;

        /**
         * @brief Fills the specified parameters
         * @param params the fitted values of the parameters
         * @param paramerrors the uncertainty on the fitted parameters
         * @param Xmat the matrix of the x^n/s^2 sums
         * @param Smat the covariance matrix
         * @param det the determinant of Xmat
         * @return true if the fit is valid (i.e. if a unique solution exists)
         */
        virtual bool FillResults(FitParameters_t& params,
                                 FitParameters_t& paramerrors,
                                 FitMatrix_t& Xmat,
                                 Data_t& det,
                                 FitMatrix_t& Smat) const = 0;

        /**
         * @brief Fills the specified parameters
         * @param params the fitted values of the parameters
         * @param paramerrors the uncertainty on the fitted parameters
         * @return true if the fit is valid (i.e. if a unique solution exists)
         */
        virtual bool FillResults(FitParameters_t& params, FitParameters_t& paramerrors) const = 0;

        /**
         * @brief Evaluates the fitted function at the specified point
         * @param x the point where to evaluate the fit function
         * @return the value of the fit function
         *
         * No check is performed whether the fit is valid.
         */
        virtual Data_t Evaluate(Data_t x) const = 0;

        /// Evaluates the fitted function; alias of Evaluate()
        Data_t operator()(Data_t x) const { return Evaluate(x); }

        /// Returns the square of the specified data value
        static constexpr Data_t sqr(Data_t v) { return v * v; }

        /// Returns the cube of the specified data value
        static constexpr Data_t cube(Data_t v) { return v * v * v; }

      protected:
        /// Computes the determinant of a matrix
        virtual Data_t Determinant(FitMatrix_t const& mat) const
        {
          return MatrixOps::Determinant(mat);
        }

        /// Computes the inverse of a matrix (using provided determinant)
        virtual FitMatrix_t InvertMatrix(FitMatrix_t const& mat, Data_t det) const
        {
          return MatrixOps::InvertSymmetricMatrix(mat, det);
        }

        /// Computes the inverse of a matrix
        virtual FitMatrix_t InvertMatrix(FitMatrix_t const& mat) const
        {
          return MatrixOps::InvertSymmetricMatrix(mat);
        }

        /// Computes the product of a FitMatrix_t and a FitParameters_t
        virtual FitParameters_t MatrixProduct(FitMatrix_t const& mat,
                                              FitParameters_t const& vec) const
        {
          return MatrixOps::MatrixVectorProduct(mat, vec);
        }

      }; // class SimpleFitterInterface<>

      /// Base class providing virtual fitting interface for polynomial fitters
      template <typename T, unsigned int D>
      class SimplePolyFitterBase : public SimpleFitterInterface<T, D + 1>,
                                   public SimplePolyFitterDataBase<T, D> {
        using Interface_t = SimpleFitterInterface<T, D + 1>; ///< interface
        using Base_t = SimplePolyFitterDataBase<T, D>;       ///< class storing input

      public:
        using Base_t::sqr;

        /// Degree of the fit
        static constexpr unsigned int Degree = Base_t::Degree;

        /// Number of parameters in the fit
        static constexpr unsigned int NParams = Interface_t::NParams;

        using Data_t = typename Base_t::Data_t; ///< type of the data

        /// type of set of fit parameters
        using FitParameters_t = typename Interface_t::FitParameters_t;

        /// type of matrix for covariance (a std::array)
        using FitMatrix_t = typename Interface_t::FitMatrix_t;

        /**
         * @brief Returns if the fit has valid results
         * @return if the fit has valid results
         *
         * The fit has no valid results if:
         * 1. insufficient data has been add()ed (no more than the fit Degree)
         * 2. if input points are vertically aligned
         *
         * Note that checking point 2 is expensive in terms of time.
         */
        virtual bool isValid() const override;

        /// @{
        /// @name Fitting

        /**
         * @brief Computes and returns all the parameters of the fit result
         * @return the full set of parameters of the fit
         * @throws std::range_error if there is no unique solution
         */
        virtual FitParameters_t FitParameters() const override;

        /**
         * @brief Computes and returns all the parameter errors of the fit result
         * @return the full set of parameter errors of the fit
         * @throws std::range_error if there is no unique solution
         */
        virtual FitParameters_t FitParameterErrors() const override;

        /**
         * @brief Computes and returns all the covariance matrix of the fit result
         * @return the the covariance matrix of the fit
         * @throws std::range_error if there is no unique solution
         *
         * The matrix is symmetric, and stored in a linear way (say, first row,
         * then second row).
         */
        virtual FitMatrix_t FitParameterCovariance() const override;

        /**
         * @brief Returns the parameter n of the fit result
         * @param n degree of the parameter; must be no larger than Degree
         * @return the parameter of the fit, in y/x^n units
         * @throws std::range_error if there is no unique solution
         */
        virtual Data_t FitParameter(unsigned int n) const override;

        /**
         * @brief Returns the error on parameter n of the fit result
         * @param n degree of the parameter; must be no larger than Degree
         * @return the error on the parameter of the fit, in y/x^n units
         * @throws std::range_error if there is no unique solution
         */
        virtual Data_t FitParameterError(unsigned int n) const override;

        /**
         * @brief Returns the @f$ \chi^{2} @f$ of the fit
         * @return the @f$ \chi^{2} @f$ of the fit (not divided by NDF())
         */
        virtual Data_t ChiSquare() const override;

        /**
         * @brief Returns the degrees of freedom in the determination of the fit
         * @return the degrees of freedom in the determination of the fit
         *
         * The return value may be 0 or negative if insufficient points have been
         * added.
         */
        virtual int NDF() const override { return Base_t::N() - NParams; }

        /// @}

        /**
         * @brief Fills the specified parameters
         * @param params the fitted values of the parameters
         * @param Xmat the matrix of the x^n/s^2 sums
         * @param Smat the covariance matrix
         * @param det the determinant of Xmat
         * @return true if the fit is valid (i.e. if a unique solution exists)
         */
        bool FillResults(FitParameters_t& params,
                         FitMatrix_t& Xmat,
                         Data_t& det,
                         FitMatrix_t& Smat) const override;

        /**
         * @brief Fills the specified parameters
         * @param params the fitted values of the parameters
         * @param paramerrors the uncertainty on the fitted parameters
         * @param Xmat the matrix of the x^n/s^2 sums
         * @param Smat the covariance matrix
         * @param det the determinant of Xmat
         * @return true if the fit is valid (i.e. if a unique solution exists)
         */
        bool FillResults(FitParameters_t& params,
                         FitParameters_t& paramerrors,
                         FitMatrix_t& Xmat,
                         Data_t& det,
                         FitMatrix_t& Smat) const override;

        /**
         * @brief Fills the specified parameters
         * @param params the fitted values of the parameters
         * @param paramerrors the uncertainty on the fitted parameters
         * @return true if the fit is valid (i.e. if a unique solution exists)
         */
        bool FillResults(FitParameters_t& params, FitParameters_t& paramerrors) const override;

        /**
         * @brief Evaluates the fitted function at the specified point
         * @param x the point where to evaluate the fit function
         * @return the value of the fit function
         *
         * No check is performed whether the fit is valid.
         */
        virtual Data_t Evaluate(Data_t x) const override;

        /// Extracts parameter errors from diagonal of the covarriance matrix
        static FitParameters_t ExtractParameterErrors(FitMatrix_t const& Smat);

      protected:
        /// Fills and returns the matrix of x^n sum coefficients ( { x^(i+j) } )
        virtual FitMatrix_t MakeMatrixX() const;

        /// Fills and returns the matrix (vector) of x^n y sum coefficients
        virtual FitParameters_t MakeMatrixY() const;

        /// Computes and returns all the parameter errors of the fit result
        virtual FitParameters_t FitParameterErrors(FitMatrix_t const& Smat) const
        {
          return ExtractParameterErrors(Smat);
        }

        //@{
        /// Returns the fitted parameters using the provided information
        virtual FitParameters_t FitParameters(FitMatrix_t const& Xmat) const;

        virtual FitParameters_t FitParameters(FitMatrix_t const& Smat, Data_t /* det */) const;
        //@}

        //@{
        /// Computes a single fit parameter using the given information
        virtual Data_t Param(unsigned int n, FitMatrix_t const& Xmat) const;
        virtual Data_t Param(unsigned int n, FitMatrix_t const& Xmat, Data_t detXmat) const;
        //@}

        // import the declaration from the interface
        using Interface_t::Determinant;
        using Interface_t::InvertMatrix;
        using Interface_t::MatrixProduct;

      }; // class SimplePolyFitterBase<>

    } // namespace details

    /** ************************************************************************
     * @brief Performs a linear regression of data
     * @tparam T type of the quantities
     * @tparam W type of the weight (as T by default)
     *
     * The linear regression connects measurements
     * @f$ ( y_{i} \pm \sigma_{y,i} ) @f$ with a parameter @f$ ( x_{i} ) @f$
     * not affected by uncertainty.
     * The returned parameters describe a straight line @f$ y = a x + b @f$
     * obtained by minimization of
     * @f$ \chi^{2} = \sum_{i} \frac{ \left(y_{i} - a x_{i} - b \right)^{2} }{ \sigma^{2}_{y,i} }@f$
     *
     * This saves having to link to ROOT for the simplest cases.
     *
     * This simple linear fitter does not store any result: each time a result
     * is requested, it is computed anew.
     * In particular that is true also for ChiSquare(), that requires the full
     * parameters set and therefore reruns the full fit (FitParameters())
     * and for the covariance matrix of the parameters.
     */
    template <typename T>
    class LinearFit : public details::SimplePolyFitterBase<T, 1U> {
      using Base_t = details::SimplePolyFitterBase<T, 1U>;

    public:
      using Base_t::Degree;
      using Base_t::NParams;
      using typename Base_t::Data_t;

      /// type of measurement without uncertainty
      using typename Base_t::Measurement_t;

      /// type of measurement with uncertainty
      using typename Base_t::MeasurementAndUncertainty_t;

      using Base_t::sqr;

      /// type of set of fit parameters
      using FitParameters_t = std::array<Data_t, NParams>;
      using FitMatrix_t = std::array<Data_t, sqr(NParams)>;

      // default constructor, destructor and all the rest

      /**
       * @brief Returns the intercept of the fit
       * @return the intercept of the fit, in y units
       * @throws std::range_error if there is no unique solution
       */
      Data_t Intercept() const { return this->FitParameter(0); }

      /**
       * @brief Returns the slope of the fit
       * @return the slope of the fit, in y/x units
       * @throws std::range_error if there is no unique solution
       */
      Data_t Slope() const { return this->FitParameter(1); }

      /**
       * @brief Returns the error on intercept of the fit
       * @return the error on intercept of the fit, in y units
       * @throws std::range_error if there is no unique solution
       */
      Data_t InterceptError() const { return this->FitParameterError(0); }

      /**
       * @brief Returns the error in slope of the fit
       * @return the error on slope of the fit, in y/x units
       * @throws std::range_error if there is no unique solution
       */
      Data_t SlopeError() const { return this->FitParameterError(1); }

      /**
       * @brief Returns the covariance between intercept and slope of the fit
       * @return the covariance between intercept and slope of the fit, in y^2 units
       * @throws std::range_error if there is no unique solution
       */
      Data_t InterceptSlopeCovariance() const
      {
        return this->FitParameterCovariance()[0 * NParams + 1];
      }

      /**
       * @brief Returns the @f$ \chi^{2} @f$ of the fit
       * @return the @f$ \chi^{2} @f$ of the fit (not divided by NDF())
       */
      virtual Data_t ChiSquare() const override;

    protected:
      //@{
      /// Aliases
      Data_t I() const { return Base_t::stats.template XN<0>(); }
      Data_t X() const { return Base_t::stats.template XN<1>(); }
      Data_t X2() const { return Base_t::stats.template XN<2>(); }
      Data_t Y() const { return Base_t::stats.template XNY<0>(); }
      Data_t XY() const { return Base_t::stats.template XNY<1>(); }
      Data_t Y2() const { return Base_t::stats.Y2(); }
      //@}

    }; // class LinearFit<>

    /** ************************************************************************
     * @brief Performs a second-degree fit of data
     * @tparam T type of the quantities
     * @tparam W type of the weight (as T by default)
     *
     * The quadratic fit connects measurements
     * @f$ ( y_{i} \pm \sigma_{y,i} ) @f$ with a parameter @f$ ( x_{i} ) @f$
     * not affected by uncertainty.
     * The returned parameters describe a quadratic curve
     * @f$ f(x) = a_{0} + a_{1} x + a_{2} x^{2} @f$ obtained by minimization of
     * @f$ \chi^{2} = \sum_{i} \frac{ \left(y_{i} - f(x_{i}) \right)^{2} }{ \sigma^{2}_{y,i} }@f$
     *
     * This saves having to link to ROOT for the simplest cases.
     *
     * This simple quadratic fitter does not store any result: each time a
     * result is requested, it is computed anew.
     * In particular that is true also for ChiSquare(), that requires the full
     * parameters set and therefore reruns the full fit (FitParameters())
     * and for the covariance matrix of the parameters.
     */
    template <typename T>
    class QuadraticFit : public details::SimplePolyFitterBase<T, 2U> {
      using Base_t = details::SimplePolyFitterBase<T, 2U>;

    public:
      using Base_t::Degree;
      using Base_t::NParams;
      using typename Base_t::Data_t;

      /// type of measurement without uncertainty
      using typename Base_t::Measurement_t;

      /// type of measurement with uncertainty
      using typename Base_t::MeasurementAndUncertainty_t;

      using Base_t::sqr;

      /// type of set of fit parameters
      using FitParameters_t = std::array<Data_t, NParams>;
      using FitMatrix_t = std::array<Data_t, sqr(NParams)>;

      // default constructor, destructor and all the rest

      /**
       * @brief Returns the @f$ \chi^{2} @f$ of the fit
       * @return the @f$ \chi^{2} @f$ of the fit (not divided by NDF())
       */
      virtual Data_t ChiSquare() const override;

    protected:
      //@{
      /// Aliases
      Data_t I() const { return Base_t::stats.template XN<0>(); }
      Data_t X() const { return Base_t::stats.template XN<1>(); }
      Data_t X2() const { return Base_t::stats.template XN<2>(); }
      Data_t X3() const { return Base_t::stats.template XN<3>(); }
      Data_t X4() const { return Base_t::stats.template XN<4>(); }
      Data_t Y() const { return Base_t::stats.template XNY<0>(); }
      Data_t XY() const { return Base_t::stats.template XNY<1>(); }
      Data_t X2Y() const { return Base_t::stats.template XNY<2>(); }
      Data_t Y2() const { return Base_t::stats.Y2(); }
      //@}

    }; // class QuadraticFit<>

    /** **********************************************************************
     * @brief "Fast" Gaussian fit
     * @tparam T data type
     *
     * This class performs a Gaussian fit on demand.
     * This fit translates the data to its logarithm and then internally
     * performs a quadratic fit.
     * Note that as a consequence this fitter does not accept negative values
     * for the y variable.
     * Negative values in input will be completely ignored.
     *
     * Methods that do not change functionality respect to the base class are
     * not documented here -- see the base class(es) documentation
     * (mostly SimplePolyFitterBase).
     */
    template <typename T>
    class GaussianFit : public details::SimpleFitterInterface<T, 3> {
      using Base_t = details::SimpleFitterInterface<T, 3>; ///< base class
      using Fitter_t = QuadraticFit<T>;                    ///< base class

    public:
      /// Number of parameters in the fit
      static constexpr unsigned int NParams = Base_t::NParams;

      using Data_t = typename Base_t::Data_t; ///< type of the data

      /// type of measurement without uncertainty
      using Measurement_t = typename Fitter_t::Measurement_t;

      /// type of measurement with uncertainty
      using MeasurementAndUncertainty_t = typename Fitter_t::MeasurementAndUncertainty_t;

      using FitParameters_t = typename Fitter_t::FitParameters_t;
      using FitMatrix_t = typename Fitter_t::FitMatrix_t;

      using Base_t::cube;
      using Base_t::sqr;

      // default constructor, destructor and all the rest

      /// @{
      /// @name Add elements
      /// @see FitDataCollector

      bool add(Data_t x, Data_t y, Data_t sy = Data_t(1.0));

      bool add(Measurement_t value, Data_t sy = Data_t(1.0))
      {
        return add(std::get<0>(value), std::get<1>(value), sy);
      }

      bool add(MeasurementAndUncertainty_t value)
      {
        return add(std::get<0>(value), std::get<1>(value), std::get<2>(value));
      }

      template <typename Iter>
      void add_without_uncertainty(Iter begin, Iter end)
      {
        add_without_uncertainty(begin, end, identity());
      }

      template <typename Iter, typename Pred>
      void add_without_uncertainty(Iter begin, Iter end, Pred extractor);

      template <typename Cont, typename Pred>
      void add_without_uncertainty(Cont cont, Pred extractor)
      {
        add_without_uncertainty(std::begin(cont), std::end(cont), extractor);
      }

      template <typename Cont>
      void add_without_uncertainty(Cont cont)
      {
        add_without_uncertainty(std::begin(cont), std::end(cont));
      }

      template <typename VIter, typename UIter, typename VPred, typename UPred = identity>
      unsigned int add_with_uncertainty(VIter begin_value,
                                        VIter end_value,
                                        UIter begin_uncertainty,
                                        VPred value_extractor,
                                        UPred uncertainty_extractor = UPred());

      template <typename Iter>
      unsigned int add_with_uncertainty(Iter begin, Iter end);

      template <typename Cont>
      unsigned int add_with_uncertainty(Cont cont)
      {
        return add_with_uncertainty(std::begin(cont), std::end(cont));
      }

      /// Clears all the input statistics
      void clear() { fitter.clear(); }

      /// Returns the number of (valid) points added
      int N() const { return fitter.N(); }

      /// Prints the collected statistics into a stream
      template <typename Stream>
      void PrintStats(Stream& out) const
      {
        fitter.PrintStats(out);
      }

      ///@}

      /// @{
      /// @name Fitting
      /**
       * @brief Returns if the fit has valid results
       * @return if the fit has valid results
       *
       * The fit has no valid results if:
       * 1. insufficient data has been add()ed (no more than the fit Degree)
       * 2. if input points are vertically aligned
       *
       * Note that checking point 2 is expensive in terms of time.
       */
      virtual bool isValid() const override { return fitter.isValid(); }

      /// @{
      /// @name Fitting

      /**
       * @brief Computes and returns all the parameters of the fit result
       * @return the full set of parameters of the fit
       * @throws std::runtime_error if there is no unique solution
       */
      virtual FitParameters_t FitParameters() const override;

      /**
       * @brief Computes and returns all the parameter errors of the fit result
       * @return the full set of parameter errors of the fit
       * @throws std::runtime_error if there is no unique solution
       */
      virtual FitParameters_t FitParameterErrors() const override;

      /**
       * @brief Computes and returns all the covariance matrix of the fit result
       * @return the the covariance matrix of the fit
       * @throws std::runtime_error if there is no unique solution
       *
       * Not supported.
       * It's fairly too complicate to fill the whole matrix.
       * Doable, on request.
       */
      virtual FitMatrix_t FitParameterCovariance() const override;

      /**
       * @brief Returns the @f$ \chi^{2} @f$ of the original fit
       * @return the @f$ \chi^{2} @f$ of the original fit (not divided by NDF())
       *
       * This is not defined in the space of the Gaussian, but in the space of
       * the internal quadratic fit. Where one is minimum, the other also is,
       * but the actual value is different.
       */
      virtual Data_t ChiSquare() const override { return fitter.ChiSquare(); }

      /**
       * @brief Returns the degrees of freedom in the determination of the fit
       * @return the degrees of freedom in the determination of the fit
       *
       * The return value may be 0 or negative if insufficient points have been
       * added.
       */
      virtual int NDF() const override { return fitter.NDF(); }

      /// @}

      /**
       * @brief Fills the specified parameters
       * @param params the fitted values of the parameters
       * @param Xmat the matrix of the x^n/s^2 sums
       * @param Smat the covariance matrix
       * @param det the determinant of Xmat
       * @return true if the fit is valid (i.e. if a unique solution exists)
       *
       * Unsupported.
       */

      virtual bool FillResults(FitParameters_t& params,
                               FitMatrix_t& Xmat,
                               Data_t& det,
                               FitMatrix_t& Smat) const override;
      /**
       * @brief Fills the specified parameters
       * @param params the fitted values of the parameters
       * @param paramerrors the uncertainty on the fitted parameters
       * @param Xmat the matrix of the x^n/s^2 sums
       * @param Smat the covariance matrix
       * @param det the determinant of Xmat
       * @return true if the fit is valid (i.e. if a unique solution exists)
       *
       * Unsupported.
       */
      virtual bool FillResults(FitParameters_t& params,
                               FitParameters_t& paramerrors,
                               FitMatrix_t& Xmat,
                               Data_t& det,
                               FitMatrix_t& Smat) const override;

      /**
       * @brief Fills the specified parameters
       * @param params the fitted values of the parameters
       * @param paramerrors the uncertainty on the fitted parameters
       * @return true if the fit is valid (i.e. if a unique solution exists)
       *
       * Only the version returning the parameters and errors is supported.
       */
      virtual bool FillResults(FitParameters_t& params,
                               FitParameters_t& paramerrors) const override;

      /**
       * @brief Evaluates the fitted function at the specified point
       * @param x the point where to evaluate the fit function
       * @return the value of the fit function
       *
       * No check is performed whether the fit is valid.
       */
      virtual Data_t Evaluate(Data_t x) const override
      {
        return Evaluate(x, FitParameters().data());
      }

      /// Returns the internal fitter (mostly for debugging)
      virtual Fitter_t const& Fitter() const { return fitter; }

      /**
       * @brief Evaluates a Gaussian with given parameters at one point
       * @param x the point where to evaluate the fit function
       * @param params Gaussian parameters: amplitude, mean, sigma
       * @return the Gaussian function evaluated at x
       */
      static Data_t Evaluate(Data_t x, Data_t const* params);

    protected:
      Fitter_t fitter; ///< the actual fitter and data holder

      /// @{
      /// @name Mumbo-jumbo to convert the values for a quadratic fit

      ///< type of value and error
      struct Value_t : public std::tuple<Data_t, Data_t> {
        using Base_t = std::tuple<Data_t, Data_t>;

        Value_t(Data_t v, Data_t e) : Base_t(v, e) {}
        Value_t(MeasurementAndUncertainty_t meas) : Base_t(std::get<1>(meas), std::get<2>(meas)) {}

        constexpr Data_t value() const { return std::get<0>(*this); }
        constexpr Data_t error() const { return std::get<1>(*this); }

        Data_t& value() { return std::get<0>(*this); }
        Data_t& error() { return std::get<1>(*this); }

      }; // Value_t

      /// Converts a value into a proper input for the quadratic fit;
      /// does not accept 0 or negative values!
      static Data_t EncodeValue(Data_t value) { return std::log(value); }

      /// Converts a value from the quadratic fit into a proper value
      static Data_t DecodeValue(Data_t value) { return std::exp(value); }

      /// Converts a value and error into a proper input for the quadratic fit
      static Value_t EncodeValue(Data_t value, Data_t error)
      {
        return {std::log(value), error / std::abs(value)};
      }

      /// Converts a value and error into a proper input for the quadratic fit
      static Value_t EncodeValue(Value_t const& value)
      {
        return EncodeValue(value.value(), value.error());
      }

      /// Converts a value from the quadratic fit into a proper value
      static Value_t DecodeValue(Value_t const& value)
      {
        const Data_t v = std::exp(value.value());
        return {v, v * value.error()};
      } // DecodeValue()

      /// Converts a value and error into a proper input for the quadratic fit
      static Measurement_t EncodeValue(Measurement_t const& meas)
      {
        return Measurement_t(std::get<0>(meas), EncodeValue(std::get<1>(meas)));
      }

      /// Converts a value and error into a proper input for the quadratic fit
      static MeasurementAndUncertainty_t EncodeValue(MeasurementAndUncertainty_t const& meas)
      {
        Value_t value = EncodeValue(Value_t(meas));
        return {std::get<0>(meas), value.value(), value.error()};
      } // EncodeValue(MeasurementAndUncertainty_t)

      /// Converts a value and error into a proper input for the quadratic fit
      static MeasurementAndUncertainty_t EncodeValue(Measurement_t const& meas, Data_t error)
      {
        Value_t value = EncodeValue(Value_t(std::get<1>(meas), error));
        return {std::get<0>(meas), value.value(), value.error()};
      } // EncodeValue(Measurement_t, Data_t)

      /// Wrapper to encode a MeasurementAndUncertainty_t from a value
      /// and a error extractor
      template <typename VPred, typename UPred = void>
      struct EncodeExtractor {
        EncodeExtractor(VPred& vpred, UPred& upred) : value_extractor(vpred), error_extractor(upred)
        {}

        // extractor takes whatever dereferencing Iter gives and returns
        // Measurement_t or MeasurementAndUncertainty_t,
        // the one the extractor returns
        template <typename Elem>
        auto operator()(Elem elem)
        {
          // use explicit casts to make sure we know what we are doing
          return EncodeValue(static_cast<Measurement_t&&>(value_extractor(elem)),
                             static_cast<Data_t&&>(error_extractor(elem)));
        } // operator()

        VPred& value_extractor; ///< value extractor
        UPred& error_extractor; ///< uncertainty extractor
      };                        // struct EncodeExtractor

      /// Wrapper to encode a Measurement_t or MeasurementAndUncertainty_t
      /// from a extractor
      template <typename Pred>
      struct EncodeExtractor<Pred, void> {
        EncodeExtractor(Pred& pred) : extractor(pred) {}

        // extractor takes whatever dereferencing Iter gives and returns
        // Measurement_t or MeasurementAndUncertainty_t,
        // the one the extractor returns;
        // as usual with enable_if, I am not sure it makes sense
        template <typename Elem, typename = std::enable_if<std::is_const<Pred>::value, void>>
        auto operator()(Elem elem) const
        {
          return EncodeValue(extractor(elem));
        }

        template <typename Elem, typename = std::enable_if<!std::is_const<Pred>::value, void>>
        auto operator()(Elem elem)
        {
          return EncodeValue(extractor(elem));
        }

        Pred& extractor;
      }; // struct EncodeExtractor<Pred>

      template <typename Pred>
      static EncodeExtractor<Pred> Encoder(Pred& pred)
      {
        return {pred};
      }

      template <typename VPred, typename UPred>
      static EncodeExtractor<VPred, UPred> Encoder(VPred& vpred, UPred& upred)
      {
        return {vpred, upred};
      }

      /// @}

      /**
       * @brief Converts the specified quadratic fit parameters into Gaussian
       * @param qpars the quadratic fit parameters
       * @return Gaussian function parameters
       */
      static FitParameters_t ConvertParameters(FitParameters_t const& qpars);

      /**
       * @brief Converts the specified quadratic fit parameters and errors
       * @param qpars the quadratic fit parameters
       * @param qparerrmat the quadratic fit parameter error matrix
       * @param params the Gaussian fit parameters
       * @param paramerrors the Gaussian fit parameter errors
       */
      static void ConvertParametersAndErrors(FitParameters_t const& qpars,
                                             FitMatrix_t const& qparerrmat,
                                             FitParameters_t& params,
                                             FitParameters_t& paramerrors);

      /**
       * @brief Converts the specified quadratic fit parameters and errors
       * @param qpars the quadratic fit parameters
       * @param qparerrmat the quadratic fit parameter error matrix
       * @param params the Gaussian fit parameters
       * @param paramvariances the Gaussian fit parameter variance
       */
      static void ConvertParametersAndVariances(FitParameters_t const& qpars,
                                                FitMatrix_t const& qparerrmat,
                                                FitParameters_t& params,
                                                FitParameters_t& paramvariances);

      /**
       * @brief Converts the specified quadratic fit parameters and errors
       * @param qpars the quadratic fit parameters
       * @param qparerrmat the quadratic fit parameter error matrix
       * @param params the Gaussian fit parameters
       * @param Smat the covariance matrix of the Gaussian fit parameters
       */
      static void ConvertParametersAndErrorMatrix(FitParameters_t const& qpars,
                                                  FitMatrix_t const& qparerrmat,
                                                  FitParameters_t& params,
                                                  FitMatrix_t& Smat);

      /**
       * @brief Returns whether the specified parameters represent a valid fit
       * @param params Gaussian parameters
       * @param qpars quadratic fit parameters
       * @return whether specified parameters represent a valid Gaussian fit
       */
      static bool isValid(FitParameters_t const& params, FitParameters_t const& qpars);

      static void ThrowNotImplemented [[noreturn]] (std::string method)
      {
        throw std::logic_error("Method " + method + "() not implemented");
      }

    }; // class GaussianFit<>

  } // namespace util
} // namespace lar

//==============================================================================
//=== template implementation
//==============================================================================

//******************************************************************************
//***  FitDataCollector<>
//***

template <typename T, unsigned int D>
bool lar::util::details::FitDataCollector<T, D>::add(Data_t x_value,
                                                     Data_t y_value,
                                                     Data_t sy /* = Data_t(1.0) */)
{
  Data_t w = UncertaintyToWeight(sy);
  if (!std::isnormal(w)) return false;
  // the x section has a 1/s^2 weight; we track that weight separately
  s2.add(w);
  x.add(x_value, w);
  // we treat the y section as if it were a x section with a y/s^2 weight;
  // we track that weight separately
  Data_t yw = y_value * w;
  y.add(yw);
  y2.add(sqr(y_value), w); // used only for chi^2
  xy.add(x_value, yw);

  return true; // we did add the value
} // FitDataCollector<>::add()

template <typename T, unsigned int D>
template <typename Iter, typename Pred>
void lar::util::details::FitDataCollector<T, D>::add_without_uncertainty(Iter begin,
                                                                         Iter end,
                                                                         Pred extractor)
{
  std::for_each(begin, end, [this, extractor](auto item) { this->add(extractor(item)); });
} // FitDataCollector<>::add_without_uncertainty(Iter, Pred)

template <typename T, unsigned int D>
template <typename VIter, typename UIter, typename VPred, typename UPred>
unsigned int lar::util::details::FitDataCollector<T, D>::add_with_uncertainty(
  VIter begin_value,
  VIter end_value,
  UIter begin_uncertainty,
  VPred value_extractor,
  UPred uncertainty_extractor /* = UPred() */
)
{
  unsigned int n = 0;
  while (begin_value != end_value) {
    if (add(value_extractor(*begin_value), uncertainty_extractor(*begin_uncertainty))) ++n;
    ++begin_value;
    ++begin_uncertainty;
  } // while
  return n;
} // FitDataCollector<>::add_with_uncertainty(VIter, VIter, UIter, VPred, UPred)

template <typename T, unsigned int D>
template <typename Iter>
unsigned int lar::util::details::FitDataCollector<T, D>::add_with_uncertainty(Iter begin, Iter end)
{
  unsigned int old_n = N();
  std::for_each(begin, end, [this](auto p) { this->add(p); });
  return N() - old_n;
} // FitDataCollector<>::add_with_uncertainty(Iter, Iter)

template <typename T, unsigned int D>
inline void lar::util::details::FitDataCollector<T, D>::clear()
{
  s2.clear();
  x.clear();
  y.clear();
  y2.clear();
  xy.clear();
} // FitDataCollector<>::clear()

template <typename T, unsigned int D>
template <typename Stream>
void lar::util::details::FitDataCollector<T, D>::Print(Stream& out) const
{

  out << "Sums  1/s^2=" << s2.Weights() << "\n      x/s^2=" << x.template SumN<1>();
  for (unsigned int degree = 2; degree <= x.Power; ++degree)
    out << "\n    x^" << degree << "/s^2=" << x.Sum(degree);
  out << "\n      y/s^2=" << y.Weights() << "\n    y^2/s^2=" << y2.Sum();
  if (xy.Power >= 1) out << "\n     xy/s^2=" << xy.template SumN<1>();
  for (unsigned int degree = 2; degree <= xy.Power; ++degree)
    out << "\n   x^" << degree << "y/s^2=" << xy.Sum(degree);
  out << std::endl;
} // FitDataCollector<>::Print()

//******************************************************************************
//***  SimplePolyFitterBase<>
//***
template <typename T, unsigned int D>
inline bool lar::util::details::SimplePolyFitterBase<T, D>::isValid() const
{
  return (Base_t::N() > (int)Degree) && std::isnormal(Determinant(MakeMatrixX()));
} // SimplePolyFitterBase<>::isValid()

template <typename T, unsigned int D>
inline auto lar::util::details::SimplePolyFitterBase<T, D>::FitParameter(unsigned int n) const
  -> Data_t
{
  return Param(n, MakeMatrixX());
} // SimplePolyFitterBase<>::FitParameter(unsigned int)

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::FitParameterError(unsigned int n) const
  -> Data_t
{
  if (n > Degree) return Data_t(0); // no parameter, no error
  return std::sqrt(FitParameterCovariance()[n * (NParams + 1)]);
} // SimplePolyFitterBase<>::FitParameterError()

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::FitParameters() const -> FitParameters_t
{
  FitMatrix_t Xmat = MakeMatrixX();
  FitParameters_t fit_params;
  for (unsigned int iParam = 0; iParam < NParams; ++iParam)
    fit_params[iParam] = Param(iParam, Xmat);
  return fit_params;
} // SimplePolyFitterBase<>::FitParameters()

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::FitParameterErrors() const -> FitParameters_t
{
  return FitParameterErrors(FitParameterCovariance());
} // SimplePolyFitterBase<>::FitParameterErrors()

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::FitParameterCovariance() const -> FitMatrix_t
{
  FitMatrix_t Xmat = MakeMatrixX();
  Data_t det = Determinant(Xmat);
  if (!std::isnormal(det)) {
    throw std::range_error(
      "SimplePolyFitterBase::FitParameterCovariance(): determinant 0 while fitting");
  }
  return InvertMatrix(Xmat, det);
} // SimplePolyFitterBase<>::FitParameterCovariance()

template <typename T, unsigned int D>
bool lar::util::details::SimplePolyFitterBase<T, D>::FillResults(FitParameters_t& params,
                                                                 FitMatrix_t& Xmat,
                                                                 Data_t& det,
                                                                 FitMatrix_t& Smat) const
{

  Xmat = MakeMatrixX();
  det = Determinant(Xmat);
  if (!std::isnormal(det)) {
    Smat.fill(Data_t(0));
    params.fill(Data_t(0));
    return false;
  }
  Smat = InvertMatrix(Xmat, det);
  params = FitParameters(Smat, det);
  return true;
} // SimplePolyFitterBase<>::FillResults(params, matrices, determinant)

template <typename T, unsigned int D>
bool lar::util::details::SimplePolyFitterBase<T, D>::FillResults(FitParameters_t& params,
                                                                 FitParameters_t& paramerrors,
                                                                 FitMatrix_t& Xmat,
                                                                 Data_t& det,
                                                                 FitMatrix_t& Smat) const
{

  if (!this->FillResults(params, Xmat, det, Smat)) return false;
  paramerrors = ExtractParameterErrors(Smat);
  return true;
} // SimplePolyFitterBase<>::FillResults(params, errors, matrices, determinant)

template <typename T, unsigned int D>
bool lar::util::details::SimplePolyFitterBase<T, D>::FillResults(FitParameters_t& params,
                                                                 FitParameters_t& paramerrors) const
{
  // to compute the parameters, we need all the stuff;
  // we just keep it local and discard it in the end. Such a waste.
  FitMatrix_t Xmat, Smat;
  Data_t det;
  return FillResults(params, paramerrors, Xmat, det, Smat);
} // SimplePolyFitterBase<>::FillResults(params, errors)

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::Evaluate(Data_t x) const -> Data_t
{
  FitParameters_t params = FitParameters();
  unsigned int iParam = NParams - 1; // point to last parameter (highest degree)
  Data_t v = params[iParam];
  while (iParam > 0)
    v = v * x + params[--iParam];
  return v;
} // SimplePolyFitterBase<>::Evaluate()

// --- protected methods follow ---
template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::MakeMatrixX() const -> FitMatrix_t
{
  FitMatrix_t Xmat;
  for (unsigned int i = 0; i < NParams; ++i) {   // row
    for (unsigned int j = i; j < NParams; ++j) { // column
      Xmat[j * NParams + i] = Xmat[i * NParams + j] = Base_t::XN(i + j);
    } // for j
  }   // for i
  return Xmat;
} // SimplePolyFitterBase<>::MakeMatrixX()

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::MakeMatrixY() const -> FitParameters_t
{
  FitParameters_t Ymat;
  for (unsigned int i = 0; i < NParams; ++i)
    Ymat[i] = Base_t::XNY(i);
  return Ymat;
} // SimplePolyFitterBase<>::MakeMatrixY()

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::FitParameters(FitMatrix_t const& Xmat) const
  -> FitParameters_t
{
  FitParameters_t fit_params;
  for (unsigned int iParam = 0; iParam < NParams; ++iParam)
    fit_params[iParam] = Param(iParam, Xmat);
  return fit_params;
} // SimplePolyFitterBase<>::FitParameters(FitMatrix_t)

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::FitParameters(FitMatrix_t const& Smat,
                                                                   Data_t /* det */) const
  -> FitParameters_t
{
  return MatrixProduct(Smat, MakeMatrixY());
} // SimplePolyFitterBase<>::FitParameters(FitMatrix_t)

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::Param(unsigned int n,
                                                           FitMatrix_t const& Xmat) const -> Data_t
{
  if (n > Degree) return Data_t(0); // no such a degree, its coefficient is 0

  Data_t detXmat = Determinant(Xmat);
  if (!std::isnormal(detXmat)) {
    throw std::range_error("SimplePolyFitterBase::Param(): Determinant 0 while fitting");
  }
  return Param(n, Xmat, detXmat);
} // SimplePolyFitterBase<>::Param(unsigned int, FitMatrix_t)

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::ExtractParameterErrors(FitMatrix_t const& Smat)
  -> FitParameters_t
{
  FitParameters_t fit_errors;
  for (unsigned int iParam = 0; iParam <= Degree; ++iParam)
    fit_errors[iParam] = std::sqrt(Smat[iParam * (NParams + 1)]);
  return fit_errors;
} // SimplePolyFitterBase<>::FitParameterErrors(FitMatrix_t)

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::Param(unsigned int n,
                                                           FitMatrix_t const& Xmat,
                                                           Data_t detXmat) const -> Data_t
{
  if (n > Degree) return Data_t(0); // no such a degree, its coefficient is 0
  // XYmat is as Xmat...
  FitMatrix_t XYmat(Xmat);
  // ... except that the N-th column is replaced with { sum x^i y }
  for (unsigned int i = 0; i < NParams; ++i)
    XYmat[i * NParams + n] = Base_t::XNY(i);

  return Determinant(XYmat) / detXmat;
} // SimplePolyFitterBase<>::Param(unsigned int, FitMatrix_t, Data_t)

template <typename T, unsigned int D>
auto lar::util::details::SimplePolyFitterBase<T, D>::ChiSquare() const -> Data_t
{
  // the generic implementation of ChiSquare from sums is complex enough that
  // I freaked out
  throw std::logic_error("SimplePolyFitterBase::ChiSquare() not implemented for generic fit");
} // SimplePolyFitterBase<>::ChiSquare()

//******************************************************************************
//***  LinearFit<>
//***

template <typename T>
auto lar::util::LinearFit<T>::ChiSquare() const -> Data_t
{
  FitParameters_t fit_params = this->FitParameters();
  const Data_t b = fit_params[0];
  const Data_t a = fit_params[1];
  return Y2() + sqr(a) * X2() + sqr(b) * I() + Data_t(2) * (a * b * X2() - a * XY() - b * Y());
} // LinearFit<T>::ChiSquare()

//******************************************************************************
//***  QuadraticFit<>
//***

template <typename T>
auto lar::util::QuadraticFit<T>::ChiSquare() const -> Data_t
{
  FitParameters_t a = this->FitParameters();
  return Y2() - Data_t(2) * (a[0] * Y() + a[1] * XY() + a[2] * X2Y()) + sqr(a[0]) * I() +
         Data_t(2) * a[0] * (a[1] * X() + a[2] * X2()) + sqr(a[1]) * X2() +
         Data_t(2) * a[1] * (a[2] * X3()) + sqr(a[2]) * X4();
} // QuadraticFit<T>::ChiSquare()

//******************************************************************************
//***  GaussianFit<>
//***

//
// data interface
//
template <typename T>
bool lar::util::GaussianFit<T>::add(Data_t x, Data_t y, Data_t sy /* = Data_t(1.0) */)
{
  if (y <= Data_t(0)) return false; // ignore the non-positive values
  Value_t value = EncodeValue(Value_t(y, sy));
  return fitter.add(x, value.value(), value.error());
} // GaussianFit<T>::add(Data_t, Data_t, Data_t)

template <typename T>
template <typename Iter, typename Pred>
void lar::util::GaussianFit<T>::add_without_uncertainty(Iter begin, Iter end, Pred extractor)
{
  return fitter.add_without_uncertainty(begin, end, Encoder(extractor));
} // GaussianFit<>::add_without_uncertainty(Iter, Iter, Pred)

template <typename T>
template <typename VIter, typename UIter, typename VPred, typename UPred>
unsigned int lar::util::GaussianFit<T>::add_with_uncertainty(
  VIter begin_value,
  VIter end_value,
  UIter begin_uncertainty,
  VPred value_extractor,
  UPred uncertainty_extractor /* = UPred() */
)
{
  return add_with_uncertainty(
    begin_value, end_value, begin_uncertainty, Encoder(value_extractor, uncertainty_extractor));
} // GaussianFit<T>::add_with_uncertainty()

template <typename T>
template <typename Iter>
unsigned int lar::util::GaussianFit<T>::add_with_uncertainty(Iter begin, Iter end)
{
  unsigned int old_n = N();
  std::for_each(begin, end, [this](auto p) { this->add(p); });
  return N() - old_n;
} // GaussianFit<T>::add_with_uncertainty()

//
// fitting interface
//
template <typename T>
auto lar::util::GaussianFit<T>::FitParameters() const -> FitParameters_t
{
  return ConvertParameters(fitter.FitParameters());
} // GaussianFit<>::FitParameters()

template <typename T>
auto lar::util::GaussianFit<T>::FitParameterErrors() const -> FitParameters_t
{
  FitParameters_t qpars, qparerrors;
  if (!FillResults(qpars, qparerrors)) {
    throw std::runtime_error("GaussianFit::FitParameterErrors() yielded invalid results");
  }
  return qparerrors;
} // GaussianFit<>::FitParameterErrors()

template <typename T>
auto lar::util::GaussianFit<T>::FitParameterCovariance() const -> FitMatrix_t
{
  // we need to go through the whole chain to get the error matrix
  FitParameters_t params;
  FitMatrix_t Xmat;
  Data_t det;
  FitMatrix_t Smat;
  if (!FillResults(params, Xmat, det, Smat)) {
    throw std::runtime_error("GaussianFit::FitParameterCovariance() yielded invalid results");
  }
  return Smat;
} // SimplePolyFitterBase<>::FitParameterCovariance()

template <typename T>
bool lar::util::GaussianFit<T>::FillResults(FitParameters_t& params,
                                            FitParameters_t& paramerrors) const
{
  FitParameters_t qpars;
  FitMatrix_t qparerrmat;
  FitMatrix_t Xmat; // not used
  Data_t det;       // not used
  if (!fitter.FillResults(qpars, Xmat, det, qparerrmat)) return false;
  ConvertParametersAndErrors(qpars, qparerrmat, params, paramerrors);
  return isValid(params, qpars);
} // GaussianFit<>::FillResults()

template <typename T>
bool lar::util::GaussianFit<T>::FillResults(FitParameters_t& params,
                                            FitMatrix_t& Xmat,
                                            Data_t& det,
                                            FitMatrix_t& Smat) const
{
  FitParameters_t qpars;
  FitMatrix_t qparerrmat;
  if (!fitter.FillResults(qpars, Xmat, det, qparerrmat)) return false;
  ConvertParametersAndErrorMatrix(qpars, qparerrmat, params, Smat);
  return isValid(params, qpars);
} // GaussianFit::FillResults()

template <typename T>
bool lar::util::GaussianFit<T>::FillResults(FitParameters_t& params,
                                            FitParameters_t& paramerrors,
                                            FitMatrix_t& Xmat,
                                            Data_t& det,
                                            FitMatrix_t& Smat) const
{
  if (!FillResults(params, Xmat, det, Smat)) return false;
  paramerrors = fitter.ExtractParameterErrors(Smat);
  return true;
} // GaussianFit::FillResults()

template <typename T>
auto lar::util::GaussianFit<T>::Evaluate(Data_t x, Data_t const* params) -> Data_t
{
  Data_t z = (x - params[1]) / params[2];
  return params[0] * std::exp(-0.5 * sqr(z));
} // GaussianFit<>::Evaluate()

template <typename T>
auto lar::util::GaussianFit<T>::ConvertParameters(FitParameters_t const& qpars) -> FitParameters_t
{
  FitParameters_t params;

  Data_t sigma2 = -0.5 / qpars[2]; // sigma^2 = -1 / (2 a2)
  params[2] = std::sqrt(sigma2);   // sigma

  params[1] = sigma2 * qpars[1]; // mean = sigma2 a1

  params[0] = std::exp(qpars[0] - 0.25 * sqr(qpars[1]) / qpars[2]);

  return params;
} // GaussianFit<>::ConvertParameters()

template <typename T>
void lar::util::GaussianFit<T>::ConvertParametersAndVariances(FitParameters_t const& qpars,
                                                              FitMatrix_t const& qparerrmat,
                                                              FitParameters_t& params,
                                                              FitParameters_t& paramvariances)
{
  params = ConvertParameters(qpars);

  FitParameters_t const& a = qpars;
  Data_t const& A = params[0];
  Data_t const& mu = params[1];
  Data_t const& sigma = params[2];

  // error on sigma
  paramvariances[2] = qparerrmat[3 * 2 + 2] / sqr(cube(sigma));

  // error on mu
  paramvariances[1] =
    sqr(mu * (+qparerrmat[3 * 1 + 1] / sqr(a[1]) - 2. * qparerrmat[3 * 2 + 1] / (a[1] * a[2]) +
              qparerrmat[3 * 2 + 2] / sqr(a[2])));

  // error on A
  paramvariances[0] =
    sqr(A * (+qparerrmat[3 * 0 + 0] + 2. * qparerrmat[3 * 0 + 1] * mu +
             (qparerrmat[3 * 1 + 1] + 2. * qparerrmat[3 * 0 + 2]) * sqr(mu) +
             2. * qparerrmat[3 * 1 + 2] * cube(mu) + qparerrmat[3 * 2 + 2] * sqr(sqr(mu))));

} // GaussianFit<>::ConvertParametersAndVariances()

template <typename T>
void lar::util::GaussianFit<T>::ConvertParametersAndErrors(FitParameters_t const& qpars,
                                                           FitMatrix_t const& qparerrmat,
                                                           FitParameters_t& params,
                                                           FitParameters_t& paramerrors)
{
  ConvertParametersAndVariances(qpars, qparerrmat, params, paramerrors);
  // paramerrors actually stores the square of the error; fix it:
  for (Data_t& paramerror : paramerrors)
    paramerror = std::sqrt(paramerror);
} // GaussianFit<>::ConvertParametersAndErrors()

template <typename T>
void lar::util::GaussianFit<T>::ConvertParametersAndErrorMatrix(FitParameters_t const& qpars,
                                                                FitMatrix_t const& qparerrmat,
                                                                FitParameters_t& params,
                                                                FitMatrix_t& Smat)
{
  FitParameters_t paramvariances;
  ConvertParametersAndVariances(qpars, qparerrmat, params, paramvariances);

  // let's call things with their names
  FitParameters_t const& a = qpars;
  Data_t const& A = params[0];
  Data_t const& mu = params[1];
  Data_t const& sigma = params[2];

  // variance on sigma
  Smat[3 * 2 + 2] = paramvariances[2];

  // variance on mu
  Smat[3 * 1 + 1] = paramvariances[1];

  // variance on A
  Smat[3 * 0 + 0] = paramvariances[0];

  // covariance on sigma and mu
  Smat[3 * 1 + 2] = Smat[3 * 2 + 1] =
    (qparerrmat[3 * 1 + 2] + 2 * mu * qparerrmat[3 * 2 + 2]) / sigma;

  // this is the sum of the derivatives of A vs. all a parameters, each one
  // multiplied by the covariance of that parameter with a2
  const Data_t dA_dak_cov_aka2 =
    A * (qparerrmat[3 * 0 + 2] + qparerrmat[3 * 1 + 2] * mu + qparerrmat[3 * 2 + 2] * sqr(mu));
  // covariance on A and sigma
  Smat[3 * 0 + 2] = Smat[3 * 2 + 0] = dA_dak_cov_aka2 / cube(sigma);

  // this other is the same as dA_dak_cov_aka2, but for a1
  const Data_t dA_dak_cov_aka1 =
    A * (qparerrmat[3 * 0 + 1] + qparerrmat[3 * 1 + 1] * mu + qparerrmat[3 * 2 + 1] * sqr(mu));

  // covariance on A and mu
  Smat[3 * 0 + 1] = Smat[3 * 1 + 0] = mu * (dA_dak_cov_aka1 / a[1] - dA_dak_cov_aka2 / a[2]);

} // GaussianFit<>::ConvertParametersAndErrors()

template <typename T>
bool lar::util::GaussianFit<T>::isValid(FitParameters_t const& params, FitParameters_t const& qpars)
{
  return (qpars[2] < Data_t(0)) && (params[0] >= Data_t(0));
} // GaussianFit<>::isValid(FitParameters_t)

//******************************************************************************

#endif // SIMPLEFITS_H
