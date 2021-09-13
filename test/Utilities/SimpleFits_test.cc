/**
 * @file    SimpleFits_test.cc
 * @brief   Tests the classes in SimpleFits.h
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    20141229
 * @version 1.0
 * @see     SimpleFits.h
 *
 * See http://www.boost.org/libs/test for the Boost test library home page.
 *
 * Timing:
 * not given yet
 */

// define this to see some fit information
#define SIMPLEFITS_TEST_DEBUG

// C/C++ standard libraries
#include <cmath>
#include <tuple>
#include <array>
#include <stdexcept> // std::range_error
#include <iterator> // std::ostream_iterator
#include <iostream>

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
#define BOOST_TEST_MODULE ( SimpleFits_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()
#include <boost/test/tools/floating_point_comparison.hpp> // BOOST_CHECK_CLOSE()

// LArSoft libraries
#include "lardata/Utilities/SimpleFits.h"


//==============================================================================
//=== Test code
//===

//------------------------------------------------------------------------------
//--- Fit tests
//---


template <typename Array>
void PrintMatrix
  (std::ostream& out, Array const& m, std::string name = "Matrix")
{
  const size_t Dim = size_t(std::sqrt(m.size()));

  out << name << " " << Dim << "x" << Dim << ":";
  for (size_t r = 0; r < Dim; ++r) {
    out << "\n |";
    for (size_t c = 0; c < Dim; ++c)
      out << " " << m[r * Dim + c];
    out << " |";
  } // for
  out << std::endl;
} // PrintMatrix()


template <typename Array>
void PrintVector
  (std::ostream& out, Array const& v, std::string name = "Vector")
{
  using Data_t = typename Array::value_type;
  const size_t Dim = v.size();

  out << name << " [" << Dim << "]: { ";
  std::copy(v.begin(), v.end(), std::ostream_iterator<Data_t>(std::cout, " "));
  out << "}" << std::endl;
} // PrintVector()


template <typename Fitter>
void PrintFitterInfo(Fitter const& fitter) {
#ifdef SIMPLEFITS_TEST_DEBUG
  typename Fitter::FitParameters_t params;
  typename Fitter::FitMatrix_t Xmat, Smat;
  typename Fitter::Data_t det;

  bool res = fitter.FillResults(params, Xmat, det, Smat);
  fitter.PrintStats(std::cout);
  PrintMatrix(std::cout, Xmat, "X matrix");
  std::cout << "det(X) = " << det << std::endl;
  PrintMatrix(std::cout, Smat, "S matrix");
  PrintVector(std::cout, params, "Fit parameters");
  if (!res) std::cout << "The fit results are marked as invalid!" << std::endl;
#else // !SIMPLEFITS_TEST_DEBUG
  fitter.isValid(); // just to avoid compiler warnings
#endif // SIMPLEFITS_TEST_DEBUG
} // PrintFitterInfo()


template <typename T>
void CheckLinearFit(
  lar::util::LinearFit<T> const& fitter,
  int n,
  T intercept,
  T slope,
  T intercept_error,
  T slope_error,
  T intercept_slope_covariance,
  T chisq,
  int NDF
) {

  BOOST_CHECK_EQUAL(fitter.N(), n);
  if (n == 0) {
    BOOST_CHECK(!fitter.isValid());
    BOOST_CHECK_THROW(fitter.Slope(), std::range_error);
    BOOST_CHECK_THROW(fitter.Intercept(), std::range_error);
    BOOST_CHECK_THROW(fitter.SlopeError(), std::range_error);
    BOOST_CHECK_THROW(fitter.InterceptError(), std::range_error);
    BOOST_CHECK_THROW(fitter.InterceptSlopeCovariance(), std::range_error);
    BOOST_CHECK_THROW(fitter.ChiSquare(), std::range_error);
    BOOST_CHECK_EQUAL(fitter.NDF(), -2);
  }
  else {
    BOOST_CHECK(fitter.isValid());

    PrintFitterInfo(fitter);

    BOOST_CHECK_CLOSE(double(fitter.Intercept()), double(intercept), 0.1);
    BOOST_CHECK_CLOSE(double(fitter.Slope()), double(slope), 0.1);
    BOOST_CHECK_CLOSE
      (double(fitter.InterceptError()), double(intercept_error), 0.1);
    BOOST_CHECK_CLOSE(double(fitter.SlopeError()), double(slope_error), 0.1);
    BOOST_CHECK_CLOSE(double(fitter.InterceptSlopeCovariance()),
      double(intercept_slope_covariance), 0.1);
    if (double(chisq) == 0.)
      BOOST_CHECK_SMALL(double(fitter.ChiSquare()), 1e-5);
    else
      BOOST_CHECK_CLOSE(double(fitter.ChiSquare()), double(chisq), 0.1);
    BOOST_CHECK_EQUAL(fitter.NDF(), NDF);
  }

} // CheckLinearFit<>()


template <typename T>
void CheckQuadraticFit(
  lar::util::QuadraticFit<T> const& fitter,
  int n,
  std::array<T, 3> const& solution,
  std::array<T, 3> const& error2,
  T chisq,
  int NDF
) {

  BOOST_CHECK_EQUAL(fitter.N(), n);
  if (n == 0) {
    BOOST_CHECK(!fitter.isValid());
    BOOST_CHECK_THROW(fitter.FitParameter(0), std::range_error);
    BOOST_CHECK_THROW(fitter.FitParameter(1), std::range_error);
    BOOST_CHECK_THROW(fitter.FitParameter(2), std::range_error);
    BOOST_CHECK_THROW(fitter.FitParameterError(0), std::range_error);
    BOOST_CHECK_THROW(fitter.FitParameterError(1), std::range_error);
    BOOST_CHECK_THROW(fitter.FitParameterError(2), std::range_error);
    BOOST_CHECK_THROW(fitter.ChiSquare(), std::range_error);
    BOOST_CHECK_EQUAL(fitter.NDF(), -3);
  }
  else {
    BOOST_CHECK(fitter.isValid());

    PrintFitterInfo(fitter);

    // tolerance: 0.1%
    BOOST_CHECK_CLOSE(double(fitter.FitParameter(0)), double(solution[0]), 0.1);
    BOOST_CHECK_CLOSE(double(fitter.FitParameter(1)), double(solution[1]), 0.1);
    BOOST_CHECK_CLOSE(double(fitter.FitParameter(2)), double(solution[2]), 0.1);
    BOOST_CHECK_CLOSE
      (double(fitter.FitParameterError(0)), std::sqrt(double(error2[0])), 0.1);
    BOOST_CHECK_CLOSE
      (double(fitter.FitParameterError(1)), std::sqrt(double(error2[1])), 0.1);
    BOOST_CHECK_CLOSE
      (double(fitter.FitParameterError(2)), std::sqrt(double(error2[2])), 0.1);
    if (double(chisq) == 0.)
      BOOST_CHECK_SMALL(double(fitter.ChiSquare()), 1e-5);
    else
      BOOST_CHECK_CLOSE(double(fitter.ChiSquare()), double(chisq), 0.1);
    BOOST_CHECK_EQUAL(fitter.NDF(), NDF);
  }

} // CheckQuadraticFit<>()


template <typename T>
void CheckGaussianFit(
  lar::util::GaussianFit<T> const& fitter,
  int n,
  std::array<T, 3> const& solution,
  std::array<T, 3> const& error2,
  T chisq,
  int NDF
) {

  BOOST_CHECK_EQUAL(fitter.N(), n);
  if (n == 0) {
    BOOST_CHECK(!fitter.isValid());
    BOOST_CHECK_THROW(fitter.FitParameters(), std::runtime_error);
    BOOST_CHECK_THROW(fitter.FitParameterErrors(), std::runtime_error);
    BOOST_CHECK_THROW(fitter.ChiSquare(), std::runtime_error);
    BOOST_CHECK_EQUAL(fitter.NDF(), -3);
  }
  else {
    BOOST_CHECK(fitter.isValid());

    using FitParameters_t = typename lar::util::GaussianFit<T>::FitParameters_t;

    // I am disabling the check on the parameter errors since I have no
    // cross check available

    PrintFitterInfo(fitter.Fitter());
    PrintFitterInfo(fitter);
    FitParameters_t  params = fitter.FitParameters();
    /* FitParameters_t perrors = */ fitter.FitParameterErrors();

    // tolerance: 0.1%
    BOOST_CHECK_CLOSE(double(params[0]), double(solution[0]), 0.1);
    BOOST_CHECK_CLOSE(double(params[1]), double(solution[1]), 0.1);
    BOOST_CHECK_CLOSE(double(params[2]), double(solution[2]), 0.1);
  //  BOOST_CHECK_CLOSE(double(perrors[0]), std::sqrt(double(error2[0])), 0.1);
  //  BOOST_CHECK_CLOSE(double(perrors[1]), std::sqrt(double(error2[1])), 0.1);
  //  BOOST_CHECK_CLOSE(double(perrors[2]), std::sqrt(double(error2[2])), 0.1);
    if (double(chisq) == 0.)
      BOOST_CHECK_SMALL(double(fitter.ChiSquare()), 1e-5);
    else
      BOOST_CHECK_CLOSE(double(fitter.ChiSquare()), double(chisq), 0.1);
    BOOST_CHECK_EQUAL(fitter.NDF(), NDF);
  }

} // CheckGaussianFit<>()



/**
 * @brief Tests LinearFit object with a known input
 */
template <typename T>
void LinearFitTest() {

  using Data_t = T;

  using PerfectItem_t = std::pair<Data_t, Data_t>;
  using UncertainItem_t = std::tuple<Data_t, Data_t, Data_t>;

  using PerfectData_t = std::vector<PerfectItem_t>;
  using UncertainData_t = std::vector<UncertainItem_t>;

  // prepare input data
  PerfectData_t perfect_data{
    { Data_t(-4),  Data_t( 8) },
    { Data_t( 0),  Data_t( 0) },
    { Data_t( 4),  Data_t(-8) }
    };

  const int      n          =         3;
  const Data_t   intercept  = Data_t( 0);
  const Data_t   slope      = Data_t(-2);
  const Data_t   perf_chisq = Data_t( 0);
  const Data_t   perf_intercept_error     = std::sqrt(Data_t(32)/Data_t(96));
  const Data_t   perf_slope_error         = std::sqrt(Data_t( 3)/Data_t(96));
  const Data_t   perf_intercept_slope_cov = - Data_t(0)/Data_t(96);
  const int      perf_DoF   =         1;

  UncertainData_t uncertain_data({
    UncertainItem_t{ Data_t(-4), Data_t( 8), Data_t(1) },
    UncertainItem_t{ Data_t( 0), Data_t( 0), Data_t(2) },
    UncertainItem_t{ Data_t( 4), Data_t(-8), Data_t(2) }
    });

  const Data_t   unc_chisq               = Data_t( 0);
  const Data_t   unc_intercept_error     = std::sqrt(Data_t(20)/Data_t(21));
  const Data_t   unc_slope_error         = std::sqrt(Data_t(1.5)/Data_t(21));
  const Data_t   unc_intercept_slope_cov = - Data_t(-3)/Data_t(21);
  const int      unc_DoF                 =         1;

  //
  // part I: construction
  //
  lar::util::LinearFit<Data_t> fitter;

  // check that everything is 0 or NaN-like
  CheckLinearFit<Data_t>(fitter, 0, 0., 0., 0., 0., 0., 0., 0);

  //
  // part II: add elements one by one
  //
  // the data is the same as uncertain_data, just inserted one by one
  // and exercising both uncertain and certain addition;
  // this part deliberately ignores directly interfaces adding pairs and tuples
  for (auto const& data: uncertain_data) {
    if (std::get<2>(data) == Data_t(1))
      fitter.add(std::get<0>(data), std::get<1>(data));
    else
      fitter.add(std::get<0>(data), std::get<1>(data), std::get<2>(data));
  } // for

  // by construction of the input, the statistics for X and Y are the same
  CheckLinearFit<Data_t>(fitter, n,
    intercept, slope,
    unc_intercept_error, unc_slope_error, unc_intercept_slope_cov,
    unc_chisq, unc_DoF
    );


  //
  // part III: add elements without uncertainty by bulk
  //

  // - III.1: clear the fitter
  fitter.clear();
  CheckLinearFit<Data_t>(fitter, 0, 0., 0., 0., 0., 0., 0., 0);

  // - III.2: fill by iterators
  fitter.add_without_uncertainty
    (std::begin(perfect_data), std::end(perfect_data));
  CheckLinearFit<Data_t>(fitter, n,
    intercept, slope,
    perf_intercept_error, perf_slope_error, perf_intercept_slope_cov,
    perf_chisq, perf_DoF
    );

  // - III.3: fill by container
  fitter.clear();
  fitter.add_without_uncertainty(perfect_data);
  CheckLinearFit<Data_t>(fitter, n,
    intercept, slope,
    perf_intercept_error, perf_slope_error, perf_intercept_slope_cov,
    perf_chisq, perf_DoF
    );

  // - III.4: fill by iterators and extractor
  fitter.clear();
  fitter.add_without_uncertainty(
    uncertain_data.begin(), uncertain_data.end(),
    [](UncertainItem_t const& d)
      { return PerfectItem_t{ std::get<0>(d), std::get<1>(d) }; }
    );
  CheckLinearFit<Data_t>(fitter, n,
    intercept, slope,
    perf_intercept_error, perf_slope_error, perf_intercept_slope_cov,
    perf_chisq, perf_DoF
    );

  // - III.5: fill by container and extractor
  fitter.clear();
  fitter.add_without_uncertainty(uncertain_data,
    [](UncertainItem_t const& d)
      { return PerfectItem_t{ std::get<0>(d), std::get<1>(d) }; }
    );
  CheckLinearFit<Data_t>(fitter, n,
    intercept, slope,
    perf_intercept_error, perf_slope_error, perf_intercept_slope_cov,
    perf_chisq, perf_DoF
    );


  //
  // part IV: add elements with uncertainty by bulk
  //

  // - IV.1: fill by iterators
  fitter.clear();
  fitter.add_with_uncertainty(uncertain_data.begin(), uncertain_data.end());
  CheckLinearFit<Data_t>(fitter, n,
    intercept, slope,
    unc_intercept_error, unc_slope_error, unc_intercept_slope_cov,
    unc_chisq, unc_DoF
    );

  // - IV.2: fill by container
  fitter.clear();
  fitter.add_with_uncertainty(uncertain_data);
  CheckLinearFit<Data_t>(fitter, n,
    intercept, slope,
    unc_intercept_error, unc_slope_error, unc_intercept_slope_cov,
    unc_chisq, unc_DoF
    );

} // LinearFitTest()


/** ****************************************************************************
 * @brief Tests QuadraticFit object with a known input
 */
template <typename T>
void QuadraticFitTest() {

  using Data_t = T;

  using PerfectItem_t = std::pair<Data_t, Data_t>;
  using UncertainItem_t = std::tuple<Data_t, Data_t, Data_t>;

  using PerfectData_t = std::vector<PerfectItem_t>;
  using UncertainData_t = std::vector<UncertainItem_t>;

  // prepare input data
  PerfectData_t perfect_data{
    { Data_t(-4),  Data_t( 9) },
    { Data_t( 0),  Data_t(-1) },
    { Data_t( 4),  Data_t( 5) },
    { Data_t( 6),  Data_t(14) }
    };

  const int                   n            =         4;
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  const std::array<Data_t, 3> solution     = {{ -1.0, -0.5, 0.5 }};
  const std::array<Data_t, 3> perf_errors2 = {{ 149./199., 163./6368., 59./25472. }};
  const Data_t                perf_chisq   = Data_t( 0);
  const int                   perf_DoF     =         1;

  UncertainData_t uncertain_data({
    UncertainItem_t{ Data_t(-4), Data_t( 9), Data_t(2) },
    UncertainItem_t{ Data_t( 0), Data_t(-1), Data_t(1) },
    UncertainItem_t{ Data_t( 4), Data_t( 5), Data_t(1) },
    UncertainItem_t{ Data_t( 6), Data_t(14), Data_t(2) }
    });

  const Data_t                unc_chisq  = Data_t( 0);
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  const std::array<Data_t, 3> unc_errors2 = {{ 517./617., 769./9872., 209./39488. }};
  const int                   unc_DoF    =         1;

  //
  // part I: construction
  //
  lar::util::QuadraticFit<Data_t> fitter;

  std::array<Data_t, 3> empty_params;
  empty_params.fill(0);

  // check that everything is 0 or NaN-like
  CheckQuadraticFit<Data_t>(fitter, 0, empty_params, empty_params, 0., 0);

  //
  // part II: add elements one by one
  //
  // the data is the same as uncertain_data, just inserted one by one
  // and exercising both uncertain and certain addition;
  // this part deliberately ignores directly interfaces adding pairs and tuples
  for (auto const& data: uncertain_data) {
    if (std::get<2>(data) == Data_t(1))
      fitter.add(std::get<0>(data), std::get<1>(data));
    else
      fitter.add(std::get<0>(data), std::get<1>(data), std::get<2>(data));
  } // for

  // by construction of the input, the statistics for X and Y are the same
  CheckQuadraticFit<Data_t>
    (fitter, n, solution, unc_errors2, unc_chisq, unc_DoF);


  //
  // part III: add elements without uncertainty by bulk
  //

  // - III.1: clear the fitter
  fitter.clear();
  CheckQuadraticFit<Data_t>(fitter, 0, empty_params, empty_params, 0., 0);

  // - III.2: fill by iterators
  fitter.add_without_uncertainty
    (std::begin(perfect_data), std::end(perfect_data));
  CheckQuadraticFit<Data_t>
    (fitter, n, solution, perf_errors2, perf_chisq, perf_DoF);

  // - III.3: fill by container
  fitter.clear();
  fitter.add_without_uncertainty(perfect_data);
  CheckQuadraticFit<Data_t>
    (fitter, n, solution, perf_errors2, perf_chisq, perf_DoF);

  // - III.4: fill by iterators and extractor
  fitter.clear();
  fitter.add_without_uncertainty(
    uncertain_data.begin(), uncertain_data.end(),
    [](UncertainItem_t const& d)
      { return PerfectItem_t{ std::get<0>(d), std::get<1>(d) }; }
    );
  CheckQuadraticFit<Data_t>
    (fitter, n, solution, perf_errors2, perf_chisq, perf_DoF);

  // - III.5: fill by container and extractor
  fitter.clear();
  fitter.add_without_uncertainty(uncertain_data,
    [](UncertainItem_t const& d)
      { return PerfectItem_t{ std::get<0>(d), std::get<1>(d) }; }
    );
  CheckQuadraticFit<Data_t>
    (fitter, n, solution, perf_errors2, perf_chisq, perf_DoF);


  //
  // part IV: add elements with uncertainty by bulk
  //

  // - IV.1: fill by iterators
  fitter.clear();
  fitter.add_with_uncertainty(uncertain_data.begin(), uncertain_data.end());
  CheckQuadraticFit<Data_t>
    (fitter, n, solution, unc_errors2, unc_chisq, unc_DoF);

  // - IV.2: fill by container
  fitter.clear();
  fitter.add_with_uncertainty(uncertain_data);
  CheckQuadraticFit<Data_t>
    (fitter, n, solution, unc_errors2, unc_chisq, unc_DoF);

} // QuadraticFitTest()


/** ****************************************************************************
 * @brief Tests GausssianFit object with a known input
 */

template <typename T>
T gaus(T x, T amplitude, T mean, T sigma) {
  const T z = (x - mean) / sigma;
  return amplitude * std::exp(-0.5*z*z);
} // gaus()


template <typename T>
void GaussianFitTest() {

  using Data_t = T;

  using PerfectItem_t = std::pair<Data_t, Data_t>;
  using UncertainItem_t = std::tuple<Data_t, Data_t, Data_t>;

  using PerfectData_t = std::vector<PerfectItem_t>;
  using UncertainData_t = std::vector<UncertainItem_t>;

  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  const std::array<Data_t, 3> solution = {{ 5.0, 1.0, 2.0 }};

  // prepare input data
  PerfectData_t perfect_data{
    { Data_t(-1), gaus(Data_t(-1), solution[0], solution[1], solution[2]) },
    { Data_t( 0), gaus(Data_t( 0), solution[0], solution[1], solution[2]) },
    { Data_t(+1), gaus(Data_t(+1), solution[0], solution[1], solution[2]) },
    { Data_t(+3), gaus(Data_t(+3), solution[0], solution[1], solution[2]) }
    };

  const int                   n            =         4;
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  const std::array<Data_t, 3> perf_errors2 = {{ 0., 0., 0. }};
  const Data_t                perf_chisq   = Data_t( 0);
  const int                   perf_DoF     =         1;

  UncertainData_t uncertain_data({
    UncertainItem_t{ Data_t(-1), gaus(Data_t(-1), solution[0], solution[1], solution[2]), Data_t(2)  },
    UncertainItem_t{ Data_t( 0), gaus(Data_t( 0), solution[0], solution[1], solution[2]), Data_t(1)  },
    UncertainItem_t{ Data_t(+1), gaus(Data_t(+1), solution[0], solution[1], solution[2]), Data_t(1)  },
    UncertainItem_t{ Data_t(+3), gaus(Data_t(+3), solution[0], solution[1], solution[2]), Data_t(2)  }
    });

  const Data_t                unc_chisq  = Data_t( 0);
  // BUG the double brace syntax is required to work around clang bug 21629
  // (https://bugs.llvm.org/show_bug.cgi?id=21629)
  const std::array<Data_t, 3> unc_errors2 = {{ 0., 0., 0. }};
  const int                   unc_DoF    =         1;

  //
  // part I: construction
  //
  lar::util::GaussianFit<Data_t> fitter;

  std::array<Data_t, 3> empty_params;
  empty_params.fill(0);

  // check that everything is 0 or NaN-like
  CheckGaussianFit<Data_t>(fitter, 0, empty_params, empty_params, 0., 0);

  //
  // part II: add elements one by one
  //
  // the data is the same as uncertain_data, just inserted one by one
  // and exercising both uncertain and certain addition;
  // this part deliberately ignores directly interfaces adding pairs and tuples
  for (auto const& data: uncertain_data) {
    if (std::get<2>(data) == Data_t(1))
      fitter.add(std::get<0>(data), std::get<1>(data));
    else
      fitter.add(std::get<0>(data), std::get<1>(data), std::get<2>(data));
  } // for

  // by construction of the input, the statistics for X and Y are the same
  CheckGaussianFit<Data_t>
    (fitter, n, solution, unc_errors2, unc_chisq, unc_DoF);


  //
  // part III: add elements without uncertainty by bulk
  //

  // - III.1: clear the fitter
  fitter.clear();
  CheckGaussianFit<Data_t>(fitter, 0, empty_params, empty_params, 0., 0);

  // - III.2: fill by iterators
  fitter.add_without_uncertainty
    (std::begin(perfect_data), std::end(perfect_data));
  CheckGaussianFit<Data_t>
    (fitter, n, solution, perf_errors2, perf_chisq, perf_DoF);

  // - III.3: fill by container
  fitter.clear();
  fitter.add_without_uncertainty(perfect_data);
  CheckGaussianFit<Data_t>
    (fitter, n, solution, perf_errors2, perf_chisq, perf_DoF);

  // - III.4: fill by iterators and extractor
  fitter.clear();
  fitter.add_without_uncertainty(
    uncertain_data.begin(), uncertain_data.end(),
    [](UncertainItem_t const& d)
      { return PerfectItem_t{ std::get<0>(d), std::get<1>(d) }; }
    );
  CheckGaussianFit<Data_t>
    (fitter, n, solution, perf_errors2, perf_chisq, perf_DoF);

  // - III.5: fill by container and extractor
  fitter.clear();
  fitter.add_without_uncertainty(uncertain_data,
    [](UncertainItem_t const& d)
      { return PerfectItem_t{ std::get<0>(d), std::get<1>(d) }; }
    );
  CheckGaussianFit<Data_t>
    (fitter, n, solution, perf_errors2, perf_chisq, perf_DoF);


  //
  // part IV: add elements with uncertainty by bulk
  //

  // - IV.1: fill by iterators
  fitter.clear();
  fitter.add_with_uncertainty(uncertain_data.begin(), uncertain_data.end());
  CheckGaussianFit<Data_t>
    (fitter, n, solution, unc_errors2, unc_chisq, unc_DoF);

  // - IV.2: fill by container
  fitter.clear();
  fitter.add_with_uncertainty(uncertain_data);
  CheckGaussianFit<Data_t>
    (fitter, n, solution, unc_errors2, unc_chisq, unc_DoF);

} // QuadraticFitTest()


//------------------------------------------------------------------------------
//--- registration of tests
//
// Boost needs now to know which tests we want to run.
// Tests are "automatically" registered, hence the BOOST_AUTO_TEST_CASE()
// macro name. The argument is a name for the test; each test will have a
// number of checks and it will fail if any of them does.
//

//
// LinearFit tests
//
BOOST_AUTO_TEST_CASE(LinearFitRealTest) {
  LinearFitTest<double>();
}

//
// QuadraticFit tests
//
BOOST_AUTO_TEST_CASE(QuadraticFitRealTest) {
  QuadraticFitTest<double>();
}

//
// GaussianFit tests
//
BOOST_AUTO_TEST_CASE(GaussianFitRealTest) {
  GaussianFitTest<double>();
}

