////////////////////////////////////////////////////////////////////////
///
/// \file   SignalShaping.h
///
/// \brief  Generic class for shaping signals on wires.
///
/// \author H. Greenlee
///
/// This is a generic class for shaping signals on wires during simulation
/// (convolution) and reconstruction (deconvolution).
///
/// This class acts as a repository for a consistent set of convolution
/// and deconvolution kernels.  It also supplies an interface for
/// convoluting either type of kernel with a time series signal.  All
/// FFT type calculations are done using LArFFT service.
///
/// This class has only a default constructor.  Configuration must be done
/// externally by calling configuration methods.  The proper method for
/// configuring this class is as follows.
///
/// 1.  Add one or more response functions using method AddReponseFunction.
/// 2.  Optionally call methods SetPeakResponseTime or ShiftResponseTime.
/// 3.  Add one or more filter functions using method AddFilterFunction.
/// 4.  Call method CalculateDeconvKernel once.
///
/// After the deconvolution kernel is calculated, the configuration is locked.
///
/// Notes on time and frequency series functions
/// ---------------------------------------------
///
/// Times and frequencies are measured in units of ticks and cycles/tick.
///
/// Time series are represented as vector<double> of length N, representing
/// sampled times on interval [0,N) ticks.   (N = LArFFT::FFTSize().)
///
/// Frequency series are represented as vector<TComplex> of length (N/2+1),
/// representing sampled frequencies on interval [0, 1/2] cycles/tick.
/// Negative frequencies (not stored) are complex conjugate of
/// corresponding positive frequency.
///
/// Update notes
/// -------------
///
/// * Yun-Tse Tsai (yuntse@slac.stanford.edu), July 17th, 2014<br/>
///     Modify
///     `void AddResponseFunction(const std::vector<double>& resp);`
///     to
///     `void AddResponseFunction(const std::vector<double>& resp, bool ResetResponse = false );`
///     If you want to reset your response, `fResponse` in this object, you can
///     do
///     `AddResponseFunction( yourResponse, true )`
///     The other part involving `AddResponseFunction` shouldn't be affected.
/// * X. Qian 2015/01/06 <br/>
///     Add the time offset variable<br/>
///     Need to add the set and extraction code
////////////////////////////////////////////////////////////////////////

#ifndef SIGNALSHAPING_H
#define SIGNALSHAPING_H

#include "TComplex.h"
#include <vector>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"

namespace util {

  class SignalShaping {
  public:
    // Constructor, destructor.
    SignalShaping();
    virtual ~SignalShaping();

    // Accessors.
    const std::vector<double>& Response() const { return fResponse; }
    const std::vector<double>& Response_save() const { return fResponse_save; }
    const std::vector<TComplex>& ConvKernel() const { return fConvKernel; }
    const std::vector<TComplex>& Filter() const { return fFilter; }
    const std::vector<TComplex>& DeconvKernel() const { return fDeconvKernel; }
    /* const int GetTimeOffset() const {return fTimeOffset;} */

    // Signal shaping methods.

    // Convolute a time series with convolution kernel.
    template <class T>
    void Convolute(std::vector<T>& func) const;

    // Convolute a time series with deconvolution kernel.
    template <class T>
    void Deconvolute(std::vector<T>& func) const;

    // Configuration methods.

    // Only reset deconvolution
    //void ResetDecon();
    // Reset this class to default-constructed state.
    void Reset();

    void save_response()
    {
      fResponse_save.clear();
      fResponse_save = fResponse;
    }
    void set_normflag(bool flag) { fNorm = flag; }

    // Add a time domain response function.
    // Updates overall response function and convolution kernel.
    void AddResponseFunction(const std::vector<double>& resp, bool ResetResponse = false);

    /* //X. Qian, set time offset */
    /* void SetTimeOffset(const int time){fTimeOffset = time;} */

    // Shift response function in time.
    // Updates overall response function and convolution kernel.
    void ShiftResponseTime(double ticks);
    void SetPeakResponseTime(double tick);

    // Add a filter function.
    void AddFilterFunction(const std::vector<TComplex>& filt);

    //Add DeconvKernel Polarity switch to decide how to normalize
    //deconvoluted signal w.r.t. RawDigits. If +1 then normalize
    //to Max ADC, if -1 to Min ADC
    void SetDeconvKernelPolarity(int pol);

    // Test and lock the current response function.
    // Does not lock filter configuration.
    void LockResponse() const;

    // Calculate deconvolution kernel using current convolution kernel
    // and filter function.
    // Fully locks configuration.
    void CalculateDeconvKernel() const;

  private:
    // Attributes.
    // unused double fMinConvKernelFrac;  ///< minimum value of convKernel/peak for deconvolution

    // Lock flags.
    mutable bool fResponseLocked;
    mutable bool fFilterLocked;

    // Overall response.
    std::vector<double> fResponse;
    std::vector<double> fResponse_save;

    // Convolution kernel (fourier transform of response function).
    std::vector<TComplex> fConvKernel;

    // Overall filter function.
    std::vector<TComplex> fFilter;

    // Deconvolution kernel (= fFilter / fConvKernel).
    mutable std::vector<TComplex> fDeconvKernel;

    // Deconvolution Kernel Polarity Flag
    // Set to +1 if deconv signal should be deconv to + ADC count
    // Set to -1 if one wants to normalize to - ADC count
    int fDeconvKernelPolarity;

    // Xin added */
    bool fNorm;
  };

}

//----------------------------------------------------------------------
// Convolute a time series with current response.
template <class T>
inline void util::SignalShaping::Convolute(std::vector<T>& func) const
{
  // Make sure response configuration is locked.
  if (!fResponseLocked) LockResponse();

  art::ServiceHandle<util::LArFFT> fft;

  // Make sure that time series has the correct size.
  if (int const n = func.size(); n != fft->FFTSize())
    throw cet::exception("SignalShaping") << "Bad time series size = " << n << "\n";

  fft->Convolute(func, const_cast<std::vector<TComplex>&>(fConvKernel));
}

//----------------------------------------------------------------------
// Convolute a time series with deconvolution kernel.
template <class T>
inline void util::SignalShaping::Deconvolute(std::vector<T>& func) const
{
  // Make sure deconvolution kernel is configured.
  if (!fFilterLocked) CalculateDeconvKernel();

  art::ServiceHandle<util::LArFFT> fft;

  // Make sure that time series has the correct size.
  if (int const n = func.size(); n != fft->FFTSize())
    throw cet::exception("SignalShaping") << "Bad time series size = " << n << "\n";

  fft->Convolute(func, fDeconvKernel);
}

#endif
