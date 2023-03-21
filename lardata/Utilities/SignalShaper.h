#ifndef SIGNALSHAPER_H
#define SIGNALSHAPER_H

#include <complex>
#include <vector>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFTW.h"
#include "lardata/Utilities/LArFFTWPlan.h"

namespace util {

  class SignalShaper {
  public:
    // Constructor, destructor.
    SignalShaper(int fftsize, std::string fftopt);
    virtual ~SignalShaper();

    // Accessors.
    const std::vector<double>& Response() const { return fResponse; }
    const std::vector<double>& Response_save() const { return fResponse_save; }
    const std::vector<std::complex<double>>& ConvKernel() const { return fConvKernel; }
    const std::vector<std::complex<double>>& Filter() const { return fFilter; }
    const std::vector<std::complex<double>>& DeconvKernel() const { return fDeconvKernel; }

    // Signal shaping methods.

    // Convolute a time series with convolution kernel.
    template <class T>
    void Convolute(std::vector<T>& func) const;

    // Convolute a time series with deconvolution kernel.
    template <class T>
    void Deconvolute(std::vector<T>& func) const;

    // Configuration methods.

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

    // Shift response function in time.
    // Updates overall response function and convolution kernel.
    void ShiftResponseTime(double ticks);
    void SetPeakResponseTime(double tick);

    // Add a filter function.
    void AddFilterFunction(const std::vector<std::complex<double>>& filt);

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
    std::vector<std::complex<double>> fConvKernel;

    // Overall filter function.
    std::vector<std::complex<double>> fFilter;

    // Deconvolution kernel (= fFilter / fConvKernel).
    mutable std::vector<std::complex<double>> fDeconvKernel;

    // Deconvolution Kernel Polarity Flag
    // Set to +1 if deconv signal should be deconv to + ADC count
    // Set to -1 if one wants to normalize to - ADC count
    int fDeconvKernelPolarity;

    // Xin added */
    bool fNorm;

    int fFFTSize;
    const void* fPlan;
    const void* rPlan;
    std::unique_ptr<util::LArFFTWPlan> fFFTPlan;
    std::unique_ptr<util::LArFFTW> fFFT;
  };

}

#endif
