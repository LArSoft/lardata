#ifndef LARFFTW_H
#define LARFFTW_H

// C/C++ standard libraries
#include <string>
#include <vector>
#include <complex>
#include <algorithm>

#include "fftw3.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/coded_exception.h"
#include "lardata/Utilities/MarqFitAlg.h"

namespace util {

class LArFFTW {

  public:

    using FloatVector = std::vector<float>;
    using DoubleVector = std::vector<double>;
    using ComplexVector = std::vector<std::complex<double>>;

    LArFFTW(int transformSize, const void* fplan, const void* rplan, int fitbins);
    ~LArFFTW();

    template <class T> void DoFFT(std::vector<T>& input);
    template <class T> void DoFFT(std::vector<T>& input, ComplexVector& output);
    template <class T> void DoInvFFT(std::vector<T>& output);
    template <class T> void DoInvFFT(ComplexVector& input, std::vector<T>& output);

    // ... Do convolution calculation (for simulation).
    template <class T> void Convolute(std::vector<T>& func, const ComplexVector& kern);
    template <class T> void Convolute(std::vector<T>& func, std::vector<T>& resp);

    // ... Do deconvolution calculation (for reconstruction).
    template <class T> void Deconvolute(std::vector<T>& func, const ComplexVector& kern);
    template <class T> void Deconvolute(std::vector<T>& func, std::vector<T>& resp);

    // ... Do correlation
    template <class T> void Correlate(std::vector<T>& func, const ComplexVector& kern);
    template <class T> void Correlate(std::vector<T>& func, std::vector<T>& resp);

    void ShiftData(ComplexVector & input, double shift);
    template <class T> void ShiftData(std::vector<T> & input, double shift);

    template <class T> void AlignedSum(std::vector<T> & input, std::vector<T> &output,
                                       bool add = true);
    template <class T> T PeakCorrelation(std::vector<T> &shape1,std::vector<T> &shape2);

  private:

    ComplexVector fKern;	// transformed response function
    ComplexVector fCompTemp;	// temporary complex data
    std::vector<float> fConvHist;	// Fit data histogram
    int fSize;			// size of transform
    int fFreqSize;		// size of frequency space
    void *fIn;
    void *fOut;
    const void *fPlan;
    void *rIn;
    void *rOut;
    const void *rPlan;
    int fFitBins;		// Bins used for peak fit

    gshf::MarqFitAlg* fMarqFitAlg;
};

}  // end namespace util

// -----------------------------------------------------------------------------
// ~~~~ Do Forward Fourier Transform - DoFFT( REAL In )
// -----------------------------------------------------------------------------
template <class T> inline void util::LArFFTW::DoFFT(std::vector<T> & input)
{
  // ..set point
  for(size_t p = 0; p < input.size(); ++p){
    ((double *)fIn)[p] = input[p];
  }

  // ..transform (using the New-array Execute Functions)
  fftw_execute_dft_r2c((fftw_plan)fPlan,(double*)fIn,(fftw_complex*)fOut);

  return;
}

// -----------------------------------------------------------------------------
// ~~~~ Do Forward Fourier Transform - DoFFT( REAL In, COMPLEX Out )
// -----------------------------------------------------------------------------
template <class T> inline void util::LArFFTW::DoFFT(std::vector<T> & input, ComplexVector& output)
{
  // ..set point
  for(size_t p = 0; p < input.size(); ++p){
    ((double *)fIn)[p] = input[p];
  }

  // ..transform (using the New-array Execute Functions)
  fftw_execute_dft_r2c((fftw_plan)fPlan,(double*)fIn,(fftw_complex*)fOut);

  for(int i = 0; i < fFreqSize; ++i){
    output[i].real(((fftw_complex*)fOut)[i][0]);
    output[i].imag(((fftw_complex*)fOut)[i][1]);
  }

  return;
}

// -----------------------------------------------------------------------------
// ~~~~ Do Inverse Fourier Transform - DoInvFFT( REAL Out )
// -----------------------------------------------------------------------------
template <class T> inline void util::LArFFTW::DoInvFFT(std::vector<T> & output)
{
  // ..transform (using the New-array Execute Functions)
  fftw_execute_dft_c2r((fftw_plan)rPlan,(fftw_complex*)rIn,(double*)rOut);

  // ..get point real
  double factor = 1.0/(double) fSize;
  const double * array =  (const double*)(rOut);
  for(int i = 0; i < fSize; ++i){
    output[i] = factor*array[i];
  }

  return;
}

// -----------------------------------------------------------------------------
// ~~~~ Do Inverse Fourier Transform - DoInvFFT( COMPLEX In, REAL Out )
// -----------------------------------------------------------------------------
template <class T> inline void util::LArFFTW::DoInvFFT(ComplexVector& input, std::vector<T> & output)
{
  // ..set point complex
  for(int i = 0; i < fFreqSize; ++i){
    ((fftw_complex*)rIn)[i][0] = input[i].real();
    ((fftw_complex*)rIn)[i][1] = input[i].imag();
  }

  // ..transform (using the New-array Execute Functions)
  fftw_execute_dft_c2r((fftw_plan)rPlan,(fftw_complex*)rIn,(double*)rOut);

  // ..get point real
  double factor = 1.0/(double) fSize;
  const double * array =  (const double*)(rOut);
  for(int i = 0; i < fSize; ++i){
    output[i] = factor*array[i];
  }

  return;
}

// -----------------------------------------------------------------------------
// ~~~~ Do Convolution: using transformed response function
// -----------------------------------------------------------------------------
template <class T>
inline void util::LArFFTW::Convolute(std::vector<T>& func,
                                          const ComplexVector& kern){

  // ... Make sure that time series and kernel have the correct size.
  int n = func.size();
  if(n != fSize){
    throw cet::exception("LArFFTW") << "Bad time series size = " << n << "\n";
  }
  n = kern.size();
  if(n != fFreqSize){
    throw cet::exception("LArFFTW") << "Bad kernel size = " << n << "\n";
  }

  DoFFT(func);

  // ..perform the convolution
  for(int i = 0; i < fFreqSize; ++i){
    double re = ((fftw_complex*)fOut)[i][0];
    double im = ((fftw_complex*)fOut)[i][1];
    ((fftw_complex*)rIn)[i][0] = re*kern[i].real()-im*kern[i].imag();
    ((fftw_complex*)rIn)[i][1] = re*kern[i].imag()+im*kern[i].real();
  }

  DoInvFFT(func);
}

// -----------------------------------------------------------------------------
// ~~~~ Do Convolution: using all time-domain information
// -----------------------------------------------------------------------------
template <class T>
inline void util::LArFFTW::Convolute(std::vector<T>& func1,
                                          std::vector<T>& func2){

  // ... Make sure that time series has the correct size.
  int n = func1.size();
  if(n != fSize){
    throw cet::exception("LArFFTW") << "Bad 1st time series size = " << n << "\n";
  }
  n = func2.size();
  if(n != fSize){
    throw cet::exception("LArFFTW") << "Bad 2nd time series size = " << n << "\n";
  }

  DoFFT(func2);
  for(int i = 0; i < fFreqSize; ++i){
    fKern[i].real(((fftw_complex*)fOut)[i][0]);
    fKern[i].imag(((fftw_complex*)fOut)[i][1]);
  }
  DoFFT(func1);

  // ..perform the convolution
  for(int i = 0; i < fFreqSize; ++i){
    double re = ((fftw_complex*)fOut)[i][0];
    double im = ((fftw_complex*)fOut)[i][1];
    ((fftw_complex*)rIn)[i][0] = re*fKern[i].real()-im*fKern[i].imag();
    ((fftw_complex*)rIn)[i][1] = re*fKern[i].imag()+im*fKern[i].real();
  }

  DoInvFFT(func1);
}

// -----------------------------------------------------------------------------
// ~~~~ Do Deconvolution: using transformed response function
// -----------------------------------------------------------------------------
template <class T>
inline void util::LArFFTW::Deconvolute(std::vector<T>& func,
                                            const ComplexVector& kern){

  // ... Make sure that time series and kernel have the correct size.
  int n = func.size();
  if(n != fSize){
    throw cet::exception("LArFFTW") << "Bad time series size = " << n << "\n";
  }
  n = kern.size();
  if(n != fFreqSize){
    throw cet::exception("LArFFTW") << "Bad kernel size = " << n << "\n";
  }

  DoFFT(func);

  // ..perform the deconvolution
  double a,b,c,d,e;
  for(int i = 0; i < fFreqSize; ++i){
    a = ((fftw_complex*)fOut)[i][0];
    b = ((fftw_complex*)fOut)[i][1];
    c = kern[i].real();
    d = kern[i].imag();
    e = 1./(c*c+d*d);
    ((fftw_complex*)rIn)[i][0] = (a*c+b*d)*e;
    ((fftw_complex*)rIn)[i][1] = (b*c-a*d)*e;
  }

  DoInvFFT(func);
}

// -----------------------------------------------------------------------------
// ~~~~ Do Deconvolution: using all time domain information
// -----------------------------------------------------------------------------
template <class T>
inline void util::LArFFTW::Deconvolute(std::vector<T>& func,
                                            std::vector<T>& resp){

  // ... Make sure that time series has the correct size.
  int n = func.size();
  if(n != fSize){
    throw cet::exception("LArFFTW") << "Bad 1st time series size = " << n << "\n";
  }
  n = resp.size();
  if(n != fSize){
    throw cet::exception("LArFFTW") << "Bad 2nd time series size = " << n << "\n";
  }

  DoFFT(resp);
  for(int i = 0; i < fFreqSize; ++i){
    fKern[i].real(((fftw_complex*)fOut)[i][0]);
    fKern[i].imag(((fftw_complex*)fOut)[i][1]);
  }
  DoFFT(func);

  // ..perform the deconvolution
  double a,b,c,d,e;
  for(int i = 0; i < fFreqSize; ++i){
    a = ((fftw_complex*)fOut)[i][0];
    b = ((fftw_complex*)fOut)[i][1];
    c = fKern[i].real();
    d = fKern[i].imag();
    e = 1./(c*c+d*d);
    ((fftw_complex*)rIn)[i][0] = (a*c+b*d)*e;
    ((fftw_complex*)rIn)[i][1] = (b*c-a*d)*e;
  }

  DoInvFFT(func);

}

// -----------------------------------------------------------------------------
// ~~~~ Do Deconvolution: using transformed response function
// -----------------------------------------------------------------------------
template <class T>
inline void util::LArFFTW::Correlate(std::vector<T>& func,
                                          const ComplexVector& kern){

  // ... Make sure that time series and kernel have the correct size.
  int n = func.size();
  if(n != fSize){
    throw cet::exception("LArFFTW") << "Bad time series size = " << n << "\n";
  }
  n = kern.size();
  if(n != fFreqSize){
    throw cet::exception("LArFFTW") << "Bad kernel size = " << n << "\n";
  }

  DoFFT(func);

  // ..perform the correlation
  for(int i = 0; i < fFreqSize; ++i){
    double re = ((fftw_complex*)fOut)[i][0];
    double im = ((fftw_complex*)fOut)[i][1];
    ((fftw_complex*)rIn)[i][0] =  re*kern[i].real()+im*kern[i].imag();
    ((fftw_complex*)rIn)[i][1] = -re*kern[i].imag()+im*kern[i].real();
  }

  DoInvFFT(func);

}

// -----------------------------------------------------------------------------
// ~~~~ Do Correlation: using all time domain information
// -----------------------------------------------------------------------------
template <class T>
inline void util::LArFFTW::Correlate(std::vector<T>& func1,
                                          std::vector<T>& func2){

  // ... Make sure that time series has the correct size.
  int n = func1.size();
  if(n != fSize){
    throw cet::exception("LArFFTW") << "Bad 1st time series size = " << n << "\n";
  }
  n = func2.size();
  if(n != fSize){
    throw cet::exception("LArFFTW") << "Bad 2nd time series size = " << n << "\n";
  }

  DoFFT(func2);
  for(int i = 0; i < fFreqSize; ++i){
    fKern[i].real(((fftw_complex*)fOut)[i][0]);
    fKern[i].imag(((fftw_complex*)fOut)[i][1]);
  }
  DoFFT(func1);

  // ..perform the correlation
  for(int i = 0; i < fFreqSize; ++i){
    double re = ((fftw_complex*)fOut)[i][0];
    double im = ((fftw_complex*)fOut)[i][1];
    ((fftw_complex*)rIn)[i][0] =  re*fKern[i].real()+im*fKern[i].imag();
    ((fftw_complex*)rIn)[i][1] = -re*fKern[i].imag()+im*fKern[i].real();
  }

  DoInvFFT(func1);

}

// -----------------------------------------------------------------------------
// ~~~~ Shifts real vectors using above ShiftData function
// -----------------------------------------------------------------------------
template <class T>
inline void util::LArFFTW::ShiftData(std::vector<T> & input, double shift)
{
  DoFFT(input,fCompTemp);
  ShiftData(fCompTemp,shift);
  DoInvFFT(fCompTemp,input);

  return;
}

// -----------------------------------------------------------------------------
// ~~~~ Scheme for adding two signals which have an arbitrary relative
//      translation.  Shape1 is translated over shape2 and is replaced with the
//      sum, or the translated result if add = false
// -----------------------------------------------------------------------------
template <class T> inline void util::LArFFTW::AlignedSum(std::vector<T> & shape1,
							std::vector<T> & shape2,
							bool add)
{
  double shift = PeakCorrelation(shape1,shape2);

  ShiftData(shape1,shift);

  if(add)for(int i = 0; i < fSize; i++) shape1[i]+=shape2[i];

  return;
}

// -----------------------------------------------------------------------------
// ~~~~ Returns the length of the translation at which the correlation
//      of 2 signals is maximal.
// -----------------------------------------------------------------------------
template <class T> inline T util::LArFFTW::PeakCorrelation(std::vector<T> & shape1,
                                                                std::vector<T> & shape2)
{
  float chiSqr = std::numeric_limits<float>::max();
  float dchiSqr = std::numeric_limits<float>::max();
  const float chiCut   = 1e-3;
  float lambda  = 0.001;	// Marquardt damping parameter
  std::vector<float> p;

  std::vector<T> holder = shape1;
  Correlate(holder,shape2);

  int	maxT   = max_element(holder.begin(), holder.end())-holder.begin();
  float startT = maxT-fFitBins/2;
  int	offset = 0;

  for(int i = 0; i < fFitBins; i++) {
    if(startT+i < 0) offset=fSize;
    else if(startT+i > fSize) offset=-fSize;
    else offset = 0;
    if(holder[i+startT+offset]<=0.) {
      fConvHist[i]=0.;
    } else {
      fConvHist[i]=holder[i+startT+offset];
    }
  }

  p[0] = *max_element(fConvHist.begin(), fConvHist.end());
  p[1] = fFitBins/2;
  p[2] = fFitBins/2;
  float p1 = p[1];	// save initial p[1] guess

  int fitResult{-1};
  int trial=0;
  lambda=-1.;		// initialize lambda on first call
  do{
    fitResult=fMarqFitAlg->mrqdtfit(lambda, &p[0], &fConvHist[0], 3, fFitBins, chiSqr, dchiSqr);
    trial++;
    if(fitResult){
        mf::LogWarning("LArFFTW") << "Peak Correlation Fitting failed";
	break;
    }else if (trial>100){
  	break;
    }
  }
  while (fabs(dchiSqr) >= chiCut);
  if (!fitResult)p1=p[1]; // if fit succeeded, use fit result

  return p1 + 0.5 + startT;
}
#endif
