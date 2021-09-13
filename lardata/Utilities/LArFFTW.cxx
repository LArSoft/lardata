#include "lardata/Utilities/LArFFTW.h"

using std::string;

util::LArFFTW::LArFFTW(int transformSize, const void* fplan, const void* rplan, int fitbins)
  : fSize    (transformSize)
  , fPlan    (fplan)
  , rPlan    (rplan)
  , fFitBins (fitbins)
{

  fFreqSize = fSize/2+1;

  // ... Real-Complex
  fIn = fftw_malloc(sizeof(double)*fSize);
  fOut= fftw_malloc(sizeof(fftw_complex)*fFreqSize);

  // ... Complex-Real
  rIn = fftw_malloc(sizeof(fftw_complex)*fFreqSize);
  rOut= fftw_malloc(sizeof(double)*fSize);

  // ... allocate other data vectors
  fCompTemp.resize(fFreqSize);
  fKern.resize(fFreqSize);
  fConvHist.resize(fFitBins);
}

util::LArFFTW::~LArFFTW()
{
  fPlan = 0;
  fftw_free(fIn);
  fIn = 0;
  fftw_free((fftw_complex*)fOut);
  fOut = 0;

  rPlan = 0;
  fftw_free((fftw_complex*)rIn);
  rIn = 0;
  fftw_free(rOut);
  rOut = 0;
}

// According to the Fourier transform identity
// f(x-a) = Inverse Transform(exp(-2*Pi*i*a*w)F(w))
// -----------------------------------------------------------------------------
void util::LArFFTW::ShiftData(ComplexVector & input, double shift)
{
  double factor = -2.0*std::acos(-1)*shift/(double)fSize;

  for(int i = 0; i < fFreqSize; i++){
    input[i] *= std::exp(std::complex<double>(0,factor*(double)i));
  }

  return;
}
