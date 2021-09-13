#include "lardata/Utilities/LArFFTWPlan.h"

using std::string;
std::mutex util::LArFFTWPlan::mutex_;

util::LArFFTWPlan::LArFFTWPlan(int transformSize, const std::string &option)
  : fSize    (transformSize)
  , fOption  (option){

  std::lock_guard<std::mutex> lock(mutex_);

  fFreqSize = fSize/2+1;
  fN = new int[1];
  fN[0] = fSize;

  fIn = fftw_malloc(sizeof(double)*fSize);
  fOut= fftw_malloc(sizeof(fftw_complex)*fFreqSize);
  fPlan = (void*)fftw_plan_dft_r2c(1, fN, (double*)fIn, (fftw_complex*)fOut, MapFFTWOption());

  rIn = fftw_malloc(sizeof(fftw_complex)*fFreqSize);
  rOut= fftw_malloc(sizeof(double)*fSize);
  rPlan = (void*)fftw_plan_dft_c2r(1, fN, (fftw_complex*)rIn, (double*)rOut, MapFFTWOption());
}

util::LArFFTWPlan::~LArFFTWPlan()
{
  fftw_destroy_plan((fftw_plan)fPlan);
  fPlan = 0;
  fftw_free(fIn);
  fIn = 0;
  fftw_free((fftw_complex*)fOut);
  fOut = 0;

  fftw_destroy_plan((fftw_plan)rPlan);
  rPlan = 0;
  fftw_free((fftw_complex*)rIn);
  rIn = 0;
  fftw_free(rOut);
  rOut = 0;

  delete [] fN;
  fN = 0;
}

unsigned int util::LArFFTWPlan::MapFFTWOption()
{
  std::transform(fOption.begin(), fOption.end(),fOption.begin(), ::toupper);
  if (fOption.find("ES")!=string::npos)
     return FFTW_ESTIMATE;
  if (fOption.find("M")!=string::npos)
     return FFTW_MEASURE;
  if (fOption.find("P")!=string::npos)
     return FFTW_PATIENT;
  if (fOption.find("EX")!=string::npos)
     return FFTW_EXHAUSTIVE;
  return FFTW_ESTIMATE;
}
