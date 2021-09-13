#ifndef LARFFTWPLAN_H
#define LARFFTWPLAN_H

// C/C++ standard libraries
#include <string>
#include <algorithm>
#include <mutex>

#include "fftw3.h"

namespace util {

class LArFFTWPlan {

  public:
    LArFFTWPlan(int transformSize, const std::string &option);
    ~LArFFTWPlan();
    void *fPlan;
    void *rPlan;
    void *fIn;
    void *fOut;
    void *rIn;
    void *rOut;

  private:
    static std::mutex mutex_;
    int fSize;		// size of transform
    int fFreqSize;	// size of frequency space
    int *fN;
    std::string fOption;	// FFTW setting

    unsigned int MapFFTWOption();

};

}  // end namespace util

#endif
