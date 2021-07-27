////////////////////////////////////////////////////////////////////////
//
// \file LArFFT_plugin
//
//  This class simplifies implementation of Fourier transforms.
//  Because all data inputs and outputs are purely real,  the
//  transforms implemented in this way get a substantial performance
//  increase ~2x.
//
// \author pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/LArFFT.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

//-----------------------------------------------
util::LArFFT::LArFFT(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  : fSize(pset.get<int>("FFTSize", 0))
  , fOption(pset.get<std::string>("FFTOption"))
  , fFitBins(pset.get<int>("FitBins"))
{
  // Default to the readout window size if the user didn't input
  // a specific size
  if (fSize <= 0) {
    // Creating a service handle to DetectorPropertiesService not only
    // creates the service if it doesn't exist, it also guarantees
    // that its callbacks are invoked before any of LArFFT's callbacks
    // are invoked.
    fSize = art::ServiceHandle<detinfo::DetectorPropertiesService const>()
              ->DataForJob()
              .ReadOutWindowSize();
    reg.sPreBeginRun.watch(this, &util::LArFFT::resetSizePerRun);
  }
  InitializeFFT();
}

//-----------------------------------------------
void
util::LArFFT::resetSizePerRun(art::Run const&)
{
  fSize = art::ServiceHandle<detinfo::DetectorPropertiesService const>()
            ->DataForJob()
            .ReadOutWindowSize();
  ReinitializeFFT(fSize, fOption, fFitBins);
}

//-----------------------------------------------
void
util::LArFFT::InitializeFFT()
{
  int i;
  for (i = 1; i < fSize; i *= 2) {}
  fSize = i;
  fFreqSize = fSize / 2 + 1;

  // allocate and setup Transform objects
  fFFT = new TFFTRealComplex(fSize, false);
  fInverseFFT = new TFFTComplexReal(fSize, false);

  int dummy[1] = {0};
  // appears to be dummy argument from root page
  fFFT->Init(fOption.c_str(), -1, dummy);
  fInverseFFT->Init(fOption.c_str(), 1, dummy);

  fPeakFit = new TF1("fPeakFit", "gaus"); //allocate function used for peak fitting
  fConvHist = new TH1D("fConvHist",
                       "Convolution Peak Data",
                       fFitBins,
                       0,
                       fFitBins); //allocate histogram for peak fitting
  //allocate other data vectors
  fCompTemp.resize(fFreqSize);
  fKern.resize(fFreqSize);
}

//------------------------------------------------
util::LArFFT::~LArFFT()
{
  delete fFFT;
  delete fInverseFFT;
  delete fPeakFit;
  delete fConvHist;
}

//------------------------------------------------
void
util::LArFFT::ReinitializeFFT(int size, std::string option, int fitbins)
{
  //delete these, which will be remade
  delete fFFT;
  delete fInverseFFT;
  delete fPeakFit;
  delete fConvHist;

  //set members
  fSize = size;
  fOption = option;
  fFitBins = fitbins;

  //now initialize
  InitializeFFT();
}

//-------------------------------------------------
// For the sake of efficiency, as all transforms should
// be of the same size, all functions expect vectors
// to be of the desired size.
// Note that because the transforms are real to real
// there is a redundancy in the information in the
// complex part in the positive and negative
// components of the FFT, thus the size of the
// frequency vectors are input_fFreqSize
// --see the FFTW3 or Root docmentation for details

//According to the Fourier transform identity
//f(x-a) = Inverse Transform(exp(-2*Pi*i*a*w)F(w))
//--------------------------------------------------
void
util::LArFFT::ShiftData(std::vector<TComplex>& input, double shift)
{
  double factor = -2.0 * TMath::Pi() * shift / (double)fSize;

  for (int i = 0; i < fFreqSize; i++)
    input[i] *= TComplex::Exp(TComplex(0, factor * (double)i));

  return;
}

DEFINE_ART_SERVICE(util::LArFFT)
