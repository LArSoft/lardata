////////////////////////////////////////////////////////////////////////
/// \file LArFFT.h
///
/// Utility FFT functions
///
/// \author  Brian Page
////////////////////////////////////////////////////////////////////////
#ifndef LARFFT_H
#define LARFFT_H

#include "TComplex.h"
#include "TFFTRealComplex.h"
#include "TFFTComplexReal.h"
#include "TF1.h"
#include "TH1D.h"
#include <vector>
#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

///General LArSoft Utilities
namespace util{
    class LArFFT {
    public:
      LArFFT(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~LArFFT();


      template <class T> void         DoFFT(std::vector<T> & input,
					    std::vector<TComplex> & output);

      template <class T> void         DoInvFFT(std::vector<TComplex> & input,
					       std::vector<T> & output);

      template <class T> void        Deconvolute(std::vector<T> & input,
						 std::vector<T> & respFunc);

      template <class T> void        Deconvolute(std::vector<T> & input,
						 std::vector<TComplex> & kern);

      template <class T> void        Convolute(std::vector<T> & input,
					       std::vector<T> & respFunc);

      template <class T> void        Convolute(std::vector<T> & input,
					       std::vector<TComplex> & kern);

      template <class T> void        Correlate(std::vector<T> & input,
					       std::vector<T> & respFunc);

      template <class T> void        Correlate(std::vector<T> & input,
					       std::vector<TComplex> & kern);

      template <class T> void        AlignedSum(std::vector<T> & input,
						std::vector<T> &output,
						bool add = true);

                         void        ShiftData(std::vector<TComplex> & input,
					       double shift);

      template <class T> void        ShiftData(std::vector<T> & input,
					       double shift);

      template <class T> T           PeakCorrelation(std::vector<T> &shape1,
						     std::vector<T> &shape2);

      int   FFTSize()          const { return fSize; }
      std::string FFTOptions() const { return fOption; }
      int FFTFitBins()         const { return fFitBins; }

      void ReinitializeFFT(int, std::string, int);

	private:

      int                    fSize;       //size of transform
      int                    fFreqSize;   //size of frequency space
      std::string            fOption;     //FFTW setting
      int                    fFitBins;    //Bins used for peak fit
      TF1                   *fPeakFit;    //Gaussian peak function
      TH1D                  *fConvHist;   //Fit data histogram
      std::vector<TComplex>  fCompTemp;   //temporary complex data
      std::vector<TComplex>  fKern;       //transformed response function

      TFFTRealComplex       *fFFT;        ///< object to do FFT
      TFFTComplexReal       *fInverseFFT; ///< object to do Inverse FF

      void InitializeFFT();
      void resetSizePerRun(art::Run const&);

    }; // class LArFFT

} //namespace util

// "Forward" Fourier Transform
//--------------------------------------------------------
template <class T> inline void util::LArFFT::DoFFT(std::vector<T> & input,
						   std::vector<TComplex> & output)
{
  double real      = 0.;  //real value holder
  double imaginary = 0.;  //imaginary value hold

  // set the points
  for(size_t p = 0; p < input.size(); ++p)
    fFFT->SetPoint(p, input[p]);

  fFFT->Transform();

  for(int i = 0; i < fFreqSize; ++i){
    fFFT->GetPointComplex(i, real, imaginary);
    output[i]=TComplex(real, imaginary);
  }

  return;
}

//Inverse Fourier Transform
//-------------------------------------------------
template <class T> inline void util::LArFFT::DoInvFFT(std::vector<TComplex> & input,
						      std::vector<T> & output)
{
  for(int i = 0; i < fFreqSize; ++i)
    fInverseFFT->SetPointComplex(i, input[i]);

  fInverseFFT->Transform();
  double factor = 1.0/(double) fSize;

  for(int i = 0; i < fSize; ++i)
    output[i] = factor*fInverseFFT->GetPointReal(i,false);

  return;
}

//Deconvolution scheme taking all time-domain
//information
//--------------------------------------------------
template <class T> inline void util::LArFFT::Deconvolute(std::vector<T> & input,
							 std::vector<T> & respFunction)
{
  DoFFT(respFunction, fKern);
  DoFFT(input, fCompTemp);

  for(int i = 0; i < fFreqSize; i++)
    fCompTemp[i]/=fKern[i];

  DoInvFFT(fCompTemp, input);

  return;
}

//Deconvolution scheme using an already transformed
//response function
//saves cpu time if same response function is used
//for many consecutive transforms
//--------------------------------------------------
template <class T> inline void util::LArFFT::Deconvolute(std::vector<T> & input,
							 std::vector<TComplex> & kern)
{
  DoFFT(input, fCompTemp);

  for(int i = 0; i < fFreqSize; i++)
    fCompTemp[i]/=kern[i];

  DoInvFFT(fCompTemp, input);

  return;
}

//Convolution scheme taking all time-domain
//information
//--------------------------------------------------
template <class T> inline void util::LArFFT::Convolute(std::vector<T> & shape1,
						       std::vector<T> & shape2)
{
  DoFFT(shape1, fKern);
  DoFFT(shape2, fCompTemp);

  for(int i = 0; i < fFreqSize; i++)
    fCompTemp[i]*=fKern[i];

  DoInvFFT(fCompTemp, shape1);

  return;
}

//Convolution scheme using an already transformed
//response function
//saves cpu time if same response function is used
//for many consecutive transforms
//--------------------------------------------------
template <class T> inline void util::LArFFT::Convolute(std::vector<T> & input,
						       std::vector<TComplex> & kern)
{
  DoFFT(input, fCompTemp);

  for(int i = 0; i < fFreqSize; i++)
    fCompTemp[i]*=kern[i];

  DoInvFFT(fCompTemp, input);

  return;
}

//Correlation taking all time domain data
//--------------------------------------------------
template <class T> inline void util::LArFFT::Correlate(std::vector<T> & shape1,
						       std::vector<T> & shape2)
{
  DoFFT(shape1, fKern);
  DoFFT(shape2, fCompTemp);

  for(int i = 0; i < fFreqSize; i++)
    fCompTemp[i]*=TComplex::Conjugate(fKern[i]);

  DoInvFFT(fCompTemp, shape1);

  return;
}

//Convolution scheme using an already transformed
//response function
//saves cpu time if same response function is used
//for many consecutive transforms
//--------------------------------------------------
template <class T> inline void util::LArFFT::Correlate(std::vector<T> & input,
						       std::vector<TComplex> & kern)
{
  DoFFT(input, fCompTemp);

  for(int i = 0; i < fFreqSize; i++)
    fCompTemp[i]*=TComplex::Conjugate(kern[i]);

  DoInvFFT(fCompTemp, input);

  return;
}

//Scheme for adding two signals which have an arbitrary
//relative translation.  Shape1 is translated over shape2
//and is replaced with the sum, or the translated result
//if add = false
//--------------------------------------------------
template <class T> inline void util::LArFFT::AlignedSum(std::vector<T> & shape1,
							std::vector<T> & shape2,
							bool add)
{
  double shift = PeakCorrelation(shape1,shape2);

  ShiftData(shape1,shift);

  if(add)for(int i = 0; i < fSize; i++) shape1[i]+=shape2[i];

  return;
}

//Shifts real vectors using above function
//--------------------------------------------------
template <class T> inline void util::LArFFT::ShiftData(std::vector<T> & input,
						       double shift)
{
  DoFFT(input,fCompTemp);
  ShiftData(fCompTemp,shift);
  DoInvFFT(fCompTemp,input);

  return;
}

//Returns the length of the translation at which the correlation
//of 2 signals is maximal.
//--------------------------------------------------
template <class T> inline T util::LArFFT::PeakCorrelation(std::vector<T> & shape1,
							  std::vector<T> & shape2)
{
  fConvHist->Reset("ICE");
  std::vector<T> holder = shape1;
  Correlate(holder,shape2);

  int   maxT   = max_element(holder.begin(), holder.end())-holder.begin();
  float startT = maxT-fFitBins/2;
  int   offset = 0;

  for(int i = 0; i < fFitBins; i++) {
    if(startT+i < 0) offset=fSize;
    else if(startT+i > fSize) offset=-fSize;
    else offset = 0;
    fConvHist->Fill(i,holder[i+startT+offset]);
  }

  fPeakFit->SetParameters(fConvHist->GetMaximum(),fFitBins/2,fFitBins/2);
  fConvHist->Fit(fPeakFit,"QWNR","",0,fFitBins);
  return fPeakFit->GetParameter(1)+startT;
}

DECLARE_ART_SERVICE(util::LArFFT, LEGACY)
#endif // LARFFT_H
