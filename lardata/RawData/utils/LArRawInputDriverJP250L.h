////////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriverJP250L.h
/// \brief Source to convert JP250L files to LArSoft files
///
/// \author  eito@post.kek.jp, brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include <string>

namespace art {
  class EventPrincipal;
  class FileBlock;
  class ProductRegistryHelper;
  class RunPrincipal;
  class SourceHelper;
  class SubRunPrincipal;
}
namespace fhicl {
  class ParameterSet;
}

class TTree;

///Conversion of binary data to root files
namespace lris {
  class LArRawInputDriverJP250L;
}

class lris::LArRawInputDriverJP250L {
  /// Class to fill the constraints on a template argument to the class,
  /// FileReaderSource
public:
  // Required constructor
  LArRawInputDriverJP250L(fhicl::ParameterSet const& pset,
                          art::ProductRegistryHelper& helper,
                          art::SourceHelper const& pm);

  // Required by FileReaderSource:
  void closeCurrentFile();
  void readFile(std::string const& name, art::FileBlock*& fb);
  bool readNext(art::RunPrincipal* const& inR,
                art::SubRunPrincipal* const& inSR,
                art::RunPrincipal*& outR,
                art::SubRunPrincipal*& outSR,
                art::EventPrincipal*& outE);

private:
  // --- data members:
  art::SourceHelper const& principalMaker_;

  // added by E.Iwai
  TTree* m_eventTree;         ///< TTree containing information from each trigger
  unsigned int m_nEvent;      ///< number of triggers in the TTree
  unsigned int m_current;     ///< current entry in the TTree
  unsigned short m_runID;     ///< run ID, has to start from 1
  unsigned int m_unixtime;    ///< unix timestamp of the start of the run
  unsigned short m_nChannels; ///< number of channels in the detector
  unsigned short m_nSamples;  ///< number of time samples per channel
  unsigned short* m_data;     ///< the ADC of each time sample for each channel

}; // LArRawInputDriverJP250L
