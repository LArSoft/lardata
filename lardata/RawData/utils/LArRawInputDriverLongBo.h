////////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriver.h
/// \brief Source to convert raw binary files to root files for Long Bo TPC
///
/// \author  brebel@fnal.gov, soderber@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "canvas/Persistency/Provenance/SubRunID.h"

#include <string>
#include <vector>

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

///Conversion of binary data to root files
namespace lris {
  class LArRawInputDriverLongBo;
}

class lris::LArRawInputDriverLongBo {
  /// Class to fill the constraints on a template argument to the class,
  /// FileReaderSource
public:
  // Required constructor
  LArRawInputDriverLongBo(fhicl::ParameterSet const& pset,
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
  typedef std::vector<std::string> stringvec_t;

  art::SourceHelper const& principalMaker_;
  std::string currentDir_;
  stringvec_t inputfiles_;
  stringvec_t::const_iterator nextfile_;
  stringvec_t::const_iterator filesdone_;
  art::SubRunID currentSubRunID_;
}; // LArRawInputDriverLongBo
