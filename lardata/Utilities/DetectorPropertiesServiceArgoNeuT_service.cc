////////////////////////////////////////////////////////////////////////
//
//  \file DetectorPropertiesServiceArgoNeuT_service.cc
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/Utilities/DatabaseUtil.h"
#include "lardata/Utilities/DetectorPropertiesServiceArgoNeuT.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Art includes
#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "canvas/Utilities/Exception.h"

#include "TFile.h"
#include "TTree.h"

namespace util {

  //--------------------------------------------------------------------
  DetectorPropertiesServiceArgoNeuT::DetectorPropertiesServiceArgoNeuT(
    fhicl::ParameterSet const& pset,
    art::ActivityRegistry& reg)
    : fDetProp{pset}
  {
    fInheritNumberTimeSamples = pset.get<bool>("InheritNumberTimeSamples", false);

    // Save the parameter set.
    fPS = pset;

    if (pset.has_key("InheritTriggerOffset"))
      throw cet::exception(__FUNCTION__)
        << "InheritTriggerOffset is a deprecated fcl parameter for "
           "DetectorPropertiesServiceArgoNeuT!";

    // Register for callbacks.
    reg.sPostOpenFile.watch(this, &DetectorPropertiesServiceArgoNeuT::postOpenFile);
  }

  //--------------------------------------------------------------------
  //  Callback called after input file is opened.

  void
  DetectorPropertiesServiceArgoNeuT::postOpenFile(const std::string& filename)
  {
    // Use this method to figure out whether to inherit configuration
    // parameters from previous jobs.
    //
    // There is no way currently to correlate parameter sets saved in
    // sqlite RootFileDB with process history (from MetaData tree).
    // Therefore, we use the approach of scanning every historical
    // parameter set in RootFileDB, and finding all parameter sets
    // that appear to be DetectorPropertiesServiceArgoNeuT configurations.  If
    // all historical parameter sets are in agreement about the value of an
    // inherited parameter, then we accept the historical value, print a
    // message, and override the configuration parameter.  In cases where the
    // historical configurations are not in agreement about the value of an
    // inherited parameter, we ignore any historical parameter values that are
    // the same as the current configured value of the parameter (that is, we
    // resolve the conflict in favor of parameters values that are different
    // than the current configuration).  If two or more historical values differ
    // from the current configuration, throw an exception. Note that it is
    // possible to give precendence to the current configuration by disabling
    // inheritance for that configuration parameter.

    // Don't do anything if no parameters are supposed to be inherited.

    if (!fInheritNumberTimeSamples) return;

    // The only way to access art service metadata from the input file
    // is to open it as a separate TFile object.  Do that now.

    if (empty(filename)) return;

    std::unique_ptr<TFile> file{TFile::Open(filename.c_str(), "READ")};
    if (file && !file->IsZombie() && file->IsOpen()) {

      // Open the sqlite datatabase.

      art::SQLite3Wrapper sqliteDB(file.get(), "RootFileDB");

      // Loop over all stored ParameterSets.

      unsigned int iNumberTimeSamples = 0; // Combined value of NumberTimeSamples.
      unsigned int nNumberTimeSamples = 0; // Number of NumberTimeSamples parameters seen.

      sqlite3_stmt* stmt = nullptr;
      sqlite3_prepare_v2(sqliteDB, "SELECT PSetBlob from ParameterSets;", -1, &stmt, nullptr);
      while (sqlite3_step(stmt) == SQLITE_ROW) {
        fhicl::ParameterSet ps;
        ps = fhicl::ParameterSet::make(reinterpret_cast<char const*>(sqlite3_column_text(stmt, 0)));
        // Is this a DetectorPropertiesServiceArgoNeuT parameter set?

        if (isDetectorPropertiesServiceArgoNeuT(ps)) {
          // Check NumberTimeSamples

          unsigned int newNumberTimeSamples = ps.get<unsigned int>("NumberTimeSamples");

          // Ignore parameter values that match the current configuration.

          if (newNumberTimeSamples != fPS.get<unsigned int>("NumberTimeSamples")) {
            if (nNumberTimeSamples == 0)
              iNumberTimeSamples = newNumberTimeSamples;
            else if (newNumberTimeSamples != iNumberTimeSamples) {
              throw cet::exception("DetectorPropertiesServiceArgoNeuT")
                << "Historical values of NumberTimeSamples do not agree: " << iNumberTimeSamples
                << " " << newNumberTimeSamples << "\n";
            }
            ++nNumberTimeSamples;
          }
        }
      }

      // Done looping over parameter sets.
      // Now decide which parameters we will actually override.

      if (nNumberTimeSamples != 0 && iNumberTimeSamples != fDetProp.NumberTimeSamples()) {
        mf::LogInfo("DetectorPropertiesServiceArgoNeuT")
          << "Overriding configuration parameter NumberTimeSamples using historical value.\n"
          << "  Configured value:        " << fDetProp.NumberTimeSamples() << "\n"
          << "  Historical (used) value: " << iNumberTimeSamples << "\n";
        fDetProp.SetNumberTimeSamples(iNumberTimeSamples);
      }
    }
  }

  //--------------------------------------------------------------------
  //  Determine whether a parameter set is a DetectorPropertiesServiceArgoNeuT
  //  configuration.

  bool
  DetectorPropertiesServiceArgoNeuT::isDetectorPropertiesServiceArgoNeuT(
    const fhicl::ParameterSet& ps)
  {
    // This method uses heuristics to determine whether the parameter
    // set passed as argument is a DetectorPropertiesServiceArgoNeuT
    // configuration parameter set.

    std::string s;
    double d;
    int i;
    unsigned int u;

    return !ps.get_if_present("module_label", s) && ps.get_if_present("TriggerOffset", i) &&
           ps.get_if_present("SamplingRate", d) && ps.get_if_present("NumberTimeSamples", u) &&
           ps.get_if_present("ReadOutWindowSize", u);
  }

} // namespace

DEFINE_ART_SERVICE_INTERFACE_IMPL(util::DetectorPropertiesServiceArgoNeuT,
                                  detinfo::DetectorPropertiesService)
