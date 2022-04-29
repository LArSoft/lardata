#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
// vim: set sw=2 expandtab :

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandardDataFor.h"

#include "TFile.h"
// Must precede include of art_root_io/detail/readMetadata.h
#include "TTree.h"

#include "art_root_io/Inputfwd.h"
#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "art_root_io/detail/readMetadata.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include "canvas/Persistency/Provenance/FileFormatVersion.h"
#include "canvas/Persistency/Provenance/ParameterSetMap.h"
#include "canvas/Persistency/Provenance/rootNames.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include <bitset>
#include <memory>
#include <string>
#include <vector>

using namespace std;

namespace detinfo {

  DetectorClocksServiceStandard::DetectorClocksServiceStandard(fhicl::ParameterSet const& pset,
                                                               art::ActivityRegistry& reg)
    : fClocks{pset}, fInheritClockConfig{pset.get<bool>("InheritClockConfig")}
  {
    reg.sPostOpenFile.watch(this, &DetectorClocksServiceStandard::postOpenFile);
    reg.sPreBeginRun.watch(this, &DetectorClocksServiceStandard::preBeginRun);
  }

  void
  DetectorClocksServiceStandard::preBeginRun(art::Run const& run)
  {
    // This callback probably is not necessary.
    fClocks.ApplyParams();
  }

  void
  DetectorClocksServiceStandard::postOpenFile(string const& filename)
  {
    if (!fInheritClockConfig) { return; }
    if (filename.empty()) { return; }
    std::unique_ptr<TFile> file{TFile::Open(filename.c_str(), "READ")};
    if (!file || file->IsZombie() || !file->IsOpen()) { return; }
    std::unique_ptr<TTree> metaDataTree{
      file->Get<TTree>(art::rootNames::metaDataTreeName().c_str())};
    if (metaDataTree == nullptr) {
      throw cet::exception("DetectorClocksServiceStandard",
                           "Input file does not contain a metadata tree!");
    }
    auto const fileFormatVersion =
      art::detail::readMetadata<art::FileFormatVersion>(metaDataTree.get());
    fhicl::ParameterSet ps;
    vector<string> const cfgName(fClocks.ConfigNames());
    vector<double> const cfgValue(fClocks.ConfigValues());
    bitset<kConfigTypeMax> config_set;
    vector<double> config_value(kConfigTypeMax, 0);

    auto count_configuration_changes =
      [&cfgName, &config_set, &config_value](fhicl::ParameterSet const& ps) {
        for (size_t i = 0; i < kConfigTypeMax; ++i) {
          auto const value_from_file = ps.get<double>(cfgName[i]);
          if (not config_set[i]) {
            config_value[i] = value_from_file;
            config_set[i] = true;
          }
          else if (config_value[i] != value_from_file) {
            throw cet::exception("DetectorClocksServiceStandard")
              << "Found historical value disagreement for " << cfgName[i] << " ... "
              << config_value[i] << " != " << value_from_file;
          }
        }
      };

    if (fileFormatVersion.value_ < 5) {
      art::ParameterSetMap psetMap;
      if (!art::detail::readMetadata(metaDataTree.get(), psetMap)) {
        throw cet::exception("DetectorClocksServiceStandard",
                             "Could not read ParameterSetMap from metadata tree!");
      }

      for (auto const& psEntry : psetMap) {
        fhicl::ParameterSet ps;
        ps = fhicl::ParameterSet::make(psEntry.second.pset_);
        if (!fClocks.IsRightConfig(ps)) { continue; }

        count_configuration_changes(ps);
      }
    }
    else {
      art::SQLite3Wrapper sqliteDB(file.get(), "RootFileDB");
      sqlite3_stmt* stmt{nullptr};
      sqlite3_prepare_v2(sqliteDB, "SELECT PSetBlob from ParameterSets;", -1, &stmt, nullptr);
      while (sqlite3_step(stmt) == SQLITE_ROW) {
        fhicl::ParameterSet ps;
        ps = fhicl::ParameterSet::make(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0)));
        if (!fClocks.IsRightConfig(ps)) { continue; }

        count_configuration_changes(ps);
      }
      sqlite3_finalize(stmt);
    }

    for (size_t i = 0; i < kConfigTypeMax; ++i) {
      if (not config_set[i]) continue;
      if (cfgValue[i] == config_value[i]) continue;

      cout << "Overriding configuration parameter " << cfgName[i] << " ... " << cfgValue[i]
           << " (fcl) => " << config_value[i] << " (data file)" << endl;
      fClocks.SetConfigValue(i, config_value[i]);
    }
    fClocks.ApplyParams();
  } // DetectorClocksServiceStandard::postOpenFile()
  
  
  DetectorClocksData DetectorClocksServiceStandard::DataFor
    (art::Event const& e) const
    { return detinfo::detectorClocksStandardDataFor(fClocks, e); }
  

} // namespace detinfo

