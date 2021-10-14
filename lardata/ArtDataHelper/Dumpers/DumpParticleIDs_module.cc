/**
 * @file   DumpParticleIDs_module.cc
 * @brief  Dump ParticleID objects.
 * @author H. Greenlee
 * @date   Oct. 14, 2021
 */

// C//C++ standard libraries
#include <string>

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Comment.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Utilities/InputTag.h"

// ... plus see below ...

namespace pid {

  /**
   * @brief Prints the content of all the partidle IDs on screen
   *
   * This analyser prints the content of all the particle IDs into the
   * LogInfo/LogVerbatim stream.
   *
   * Configuration parameters
   * =========================
   *
   * - *ParticleIDModuleLabel* (string): label of the producer used to create the
   *   anab::ParticleID collection
   * - *OutputCategory* (string, default: "DumpParticleIDs"): the category
   *   used for the output (useful for filtering)
   *
   */
  class DumpParticleIDs: public art::EDAnalyzer {
      public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> ParticleIDModuleLabel{
        Name("ParticleIDModuleLabel"),
        Comment("tag of the producer used to create the anab::ParticleID collection")
        };

      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("the messagefacility category used for the output"),
        "DumpParticleIDs"
        };

    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;


    /// Default constructor
    explicit DumpParticleIDs(Parameters const& config);

    /// Does the printing
    void analyze (const art::Event& evt);

      private:

    art::InputTag fParticleIDsModuleLabel; ///< name of module that produced the pids
    std::string fOutputCategory;    ///< category for LogInfo output

  }; // class DumpParticleIDs

} // namespace pid


//------------------------------------------------------------------------------
//---  module implementation
//---
// C//C++ standard libraries
#include <memory> // std::unique_ptr<>

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Framework/Principal/Handle.h"

// LArSoft includes
#include "lardataobj/AnalysisBase/ParticleID.h"


namespace pid {

  //-------------------------------------------------
  DumpParticleIDs::DumpParticleIDs(Parameters const& config)
    : EDAnalyzer              (config)
    , fParticleIDsModuleLabel (config().ParticleIDModuleLabel())
    , fOutputCategory         (config().OutputCategory())
    {}


  //-------------------------------------------------
  void DumpParticleIDs::analyze(const art::Event& evt) {

    // fetch the data to be dumped on screen
    auto ParticleIDs = evt.getValidHandle<std::vector<anab::ParticleID>>(fParticleIDsModuleLabel);

    mf::LogInfo(fOutputCategory)
      << "The event contains " << ParticleIDs->size() << " '"
      << fParticleIDsModuleLabel.encode() << "' particle IDs";

    unsigned int ipid = 0;
    for (const anab::ParticleID& pid: *ParticleIDs) {

      // print a header for the cluster
      mf::LogVerbatim log(fOutputCategory);
      log << "ParticleID #" << ipid << '\n';
      if(pid.PlaneID()) {
	log << "  Plane ID    = " << pid.PlaneID() << '\n'
	    << "  NDF         = " << pid.Ndf() << '\n'
	    << "  MinChi2     = " << pid.MinChi2() << '\n'
	    << "  DeltaChi2   = " << pid.DeltaChi2() << '\n'
	    << "  Chi2Proton  = " << pid.Chi2Proton() << '\n'
	    << "  Chi2Kaon    = " << pid.Chi2Kaon() << '\n'
	    << "  Chi2Muon    = " << pid.Chi2Muon() << '\n'
	    << "  MissingE    = " << pid.MissingE() << '\n'
	    << "  MissingEavg = " << pid.MissingEavg() << '\n'
	    << "  PIDA        = " << pid.PIDA() << std::endl;
      }
      else
	log << "  Invalid" << '\n';

      ++ipid;
    } // for pids

  } // DumpParticleIDs::analyze()

  DEFINE_ART_MODULE(DumpParticleIDs)

} // namespace pid
