/**
 * @file   DumpOpHits_module.cc
 * @brief  Dumps on screen the content of the reconstructed optical hits.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 10, 2023
 */

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "lardataobj/RecoBase/OpHit.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C//C++ standard libraries
#include <ostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
namespace ophit {
  class DumpOpHits;
}

/**
 * @brief Prints the content of all the hits on screen
 *
 * This analyser prints the content of all the hits into the
 * `LogInfo`/`LogVerbatim` stream.
 *
 * Configuration parameters
 * =========================
 *
 * * **OpHitModuleLabel** (input tag): tag of the `recob::OpHit` collection data
 *    product to be dumped.
 * * **OutputCategory** (string, default: `"DumpOpHits"`): the category
 *   used for the output (useful for filtering).
 *
 */
class ophit::DumpOpHits : public art::EDAnalyzer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> OpHitModuleLabel{
      Name{"OpHitModuleLabel"},
      Comment{"tag of the recob::OpHit collection to be dumped"}};

    fhicl::Atom<std::string> OutputCategory{
      Name("OutputCategory"),
      Comment("the messagefacility category used for the output"),
      "DumpOpHits"};

  }; // Config

  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Default constructor
  explicit DumpOpHits(Parameters const& config);

  /// Does the printing
  void analyze(const art::Event& evt) override;

private:
  art::InputTag const fOpHitModuleTag; ///< Optical hit data product tag.
  std::string const fOutputCategory;   ///< Category for `mf::LogInfo` output.

}; // class ophit::DumpOpHits

// -----------------------------------------------------------------------------
// ---  implementation
// -----------------------------------------------------------------------------
ophit::DumpOpHits::DumpOpHits(Parameters const& config)
  : art::EDAnalyzer(config)
  , fOpHitModuleTag(config().OpHitModuleLabel())
  , fOutputCategory(config().OutputCategory())
{
  consumes<std::vector<recob::OpHit>>(fOpHitModuleTag);
}

//------------------------------------------------------------------------------
void ophit::DumpOpHits::analyze(art::Event const& event)
{

  // fetch the data to be dumped on screen
  auto const& OpHits = event.getProduct<std::vector<recob::OpHit>>(fOpHitModuleTag);

  mf::LogVerbatim(fOutputCategory) << "Event " << event.id() << " contains " << OpHits.size()
                                   << " '" << fOpHitModuleTag.encode() << "' optical hits.";

  for (auto const& [iHit, hit] : util::enumerate(OpHits))
    mf::LogVerbatim(fOutputCategory) << "OpHit #" << iHit << ": " << hit;

} // ophit::DumpOpHits::analyze()

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(ophit::DumpOpHits)

//------------------------------------------------------------------------------
