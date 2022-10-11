/**
 * @file   DumpChargedSpacePoints_module.cc
 * @brief  Dumps on screen the content of space points and associated charge.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 22nd, 2017
 */

// LArSoft includes
#include "lardata/RecoBaseProxy/ChargedSpacePoints.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

// support libraries
#include "fhiclcpp/types/Atom.h" // also pulls in fhicl::Name and fhicl::Comment
#include "messagefacility/MessageLogger/MessageLogger.h"

// C//C++ standard libraries
#include <string>

namespace recob {

  /**
   * @brief Prints the content of all the space points and charge on screen.
   *
   * This analyser prints the content of all the space points into the
   * LogInfo/LogVerbatim stream.
   *
   * The space point and charge data products must fulfil the requirements of
   * the `proxy::ChargedSpacePoints` proxy.
   *
   * Configuration parameters
   * =========================
   *
   * - *SpacePointLabel* (`art::InputTag`, mandatory): label of the
   *   producer used to create the `recob::SpacePoint` _and_
   *   `recob::PointCharge` collections to be dumped
   * - *OutputCategory* (string, default: "DumpChargedSpacePoints"): the
   *   category used for the output (useful for filtering)
   *
   */
  class DumpChargedSpacePoints : public art::EDAnalyzer {
  public:
    /// Configuration parameters
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> SpacePointTag{
        Name("SpacePointLabel"),
        Comment("label of the producer used to create"
                " the recob::SpacePoint collection to be dumped")};
      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("the category used for the output (useful for filtering)"),
        "DumpChargedSpacePoints" /* default value */
      };

    }; // struct Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// Constructor.
    explicit DumpChargedSpacePoints(Parameters const& config);

    /// Does the printing.
    virtual void analyze(art::Event const& event) override;

  private:
    art::InputTag fInputTag;     ///< Input tag of the SpacePoint product.
    std::string fOutputCategory; ///< Category for LogInfo output.

  }; // class DumpChargedSpacePoints

} // namespace recob

//==============================================================================
//===  Implementation section
//==============================================================================

//----------------------------------------------------------------------------
recob::DumpChargedSpacePoints::DumpChargedSpacePoints(art::EDAnalyzer::Table<Config> const& config)
  : EDAnalyzer(config)
  , fInputTag(config().SpacePointTag())
  , fOutputCategory(config().OutputCategory())
{}

//----------------------------------------------------------------------------
void recob::DumpChargedSpacePoints::analyze(art::Event const& event)
{

  //
  // collect all the available information
  //
  // fetch the data to be dumped on screen
  auto const& points = proxy::getChargedSpacePoints(event, fInputTag);

  size_t const nPoints = points.size();
  mf::LogVerbatim log(fOutputCategory);
  log << "The event " << event.id() << " contains " << nPoints << " space points from '"
      << fInputTag.encode() << "'";

  for (auto const& point : points) {

    log << "\n [#" << point.index() << "] " << point.point() << " " << point.charge();

  } // for

  log << "\n"; // two empty lines

} // DumpChargedSpacePoints::analyze()

//----------------------------------------------------------------------------
DEFINE_ART_MODULE(recob::DumpChargedSpacePoints)

//----------------------------------------------------------------------------
