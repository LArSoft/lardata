/**
 * @file   DumpHits_module.cc
 * @brief  Dumps on screen the content of the hits
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   Match 9th, 2015
 */

// C//C++ standard libraries
#include <string>

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Utilities/InputTag.h"

// ... plus see below ...

namespace hit {

  /**
   * @brief Prints the content of all the hits on screen
   *
   * This analyser prints the content of all the hits into the
   * LogInfo/LogVerbatim stream.
   *
   * Configuration parameters
   * =========================
   *
   * - *HitModuleLabel* (string): label of the producer used to create the
   *   recob::Hit collection
   * - *OutputCategory* (string, default: "DumpHits"): the category
   *   used for the output (useful for filtering)
   * - *CheckWireAssociation* (string, default: false): if set, verifies
   *   that the associated wire are on the same channel as the hit
   * - *CheckRawDigitAssociation* (string, default: false): if set, verifies
   *   that the associated raw digits are on the same channel as the hit
   *
   */
  class DumpHits : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> HitModuleLabel{
        Name("HitModuleLabel"),
        Comment("tag of the producer used to create the recob::Hit collection")};

      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("the messagefacility category used for the output"),
        "DumpHits"};

      fhicl::Atom<bool> CheckRawDigitAssociation{
        Name("CheckRawDigitAssociation"),
        Comment("verify the associated raw digits are on the same channel as the hit"),
        false}; // CheckRawDigitAssociation

      fhicl::Atom<bool> CheckWireAssociation{
        Name("CheckWireAssociation"),
        Comment("verify the associated wire is on the same channel as the hit"),
        false}; // CheckWireAssociation

    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// Default constructor
    explicit DumpHits(Parameters const& config);

    /// Does the printing
    void analyze(const art::Event& evt);

  private:
    art::InputTag fHitsModuleLabel; ///< name of module that produced the hits
    std::string fOutputCategory;    ///< category for LogInfo output
    bool bCheckRawDigits;           ///< check associations with raw digits
    bool bCheckWires;               ///< check associations with wires

  }; // class DumpHits

} // namespace hit

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
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

namespace hit {

  //-------------------------------------------------
  DumpHits::DumpHits(Parameters const& config)
    : EDAnalyzer(config)
    , fHitsModuleLabel(config().HitModuleLabel())
    , fOutputCategory(config().OutputCategory())
    , bCheckRawDigits(config().CheckRawDigitAssociation())
    , bCheckWires(config().CheckWireAssociation())
  {}

  //-------------------------------------------------
  void DumpHits::analyze(const art::Event& evt)
  {

    // fetch the data to be dumped on screen
    auto Hits = evt.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);

    mf::LogInfo(fOutputCategory) << "The event contains " << Hits->size() << " '"
                                 << fHitsModuleLabel.encode() << "' hits";

    std::unique_ptr<art::FindOne<raw::RawDigit>> HitToRawDigit;
    if (bCheckRawDigits) {
      HitToRawDigit.reset(new art::FindOne<raw::RawDigit>(Hits, evt, fHitsModuleLabel));
      if (!HitToRawDigit->isValid()) {
        throw art::Exception(art::errors::ProductNotFound)
          << "DumpHits: can't find associations between raw digits and hits from '"
          << fHitsModuleLabel << "'";
      }
    } // if check raw digits

    std::unique_ptr<art::FindOne<recob::Wire>> HitToWire;
    if (bCheckWires) {
      HitToWire.reset(new art::FindOne<recob::Wire>(Hits, evt, fHitsModuleLabel));
      if (!HitToWire->isValid()) {
        throw art::Exception(art::errors::ProductNotFound)
          << "DumpHits: can't find associations between wires and hits from '" << fHitsModuleLabel
          << "'";
      }
    } // if check wires

    unsigned int iHit = 0;
    for (const recob::Hit& hit : *Hits) {

      // print a header for the cluster
      mf::LogVerbatim(fOutputCategory) << "Hit #" << iHit << ": " << hit;

      if (HitToRawDigit) {
        raw::ChannelID_t assChannelID = HitToRawDigit->at(iHit).ref().Channel();
        if (assChannelID != hit.Channel()) {
          throw art::Exception(art::errors::DataCorruption)
            << "Hit #" << iHit << " on channel " << hit.Channel()
            << " is associated with raw digit on channel " << assChannelID << "!!";
        } // mismatch
      }   // raw digit check

      if (HitToWire) {
        raw::ChannelID_t assChannelID = HitToWire->at(iHit).ref().Channel();
        if (assChannelID != hit.Channel()) {
          throw art::Exception(art::errors::DataCorruption)
            << "Hit #" << iHit << " on channel " << hit.Channel()
            << " is associated with wire on channel " << assChannelID << "!!";
        } // mismatch
      }   // wire check

      ++iHit;
    } // for hits

  } // DumpHits::analyze()

  DEFINE_ART_MODULE(DumpHits)

} // namespace hit
