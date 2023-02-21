/**
 * @file   DumpOpFlashes_module.cc
 * @brief  Dumps on screen the content of the reconstructed optical flashes.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 20, 2023
 */

// LArSoft libraries
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/values.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C//C++ standard libraries
#include <cassert>
#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
namespace ophit {
  class DumpOpFlashes;
}

/**
 * @brief Prints the content of all the flashes on screen
 *
 * This analyser prints the content of all the flashes into the
 * `LogInfo`/`LogVerbatim` stream.
 *
 *
 * Configuration parameters
 * =========================
 *
 * * **OpFlashModuleLabel** (input tag): tag of the `recob::OpFlash` collection
 *    data product to be dumped.
 * * **PrintOpHitAssociations** (flag, default: `true`): also print a list of
 *    optical hits associated to each flash.
 * * **OutputCategory** (string, default: `"DumpOpFlashes"`): the category
 *   used for the output (useful for filtering).
 *
 *
 * Input data products
 * ====================
 *
 * * `std::vector<recob::OpFlash>` (`OpFlashModuleLabel`): the flash collection
 *     being dumped
 * * `art::Assns<recob::OpFlash, recob::OpHit>` (`OpFlashModuleLabel`): the
 *     associated optical hits. The index of each hit will be printed.
 *     In addition, if the optical hit data product is available, it will be
 *     fully read and channel, start time and p.e. of the hits will be printed.
 *
 */
class ophit::DumpOpFlashes : public art::EDAnalyzer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> OpFlashModuleLabel{
      Name{"OpFlashModuleLabel"},
      Comment{"tag of the recob::OpFlash collection to be dumped"}};

    fhicl::Atom<bool> PrintOpHitAssociations{
      Name("PrintOpHitAssociations"),
      Comment("also prints a list of optical hits associated to the flash"),
      true};

    fhicl::Atom<std::string> OutputCategory{
      Name("OutputCategory"),
      Comment("the messagefacility category used for the output"),
      "DumpOpFlashes"};

  }; // Config

  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Default constructor
  explicit DumpOpFlashes(Parameters const& config);

  /// Does the printing
  void analyze(const art::Event& evt) override;

private:
  art::InputTag const fOpFlashModuleTag; ///< Optical hit data product tag.
  bool const fPrintOpHitAssociations;    ///< Whether to print optical hits.
  std::string const fOutputCategory;     ///< Category for `mf::LogInfo` output.

}; // class ophit::DumpOpFlashes

// -----------------------------------------------------------------------------
// ---  implementation
// -----------------------------------------------------------------------------
namespace {

  struct ProductSourceEntry_t {
    std::size_t index = 0;                        ///< An incremental index.
    art::BranchDescription const* info = nullptr; ///< Pointer to the product information.
  };
  using ProductSourceDirectory_t = std::unordered_map<art::ProductID, ProductSourceEntry_t>;

  /**
   * @brief Returns a map from product ID to product information for all
   *        pointers on one side of an association.
   * @tparam Side which side of the association (either `0` or `1`)
   * @tparam Assns the type of association
   * @tparam Event the type of event (makes this code gallery-compatible)
   * @param assns the association to be parsed
   * @param event the event holding the information
   * @return a map from product ID to product information
   *
   * The product information is a pair:
   *  * `first`: an incremental index
   *  * `second`: a pointer to the product information (_art_ branch description)
   */
  template <std::size_t Side = 1U, typename Assns, typename Event>
  ProductSourceDirectory_t buildAssnsSourceDirectory(Assns const& assns, Event const& event)
  {

    ProductSourceDirectory_t dir;

    // extract all the unique IDs
    for (auto const& ptrs : assns)
      dir.try_emplace(std::get<Side>(ptrs).id());

    // assign the index and pointer to each entry (order is not defined)
    std::size_t index{0};
    for (auto& [id, info] : dir)
      info = {index++, event.getProductDescription(id).get()};

    return dir;
  } // buildAssnsSourceDirectory()

  //----------------------------------------------------------------------------
  void printAssociatedOpHits(std::ostream& out,
                             std::vector<art::Ptr<recob::OpHit>> const& hits,
                             ProductSourceDirectory_t const* sourceMap,
                             std::string const& indent,
                             int const itemsPerLine)
  {
    out << hits.size() << " hits associated:";
    int itemsLeft = 0;
    const bool bPrintSources = sourceMap && (sourceMap->size() >= 2);
    for (art::Ptr<recob::OpHit> const& hitPtr : hits) {
      if (itemsLeft-- <= 0) {
        itemsLeft += itemsPerLine;
        out << "\n" << indent;
      }

      out << "  [";
      if (bPrintSources) out << "S#" << sourceMap->at(hitPtr.id()).index << ";";
      out << "#" << hitPtr.key() << "]";
      if (!hitPtr.isAvailable()) continue;

      itemsLeft -= 5; // full hit information is equivalent to 6 "items"

      recob::OpHit const& hit = *hitPtr;
      out << " ch=" << hit.OpChannel();
      if (hit.HasStartTime())
        out << " time=" << hit.StartTime();
      else
        out << " peak time=" << hit.PeakTime();
      out << " p.e.=" << hit.PE();

    } // for hits
  }   // printAssociatedOpHits()

  //----------------------------------------------------------------------------
  struct dumpAssociatedOpHits {
    std::vector<art::Ptr<recob::OpHit>> const& hits;
    ProductSourceDirectory_t const* sourceMap;
  };

  std::ostream& operator<<(std::ostream& out, dumpAssociatedOpHits const& hits)
  {
    printAssociatedOpHits(out, hits.hits, hits.sourceMap, "", 18);
    return out;
  }

  //----------------------------------------------------------------------------

} // local namespace

// -----------------------------------------------------------------------------
ophit::DumpOpFlashes::DumpOpFlashes(Parameters const& config)
  : art::EDAnalyzer(config)
  , fOpFlashModuleTag(config().OpFlashModuleLabel())
  , fPrintOpHitAssociations(config().PrintOpHitAssociations())
  , fOutputCategory(config().OutputCategory())
{
  consumes<std::vector<recob::OpFlash>>(fOpFlashModuleTag);
  if (fPrintOpHitAssociations)
    consumes<art::Assns<recob::OpFlash, recob::OpHit>>(fOpFlashModuleTag);
}

//------------------------------------------------------------------------------
void ophit::DumpOpFlashes::analyze(art::Event const& event)
{

  // fetch the data to be dumped on screen
  auto const& OpFlashHandle = event.getValidHandle<std::vector<recob::OpFlash>>(fOpFlashModuleTag);

  mf::LogVerbatim(fOutputCategory)
    << "Event " << event.id() << " contains " << OpFlashHandle->size() << " '"
    << fOpFlashModuleTag.encode() << "' optical flashes.";

  std::optional const flashHits =
    fPrintOpHitAssociations ?
      std::make_optional<art::FindManyP<recob::OpHit>>(OpFlashHandle, event, fOpFlashModuleTag) :
      std::nullopt;
  ProductSourceDirectory_t const sourceMap =
    fPrintOpHitAssociations ?
      buildAssnsSourceDirectory<1>(
        event.getProduct<art::Assns<recob::OpFlash, recob::OpHit>>(fOpFlashModuleTag), event) :
      ProductSourceDirectory_t{};
  assert(!fPrintOpHitAssociations || (flashHits->size() == 0) || !sourceMap.empty());

  for (auto const& [iFlash, flash] : util::enumerate(*OpFlashHandle)) {

    mf::LogVerbatim out{fOutputCategory};
    out << "OpFlash #" << iFlash << ": " << flash;

    if (flashHits) out << "\n  " << dumpAssociatedOpHits{flashHits->at(iFlash), &sourceMap};

  } // for flashes

  if (sourceMap.size() > 1) {
    mf::LogVerbatim out{fOutputCategory};
    out << "Hits were from " << sourceMap.size() << " data products (SRC):";
    for (auto const& [index, info] : util::values(sourceMap))
      out << "  [S#" << index << "]  '" << info->inputTag().encode() << "'";
    out << ".";
  }
  else if (sourceMap.size() == 1) {
    mf::LogVerbatim{fOutputCategory} << "All hits were from the data product '"
                                     << sourceMap.begin()->second.info->inputTag().encode() << "'.";
  }

} // ophit::DumpOpFlashes::analyze()

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(ophit::DumpOpFlashes)

//------------------------------------------------------------------------------
