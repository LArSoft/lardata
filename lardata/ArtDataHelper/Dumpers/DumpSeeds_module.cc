/**
 * @file   DumpSeeds_module.cc
 * @brief  Dumps on screen the content of seeds
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 29th, 2015
 */

// LArSoft includes
#include "lardata/ArtDataHelper/Dumpers/hexfloat.h"
#include "lardataobj/RecoBase/Seed.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"

// C//C++ standard libraries
#include <string>

// ... and more in the implementation part

namespace recob {

  /**
   * @brief Prints the content of all the seeds on screen
   *
   * This analyser prints the content of all the seeds into the
   * LogInfo/LogVerbatim stream.
   *
   * Configuration parameters
   * =========================
   *
   * - *SeedModuleLabel* (art::InputTag, mandatory): label of the
   *   producer used to create the recob::Seed collection to be dumped
   * - *OutputCategory* (string, default: "DumpSeeds"): the category used
   *   for the output (useful for filtering)
   * - *PrintHexFloats* (boolean, default: `false`): print all the floating
   *   point numbers in base 16
   *
   */
  class DumpSeeds : public art::EDAnalyzer {
  public:
    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> SeedModuleLabel{
        Name("SeedModuleLabel"),
        Comment("tag of the recob::Seed collection data product to be dumped")};

      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("name of the message facility category to be used for output"),
        "DumpSeeds"};

      fhicl::Atom<bool> PrintHexFloats{Name("PrintHexFloats"),
                                       Comment("print all the floating point numbers in base 16"),
                                       false};

    }; // struct Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// Default constructor
    explicit DumpSeeds(Parameters const& config);

    /// Does the printing
    virtual void analyze(const art::Event& evt) override;

  private:
    art::InputTag fInputTag;     ///< input tag of the Seed product
    std::string fOutputCategory; ///< category for LogInfo output
    bool fPrintHexFloats;        ///< whether to print floats in base 16

  }; // class DumpSeeds

} // namespace recob

//==============================================================================
//===  Implementation section
//==============================================================================

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Seed.h"

// art libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindMany.h"

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C//C++ standard libraries
#include <array>

namespace {

  //----------------------------------------------------------------------------
  class SeedDumper {
  public:
    /// Collection of available printing style options
    struct PrintOptions_t {
      bool hexFloats = false; ///< print all floating point numbers in base 16
      std::string indent;     ///< indentation string
    };                        // PrintOptions_t

    /// Constructor; will dump seeds from the specified list.
    SeedDumper(std::vector<recob::Seed> const& seed_list) : SeedDumper(seed_list, {}) {}

    /// Constructor; will dump seeds from the specified list
    SeedDumper(std::vector<recob::Seed> const& seed_list, PrintOptions_t print_options)
      : seeds(seed_list), options(print_options)
    {}

    /// Sets the hits associated to each seed
    void SetHits(art::FindMany<recob::Hit> const* hit_query) { hits = hit_query; }

    /// Dump a seed specified by its index in the input particle list
    template <typename Stream>
    void DumpSeed(Stream&& out, size_t iSeed) const
    {
      lar::OptionalHexFloat hexfloat(options.hexFloats);
      std::string const& indentstr = options.indent;

      recob::Seed const& seed = seeds.at(iSeed);
      //
      // intro
      //
      out << "\n" << indentstr << "[#" << iSeed << "]";
      if (!seed.IsValid())
        out << " invalid!";
      else {
        std::array<double, 3> start, dir;
        seed.GetDirection(dir.data(), nullptr);
        seed.GetPoint(start.data(), nullptr);
        out << " starts at (" << hexfloat(start[0]) << "," << hexfloat(start[1]) << ","
            << hexfloat(start[2]) << ") toward (" << hexfloat(dir[0]) << "," << hexfloat(dir[1])
            << "," << hexfloat(dir[2]) << "); length: " << hexfloat(seed.GetLength()) << " cm";
      }

      //
      // hits
      //
      if (hits) {
        std::vector<recob::Hit const*> myHits = hits->at(iSeed);
        if (!myHits.empty()) {
          // we do not honour the base 16 printout requirement here, because
          // these data members are single precision and there is no printf()
          // flag taking a float as an argument;, and a promotion to double
          // would silently occurr, which we want to avoid
          out << "; " << myHits.size() << " hits:";
          for (recob::Hit const* hit : myHits) {
            out << "\n"
                << indentstr << "  on " << hit->WireID() << ", peak at tick " << hit->PeakTime()
                << ", " << hit->PeakAmplitude() << " ADC, RMS: " << hit->RMS()
                << " (channel: " << hit->Channel() << ")";
          } // for hits
        }   // if we have hits
      }     // if we have hit information

      //
      // done
      //

    } // DumpSeed()

    /// Dumps all seeds in the input list
    template <typename Stream>
    void DumpAllSeeds(Stream&& out) const
    {
      size_t const nSeeds = seeds.size();
      for (size_t iSeed = 0; iSeed < nSeeds; ++iSeed)
        DumpSeed(out, iSeed);
    } // DumpAllSeeds()

  protected:
    std::vector<recob::Seed> const& seeds; ///< input list

    PrintOptions_t options; ///< printing and formatting options

    /// Associated hits (expected same order as for seeds)
    art::FindMany<recob::Hit> const* hits = nullptr;

  }; // SeedDumper

  //----------------------------------------------------------------------------

} // local namespace

namespace recob {

  //----------------------------------------------------------------------------
  DumpSeeds::DumpSeeds(Parameters const& config)
    : EDAnalyzer(config)
    , fInputTag(config().SeedModuleLabel())
    , fOutputCategory(config().OutputCategory())
    , fPrintHexFloats(config().PrintHexFloats())
  {}

  //----------------------------------------------------------------------------
  void DumpSeeds::analyze(const art::Event& evt)
  {

    //
    // collect all the available information
    //
    // fetch the data to be dumped on screen
    auto Seeds = evt.getValidHandle<std::vector<recob::Seed>>(fInputTag);

    art::FindMany<recob::Hit> const SeedHits(Seeds, evt, fInputTag);

    size_t const nSeeds = Seeds->size();
    mf::LogVerbatim(fOutputCategory) << "Event " << evt.id() << " contains " << nSeeds
                                     << " seeds from '" << fInputTag.encode() << "'";

    // prepare the dumper
    SeedDumper::PrintOptions_t options;
    options.hexFloats = fPrintHexFloats;
    options.indent = "  ";
    SeedDumper dumper(*Seeds, options);

    if (SeedHits.isValid())
      dumper.SetHits(&SeedHits);
    else
      mf::LogWarning("DumpSeeds") << "hit information not avaialble";

    dumper.DumpAllSeeds(mf::LogVerbatim(fOutputCategory));

    mf::LogVerbatim(fOutputCategory) << "\n"; // two empty lines

  } // DumpSeeds::analyze()

  DEFINE_ART_MODULE(DumpSeeds)

} // namespace recob
