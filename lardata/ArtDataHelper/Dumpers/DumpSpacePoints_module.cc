/**
 * @file   DumpSpacePoints_module.cc
 * @brief  Dumps on screen the content of space points
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 18th, 2015
 */

// LArSoft includes
#include "lardata/ArtDataHelper/Dumpers/NewLine.h" // recob::dumper::makeNewLine()
#include "lardata/ArtDataHelper/Dumpers/SpacePointDumpers.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

// support libraries
#include "fhiclcpp/types/Atom.h" // also pulls in fhicl::Name and fhicl::Comment

// C//C++ standard libraries
#include <string>

// ... and more in the implementation part

namespace recob {

  /**
   * @brief Prints the content of all the space points on screen
   *
   * This analyser prints the content of all the space points into the
   * LogInfo/LogVerbatim stream.
   *
   * Configuration parameters
   * =========================
   *
   * - *SpacePointModuleLabel* (art::InputTag, mandatory): label of the
   *   producer used to create the recob::SpacePoint collection to be dumped
   * - *OutputCategory* (string, default: "DumpSpacePoints"): the category used
   *   for the output (useful for filtering)
   * - *PrintHexFloats* (boolean, default: `false`): print all the floating
   *   point numbers in base 16
   *
   */
  class DumpSpacePoints : public art::EDAnalyzer {
  public:
    /// Configuration parameters
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> SpacePointModuleLabel{
        Name("SpacePointModuleLabel"),
        Comment(
          "label of the producer used to create the recob::SpacePoint collection to be dumped")};
      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("the category used for the output (useful for filtering) [\"DumpSpacePoints\"]"),
        "DumpSpacePoints" /* default value */
      };
      fhicl::Atom<bool> PrintHexFloats{
        Name("PrintHexFloats"),
        Comment("print floating point numbers in base 16 [false]"),
        false /* default value */
      };

    }; // struct Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// Default constructor
    explicit DumpSpacePoints(Parameters const& config);

    /// Does the printing
    virtual void analyze(const art::Event& evt) override;

  private:
    art::InputTag fInputTag;     ///< input tag of the SpacePoint product
    std::string fOutputCategory; ///< category for LogInfo output
    bool fPrintHexFloats;        ///< whether to print floats in base 16

  }; // class DumpSpacePoints

} // namespace recob

//==============================================================================
//===  Implementation section
//==============================================================================

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

// art libraries
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindMany.h"

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// C//C++ standard libraries

namespace {

  //----------------------------------------------------------------------------
  class SpacePointDumper {
  public:
    using PrintOptions_t = recob::dumper::SpacePointPrintOptions_t;

    /// Constructor; will dump space points from the specified list
    SpacePointDumper(std::vector<recob::SpacePoint> const& point_list,
                     PrintOptions_t const& printOptions = {})
      : points(point_list), options(printOptions)
    {}

    /// Sets the hits associated to each space point
    void SetHits(art::FindMany<recob::Hit> const* hit_query) { hits = hit_query; }

    /// Dump a space point specified by its index in the input list
    template <typename Stream>
    void DumpSpacePoint(Stream&& out, size_t iPoint) const
    {
      DumpSpacePoint(std::forward<Stream>(out), iPoint, options);
    }

    /// Dump a space point specified by its index in the input list
    template <typename Stream>
    void DumpSpacePoint(Stream&& out, size_t iPoint, std::string indentstr) const
    {
      PrintOptions_t localOptions(options);
      localOptions.indent.indent = indentstr;
      DumpSpacePoint(std::forward<Stream>(out), iPoint, localOptions);
    }

    /// Dump a space point specified by its index in the input list
    template <typename Stream>
    void DumpSpacePoint(Stream&& out, size_t iPoint, PrintOptions_t const& localOptions) const
    {
      recob::SpacePoint const& point = points.at(iPoint);

      //
      // intro
      //
      auto first_nl = recob::dumper::makeNewLine(out, localOptions.indent);
      first_nl() << "[#" << iPoint << "] ";

      PrintOptions_t indentedOptions(localOptions);
      indentedOptions.indent.appendIndentation("  ");
      recob::dumper::DumpSpacePoint(std::forward<Stream>(out), point, indentedOptions);

      //
      // hits
      //
      if (hits) {
        std::vector<recob::Hit const*> myHits = hits->at(iPoint);
        if (myHits.empty()) { out << "; no associated hits"; }
        else {
          auto nl = recob::dumper::makeNewLine(out, indentedOptions.indent);
          out << "; " << myHits.size() << " hits:";
          for (recob::Hit const* hit : myHits) {
            nl() << "  on " << hit->WireID() << ", peak at tick " << hit->PeakTime() << ", "
                 << hit->PeakAmplitude() << " ADC, RMS: " << hit->RMS()
                 << " (channel: " << hit->Channel() << ")";
          } // for hits
        }   // if we have hits
      }     // if we have hit information

      //
      // done
      //

    } // DumpSpacePoints()

    /// Dumps all space points in the input list
    template <typename Stream>
    void DumpAllSpacePoints(Stream&& out, std::string indentstr = "") const
    {
      auto localOptions = options;
      localOptions.indent.appendIndentation(indentstr);
      size_t const nPoints = points.size();
      for (size_t iPoint = 0; iPoint < nPoints; ++iPoint)
        DumpSpacePoint(std::forward<Stream>(out), iPoint, localOptions);
    } // DumpAllSpacePoints()

  protected:
    std::vector<recob::SpacePoint> const& points; ///< input list
    PrintOptions_t options;                       ///< formatting and indentation options

    /// Associated hits (expected same order as for space points)
    art::FindMany<recob::Hit> const* hits = nullptr;

  }; // SpacePointDumper

  //----------------------------------------------------------------------------

} // local namespace

namespace recob {

  //----------------------------------------------------------------------------
  DumpSpacePoints::DumpSpacePoints(Parameters const& config)
    : EDAnalyzer(config)
    , fInputTag(config().SpacePointModuleLabel())
    , fOutputCategory(config().OutputCategory())
    , fPrintHexFloats(config().PrintHexFloats())
  {}

  //----------------------------------------------------------------------------
  void DumpSpacePoints::analyze(const art::Event& evt)
  {

    //
    // collect all the available information
    //
    // fetch the data to be dumped on screen
    auto SpacePoints = evt.getValidHandle<std::vector<recob::SpacePoint>>(fInputTag);

    art::FindMany<recob::Hit> const PointHits(SpacePoints, evt, fInputTag);

    size_t const nPoints = SpacePoints->size();
    mf::LogInfo(fOutputCategory) << "The event contains " << nPoints << " space points from '"
                                 << fInputTag.encode() << "'";

    // prepare the dumper
    SpacePointDumper dumper(*SpacePoints);
    if (PointHits.isValid())
      dumper.SetHits(&PointHits);
    else
      mf::LogWarning("DumpSpacePoints") << "hit information not avaialble";

    dumper.DumpAllSpacePoints(mf::LogVerbatim(fOutputCategory), "  ");

    mf::LogVerbatim(fOutputCategory) << "\n"; // two empty lines

  } // DumpSpacePoints::analyze()

  DEFINE_ART_MODULE(DumpSpacePoints)

} // namespace recob
