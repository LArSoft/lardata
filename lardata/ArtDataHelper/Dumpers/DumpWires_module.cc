/**
 * @file   DumpWires_module.cc
 * @brief  Dumps on screen the content of the wires
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 31st, 2014
 */

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataalg/Utilities/StatCollector.h" // lar::util::MinMaxCollector<>
#include "lardataobj/RecoBase/Wire.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

// support libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C//C++ standard libraries
#include <iomanip> // std::setprecision(), std::setw()
#include <ios>     // std::fixed
#include <string>

namespace {

  // Copied from geo::PlaneGeo, local so far
  std::string ViewName(geo::View_t view)
  {
    switch (view) {
    case geo::kU: return "U";
    case geo::kV: return "V";
    case geo::kZ:
      return "Z";
      //  case geo::kY:       return "Y";
      //  case geo::kX:       return "X";
    case geo::k3D: return "3D";
    case geo::kUnknown: return "?";
    default: return "<UNSUPPORTED (" + std::to_string((int)view) + ")>";
    } // switch
  }   // ViewName()

} // local namespace

namespace caldata {

  /**
   * @brief Prints the content of all the wires on screen.
   *
   * This analyser prints the content of all the wires into the
   * `LogVerbatim` stream.
   *
   * Configuration parameters
   * =========================
   *
   * - *CalWireModuleLabel* (string, default: `"caldata"`): label of the
   *   producer used to create the `recob::Wire` collection to be dumped
   * - *OutputCategory* (string, default: `"DumpWires"`): the category used
   *   for the output (useful for filtering)
   * - *DigitsPerLine* (integer, default: `20`): the dump of digits and ticks
   *   will put this many of them for each line; `0` suppresses digit printout
   */
  class DumpWires : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> CalWireModuleLabel{
        Name("CalWireModuleLabel"),
        Comment("tag of producer used to create the recob::Wire collection"),
        "caldata" /* default */
      };

      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("the messagefacility category used for the output"),
        "DumpWires" /* default */
      };

      fhicl::Atom<unsigned int> DigitsPerLine{
        Name("DigitsPerLine"),
        Comment("number of digits printed per line (0: don't print digits)"),
        20 /* default */
      };

    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// Constructor.
    explicit DumpWires(Parameters const& config);

    /// Does the printing.
    virtual void analyze(art::Event const& evt) override;

  private:
    art::InputTag fCalWireModuleLabel; ///< Input tag for wires.
    std::string fOutputCategory;       ///< Category for `LogVerbatim` output.
    unsigned int fDigitsPerLine;       ///< Ticks/digits per line in the output.

    /// Dumps a single `recob:Wire` to the specified output stream.
    template <typename Stream>
    void PrintWire(Stream&& out,
                   recob::Wire const& wire,
                   std::string indent = "  ",
                   std::string firstIndent = "  ") const;

  }; // class DumpWires

} // namespace caldata

//------------------------------------------------------------------------------
//---  Implementation
//------------------------------------------------------------------------------
caldata::DumpWires::DumpWires(Parameters const& config)
  : EDAnalyzer(config)
  , fCalWireModuleLabel(config().CalWireModuleLabel())
  , fOutputCategory(config().OutputCategory())
  , fDigitsPerLine(config().DigitsPerLine())
{}

//------------------------------------------------------------------------------
void caldata::DumpWires::analyze(art::Event const& evt)
{

  auto const& Wires = *(evt.getValidHandle<std::vector<recob::Wire>>(fCalWireModuleLabel));

  mf::LogVerbatim(fOutputCategory) << "Event " << evt.id() << " contains " << Wires.size() << " '"
                                   << fCalWireModuleLabel.encode() << "' wires";

  for (recob::Wire const& wire : Wires) {

    PrintWire(mf::LogVerbatim(fOutputCategory), wire);

  } // for wire

} // caldata::DumpWires::analyze()

//------------------------------------------------------------------------------
template <typename Stream>
void caldata::DumpWires::PrintWire(Stream&& out,
                                   recob::Wire const& wire,
                                   std::string indent /* = "  " */,
                                   std::string firstIndent /* = "  " */
                                   ) const
{

  using RegionsOfInterest_t = recob::Wire::RegionsOfInterest_t;

  RegionsOfInterest_t const& RoIs = wire.SignalROI();

  //
  // print a header for the wire
  //
  out << firstIndent << "channel #" << wire.Channel() << " on view " << ViewName(wire.View())
      << "; " << wire.NSignal() << " time ticks";
  if (wire.NSignal() != RoIs.size()) out << " [!!! EXPECTED " << RoIs.size() << "]";
  if (RoIs.n_ranges() == 0) {
    out << " with nothing in them";
    return;
  }
  out << " with " << RoIs.n_ranges() << " regions of interest:";

  //
  // print the list of regions of interest
  //
  for (RegionsOfInterest_t::datarange_t const& RoI : RoIs.get_ranges()) {
    out << "\n" << indent << "  from " << RoI.offset << " for " << RoI.size() << " ticks";
  } // for

  //
  // print the content of the wire
  //
  if (fDigitsPerLine > 0) {

    std::vector<RegionsOfInterest_t::value_type> DigitBuffer(fDigitsPerLine), LastBuffer;

    unsigned int repeat_count = 0; // additional lines like the last one
    unsigned int index = 0;
    lar::util::MinMaxCollector<RegionsOfInterest_t::value_type> Extrema;
    out << "\n" << indent << "  content of the wire (" << fDigitsPerLine << " ticks per line):";
    auto iTick = RoIs.cbegin(), tend = RoIs.cend();
    while (iTick != tend) {
      // the next line will show at most fDigitsPerLine ticks
      unsigned int line_size = std::min(fDigitsPerLine, (unsigned int)RoIs.size() - index);
      if (line_size == 0) break; // no more ticks

      // fill the new buffer (iTick will move forward)
      DigitBuffer.resize(line_size);
      auto iBuf = DigitBuffer.begin(), bend = DigitBuffer.end();
      while ((iBuf != bend) && (iTick != tend))
        Extrema.add(*(iBuf++) = *(iTick++));
      index += line_size;

      // if the new buffer is the same as the old one, just mark it
      if (DigitBuffer == LastBuffer) {
        repeat_count += 1;
        continue;
      }

      // if there are previous repeats, write that on screen
      // before the new, different line
      if (repeat_count > 0) {
        out << "\n"
            << indent << "  [ ... repeated " << repeat_count << " more times, "
            << (repeat_count * LastBuffer.size()) << " ticks ]";
        repeat_count = 0;
      }

      // dump the new line of ticks
      out << "\n" << indent << " " << std::fixed << std::setprecision(3);
      for (auto digit : DigitBuffer)
        out << std::setw(8) << digit;

      // quick way to assign DigitBuffer to LastBuffer
      // (we don't care we lose the former)
      std::swap(LastBuffer, DigitBuffer);

    } // while
    if (repeat_count > 0) {
      out << "\n" << indent << "  [ ... repeated " << repeat_count << " more times to the end ]";
    }
    if (Extrema.min() < Extrema.max()) {
      out << "\n"
          << indent << "    range of " << index << " samples: [" << Extrema.min() << ";"
          << Extrema.max() << "]";
    }
  } // if dumping the ticks

} // caldata::DumpWires::PrintWire()

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(caldata::DumpWires)

//------------------------------------------------------------------------------
