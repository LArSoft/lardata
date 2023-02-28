/**
 * @file   DumpRawDigits_module.cc
 * @brief  Dumps on screen the content of the raw digits
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 14th, 2015
 */

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "lardataalg/Utilities/StatCollector.h"          // lar::util::MinMaxCollector<>
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h" // raw::Uncompress()

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
#include <algorithm> // std::min(), std::copy_n()
#include <iomanip>   // std::setprecision(), std::setw()
#include <string>

namespace detsim {

  /**
   * @brief Prints the content of all the raw digits on screen.
   *
   * This analyser prints the content of all the raw digits into the
   * `LogVerbatim` stream.
   *
   * Configuration parameters
   * =========================
   *
   * - *DetSimModuleLabel* (string, default: `"daq"`): label of the
   *   producer used to create the `raw::RawDigits` collection
   * - *OutputCategory* (string, default: `"DumpDigits"`): the category used
   *   for the output (useful for filtering)
   * - *DigitsPerLine* (integer, default: `20`): the dump of digits and ticks
   *   will put this many of them for each line
   * - *Pedestal* (integer, default: `0`): digit values are written relative
   *   to this number
   *
   */
  class DumpRawDigits : public art::EDAnalyzer {

    /// Type to represent a digit.
    using Digit_t = raw::RawDigit::ADCvector_t::value_type;

    /// Type to represent a pedestal.
    using Pedestal_t = Digit_t;

  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> DetSimModuleLabel{
        Name("DetSimModuleLabel"),
        Comment("tag of producer used to create the raw::RawDigit collection"),
        "daq" /* default */
      };

      fhicl::Atom<std::string> OutputCategory{
        Name("OutputCategory"),
        Comment("the messagefacility category used for the output"),
        "DumpDigits" /* default */
      };

      fhicl::Atom<unsigned int> DigitsPerLine{
        Name("DigitsPerLine"),
        Comment("number of digits printed per line (0: don't print digits)"),
        20 /* default */
      };

      fhicl::Atom<Pedestal_t> Pedestal{
        Name("Pedestal"),
        Comment("digit values are written relative to this number"),
        0 /* default */
      };

    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    /// Constructor.
    explicit DumpRawDigits(Parameters const& config);

    /// Prints an introduction.
    virtual void beginJob() override;

    /// Does the printing.
    virtual void analyze(art::Event const& evt) override;

  private:
    art::InputTag fDetSimModuleLabel; ///< Tag for digits data product.
    std::string fOutputCategory;      ///< Category for `LogVerbatim` output.
    unsigned int fDigitsPerLine;      ///< Ticks/digits per line in the output.
    Pedestal_t fPedestal;             ///< ADC pedestal, will be subtracted from digits.

    /// Dumps a single `recob:Wire` to the specified output stream.
    template <typename Stream>
    void PrintRawDigit(Stream&& out,
                       raw::RawDigit const& digits,
                       std::string indent = "  ",
                       std::string firstIndent = "  ") const;

  }; // class DumpRawDigits

} // namespace detsim

//------------------------------------------------------------------------------
//---  Implementation
//------------------------------------------------------------------------------
detsim::DumpRawDigits::DumpRawDigits(Parameters const& config)
  : EDAnalyzer(config)
  , fDetSimModuleLabel(config().DetSimModuleLabel())
  , fOutputCategory(config().OutputCategory())
  , fDigitsPerLine(config().DigitsPerLine())
  , fPedestal(config().Pedestal())
{}

//------------------------------------------------------------------------------
void detsim::DumpRawDigits::beginJob()
{

  if (fPedestal != 0) {
    mf::LogVerbatim(fOutputCategory)
      << "A pedestal of " << fPedestal << " will be subtracted from all raw digits";
  } // if pedestal

} // detsim::DumpRawDigits::beginJob()

//------------------------------------------------------------------------------
void detsim::DumpRawDigits::analyze(art::Event const& evt)
{

  auto const& RawDigits = *(evt.getValidHandle<std::vector<raw::RawDigit>>(fDetSimModuleLabel));

  mf::LogVerbatim(fOutputCategory) << "Event " << evt.id() << " contains " << RawDigits.size()
                                   << " '" << fDetSimModuleLabel.encode() << "' waveforms";
  for (raw::RawDigit const& digits : RawDigits) {

    PrintRawDigit(mf::LogVerbatim(fOutputCategory), digits);

  } // for digits

} // caldata::DumpWires::analyze()

//------------------------------------------------------------------------------
template <typename Stream>
void detsim::DumpRawDigits::PrintRawDigit(Stream&& out,
                                          raw::RawDigit const& digits,
                                          std::string indent /* = "  " */,
                                          std::string firstIndent /* = "  " */
                                          ) const
{

  //
  // uncompress the digits
  //
  raw::RawDigit::ADCvector_t ADCs(digits.Samples());
  raw::Uncompress(digits.ADCs(), ADCs, digits.Compression());

  //
  // print a header for the raw digits
  //
  out << firstIndent << "  #" << digits.Channel() << ": " << ADCs.size() << " time ticks";
  if (digits.Samples() != ADCs.size()) out << " [!!! EXPECTED " << digits.Samples() << "] ";
  out << " (" << digits.NADC() << " after compression); compression type: ";
  switch (digits.Compression()) {
  case raw::kNone: out << "no compression"; break;
  case raw::kHuffman: out << "Huffman encoding"; break;
  case raw::kZeroSuppression: out << "zero suppression"; break;
  case raw::kZeroHuffman: out << "zero suppression + Huffman encoding"; break;
  case raw::kDynamicDec: out << "dynamic decimation"; break;
  default: out << "unknown (#" << ((int)digits.Compression()) << ")"; break;
  } // switch

  // print the content of the channel
  if (fDigitsPerLine > 0) {
    std::vector<Digit_t> DigitBuffer(fDigitsPerLine), LastBuffer;

    unsigned int repeat_count = 0; // additional lines like the last one
    unsigned int index = 0;
    lar::util::MinMaxCollector<Digit_t> Extrema;
    out << "\n" << indent << "content of the channel (" << fDigitsPerLine << " ticks per line):";
    auto iTick = ADCs.cbegin(), tend = ADCs.cend(); // const iterators
    while (iTick != tend) {
      // the next line will show at most fDigitsPerLine ticks
      unsigned int line_size = std::min(fDigitsPerLine, (unsigned int)ADCs.size() - index);
      if (line_size == 0) break; // no more ticks

      // fill the new buffer (iTick will move forward)
      DigitBuffer.resize(line_size);
      auto iBuf = DigitBuffer.begin(), bend = DigitBuffer.end();
      while ((iBuf != bend) && (iTick != tend))
        Extrema.add(*(iBuf++) = (*(iTick++) - fPedestal));
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
      out << "\n" << indent << " ";
      for (auto digit : DigitBuffer)
        out << " " << std::setw(4) << digit;

      // quick way to assign DigitBuffer to LastBuffer
      // (we don't care we lose the former)
      std::swap(LastBuffer, DigitBuffer);

    } // while
    if (repeat_count > 0) {
      out << "\n" << indent << "  [ ... repeated " << repeat_count << " more times to the end ]";
    }
    if (Extrema.min() < Extrema.max()) {
      out << "\n"
          << indent << "  range of " << index << " samples: [" << Extrema.min() << ";"
          << Extrema.max() << "]";
    }
  } // if dumping the ticks

} // detsim::DumpRawDigits::analyze()

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(detsim::DumpRawDigits)

//------------------------------------------------------------------------------
