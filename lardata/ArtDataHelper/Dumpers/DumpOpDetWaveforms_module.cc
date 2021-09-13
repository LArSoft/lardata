/**
 * @file   DumpOpDetWaveforms_module.cc
 * @brief  Dumps on screen the content of the raw optical detector waveforms.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 8th, 2017
 */

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/Dumpers/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/OpDetWaveform.h"

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
#include <algorithm> // std::sort()
#include <string>
#include <vector>

namespace detsim {

  /**
   * @brief Prints the content of all optical detector waveforms on screen.
   *
   * This analyser prints the content of all the raw digits into the
   * LogInfo/LogVerbatim stream.
   *
   * Configuration parameters
   * -------------------------
   *
   * - *OpDetWaveformsTag* (string, default: `daq`): input tag of the
   *   raw::OpDetWaveform collection to be dumped
   * - *OutputCategory* (string, default: `DumpOpDetWaveforms`): the category
   *   used for the output (useful for filtering)
   * - *DigitsPerLine* (integer, default: `20`): the dump of ADC readings
   *   will put this many of them for each line
   * - *Pedestal* (integer, default: `0`): ADC readings are written relative
   *   to this number
   * - *TickLabel* (string, default: `"tick"`): a preamble to each line of the
   *   dumped waveform digits, chosen among:
   *     - `"tick"`: the tick number of the waveform is printed (starts at `0`)
   *     - `"time"`: timestamp (&micro;s) of the first tick in the row
   *     - `"none"`: no preamble written at all
   *
   */
  class DumpOpDetWaveforms : public art::EDAnalyzer {
  public:
    struct Config {

      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> OpDetWaveformsTag{
        Name("OpDetWaveformsTag"),
        Comment("input tag of the raw::OpDetWaveform collection to be dumped")};

      fhicl::Atom<std::string> OutputCategory{Name("OutputCategory"),
                                              Comment("name of the category used for the output"),
                                              "DumpOpDetWaveforms"};

      fhicl::Atom<unsigned int> DigitsPerLine{
        Name("DigitsPerLine"),
        Comment("the dump of ADC readings will put this many of them for each line"),
        20U};

      fhicl::Atom<raw::ADC_Count_t> Pedestal{
        Name("Pedestal"),
        Comment("ADC readings are written relative to this number"),
        0};

      fhicl::Atom<bool> SortByChannelAndTime{
        Name("SortByChannelAndTime"),
        Comment
          ("waveforms are dumped in channel number order, and then timestamp"),
        true};

      fhicl::Atom<std::string> TickLabel{
        Name("TickLabel"),
        Comment("write an index in front of each digit dump line; choose between:"
                " \"tick\" (waveform tick number)"
                ", \"timestamp\" (electronics clock time in microseconds)"
                ", \"none\" (no tick label)"),
        "tick"};

    }; // struct Config

    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit DumpOpDetWaveforms(Parameters const& config);

    /// Does the printing.
    void analyze(const art::Event& evt);

  private:
    
    enum class sortMode {
      DataProductOrder, ///< Unsorted: same order as the input data product.
      ChannelAndTime    ///< Sort by PMT channel number, then timestamp.
    }; // sortMode
    
    /// Base functor for printing time according to tick number.
    struct TimestampLabelMaker : public dump::raw::OpDetWaveformDumper::TimeLabelMaker {
      double tickDuration; // should this me `util::quantities::microsecond`?

      /// Constructor: specify the duration of optical clock tick [&micro;s].
      TimestampLabelMaker(double tickDuration) : tickDuration(tickDuration) {}

      /// Returns the electronics time of the specified waveform tick.
      virtual std::string
      label(raw::OpDetWaveform const& waveform, unsigned int tick) const override
      {
        return std::to_string(waveform.TimeStamp() + tick * tickDuration);
      }

    }; // struct TimestampLabelMaker

    art::InputTag fOpDetWaveformsTag; ///< Input tag of data product to dump.
    std::string fOutputCategory;      ///< Category for `mf::LogInfo` output.
    unsigned int fDigitsPerLine;      ///< ADC readings per line in the output.
    raw::ADC_Count_t fPedestal;       ///< ADC pedestal (subtracted from readings).
    sortMode fSortByChannelAndTime; ///< How to sort the waveforms in the dump.
    

    /// The object used to print tick labels.
    std::unique_ptr<dump::raw::OpDetWaveformDumper::TimeLabelMaker> fTimeLabel;

    /// Returns pointers to all waveforms in a vector with channel as index.
    static std::vector<std::vector<raw::OpDetWaveform const*>> groupByChannel(
      std::vector<raw::OpDetWaveform> const& waveforms);

    /// Sorts all the waveforms in the vector by growing timestamp.
    static void sortByTimestamp(std::vector<raw::OpDetWaveform const*>& waveforms);

  }; // class DumpOpDetWaveforms

} // namespace detsim

namespace detsim {

  //-------------------------------------------------
  DumpOpDetWaveforms::DumpOpDetWaveforms(Parameters const& config)
    : EDAnalyzer(config)
    , fOpDetWaveformsTag(config().OpDetWaveformsTag())
    , fOutputCategory(config().OutputCategory())
    , fDigitsPerLine(config().DigitsPerLine())
    , fPedestal(config().Pedestal())
    , fSortByChannelAndTime(config().SortByChannelAndTime()
        ? sortMode::ChannelAndTime: sortMode::DataProductOrder
      )
  {
    std::string const tickLabelStr = config().TickLabel();
    if (tickLabelStr == "none") {
      fTimeLabel.reset(); // not very useful, but let's be explicit
    }
    else if (tickLabelStr == "tick") {
      fTimeLabel = std::make_unique<dump::raw::OpDetWaveformDumper::TickLabelMaker>();
    }
    else if (tickLabelStr == "time") {
      auto const clock_data =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
      fTimeLabel = std::make_unique<TimestampLabelMaker>(clock_data.OpticalClock().TickPeriod());
    }
    else {
      throw art::Exception(art::errors::Configuration)
        << "Invalid choice '" << tickLabelStr << "' for time label.\n";
    }

  } // DumpOpDetWaveforms::DumpOpDetWaveforms()

  //-------------------------------------------------
  void
  DumpOpDetWaveforms::analyze(const art::Event& event)
  {

    // fetch the data to be dumped on screen
    auto const& Waveforms
      = event.getProduct<std::vector<raw::OpDetWaveform>>(fOpDetWaveformsTag);

    dump::raw::OpDetWaveformDumper dump(fPedestal, fDigitsPerLine);
    dump.setTimeLabelMaker(fTimeLabel.get());

    mf::LogVerbatim(fOutputCategory) << "The event " << event.id() << " contains data for "
                                     << Waveforms.size() << " optical detector channels";
    if (fPedestal != 0) {
      mf::LogVerbatim(fOutputCategory)
        << "A pedestal of " << fPedestal << " counts will be subtracted from all ADC readings.";
    } // if pedestal

    switch (fSortByChannelAndTime) {
      case sortMode::DataProductOrder:
        {
          dump.setIndent("  ");
          for (raw::OpDetWaveform const& waveform : Waveforms) {
            mf::LogVerbatim log(fOutputCategory);
            dump(log, waveform);
          } // for waveforms on channel
        }
        break; // sortMode::DataProductOrder
      case sortMode::ChannelAndTime:
        {
          dump.setIndent("    ");

          auto groupedWaveforms = groupByChannel(Waveforms);
          for (auto& channelWaveforms : groupedWaveforms) {
            if (channelWaveforms.empty()) continue;

            sortByTimestamp(channelWaveforms);
            auto const channel = channelWaveforms.front()->ChannelNumber();

            mf::LogVerbatim(fOutputCategory)
              << "  optical detector channel #" << channel << " has "
              << channelWaveforms.size() << " waveforms:";

            for (raw::OpDetWaveform const* pWaveform : channelWaveforms) {
              mf::LogVerbatim log(fOutputCategory);
              dump(log, *pWaveform);
            } // for waveforms on channel

          } // for all channels
        }
        break; // sortMode::ChannelAndTime
    } // switch (fSortMode)

  } // DumpOpDetWaveforms::analyze()

  //----------------------------------------------------------------------------
  std::vector<std::vector<raw::OpDetWaveform const*>>
  DumpOpDetWaveforms::groupByChannel(std::vector<raw::OpDetWaveform> const& waveforms)
  {
    std::vector<std::vector<raw::OpDetWaveform const*>> groups;
    for (auto const& waveform : waveforms) {
      auto const channel = waveform.ChannelNumber();
      if (groups.size() <= channel) groups.resize(channel + 1);
      groups[channel].push_back(&waveform);
    } // for
    return groups;
  } // DumpOpDetWaveforms::groupByChannel()

  //----------------------------------------------------------------------------
  void
  DumpOpDetWaveforms::sortByTimestamp(std::vector<raw::OpDetWaveform const*>& waveforms)
  {

    struct ChannelSorter {

      static bool
      sort(raw::OpDetWaveform const& a, raw::OpDetWaveform const& b)
      {
        if (a.ChannelNumber() < b.ChannelNumber()) return true;
        if (a.ChannelNumber() > b.ChannelNumber()) return false;

        return (a.TimeStamp() < b.TimeStamp());
      }

      bool
      operator()(raw::OpDetWaveform const& a, raw::OpDetWaveform const& b) const
      {
        return sort(a, b);
      }

      bool
      operator()(raw::OpDetWaveform const* a, raw::OpDetWaveform const* b) const
      {
        return sort(*a, *b);
      }

    }; // struct ChannelSorter
    std::sort(waveforms.begin(), waveforms.end(), ChannelSorter());

  } // DumpOpDetWaveforms::sortByTimestamp()

  //----------------------------------------------------------------------------
  DEFINE_ART_MODULE(DumpOpDetWaveforms)

  //----------------------------------------------------------------------------

} // namespace detsim
