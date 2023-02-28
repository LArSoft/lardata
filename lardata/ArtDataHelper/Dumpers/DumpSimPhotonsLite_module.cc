/**
 * @file   DumpSimPhotonsLite_module.cc
 * @brief  Module dumping SimPhotonsLite information on screen
 * @date   March 8, 2017
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 */

// LArSoft libraries
#include "lardataobj/Simulation/SimPhotons.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <iomanip> // std::setw()
#include <numeric> // std::accumulate()
#include <string>
#include <utility> // std::forward()

namespace sim {
  class DumpSimPhotonsLite;
} // namespace sim

namespace {
  using namespace fhicl;

  /// Collection of configuration parameters for the module
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> InputPhotons{
      Name("InputPhotons"),
      Comment("data product with the SimPhotonsLite to be dumped")};

    fhicl::Atom<std::string> OutputCategory{
      Name("OutputCategory"),
      Comment("name of the output stream (managed by the message facility)"),
      "DumpSimPhotonsLite" /* default value */
    };

  }; // struct Config

} // local namespace

class sim::DumpSimPhotonsLite : public art::EDAnalyzer {
public:
  // type to enable module parameters description by art
  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Configuration-checking constructor
  explicit DumpSimPhotonsLite(Parameters const& config);

  // Plugins should not be copied or assigned.
  DumpSimPhotonsLite(DumpSimPhotonsLite const&) = delete;
  DumpSimPhotonsLite(DumpSimPhotonsLite&&) = delete;
  DumpSimPhotonsLite& operator=(DumpSimPhotonsLite const&) = delete;
  DumpSimPhotonsLite& operator=(DumpSimPhotonsLite&&) = delete;

  // Operates on the event
  void analyze(art::Event const& event) override;

  /**
   * @brief Dumps the content of specified SimPhotonsLite in the output stream.
   * @tparam Stream the type of output stream
   * @param out the output stream
   * @param photons the SimPhotonsLite to be dumped
   * @param indent base indentation string _(default: none)_
   * @param firstIndent if first output line should be indented _(default: yes)_
   *
   * The indent string is prepended to every line of output, with the possible
   * exception of the first one, in case bIndentFirst is true.
   *
   * The output starts on the current line, and the last line is *not* broken.
   */
  template <typename Stream>
  void DumpPhoton(Stream&& out,
                  sim::SimPhotonsLite const& photons,
                  std::string indent,
                  std::string firstIndent) const;

  template <typename Stream>
  void DumpPhoton(Stream&& out, sim::SimPhotonsLite const& photons, std::string indent = "") const
  {
    DumpPhoton(std::forward<Stream>(out), photons, indent, indent);
  }

private:
  art::InputTag fInputPhotons; ///< name of SimPhotons's data product
  std::string fOutputCategory; ///< name of the stream for output

}; // class sim::DumpSimPhotonsLite

//------------------------------------------------------------------------------
//---  module implementation
//---
//------------------------------------------------------------------------------
sim::DumpSimPhotonsLite::DumpSimPhotonsLite(Parameters const& config)
  : EDAnalyzer(config)
  , fInputPhotons(config().InputPhotons())
  , fOutputCategory(config().OutputCategory())
{}

//------------------------------------------------------------------------------
template <typename Stream>
void sim::DumpSimPhotonsLite::DumpPhoton(Stream&& out,
                                         sim::SimPhotonsLite const& photons,
                                         std::string indent,
                                         std::string firstIndent) const
{

  unsigned int const nPhotons =
    std::accumulate(photons.DetectedPhotons.begin(),
                    photons.DetectedPhotons.end(),
                    0U,
                    [](auto sum, auto const& entry) { return sum + entry.second; });

  out << firstIndent << "channel=" << photons.OpChannel << " has ";
  if (nPhotons) {
    out << nPhotons << " photons (format: [tick] photons):";
    constexpr unsigned int PageSize = 5;
    unsigned int pager = 0;
    for (auto const& pair : photons.DetectedPhotons) {
      if (pager-- == 0) {
        pager = PageSize - 1;
        out << "\n" << indent << " ";
      }
      out << " [" << pair.first << "] " << std::setw(6) << pair.second;
    } // for
  }
  else {
    out << "no photons";
  }

} // sim::DumpSimPhotonsLite::DumpPhoton()

//------------------------------------------------------------------------------
void sim::DumpSimPhotonsLite::analyze(art::Event const& event)
{

  // get the particles from the event
  auto const& Photons = *(event.getValidHandle<std::vector<sim::SimPhotonsLite>>(fInputPhotons));

  mf::LogVerbatim(fOutputCategory)
    << "Event " << event.id() << " : data product '" << fInputPhotons.encode() << "' contains "
    << Photons.size() << " SimPhotonsLite";

  unsigned int iChannel = 0;
  for (sim::SimPhotonsLite const& photons : Photons) {

    mf::LogVerbatim log(fOutputCategory);
    // a bit of a header
    log << "[#" << (iChannel++) << "] ";
    DumpPhoton(log, photons, "  ");

  }                                         // for
  mf::LogVerbatim(fOutputCategory) << "\n"; // just an empty line

} // sim::DumpSimPhotonsLite::analyze()

//------------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpSimPhotonsLite)

//------------------------------------------------------------------------------
