/**
 * @file   DumpSimEnergyDeposits_module.cc
 * @brief  Dumps on screen the content of the `sim::SimEnergyDeposit` objects.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 11, 2020
 */

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "lardataalg/MCDumpers/MCDumperUtils.h"        // sim::ParticleName()
#include "lardataalg/Utilities/quantities/energy.h"    // MeV
#include "lardataalg/Utilities/quantities/spacetime.h" // cm
#include "lardataobj/Simulation/SimEnergyDeposit.h"

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
#include "range/v3/view/enumerate.hpp"

// C//C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <string>

// -----------------------------------------------------------------------------
namespace sim {
  class DumpSimEnergyDeposits;
}
/**
 * @brief Prints the content of all the deposited energies on screen.
 *
 * This analyzer prints the content of all the hits into the
 * LogInfo/LogVerbatim stream.
 *
 * Configuration parameters
 * =========================
 *
 * - *EnergyDepositTag* (input tag, default: `"largeant:TPCActive"`):
 *   tag of data product containing the deposits to dump (memento: format is
 *   `"moduleLabel:instanceName"`;
 * - *OutputCategory* (string, default: "DumpSimEnergyDeposits"): the category
 *   used for the output (useful for filtering)
 *
 */
class sim::DumpSimEnergyDeposits : public art::EDAnalyzer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> EnergyDepositTag{
      Name("EnergyDepositTag"),
      Comment("tag of data product containing the `sim::SimEnergyDeposit` to dump"),
      art::InputTag{"largeant", "TPCActive"}};

    fhicl::Atom<bool> ShowLocation{Name("ShowLocation"),
                                   Comment("whether to show where the deposition took place"),
                                   true};

    fhicl::Atom<bool> ShowStep{
      Name("ShowStep"),
      Comment("whether to show start and end position of the particle step"),
      false};

    fhicl::Atom<bool> ShowEmission{
      Name("ShowEmission"),
      Comment("whether to show the number of photons and electrons generated"),
      true};

    fhicl::Atom<bool> SplitPhotons{
      Name("SplitPhotons"),
      Comment("whether to list fast- and slow-emitted photons separately"),
      [this]() { return ShowEmission(); },
      true};

    fhicl::Atom<std::string> OutputCategory{
      Name("OutputCategory"),
      Comment("the messagefacility category used for the output"),
      "DumpSimEnergyDeposits"};

  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;

  /// Constructor: reads the configuration.
  explicit DumpSimEnergyDeposits(Parameters const& config);

  /// Does the printing.
  void analyze(art::Event const& evt);

private:
  art::InputTag fEnergyDepositTag; ///< Tag for input data product.
  std::string fOutputCategory;     ///< Category for LogInfo output.

  bool bShowLocation = true; ///< Print the center of the deposition.
  bool bShowStep = true;     ///< Print the step ends.
  bool bShowEmission = true; ///< Print the photons and electrons emitted.
  bool bSplitPhotons = true; ///< Print photons by emission speed.

  template <typename Stream>
  void dumpEnergyDeposit(Stream& out, sim::SimEnergyDeposit const& dep) const;

}; // class sim::DumpSimEnergyDeposits

//------------------------------------------------------------------------------
//---  module implementation
//------------------------------------------------------------------------------
sim::DumpSimEnergyDeposits::DumpSimEnergyDeposits(Parameters const& config)
  : EDAnalyzer(config)
  , fEnergyDepositTag(config().EnergyDepositTag())
  , fOutputCategory(config().OutputCategory())
  , bShowLocation(config().ShowLocation())
  , bShowStep(config().ShowStep())
  , bShowEmission(config().ShowEmission())
  , bSplitPhotons(config().SplitPhotons())
{}

//------------------------------------------------------------------------------
void sim::DumpSimEnergyDeposits::analyze(art::Event const& event)
{

  using namespace util::quantities::energy_literals;
  using namespace util::quantities::space_literals;
  using util::quantities::centimeter;
  using util::quantities::megaelectronvolt;

  // fetch the data to be dumped on screen
  auto const& Deps = *(event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fEnergyDepositTag));

  mf::LogVerbatim(fOutputCategory) << "Event " << event.id() << " contains " << Deps.size() << " '"
                                   << fEnergyDepositTag.encode() << "' energy deposits";

  megaelectronvolt TotalE = 0_MeV;
  centimeter TotalLength{0.0};
  unsigned int TotalElectrons = 0U, TotalPhotons = 0U, TotalPhotonsFast = 0U, TotalPhotonsSlow = 0U;

  for (auto const& [iDep, dep] : Deps | ranges::views::enumerate) {

    // print a header for the cluster
    mf::LogVerbatim log(fOutputCategory);
    log << "[#" << iDep << "]  ";
    dumpEnergyDeposit(log, dep);

    // collect statistics
    TotalE += megaelectronvolt{dep.Energy()};
    TotalLength += centimeter{dep.StepLength()};
    TotalElectrons += dep.NumElectrons();
    TotalPhotons += dep.NumPhotons();
    TotalPhotonsSlow += dep.NumSPhotons();
    TotalPhotonsFast += dep.NumFPhotons();

  } // for depositions

  mf::LogVerbatim(fOutputCategory)
    << "Event " << event.id() << " energy deposits '" << fEnergyDepositTag.encode() << "' include "
    << TotalE << " worth of energy, " << TotalElectrons << " electrons and " << TotalPhotons
    << " photons (" << TotalPhotonsFast << " fast and " << TotalPhotonsSlow
    << " slow); tracked particles crossed " << TotalLength << " of space.";

} // sim::DumpSimEnergyDeposits::analyze()

// -----------------------------------------------------------------------------
template <typename Stream>
void sim::DumpSimEnergyDeposits::dumpEnergyDeposit(Stream& out,
                                                   sim::SimEnergyDeposit const& dep) const
{
  using util::quantities::centimeter;
  using util::quantities::megaelectronvolt, util::quantities::gigaelectronvolt;
  using util::quantities::nanosecond;

  auto const time{nanosecond(dep.Time())};
  auto const energy{megaelectronvolt(dep.Energy())};
  auto const length{centimeter(dep.StepLength())};

  out << "TrkID=" << dep.TrackID() << " (" << sim::ParticleName(dep.PdgCode()) << "): " << energy
      << " on " << time;
  if (bShowLocation) out << " at " << dep.MidPoint();
  if (bShowStep) out << " from " << dep.Start() << " to " << dep.End();
  out << " (step: " << length << ")";
  if (bShowEmission) {
    out << "; electrons: " << dep.NumElectrons();
    if (bSplitPhotons) {
      out << "; photons: " << dep.NumFPhotons() << " (fast), " << dep.NumSPhotons() << " (slow)";
    }
    else
      out << "; photons: " << dep.NumPhotons();
  }

} // sim::DumpSimEnergyDeposits::dumpEnergyDeposit()

// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(sim::DumpSimEnergyDeposits)

// -----------------------------------------------------------------------------
