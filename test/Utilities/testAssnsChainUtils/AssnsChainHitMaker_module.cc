/**
 * @file   AssnsChainHitMaker_module.cc
 * @brief  Test producer creating a few dummy hits.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 23, 2017
 */

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"  // raw namespace
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo namespace
#include "lardataobj/RecoBase/Hit.h"

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <memory>  // std::make_unique()
#include <utility> // std::move()

namespace lar {
  namespace test {

    // -------------------------------------------------------------------------
    /**
    * @brief Creates some dummy hits.
    *
    * The produced hits are not associated to wires or raw digits.
    *
    * Configuration parameters
    * =========================
    *
    * * *nHits* (unsigned integer, default: 100): number of hits to produce
    *
    */
    class AssnsChainHitMaker : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<unsigned int> nHits{Name("nHits"),
                                        Comment("number of dummy hits to be generated"),
                                        100};

      }; // struct Config

      using Parameters = art::EDProducer::Table<Config>;

      explicit AssnsChainHitMaker(Parameters const& config)
        : EDProducer{config}, nHits(config().nHits())
      {
        produces<std::vector<recob::Hit>>();
      }

      virtual void produce(art::Event& event) override;

    private:
      unsigned int nHits; ///< Number of hits to be generated.

    }; // AssnsChainHitMaker

    // -------------------------------------------------------------------------

  } // namespace test
} // namespace lar

// -----------------------------------------------------------------------------
void lar::test::AssnsChainHitMaker::produce(art::Event& event)
{

  auto hits = std::make_unique<std::vector<recob::Hit>>();

  for (unsigned int i = 0; i < nHits; ++i) {

    hits->emplace_back(raw::ChannelID_t(i + 1),     // channel
                       raw::TDCtick_t(10 * i),      // start_tick
                       raw::TDCtick_t(10 * i + 40), // end_tick
                       10.0 * i,                    // peak_time
                       1.0,                         // sigma_peak_time
                       0.5,                         // RMS
                       200.0 + i,                   // peak_amplitude
                       10.0,                        // sigma_peak_amplitude
                       400.0,                       // summedADC
                       400.0,                       // hit_integral
                       10.0,                        // hit_sigma_integral
                       1,                           // multiplicity
                       0,                           // local_index
                       1.0,                         // goodness_of_fit
                       37,                          // DOF
                       geo::kUnknown,               // view
                       geo::kMysteryType,           // signal_type
                       geo::WireID{0, 1, 2, i}      // wire ID
    );

  } // for

  mf::LogInfo("AssnsChainHitMaker") << "Produced " << hits->size() << " hits.";

  event.put(std::move(hits));

} // lar::test::AssnsChainHitMaker::produce()

// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(lar::test::AssnsChainHitMaker)

// -----------------------------------------------------------------------------
