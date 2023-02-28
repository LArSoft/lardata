/**
 * @file   AssnsChainShowerMaker_module.cc
 * @brief  Test producer creating a few dummy showers associated with
 *         PFParticle.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 26, 2017
 */

// LArSoft libraries
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <memory>  // std::make_unique()
#include <utility> // std::move()

namespace lar {
  namespace test {

    // -------------------------------------------------------------------------
    /**
    * @brief Creates some dummy showers and associations to PFParticle objects.
    *
    * Configuration parameters
    * =========================
    *
    * * *particles* (list of input tags): collections of the particle flow
    *     objects to be made into showers
    *
    */
    class AssnsChainShowerMaker : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Sequence<art::InputTag> particles{
          Name("particles"),
          Comment("collections of particle flow objects to be made into showers")};

      }; // struct Config

      using Parameters = art::EDProducer::Table<Config>;

      explicit AssnsChainShowerMaker(Parameters const& config)
        : EDProducer{config}, particleTags(config().particles())
      {
        produces<std::vector<recob::Shower>>();
        produces<art::Assns<recob::PFParticle, recob::Shower>>();
      }

      virtual void produce(art::Event& event) override;

    private:
      std::vector<art::InputTag> particleTags; ///< List of PFParticle tags.

      /// Returns a list of PFParticle objects to be made into showers.
      std::vector<art::Ptr<recob::PFParticle>> collectPFOs(art::Event const& event) const;

    }; // AssnsChainShowerMaker

    // -------------------------------------------------------------------------

  } // namespace test
} // namespace lar

// -----------------------------------------------------------------------------
void lar::test::AssnsChainShowerMaker::produce(art::Event& event)
{

  //
  // prepare input: merge all hits in a single collection
  //
  std::vector<art::Ptr<recob::PFParticle>> particles = collectPFOs(event);

  //
  // prepare output
  //
  auto showers = std::make_unique<std::vector<recob::Shower>>();
  auto PFOshowerAssns = std::make_unique<art::Assns<recob::PFParticle, recob::Shower>>();

  //
  // create the showers
  //
  unsigned int nShowers = particles.size();

  art::PtrMaker<recob::Shower> ptrMaker(event);

  for (unsigned int i = 0; i < nShowers; ++i) {

    //
    // generate the shower
    //
    showers->push_back(recob::Shower({0.0, 0.0, 1.0}, // dcosVtx
                                     {0.1, 0.1, 0.1}, // dcosVtxErr
                                     {0.0, 0.0, 0.0}, // xyz
                                     {1.0, 1.0, 1.0}, // xyzErr
                                     {1.0, 1.0, 1.0}, // TotalEnergy
                                     {0.1, 0.1, 0.1}, // TotalEnergyErr
                                     {2.0, 2.0, 2.0}, // dEdx
                                     {0.1, 0.1, 0.1}, // dEdxErr
                                     0,               // bestplane
                                     i,               // id
                                     1.0,             // length
                                     1.0              // openAngle
                                     ));

    //
    // generate associations
    //
    PFOshowerAssns->addSingle(particles[i], ptrMaker(i));

  } // for

  mf::LogInfo("AssnsChainShowerMaker")
    << "Created " << showers->size() << " showers from " << particles.size()
    << " particle flow objects and " << PFOshowerAssns->size() << " associations from "
    << particleTags.size() << " collections";

  event.put(std::move(showers));
  event.put(std::move(PFOshowerAssns));

} // lar::test::AssnsChainShowerMaker::produce()

// -----------------------------------------------------------------------------
std::vector<art::Ptr<recob::PFParticle>> lar::test::AssnsChainShowerMaker::collectPFOs(
  art::Event const& event) const
{

  std::vector<art::Ptr<recob::PFParticle>> allPFOs;

  for (auto const& tag : particleTags) {
    auto PFOs = event.getValidHandle<std::vector<recob::PFParticle>>(tag);

    std::size_t const nPFOs = PFOs->size();
    for (std::size_t i = 0; i < nPFOs; ++i)
      allPFOs.emplace_back(PFOs, i);

  } // for

  return allPFOs;
} // lar::test::AssnsChainShowerMaker::collectHits()

// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(lar::test::AssnsChainShowerMaker)

// -----------------------------------------------------------------------------
