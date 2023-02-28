/**
 * @file   AssnsChainPFParticleMaker_module.cc
 * @brief  Test producer creating a few dummy PFParticles associated with
 *         clusters.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 26, 2017
 */

// LArSoft libraries
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"

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
    * @brief Creates some dummy PFParticles and associations to clusters.
    *
    * Configuration parameters
    * =========================
    *
    * * *clusters* (list of input tags): collections of the clusters to be
    *     used if PFParticles
    * * *clustersPerPFO* (unsigned integer, default: 3): number of clusters
    *   combined into each PFParticle
    *
    */
    class AssnsChainPFParticleMaker : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Sequence<art::InputTag> clusters{Name("clusters"),
                                                Comment("collections of clusters to be combined")};

        fhicl::Atom<unsigned int> clustersPerPFO{
          Name("clustersPerPFO"),
          Comment("number of clusters combined into each PFParticle"),
          3};

      }; // struct Config

      using Parameters = art::EDProducer::Table<Config>;

      explicit AssnsChainPFParticleMaker(Parameters const& config)
        : EDProducer{config}
        , clusterTags(config().clusters())
        , nClustersPerPFO(config().clustersPerPFO())
      {
        produces<std::vector<recob::PFParticle>>();
        produces<art::Assns<recob::Cluster, recob::PFParticle>>();
      }

      virtual void produce(art::Event& event) override;

    private:
      std::vector<art::InputTag> clusterTags; ///< List of cluster tags.
      unsigned int nClustersPerPFO;           ///< Maximum number of clusters per PFO.

      /// Returns a list of clusters to be combined.
      std::vector<art::Ptr<recob::Cluster>> collectClusters(art::Event const& event) const;

    }; // AssnsChainPFParticleMaker

    // -------------------------------------------------------------------------

  } // namespace test
} // namespace lar

// -----------------------------------------------------------------------------
void lar::test::AssnsChainPFParticleMaker::produce(art::Event& event)
{

  //
  // prepare input: merge all hits in a single collection
  //
  std::vector<art::Ptr<recob::Cluster>> clusters = collectClusters(event);

  //
  // prepare output
  //
  auto PFOs = std::make_unique<std::vector<recob::PFParticle>>();
  auto clusterPFOassns = std::make_unique<art::Assns<recob::Cluster, recob::PFParticle>>();

  //
  // create the PFParticles
  //
  unsigned int nPFOs = clusters.size() / nClustersPerPFO;
  if (nPFOs * nClustersPerPFO < clusters.size()) ++nPFOs;

  art::PtrMaker<recob::PFParticle> ptrMaker(event);

  unsigned int nDaughtersPerParticle = 2;
  unsigned int nParticlesInTier = 1;
  std::size_t firstPFOinThisTier = 0;
  std::size_t firstPFOinNextTier = firstPFOinThisTier + nParticlesInTier;
  std::size_t nextDaughter = firstPFOinNextTier;
  std::vector<std::size_t> parents(nPFOs, std::size_t(recob::PFParticle::kPFParticlePrimary));
  for (unsigned int i = 0; i < nPFOs; ++i) {

    // prepare for the next tier, if needed
    if (i >= firstPFOinNextTier) {
      firstPFOinThisTier = firstPFOinNextTier;

      nParticlesInTier *= nDaughtersPerParticle;
      nDaughtersPerParticle += 1;
      firstPFOinNextTier += nParticlesInTier;
    } // if

    //
    // assign clusters to PFO
    //
    std::vector<art::Ptr<recob::Cluster>> PFOclusters;
    std::size_t iCluster = i;
    while (iCluster < clusters.size()) {
      PFOclusters.push_back(clusters[iCluster]);
      iCluster += nPFOs;
    } // while

    //
    // generate the PFParticle
    //
    std::size_t const endDaughter =
      std::min(nextDaughter + nDaughtersPerParticle, (std::size_t)nPFOs);
    std::vector<std::size_t> daughters;
    daughters.reserve(endDaughter - nextDaughter);
    while (nextDaughter < endDaughter) {
      parents[nextDaughter] = i;
      daughters.push_back(nextDaughter);
      ++nextDaughter;
    } // while

    PFOs->emplace_back(11,         // pdgCode (11 = shower-like)
                       i,          // self
                       parents[i], // parent
                       std::move(daughters));

    //
    // generate associations
    //
    auto const PFOptr = ptrMaker(i); // art pointer to the new PFParticle
    for (art::Ptr<recob::Cluster> const& cluster : PFOclusters) {
      mf::LogVerbatim("AssnsChainPFParticleMaker")
        << "Associating cluster " << cluster << " with PFO " << PFOptr;
      clusterPFOassns->addSingle(cluster, PFOptr);
    }

  } // for

  mf::LogInfo("AssnsChainPFParticleMaker")
    << "Created " << PFOs->size() << " particle flow objects with about " << nClustersPerPFO
    << " clusters each from " << clusters.size() << " clusters and " << clusterPFOassns->size()
    << " associations from " << clusterTags.size() << " collections";

  event.put(std::move(PFOs));
  event.put(std::move(clusterPFOassns));

} // lar::test::AssnsChainClusterMaker::produce()

// -----------------------------------------------------------------------------
std::vector<art::Ptr<recob::Cluster>> lar::test::AssnsChainPFParticleMaker::collectClusters(
  art::Event const& event) const
{

  std::vector<art::Ptr<recob::Cluster>> allClusters;

  for (auto const& tag : clusterTags) {
    auto clusters = event.getValidHandle<std::vector<recob::Cluster>>(tag);

    std::size_t const nClusters = clusters->size();
    for (std::size_t i = 0; i < nClusters; ++i)
      allClusters.emplace_back(clusters, i);

  } // for

  return allClusters;
} // lar::test::AssnsChainPFParticleMaker::collectClusters()

// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(lar::test::AssnsChainPFParticleMaker)

// -----------------------------------------------------------------------------
