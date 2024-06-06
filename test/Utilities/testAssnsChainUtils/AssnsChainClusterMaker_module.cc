/**
 * @file   AssnsChainClusterMaker_module.cc
 * @brief  Test producer creating a few dummy clusters associated with hits.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 23, 2017
 */

// LArSoft libraries
#include "art/Persistency/Common/PtrMaker.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo namespace
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

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
#include <cmath>   // std::sqrt()
#include <memory>  // std::make_unique()
#include <numeric> // std::accumulate()
#include <utility> // std::move()

namespace lar {
  namespace test {

    // -------------------------------------------------------------------------
    /**
    * @brief Creates some dummy clusters and associations to hits.
    *
    * Configuration parameters
    * =========================
    *
    * * *hits* (list of input tags): collections of the hits to be clustered
    * * *hitsPerCluster* (unsigned integer, default: 100): number of hits
    *   associated with each cluster
    *
    */
    class AssnsChainClusterMaker : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Sequence<art::InputTag> hits{Name("hits"),
                                            Comment("collections of hits to be clustered")};

        fhicl::Atom<unsigned int> hitsPerCluster{
          Name("hitsPerCluster"),
          Comment("number of hits associated with each cluster"),
          100};

      }; // struct Config

      using Parameters = art::EDProducer::Table<Config>;

      explicit AssnsChainClusterMaker(Parameters const& config)
        : EDProducer{config}, hitTags(config().hits()), nHitsPerCluster(config().hitsPerCluster())
      {
        produces<std::vector<recob::Cluster>>();
        produces<art::Assns<recob::Hit, recob::Cluster>>();
      }

      virtual void produce(art::Event& event) override;

    private:
      std::vector<art::InputTag> hitTags; ///< List of hit tags for clustering.
      unsigned int nHitsPerCluster;       ///< Maximum number of hits per cluster.

      /// Returns a list of hits to be clustered.
      std::vector<art::Ptr<recob::Hit>> collectHits(art::Event const& event) const;

    }; // AssnsChainClusterMaker

    // -------------------------------------------------------------------------

  } // namespace test
} // namespace lar

// -----------------------------------------------------------------------------
namespace {

  template <typename T>
  inline T sqr(T v)
  {
    return v * v;
  }

} // local namespace

// -----------------------------------------------------------------------------
void lar::test::AssnsChainClusterMaker::produce(art::Event& event)
{

  //
  // prepare input: merge all hits in a single collection
  //
  std::vector<art::Ptr<recob::Hit>> hits = collectHits(event);

  //
  // prepare output
  //
  auto clusters = std::make_unique<std::vector<recob::Cluster>>();
  auto hitClusterAssns = std::make_unique<art::Assns<recob::Hit, recob::Cluster>>();

  //
  // create the clusters
  //
  unsigned int nClusters = hits.size() / nHitsPerCluster;
  if (nClusters * nHitsPerCluster < hits.size()) ++nClusters;

  art::PtrMaker<recob::Cluster> ptrMaker(event);

  for (unsigned int i = 0; i < nClusters; ++i) {

    //
    // assign hits to cluster
    //
    std::vector<art::Ptr<recob::Hit>> clusterHits;
    std::size_t iHit = i;
    while (iHit < hits.size()) {
      clusterHits.push_back(hits[iHit]);
      iHit += nClusters;
    } // while

    //
    // generate the cluster
    //
    clusters->emplace_back(
      float(clusterHits.front()->WireID().Wire), // start_wire
      1.0,                                       // sigma_start_wire
      clusterHits.front()->PeakTime(),           // start_tick
      clusterHits.front()->SigmaPeakTime(),      // sigma_start_tick
      clusterHits.front()->Integral(),           // start_charge
      0.0,                                       // start_angle
      0.0,                                       // start_opening
      float(clusterHits.back()->WireID().Wire),  // end_wire
      1.0,                                       // sigma_end_wire
      clusterHits.back()->PeakTime(),            // end_tick
      clusterHits.back()->SigmaPeakTime(),       // sigma_end_tick
      clusterHits.back()->Integral(),            // end_charge
      0.0,                                       // end_angle
      0.0,                                       // end_opening
      std::accumulate(clusterHits.cbegin(),
                      clusterHits.cend(),
                      0.0,
                      [](float sum, art::Ptr<recob::Hit> const& hit) {
                        return sum + hit->Integral();
                      }), // integral
      std::sqrt(std::accumulate(clusterHits.cbegin(),
                                clusterHits.cend(),
                                0.0,
                                [](float sum, art::Ptr<recob::Hit> const& hit) {
                                  return sum + sqr(hit->SigmaIntegral());
                                })), // integral_stddev
      std::accumulate(clusterHits.cbegin(),
                      clusterHits.cend(),
                      0.0,
                      [](float sum, art::Ptr<recob::Hit> const& hit) {
                        return sum + hit->ROISummedADC();
                      }), // summedADC
      std::sqrt(std::accumulate(clusterHits.cbegin(),
                                clusterHits.cend(),
                                0.0,
                                [](float sum, art::Ptr<recob::Hit> const& hit) {
                                  return sum + hit->ROISummedADC();
                                })),             // summedADC_stddev
      clusterHits.size(),                        // n_hits
      0.0,                                       // multiple_hit_density
      2.0,                                       // width
      recob::Cluster::ID_t(i + 1),               // ID
      clusterHits.front()->View(),               // view
      clusterHits.front()->WireID().asPlaneID(), // plane
      recob::Cluster::Sentry);

    //
    // generate associations
    //
    auto const clusterPtr = ptrMaker(i); // art pointer to the new cluster
    for (art::Ptr<recob::Hit> const& hit : clusterHits)
      hitClusterAssns->addSingle(hit, clusterPtr);

  } // for

  mf::LogInfo("AssnsChainClusterMaker")
    << "Created " << clusters->size() << " clusters with about " << nHitsPerCluster
    << " hits each from " << hits.size() << " hits and " << hitClusterAssns->size()
    << " associations from " << hitTags.size() << " collections";

  event.put(std::move(clusters));
  event.put(std::move(hitClusterAssns));

} // lar::test::AssnsChainClusterMaker::produce()

// -----------------------------------------------------------------------------
std::vector<art::Ptr<recob::Hit>> lar::test::AssnsChainClusterMaker::collectHits(
  art::Event const& event) const
{

  std::vector<art::Ptr<recob::Hit>> allHits;

  for (auto const& tag : hitTags) {
    auto hits = event.getValidHandle<std::vector<recob::Hit>>(tag);

    std::size_t const nHits = hits->size();
    for (std::size_t i = 0; i < nHits; ++i)
      allHits.emplace_back(hits, i);

  } // for

  return allHits;
} // lar::test::AssnsChainClusterMaker::collectHits()

// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(lar::test::AssnsChainClusterMaker)

// -----------------------------------------------------------------------------
