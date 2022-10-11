/**
 * @file   AssnsChainTest_module.cc
 * @brief  Module running through a predefined set of associations.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 23, 2017
 *
 * This test does not test almost anything.
 * The output should be manually checked to make sure it makes sense.
 */

// LArSoft libraries
#include "larcorealg/CoreUtils/MetaUtils.h" // std::always_true_v
#include "lardata/Utilities/FindManyInChainP.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <cassert>
#include <set>

namespace lar {
  namespace test {

    // -------------------------------------------------------------------------
    /**
    * @brief Prints all the hits associated to the specified shower.
    *
    * The hits are searched transversing associations from showers to particle
    * flow objects to clusters to hits.
    *
    * Configuration parameters
    * =========================
    *
    * * *showers* (input tag, required): the shower collection
    * * *hitsPerLine* (positive integer, default: 10): how many associated hits
    *     to print per line
    * * *nShowers* _(positive integer, mandatory)_: total number of expected
    *     showers
    * * *nParticles* _(positive integer, mandatory)_: total number of expected
    *     particle flow objects
    * * *nClusters* _(positive integer, mandatory)_: total number of expected
    *     clusters
    * * *nHits* _(positive integer, mandatory)_: total number of expected hits
    *
    */
    class AssnsChainTest : public art::EDAnalyzer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> showers{
          Name("showers"),
          Comment("label of the shower collection to be explored")};

        fhicl::Atom<unsigned int> hitsPerLine{Name("hitsPerLine"),
                                              Comment("how many associated hits to print per line"),
                                              10};

        fhicl::Atom<unsigned int> nHits{
          Name("nHits"),
          Comment("total number of expected hits"),
        };

        fhicl::Atom<unsigned int> nClusters{
          Name("nClusters"),
          Comment("total number of expected clusters"),
        };

        fhicl::Atom<unsigned int> nParticles{
          Name("nParticles"),
          Comment("total number of expected particle flow objects"),
        };

        fhicl::Atom<unsigned int> nShowers{
          Name("nShowers"),
          Comment("total number of expected showers"),
        };

      }; // struct Config

      using Parameters = art::EDAnalyzer::Table<Config>;

      explicit AssnsChainTest(Parameters const& config);

      virtual void analyze(art::Event const& event) override;

    private:
      art::InputTag showerTag;      ///< Label of the input collection of showers.
      unsigned int nObjectsPerLine; ///< Number of objects to print on one line.
      unsigned int nShowers;        ///< Total number of expected showers.
      unsigned int nPFOs;           ///< Total number of expected particles.
      unsigned int nClusters;       ///< Total number of expected clusters.
      unsigned int nHits;           ///< Total number of expected hits.

      template <typename ShowerHandle>
      void printAssociatedPFOs(art::Event const& event, ShowerHandle const& showers) const;

      template <typename ShowerHandle>
      void printAssociatedClusters(art::Event const& event, ShowerHandle const& showers) const;

      template <typename ShowerHandle>
      void printAssociatedHits(art::Event const& event, ShowerHandle const& showers) const;

    }; // AssnsChainTest

    // -------------------------------------------------------------------------

  } // namespace test
} // namespace lar

// -----------------------------------------------------------------------------
lar::test::AssnsChainTest::AssnsChainTest(Parameters const& config)
  : art::EDAnalyzer(config)
  , showerTag(config().showers())
  , nObjectsPerLine(config().hitsPerLine())
  , nShowers(config().nShowers())
  , nPFOs(config().nParticles())
  , nClusters(config().nClusters())
  , nHits(config().nHits())
{}

// -----------------------------------------------------------------------------
void lar::test::AssnsChainTest::analyze(art::Event const& event)
{

  //
  // read the input collection
  //
  auto showers = event.getValidHandle<std::vector<recob::Shower>>(showerTag);
  mf::LogVerbatim("AssnsChainTest") << event.id() << " contains " << showers->size()
                                    << " showers from '" << showerTag.encode() << "'";

  if (showers->size() != nShowers) {
    throw cet::exception("AssnsChainTest")
      << "Data product '" << showerTag.encode() << "' contains " << showers->size() << " showers, "
      << nShowers << " were expected.\n";
  }

  mf::LogVerbatim("AssnsChainTest") << "\nPrinting: shower particle";
  printAssociatedPFOs(event, showers);

  mf::LogVerbatim("AssnsChainTest") << "\nPrinting: shower clusters";
  printAssociatedClusters(event, showers);

  mf::LogVerbatim("AssnsChainTest") << "\nPrinting: shower hits";
  printAssociatedHits(event, showers);

} // lar::test::AssnsChainTest::analyze()

// -----------------------------------------------------------------------------
template <typename ShowerHandle>
void lar::test::AssnsChainTest::printAssociatedHits(art::Event const& event,
                                                    ShowerHandle const& showers) const
{

  //
  // get the associated hits
  //
  lar::FindManyInChainP<recob::Hit, recob::Cluster, recob::PFParticle> showerHits(
    showers, event, showerTag);
  assert(showerHits.size() == showers->size());

  //
  // print the associated hits (just the art pointer so far)
  //
  unsigned int nDuplicates = 0;
  std::set<art::Ptr<recob::Hit>> foundHits; // all hits found
  std::set<art::ProductID> foundHitProducts;

  unsigned int const pageSize = nObjectsPerLine; // hits per line
  mf::LogVerbatim log("AssnsChainTest");
  for (std::size_t iShower = 0; iShower < showers->size(); ++iShower) {
    auto const& hits = showerHits.at(iShower);
    log << "\n #" << iShower << ": " << hits.size() << " hits";
    if (hits.empty()) continue;

    unsigned int hitsLeft = pageSize;
    for (art::Ptr<recob::Hit> const& hit : hits) {
      if (foundHits.count(hit) == 0) {
        foundHits.insert(hit);
        foundHitProducts.insert(hit.id());
      }
      else {
        mf::LogProblem("AssnsChainTest")
          << "ERROR: Hit " << hit << " appears in more than one shower!";
        ++nDuplicates;
      }
      if (hitsLeft-- == 0) {
        hitsLeft = pageSize;
        log << "\n  ";
      }
      log << " " << hit; // just print the art pointer
    }                    // for hits

  } // for iShower

  mf::LogVerbatim("AssnsChainTest") << foundHits.size() << " hits collected for " << showers->size()
                                    << " showers ('" << showers.provenance()->inputTag().encode()
                                    << "') from " << foundHitProducts.size() << " data products:";
  for (art::ProductID const& PID : foundHitProducts) {
    art::Handle<std::vector<recob::Hit>> hits;
    if (event.get(PID, hits)) {
      mf::LogVerbatim("AssnsChainTest") << " - '" << hits.provenance()->inputTag().encode()
                                        << "' (contains " << hits->size() << " hits)";
    }
    else {
      mf::LogVerbatim("AssnsChainTest") << " - <" << PID << "> (not found!)";
    }
  } // for PIDs

  if (nDuplicates > 0) {
    throw cet::exception("AssnsChainTest")
      << "Test failed: " << nDuplicates << " hits appear in more than one shower.\n";
  }
  if (foundHits.size() != nHits) {
    throw cet::exception("AssnsChainTest")
      << "Test failed: counted " << foundHits.size() << " hits, expected " << nHits << ".\n";
  }

} // lar::test::AssnsChainTest::printAssociatedHits()

// -----------------------------------------------------------------------------
template <typename ShowerHandle>
void lar::test::AssnsChainTest::printAssociatedClusters(art::Event const& event,
                                                        ShowerHandle const& showers) const
{
  std::string const objectsDesc = "clusters";

  //
  // get the associated objects
  //
  lar::FindManyInChainP<recob::Cluster, recob::PFParticle> showerObjects(showers, event, showerTag);
  assert(showerObjects.size() == showers->size());

  using ObjectPtr_t = decltype(showerObjects)::TargetPtr_t;

  //
  // print the associated objects (just the art pointer so far)
  //
  unsigned int nDuplicates = 0;
  std::set<art::Ptr<recob::Cluster>> foundObjects; // all objects found
  std::set<art::ProductID> foundObjectProducts;

  unsigned int const pageSize = nObjectsPerLine; // objects per line
  mf::LogVerbatim log("AssnsChainTest");
  for (std::size_t iShower = 0; iShower < showers->size(); ++iShower) {
    auto const& objects = showerObjects.at(iShower);
    log << "\n #" << iShower << ": " << objects.size() << " " << objectsDesc;
    if (objects.empty()) continue;

    unsigned int objectsLeft = pageSize;
    for (ObjectPtr_t const& object : objects) {
      if (foundObjects.count(object) == 0) {
        foundObjects.insert(object);
        foundObjectProducts.insert(object.id());
      }
      else {
        mf::LogProblem("AssnsChainTest")
          << "ERROR: Cluster " << object << " appears in more than one shower!";
        ++nDuplicates;
      }
      if (objectsLeft-- == 0) {
        objectsLeft = pageSize;
        log << "\n  ";
      }
      log << " " << object; // just print the art pointer
    }                       // for objects

  } // for iShower

  mf::LogVerbatim("AssnsChainTest")
    << foundObjects.size() << " " << objectsDesc << " collected for " << showers->size()
    << " showers ('" << showers.provenance()->inputTag().encode() << "') from "
    << foundObjectProducts.size() << " data products:";
  for (art::ProductID const& PID : foundObjectProducts) {
    art::Handle<std::vector<recob::Cluster>> objects;
    if (event.get(PID, objects)) {
      mf::LogVerbatim("AssnsChainTest")
        << " - '" << objects.provenance()->inputTag().encode() << "' (contains " << objects->size()
        << " " << objectsDesc << ")";
    }
    else {
      mf::LogVerbatim("AssnsChainTest") << " - <" << PID << "> (not found!)";
    }
  } // for PIDs

  if (nDuplicates > 0) {
    throw cet::exception("AssnsChainTest")
      << "Test failed: " << nDuplicates << " clusters appear in more than one shower.\n";
  }
  if (foundObjects.size() != nClusters) {
    throw cet::exception("AssnsChainTest") << "Test failed: counted " << foundObjects.size()
                                           << " clusters, expected " << nClusters << ".\n";
  }

} // lar::test::AssnsChainTest::printAssociatedClusters()

// -----------------------------------------------------------------------------
template <typename ShowerHandle>
void lar::test::AssnsChainTest::printAssociatedPFOs(art::Event const& event,
                                                    ShowerHandle const& showers) const
{
  std::string const objectsDesc = "particle flow objects";

  //
  // get the associated objects
  //
  lar::FindManyInChainP<recob::PFParticle> showerObjects(showers, event, showerTag);
  assert(showerObjects.size() == showers->size());

  using ObjectPtr_t = decltype(showerObjects)::TargetPtr_t;

  //
  // print the associated objects (just the art pointer so far)
  //
  unsigned int nDuplicates = 0;
  std::set<art::Ptr<recob::PFParticle>> foundObjects; // all objects found
  std::set<art::ProductID> foundObjectProducts;

  unsigned int const pageSize = nObjectsPerLine; // objects per line
  mf::LogVerbatim log("AssnsChainTest");
  for (std::size_t iShower = 0; iShower < showers->size(); ++iShower) {
    auto const& objects = showerObjects.at(iShower);
    log << "\n #" << iShower << ": " << objects.size() << " " << objectsDesc;

    unsigned int objectsLeft = pageSize;
    for (ObjectPtr_t const& object : objects) {
      if (foundObjects.count(object) == 0) {
        foundObjects.insert(object);
        foundObjectProducts.insert(object.id());
      }
      else {
        mf::LogProblem("AssnsChainTest")
          << "ERROR: Particle " << object << " appears in more than one shower!";
        ++nDuplicates;
      }
      if (objectsLeft-- == 0) {
        objectsLeft = pageSize;
        log << "\n  ";
      }
      log << " " << object; // just print the art pointer
    }                       // for objects

    if (objects.size() != 1) {
      throw cet::exception("AssnsChainTest")
        << "all showers are expected to have one PFParticle associated, while #" << iShower
        << " has " << objects.size() << "\n";
    }

  } // for iShower

  mf::LogVerbatim("AssnsChainTest")
    << foundObjects.size() << " " << objectsDesc << " collected for " << showers->size()
    << " showers ('" << showers.provenance()->inputTag().encode() << "') from "
    << foundObjectProducts.size() << " data products:";
  for (art::ProductID const& PID : foundObjectProducts) {
    art::Handle<std::vector<recob::PFParticle>> objects;
    if (event.get(PID, objects)) {
      mf::LogVerbatim("AssnsChainTest")
        << " - '" << objects.provenance()->inputTag().encode() << "' (contains " << objects->size()
        << " " << objectsDesc << ")";
    }
    else {
      mf::LogVerbatim("AssnsChainTest") << " - <" << PID << "> (not found!)";
    }
  } // for PIDs

  if (nDuplicates > 0) {
    throw cet::exception("AssnsChainTest")
      << "Test failed: " << nDuplicates
      << " particle flow objects appear in more than one shower.\n";
  }
  if (foundObjects.size() != nPFOs) {
    throw cet::exception("AssnsChainTest") << "Test failed: counted " << foundObjects.size()
                                           << " particle flow objects, expected " << nPFOs << ".\n";
  }

} // lar::test::AssnsChainTest::printAssociatedPFOs()

// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(lar::test::AssnsChainTest)

// -----------------------------------------------------------------------------
