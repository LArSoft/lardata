/**
 * @file   HitCollectionCreatorTest_module.cc
 * @brief  Tests classes derived from `recob::HitAndAssociationsWriterBase`.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 7, 2017
 * @see    lardata/ArtDataHelper/HitCollector.h
 */

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"  // raw::InvalidChannelID
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::kMysteryType, ...
#include "lardata/ArtDataHelper/HitCreator.h"             // recob::HitCollectionCreator
#include "lardataobj/RecoBase/Hit.h"

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// C/C++ standard libraries
#include <string>

namespace recob {
  namespace test {

    /**
     * @brief Test module for `recob::HitCollector`.
     *
     * Currently exercises:
     * @todo
     *
     * Throws an exception on failure.
     *
     * Service requirements
     * =====================
     *
     * This module requires no service.
     *
     * Configuration parameters
     * =========================
     *
     * * *instanceName* (string, default: empty): name of the data product
     *     instance to produce
     *
     */
    class HitCollectionCreatorTest : public art::EDProducer {

    public:
      struct Config {

        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<std::string> instanceName{
          Name("instanceName"),
          Comment("name of the data product instance to produce"),
          "" /* default: empty */
        };

      }; // Config

      using Parameters = art::EDProducer::Table<Config>;

      explicit HitCollectionCreatorTest(Parameters const& config);

      virtual void produce(art::Event& event) override;

    private:
      recob::HitCollectionCreatorManager hitCollManager;

      std::string fInstanceName; ///< Instance name to be used for products.

      /// Produces a collection of hits and stores it into the event.
      void produceHits(art::Event& event, std::string instanceName);

    }; // HitCollectionCreatorTest

    DEFINE_ART_MODULE(HitCollectionCreatorTest)

  } // namespace test
} // namespace recob::

//------------------------------------------------------------------------------
//--- implementation
//---
//----------------------------------------------------------------------------
recob::test::HitCollectionCreatorTest::HitCollectionCreatorTest(Parameters const& config)
  : art::EDProducer(config)
  , hitCollManager(producesCollector(),
                   config().instanceName(),
                   false /* doWireAssns */,
                   false /* doRawDigitAssns */
                   )     // produces<>() hit collections
{}

//----------------------------------------------------------------------------
void recob::test::HitCollectionCreatorTest::produce(art::Event& event)
{
  produceHits(event, fInstanceName);
} // HitCollectionCreatorTest::produce()

//----------------------------------------------------------------------------
void recob::test::HitCollectionCreatorTest::produceHits(art::Event& event, std::string instanceName)
{

  // this object will contain al the hits until they are moved into the event;
  // while it's useful to test the creation of the associations, that is too
  // onerous for this test
  auto Hits = hitCollManager.collectionWriter(event);

  // create hits, one by one
  for (double time : {0.0, 200.0, 400.0}) {
    Hits.emplace_back(recob::Hit(raw::InvalidChannelID,       /* channel */
                                 raw::TDCtick_t(1000 + time), /* start_tick */
                                 raw::TDCtick_t(1010 + time), /* end_tick */
                                 time,                        /* peak_time */
                                 1.0,                         /* sigma_peak_time */
                                 5.0,                         /* rms */
                                 100.0,                       /* peak_amplitude */
                                 1.0,                         /* sigma_peak_amplitude */
                                 500.0,                       /* ROIsummedADC */
				 500.0,                       /* HitsummedADC */
                                 500.0,                       /* hit_integral */
                                 1.0,                         /* hit_sigma_integral */
                                 1,                           /* multiplicity */
                                 0,                           /* local_index */
                                 1.0,                         /* goodness_of_fit */
                                 7,                           /* dof */
                                 geo::kUnknown,               /* view */
                                 geo::kMysteryType,           /* signal_type */
                                 geo::WireID{}                /* wireID */
                                 ));
  } // for hits

  Hits.put_into(event);

} // recob::test::HitCollectionCreatorTest::produceHits()

//----------------------------------------------------------------------------
