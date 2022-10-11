/**
 * @file   TrackProxyTrackMaker_module.cc
 * @brief  Test producer creating a few dummy hits.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 31, 2017
 *
 * This started as a copy of `lar::test::AssnsChainHitMaker` test module.
 */

// LArSoft libraries
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/TrajectoryPointFlags.h"

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT
#include "Math/SMatrix.h" // ROOT::Math::SMatrixIdentity

// C/C++ standard libraries
#include <memory>  // std::make_unique()
#include <utility> // std::move()

namespace lar {
  namespace test {

    // -------------------------------------------------------------------------
    /**
    * @brief Creates some dummy hits.
    *
    * The produced tracks have completely dummy content.
    *
    * Configuration parameters
    * =========================
    *
    * * *hits* (input tag, mandatory): the data product to read the hits from
    * * *hitsPerTrack* (list of unsigned integers, mandatory): number of hits
    *     for each produced track. If there are hits left after all the tracks
    *     specified here have been created, an additional track with all those
    *     hits is created. If there are not enough hits, an exception is thrown.
    *
    */
    class TrackProxyTrackMaker : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> hitsTag{
          Name("hits"),
          Comment("tag of the recob::Hit data products to produce tracks with.")};

        fhicl::Sequence<unsigned int> hitsPerTrack{
          Name("hitsPerTrack"),
          Comment("number of hits per track; last takes all remaining ones.")};

      }; // struct Config

      using Parameters = art::EDProducer::Table<Config>;

      explicit TrackProxyTrackMaker(Parameters const& config)
        : EDProducer{config}, hitsTag(config().hitsTag()), hitsPerTrack(config().hitsPerTrack())
      {
        produces<std::vector<recob::TrackTrajectory>>();
        produces<art::Assns<recob::TrackTrajectory, recob::Hit>>();
        produces<std::vector<recob::Track>>();
        produces<std::vector<std::vector<recob::TrackFitHitInfo>>>();
        produces<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();
        produces<art::Assns<recob::Track, recob::TrackTrajectory>>();
      }

      virtual void produce(art::Event& event) override;

    private:
      art::InputTag hitsTag;                  ///< Input hit collection label.
      std::vector<unsigned int> hitsPerTrack; ///< Hits per produced track.

    }; // TrackProxyTrackMaker

    // -------------------------------------------------------------------------

  } // namespace test
} // namespace lar

// -----------------------------------------------------------------------------
void lar::test::TrackProxyTrackMaker::produce(art::Event& event)
{

  auto trajectories = std::make_unique<std::vector<recob::TrackTrajectory>>();
  auto tracks = std::make_unique<std::vector<recob::Track>>();
  auto trackFitInfo = std::make_unique<std::vector<std::vector<recob::TrackFitHitInfo>>>();
  auto hitTrackAssn = std::make_unique<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>();
  auto hitTrajectoryAssn = std::make_unique<art::Assns<recob::TrackTrajectory, recob::Hit>>();
  auto trackTrajectoryAssn = std::make_unique<art::Assns<recob::Track, recob::TrackTrajectory>>();

  auto hitHandle = event.getValidHandle<std::vector<recob::Hit>>(hitsTag);
  auto const& hits = *hitHandle;

  art::PtrMaker<recob::TrackTrajectory> trajectoryPtrMaker(event);
  art::PtrMaker<recob::Track> trackPtrMaker(event);
  unsigned int iTrack = 0;
  unsigned int usedHits = 0;
  while (usedHits < hits.size()) {

    // how many hits for this track:
    unsigned int const nTrackHits =
      (iTrack < hitsPerTrack.size()) ?
        std::min(hitsPerTrack[iTrack], (unsigned int)(hits.size() - usedHits)) :
        (hits.size() - usedHits);

    //
    // create the track trajectory and fit information
    //
    std::size_t const firstHit = usedHits;
    recob::TrackTrajectory::Positions_t pos;
    recob::TrackTrajectory::Momenta_t mom;
    recob::TrackTrajectory::Flags_t flags;
    std::vector<recob::TrackFitHitInfo> fitInfo;
    for (unsigned int iPoint = 0; iPoint < nTrackHits; ++iPoint) {
      //
      // fill base track information
      //
      using Mask_t = recob::TrajectoryPointFlags::flag::Mask_t;
      Mask_t pointFlags{-recob::TrajectoryPointFlags::flag::NoPoint -
                        recob::TrajectoryPointFlags::flag::HitIgnored -
                        recob::TrajectoryPointFlags::flag::Suspicious -
                        recob::TrajectoryPointFlags::flag::DetectorIssue};

      // one point out of seven has no valid position at all
      if (iPoint % 7 == 2) // make sure there are at least two valid points
        pointFlags.set(recob::TrajectoryPointFlags::flag::NoPoint);
      // one point out of five was made ignoring the hit
      if (iPoint % 5) pointFlags.set(recob::TrajectoryPointFlags::flag::HitIgnored);
      // one point out of three is suspicious
      if (iPoint % 3) pointFlags.set(recob::TrajectoryPointFlags::flag::Suspicious);
      // every other point has issues
      if (iPoint % 2) pointFlags.set(recob::TrajectoryPointFlags::flag::DetectorIssue);

      pos.emplace_back(iPoint, iPoint, iPoint);
      mom.emplace_back(2.0, 1.0, 0.0);
      flags.emplace_back(usedHits++, pointFlags);
      //
      // fill optional information
      //
      fitInfo.push_back({
        double(iPoint) * 2.5,            // aHitMeas
        double(iPoint) * 1.5,            // aHitMeasErr2
        {},                              // aTrackStatePar
        {ROOT::Math::SMatrixIdentity{}}, // aTrackStateCov
        hits[usedHits + iPoint].WireID() // aWireId
      });
    } // for

    // produce some "additional" trajectory (pretty much invalid)
    trajectories->emplace_back();

    trajectories->emplace_back(std::move(pos), std::move(mom), std::move(flags), true);

    //
    // create the trajectory-hit associations
    // (no hits for the invalid trajectory)
    //
    auto const trajPtr = trajectoryPtrMaker(trajectories->size() - 1U);
    for (std::size_t iHit = firstHit; iHit < usedHits; ++iHit)
      hitTrajectoryAssn->addSingle(trajPtr, {hitHandle, iHit});

    //
    // create the track
    //
    recob::Track track(trajectories->back(), 2112, 1.0, nTrackHits, {}, {}, iTrack);
    tracks->push_back(std::move(track));

    //
    // and the additional objects
    //
    trackFitInfo->push_back(std::move(fitInfo));

    //
    // create the track-hit associations
    //
    auto const trackPtr = trackPtrMaker(iTrack);
    for (std::size_t iHit = firstHit; iHit < usedHits; ++iHit) {

      auto const hitIndex = iHit - firstHit;
      recob::TrackHitMeta const hitInfo(hitIndex, 2.0 * hitIndex);
      hitTrackAssn->addSingle(trackPtr, {hitHandle, iHit}, hitInfo);

    } // for

    //
    // create the track-trajectory associations
    //
    trackTrajectoryAssn->addSingle(trackPtr, trajPtr);

    mf::LogVerbatim("TrackProxyTrackMaker")
      << "New track #" << tracks->back().ID() << " with " << nTrackHits << " hits";

    //
    // prepare for the next track
    //
    ++iTrack;

  } // while

  mf::LogInfo("TrackProxyTrackMaker")
    << "Produced " << tracks->size() << " tracks from " << usedHits << " hits.";

  event.put(std::move(trajectories));
  event.put(std::move(hitTrajectoryAssn));
  event.put(std::move(tracks));
  event.put(std::move(trackFitInfo));
  event.put(std::move(hitTrackAssn));
  event.put(std::move(trackTrajectoryAssn));

} // lar::test::TrackProxyTrackMaker::produce()

// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(lar::test::TrackProxyTrackMaker)

// -----------------------------------------------------------------------------
