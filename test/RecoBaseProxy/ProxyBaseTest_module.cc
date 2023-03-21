/**
 * @file   ProxyBaseTest_module.cc
 * @brief  Tests feattures of `ProxyBase.h`.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 8, 2017
 *
 */

// LArSoft libraries
#include "lardata/RecoBaseProxy/Track.h" // proxy namespace
#include "lardataalg/Utilities/StatCollector.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackFitHitInfo.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

// #include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::demangle()
// #include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"

// utility libraries
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Boost libraries
#include "boost/test/unit_test.hpp"

// C/C++ libraries
#include <algorithm> // std::for_each(), std::find()
// #include <initializer_list>
#include <cstring>     // std::strlen(), std::strcpy()
#include <memory>      // std::unique_ptr<>
#include <type_traits> // std::is_rvalue_reference<>

//------------------------------------------------------------------------------
/**
 * @brief Runs a test of `proxy::Tracks` interface.
 *
 * This module is that it uses Boost unit test library, and as such it must be
 * run with `lar_ut` instead of `lar`.
 */
class ProxyBaseTest : public art::EDAnalyzer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> tracksTag{
      Name("tracks"),
      Comment("tag of the recob::Track data products to run the test on.")};

  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit ProxyBaseTest(Parameters const& config)
    : art::EDAnalyzer(config), tracksTag(config().tracksTag())
  {}

  // Plugins should not be copied or assigned.
  ProxyBaseTest(ProxyBaseTest const&) = delete;
  ProxyBaseTest(ProxyBaseTest&&) = delete;
  ProxyBaseTest& operator=(ProxyBaseTest const&) = delete;
  ProxyBaseTest& operator=(ProxyBaseTest&&) = delete;

  virtual void analyze(art::Event const& event) override;

private:
  art::InputTag tracksTag; ///< Tag for the input tracks.

  /// An example of how to access the information via track proxy.
  void proxyUsageExample(art::Event const& event) const;

  /// Returns proxies to tracks longer than a certain length.
  auto getLongTracks(art::Event const& event, double minLength) const;

  /// Tests proxy composition.
  void testProxyComposition(art::Event const& event) const;

  /// Performs the actual test.
  void testTracks(art::Event const& event) const;

  /// Single-track processing function example.
  template <typename Track>
  void processTrack(Track const& track) const;

}; // class ProxyBaseTest

//------------------------------------------------------------------------------
namespace {

  template <typename Cont, typename Value>
  std::size_t indexOf(Cont const& cont, Value const& value)
  {
    using std::begin;
    using std::end;
    auto const cbegin = begin(cont);
    auto const cend = end(cont);
    auto const it = std::find(cbegin, cend, value);
    return (it == cend) ? std::numeric_limits<std::size_t>::max() : (it - cbegin);
  } // indexOf()

  // this is not a very good test, since it assumes that R-values are
  // necessarily different; which may be not the case
  template <typename T>
  bool areSameObject(T const& a, T const& b)
  {
    return std::addressof(a) == std::addressof(b);
  }

  template <typename T>
  std::enable_if_t<std::is_rvalue_reference<T>::value, bool> areSameObject(T&& a, T&& b)
  {
    return false;
  }

} // local namespace

//------------------------------------------------------------------------------
template <typename Track>
void ProxyBaseTest::processTrack(Track const& track) const
{

  recob::Track const& trackRef = *track;

  mf::LogVerbatim("ProxyBaseTest")
    << "[#" << track.index() << "] track " << trackRef << "  " << track->Length()
    << " cm long, with " << trackRef.NPoints() << " points and "
    << track.template get<recob::Hit>().size() << " hits:";

} // ProxyBaseTest::processTrack()

//------------------------------------------------------------------------------
namespace tag {
  struct SpecialHits {};
  struct MetadataHits {};
  struct DirectHitAssns {};
  struct DirectFitInfo {};
  struct TrackSubproxy {};
  struct FitInfoProxy {};
}

//------------------------------------------------------------------------------
void ProxyBaseTest::proxyUsageExample(art::Event const& event) const
{

  auto tracks = proxy::getCollection<std::vector<recob::Track>>(
    event, tracksTag, proxy::withAssociatedMeta<recob::Hit, recob::TrackHitMeta>());

  if (tracks.empty()) {
    mf::LogVerbatim("ProxyBaseTest") << "No tracks in '" << tracksTag.encode() << "'";
    return;
  }

  mf::LogVerbatim("ProxyBaseTest")
    << "Collection '" << tracksTag.encode() << "' contains " << tracks.size() << " tracks.";

  auto onCollection = [](recob::Hit const& hit) { return hit.SignalType() == geo::kCollection; };

  for (auto trackInfo : tracks) {

    recob::Track const& track = *trackInfo; // access to the track
    double const startTheta = track.Theta();

    double const length = trackInfo->Length(); // access to track members

    // access to associated data (returns random-access collection-like object)
    decltype(auto) hits = trackInfo.get<recob::Hit>();

    lar::util::StatCollector<float> dQds;
    double charge = 0.0;
    for (auto const& hitInfo : hits) {
      // hitInfo is equivalent to a art::Ptr<recob::Hit>
      double const hitCharge = hitInfo->Integral();
      charge += hitCharge;

      if (onCollection(*hitInfo)) {            //
        double const ds = hitInfo.data().Dx(); // access recob::TrackHitMeta
        if (ds > 0.0) dQds.add(hitCharge / ds);
      } // if on collection

    } // for hits

    mf::LogVerbatim log("ProxyBaseTest");
    log << "[#" << trackInfo.index() << "] track ID=" << track.ID() << " (" << length
        << " cm, starting with theta=" << startTheta << " rad) deposited charge=" << charge
        << " with " << hits.size() << " hits";
    if (dQds.N() > 0) {
      log << " (<dQ/ds> = " << dQds.Average() << " +/- " << dQds.RMS() << " Q/cm from " << dQds.N()
          << " hits in collection planes)";
    }

  } // for tracks

} // ProxyBaseTest::proxyUsageExample()

//------------------------------------------------------------------------------
auto ProxyBaseTest::getLongTracks(art::Event const& event, double minLength) const
{
  //
  // this code is not a particularly good practice, but it is aimed to check
  // that after the proxy collection is out of scope, elements copied from it
  // are still valid
  //
  auto tracks = proxy::getCollection<std::vector<recob::Track>>(
    event, tracksTag, proxy::withAssociated<recob::Hit>());

  std::vector<decltype(tracks)::element_proxy_t> longTracks;
  for (auto track : tracks) {
    if (track->Length() >= minLength) longTracks.push_back(track);
  } // for track
  return longTracks;

} // ProxyBaseTest::getLongTracks()

//------------------------------------------------------------------------------
void ProxyBaseTest::testProxyComposition(art::Event const& event) const
{

  auto const& expectedTracks = event.getProduct<std::vector<recob::Track>>(tracksTag);

  mf::LogInfo("ProxyBaseTest") << "Starting test on " << expectedTracks.size() << " tracks from '"
                               << tracksTag.encode() << "'";

  auto directTracks = proxy::getCollection<std::vector<recob::Track>>(
    event, tracksTag, proxy::withParallelData<std::vector<recob::TrackFitHitInfo>>());

  auto tracks = proxy::getCollection<std::vector<recob::Track>>(
    event,
    tracksTag,
    proxy::withParallelData<std::vector<recob::TrackFitHitInfo>>(),
    proxy::wrapParallelDataAs<tag::TrackSubproxy>(directTracks)
    /*
    , proxy::withCollectionProxyAs
      <std::vector<recob::TrackFitHitInfo>, tag::FitInfoProxy>
      (tracksTag, proxy::withParallelData<std::vector<recob::Track>>())
    */
  );
  BOOST_TEST(tracks.get<tag::TrackSubproxy>().data() == std::addressof(directTracks));

  static_assert(std::is_lvalue_reference<decltype(tracks.get<tag::TrackSubproxy>())>(),
                "Not reference!");

  auto const& expectedFitHitInfo =
    *(event.getValidHandle<std::vector<std::vector<recob::TrackFitHitInfo>>>(tracksTag));

  std::size_t iExpectedTrack = 0;
  for (auto trackProxy : tracks) {
    auto const& expectedTrack = expectedTracks[iExpectedTrack];
    auto const& expectedTrackFitInfo = expectedFitHitInfo[iExpectedTrack];

    auto directTrackProxy = trackProxy.get<tag::TrackSubproxy>();
    BOOST_TEST(std::addressof(*directTrackProxy) == std::addressof(expectedTrack));
    BOOST_TEST(directTrackProxy->ID() == expectedTrack.ID());
    BOOST_TEST(directTrackProxy->Length() == expectedTrack.Length());

    BOOST_TEST(std::addressof(directTrackProxy.get<std::vector<recob::TrackFitHitInfo>>()) ==
               std::addressof(expectedTrackFitInfo));
    /*
    auto fitInfoProxy = trackProxy.get<tag::FitInfoProxy>();
    BOOST_TEST
      (std::addressof(*fitInfoProxy) == std::addressof(expectedFitHitInfo));
    BOOST_TEST(
      std::addressof(fitInfoProxy.get<recob::Track>()) ==
      std::addressof(expectedTrack)
      );
    */
    ++iExpectedTrack;
  } // for

  BOOST_TEST(iExpectedTrack == expectedTracks.size());

} // ProxyBaseTest::testProxyComposition()

//------------------------------------------------------------------------------
void ProxyBaseTest::testTracks(art::Event const& event) const
{

  auto expectedTracksHandle = event.getValidHandle<std::vector<recob::Track>>(tracksTag);
  auto const& expectedTracks = *expectedTracksHandle;

  auto const& expectedTrackHitAssns =
    *(event.getValidHandle<art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta>>(tracksTag));

  mf::LogInfo("ProxyBaseTest") << "Starting test on " << expectedTracks.size() << " tracks from '"
                               << tracksTag.encode() << "'";

  art::FindManyP<recob::Hit, recob::TrackHitMeta> hitsPerTrack(
    expectedTracksHandle, event, tracksTag);

  art::FindOneP<recob::TrackTrajectory> trajectoryPerTrack(expectedTracksHandle, event, tracksTag);

  auto const& expectedTrackFitHitInfo =
    *(event.getValidHandle<std::vector<std::vector<recob::TrackFitHitInfo>>>(tracksTag));

  auto directTracks = proxy::getCollection<std::vector<recob::Track>>(
    event, tracksTag, proxy::withParallelData<std::vector<recob::TrackFitHitInfo>>());

  auto tracks = proxy::getCollection<std::vector<recob::Track>>(
    event,
    tracksTag,
    proxy::withAssociated<recob::Hit>(),
    proxy::withAssociatedAs<recob::Hit, tag::SpecialHits>(),
    proxy::withAssociatedMetaAs<recob::Hit, recob::TrackHitMeta, tag::MetadataHits>(),
    proxy::withParallelData<std::vector<recob::TrackFitHitInfo>>()
    //     , proxy::withCollectionProxyAs
    //         <std::vector<recob::Track>, tag::TrackSubproxy>
    //         (tracksTag, proxy::withParallelData<std::vector<recob::TrackFitHitInfo>>())
    ,
    proxy::wrapAssociatedAs<tag::DirectHitAssns>(expectedTrackHitAssns),
    proxy::wrapParallelDataAs<tag::DirectFitInfo>(expectedTrackFitHitInfo),
    proxy::wrapParallelDataAs<tag::TrackSubproxy>(directTracks),
    proxy::withZeroOrOne<recob::TrackTrajectory>(tracksTag));

  //
  // we try to access something we did not "register" in the proxy: space points
  //
  static_assert(!tracks.has<recob::SpacePoint>(),
                "Track proxy does NOT have space points available!!!");

  static_assert(tracks.has<std::vector<recob::TrackFitHitInfo>>(),
                "recob::TrackFitHitInfo not found!!!");

  BOOST_TEST(tracks.empty() == expectedTracks.empty());
  BOOST_TEST(tracks.size() == expectedTracks.size());

  BOOST_TEST(tracks.size() == expectedTrackFitHitInfo.size());
  decltype(auto) allFitHitInfo = tracks.get<std::vector<recob::TrackFitHitInfo>>();
  static_assert(std::is_lvalue_reference<decltype(allFitHitInfo)>(), "Copy of parallel data!");
  BOOST_TEST(allFitHitInfo.data() == std::addressof(expectedTrackFitHitInfo));

  BOOST_TEST(tracks.get<tag::DirectFitInfo>().data() == std::addressof(expectedTrackFitHitInfo));

  BOOST_TEST(directTracks.get<std::vector<recob::TrackFitHitInfo>>().data() ==
             std::addressof(expectedTrackFitHitInfo));

  BOOST_TEST(tracks.get<tag::TrackSubproxy>().data() == std::addressof(directTracks));

  auto fitHitInfoSize = std::distance(allFitHitInfo.begin(), allFitHitInfo.end());
  BOOST_TEST(fitHitInfoSize == expectedTrackFitHitInfo.size());

  std::size_t iExpectedTrack = 0;
  for (auto trackProxy : tracks) {
    auto const& expectedTrack = expectedTracks[iExpectedTrack];
    art::Ptr<recob::Track> const expectedTrackPtr{expectedTracksHandle, iExpectedTrack};
    auto const& expectedHits = hitsPerTrack.at(iExpectedTrack);
    auto const& expectedHitMeta = hitsPerTrack.data(iExpectedTrack);
    auto const& expectedFitHitInfo = expectedTrackFitHitInfo[iExpectedTrack];
    art::Ptr<recob::TrackTrajectory> const expectedTrajPtr = trajectoryPerTrack.at(iExpectedTrack);

    // proxies deliver temporary objects as elements, each time a new one
    // (although an exceedingly smart compiler might decide otherwise)
    BOOST_TEST(!areSameObject(tracks[iExpectedTrack], tracks[iExpectedTrack]));

    recob::Track const& trackRef = *trackProxy;

    auto const trackProxyCopy = trackProxy;
    BOOST_TEST(std::addressof(trackProxyCopy) != std::addressof(trackProxy));

    BOOST_TEST(std::addressof(trackRef) == std::addressof(expectedTrack));
    BOOST_TEST(std::addressof(*trackProxy) == std::addressof(expectedTrack));

    // hits
    BOOST_TEST(trackProxy.get<recob::Hit>().size() == expectedHits.size());
    for (art::Ptr<recob::Hit> const& hitPtr : trackProxy.get<recob::Hit>()) {

      // with this check we just ask the hit is there
      // (the order is not guaranteed to be the same in expected and fetched)
      BOOST_TEST(indexOf(expectedHits, hitPtr) != std::numeric_limits<std::size_t>::max());

    } // for hit

    BOOST_TEST(trackProxy.index() == iExpectedTrack);

    std::vector<recob::TrackFitHitInfo> const& fitHitInfo =
      trackProxy.get<std::vector<recob::TrackFitHitInfo>>();
    static_assert(std::is_lvalue_reference<decltype(fitHitInfo)>(),
                  "Copy of parallel data element!");
    BOOST_TEST(std::addressof(fitHitInfo) == std::addressof(expectedFitHitInfo));
    BOOST_TEST(fitHitInfo.size() == expectedFitHitInfo.size());

    BOOST_TEST(std::addressof(trackProxy.get<tag::DirectFitInfo>()) ==
               std::addressof(expectedTrackFitHitInfo[iExpectedTrack]));

    BOOST_TEST(std::addressof(trackProxyCopy.get<tag::DirectFitInfo>()) ==
               std::addressof(trackProxy.get<tag::DirectFitInfo>()));

    // subproxy elements are typically temporaries
    BOOST_TEST(std::addressof(trackProxyCopy.get<tag::TrackSubproxy>()) !=
               std::addressof(trackProxy.get<tag::TrackSubproxy>()));

    auto directTrackProxy = trackProxy.get<tag::TrackSubproxy>();
    BOOST_TEST(std::addressof(*directTrackProxy) == std::addressof(expectedTrack));
    BOOST_TEST(directTrackProxy->ID() == expectedTrack.ID());
    BOOST_TEST(directTrackProxy->Length() == expectedTrack.Length());
    BOOST_TEST(std::addressof(directTrackProxy.get<std::vector<recob::TrackFitHitInfo>>()) ==
               std::addressof(fitHitInfo));

    // "special" hits
    BOOST_TEST(trackProxy.get<tag::SpecialHits>().size() == expectedHits.size());
    for (auto const& hitPtr : trackProxy.get<tag::SpecialHits>()) {
      // hitPtr is actually not a art::Ptr, but it can be compared to one

      // the syntax of the static check is pretty horrible...
      static_assert(!std::decay_t<decltype(hitPtr)>::hasMetadata(),
                    "Expected no metadata for tag::SpecialHits");
      // easier when non-static:
      BOOST_TEST(!hitPtr.hasMetadata());

      // with this check we just ask the hit is there
      // (the order is not guaranteed to be the same in expected and fetched)
      BOOST_TEST(indexOf(expectedHits, hitPtr) != std::numeric_limits<std::size_t>::max());

    } // for special hit

    // hits with metadata
    auto const& hits = trackProxy.get<tag::MetadataHits>();
    BOOST_TEST(hits.size() == expectedHits.size());
    // - range-for loop
    unsigned int nSpecialHits = 0U;
    for (auto const& hitInfo : hits) {
      ++nSpecialHits;
      // hitPtr is actually not a art::Ptr,
      // but it can be implicitly converted into one

      // the syntax of the static check is still pretty horrible...
      static_assert(std::decay_t<decltype(hitInfo)>::hasMetadata(),
                    "Expected metadata for tag::MetadataHits");
      BOOST_TEST(hitInfo.hasMetadata());

      // conversion, as reference
      art::Ptr<recob::Hit> const& hitPtr = hitInfo;

      // with this check we just ask the hit is there
      // (the order is not guaranteed to be the same in expected and fetched)
      auto const index = indexOf(expectedHits, hitPtr);
      BOOST_TEST(index != std::numeric_limits<std::size_t>::max());

      BOOST_TEST(&(hitInfo.main()) == &expectedTrack);
      BOOST_TEST(hitInfo.mainPtr() == expectedTrackPtr);

      if (index < expectedHitMeta.size()) {
        art::Ptr<recob::Hit> const& expectedHitPtr = expectedHits.at(index);
        auto const& expectedMetadata = expectedHitMeta.at(index);

        BOOST_TEST(hitInfo.valuePtr() == hitPtr);
        BOOST_TEST(std::addressof(hitInfo.value()) == std::addressof(*hitPtr));
        BOOST_TEST(hitInfo.key() == hitPtr.key());
        BOOST_TEST(hitInfo.id() == hitPtr.id());

        if (expectedHitPtr) {
          BOOST_TEST(hitPtr);
          BOOST_TEST(hitInfo.operator->() == hitPtr);
          recob::Hit const& hit = *expectedHitPtr;
          BOOST_TEST(std::addressof(*hitInfo) == std::addressof(hit));
        }

        BOOST_TEST(hitInfo.dataPtr() == expectedMetadata);
        BOOST_TEST(&(hitInfo.data()) == expectedMetadata);

        auto hitInfoCopy = hitInfo; // copy
        BOOST_TEST(static_cast<art::Ptr<recob::Hit> const&>(hitInfo) == hitPtr);
        BOOST_TEST(&(static_cast<art::Ptr<recob::Hit> const&>(hitInfo)) == &hitPtr);

        art::Ptr<recob::Hit> hitPtrMoved = std::move(hitInfoCopy);
        BOOST_TEST(hitPtrMoved == hitPtr);

      } // if the hit is correct

    } // for special hit
    BOOST_TEST(nSpecialHits == expectedHits.size());
    // - iterator loop
    nSpecialHits = 0U;
    for (auto iHit = hits.begin(); iHit != hits.end(); ++iHit) {
      ++nSpecialHits;
      static_assert(iHit.hasMetadata(), "Expected metadata for tag::MetadataHits");

      // conversion, as reference
      art::Ptr<recob::Hit> const& hitPtr = *iHit;

      // with this check we just ask the hit is there
      // (the order is not guaranteed to be the same in expected and fetched)
      auto const index = indexOf(expectedHits, hitPtr);
      BOOST_TEST(index != std::numeric_limits<std::size_t>::max());

      BOOST_TEST(&(iHit.main()) == &expectedTrack);
      BOOST_TEST(iHit.mainPtr() == expectedTrackPtr);

      if (index < expectedHitMeta.size()) {
        art::Ptr<recob::Hit> const& expectedHitPtr = expectedHits.at(index);
        auto const* expectedMetadata = expectedHitMeta.at(index);

        BOOST_TEST(iHit.valuePtr() == hitPtr);
        BOOST_TEST(std::addressof(iHit.value()) == std::addressof(*hitPtr));

        BOOST_TEST(iHit.dataPtr() == expectedMetadata);
        BOOST_TEST(std::addressof(iHit.data()) == expectedMetadata);

        BOOST_TEST(iHit.valuePtr() == expectedHitPtr);
        BOOST_TEST(iHit.dataPtr() == expectedMetadata);

      } // if the hit is correct

    } // for special hit (iterator)
    BOOST_TEST(nSpecialHits == expectedHits.size());

    BOOST_TEST(trackProxy.get<tag::DirectHitAssns>().size() == expectedHits.size());

    // direct interface to recob::Track
    BOOST_TEST(trackProxy->NPoints() == expectedTrack.NPoints());

    // trajectory?
    BOOST_TEST(trackProxy.has<recob::TrackTrajectory>() == !expectedTrajPtr.isNull());
    if (expectedTrajPtr.isNull()) { BOOST_TEST(!(trackProxy.get<recob::TrackTrajectory>())); }
    else {
      BOOST_TEST(trackProxy.get<recob::TrackTrajectory>(), expectedTrajPtr);
    }
    ++iExpectedTrack;
  } // for
  BOOST_TEST(iExpectedTrack == expectedTracks.size());

} // ProxyBaseTest::testTracks()

//------------------------------------------------------------------------------
void ProxyBaseTest::analyze(art::Event const& event)
{

  // "test" that track proxies survive their collection (part I)
  const double minLength = 30.0;
  auto const longTracks = getLongTracks(event, minLength);

  // usage example (supposed to be educational)
  proxyUsageExample(event);

  // actual test
  testTracks(event);

  // test proxy composition
  testProxyComposition(event);

  // "test" that track proxies survive their collection (part II)
  mf::LogVerbatim("ProxyBaseTest")
    << longTracks.size() << " tracks are longer than " << minLength << " cm:";
  std::for_each(
    longTracks.begin(), longTracks.end(), [this](auto const& track) { this->processTrack(track); });

} // ProxyBaseTest::analyze()

//------------------------------------------------------------------------------

DEFINE_ART_MODULE(ProxyBaseTest)
