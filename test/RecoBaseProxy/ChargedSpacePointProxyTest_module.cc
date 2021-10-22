/**
 * @file   ChargedSpacePointProxyTest_module.cc
 * @brief  Tests `proxy::ChargedSpacePoints` proxy.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 20, 2017
 *
 */


// LArSoft libraries
#include "lardata/RecoBaseProxy/ChargedSpacePoints.h" // proxy namespace
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect namespace
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"

// framework libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

// utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Comment.h"

// Boost libraries
#include "boost/test/unit_test.hpp"

// C/C++ libraries
#include <memory> // std::addressof()
#include <type_traits>
#include <cassert>


//------------------------------------------------------------------------------
/**
 * @brief Runs a test of `proxy::ChargedSpacePoints` interface.
 *
 * This module is that it uses Boost unit test library, and as such it must be
 * run with `lar_ut` instead of `lar`.
 */
class ChargedSpacePointProxyTest: public art::EDAnalyzer {
    public:

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> pointsTag{
      Name("points"),
      Comment
        ("tag of the recob::SpacePoint and recob::PointCharge data products.")
      };

  }; // struct Config

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit ChargedSpacePointProxyTest(Parameters const& config)
    : art::EDAnalyzer(config)
    , pointsTag(config().pointsTag())
    {}

  // Plugins should not be copied or assigned.
  ChargedSpacePointProxyTest(ChargedSpacePointProxyTest const &) = delete;
  ChargedSpacePointProxyTest(ChargedSpacePointProxyTest&&) = delete;
  ChargedSpacePointProxyTest& operator=
    (ChargedSpacePointProxyTest const &) = delete;
  ChargedSpacePointProxyTest& operator= (ChargedSpacePointProxyTest&&) = delete;

  virtual void analyze(art::Event const& event) override;

    private:
  art::InputTag pointsTag; ///< Tag for the input.

  /// An example of how to access the information via proxy.
  void proxyUsageExample(art::Event const& event) const;

  /// Performs the actual test.
  void testChargedSpacePoints(art::Event const& event);


}; // class ChargedSpacePointProxyTest


//------------------------------------------------------------------------------
void ChargedSpacePointProxyTest::proxyUsageExample
  (art::Event const& event) const
{

  auto points = proxy::getChargedSpacePoints(event, pointsTag);

  if (points.empty()) {
    mf::LogVerbatim("ProxyTest")
      << "No points in '" << pointsTag.encode() << "'";
    return;
  }

  mf::LogVerbatim log("ProxyTest");
  for (auto point: points) {
    log << "\nPoint at " << point.position() << " (ID=" << point.ID()
      << ") has ";
    if (point.hasCharge()) log << "charge " << point.charge();
    else                   log << "no charge";
  } // for point

  mf::LogVerbatim("ProxyTest") << "Collection '" << pointsTag.encode()
    << "' contains " << points.size() << " points.";

} // ChargedSpacePointProxyTest::proxyUsageExample()


//------------------------------------------------------------------------------
void ChargedSpacePointProxyTest::testChargedSpacePoints
  (art::Event const& event)
{

  auto const& expectedSpacePoints
    = *(event.getValidHandle<std::vector<recob::SpacePoint>>(pointsTag));

  auto const& expectedCharges
    = *(event.getValidHandle<std::vector<recob::PointCharge>>(pointsTag));

  mf::LogInfo("ProxyTest")
    << "Starting test on " << expectedSpacePoints.size() << " points and "
    << expectedCharges.size() << " charges from '"
    << pointsTag.encode() << "'";

  // this assertion fails on invalid input (test bug)
  assert(expectedSpacePoints.size() == expectedCharges.size());

  auto points = proxy::getChargedSpacePoints(event, pointsTag);

  static_assert(points.has<recob::PointCharge>(), "recob::Charge not found!!!");

  BOOST_TEST(points.empty() == expectedSpacePoints.empty());
  BOOST_TEST(points.size() == expectedSpacePoints.size());

  decltype(auto) spacePoints = points.spacePoints();
  BOOST_TEST
    (std::addressof(spacePoints) == std::addressof(expectedSpacePoints));
  BOOST_TEST(spacePoints.size() == expectedSpacePoints.size());

  decltype(auto) charges = points.charges();
  BOOST_TEST(std::addressof(charges) == std::addressof(expectedCharges));
  BOOST_TEST(charges.size() == expectedCharges.size());

  std::size_t iExpectedPoint = 0;
  for (auto pointProxy: points) {
    auto const& expectedSpacePoint = expectedSpacePoints[iExpectedPoint];
    auto const& expectedChargeInfo = expectedCharges[iExpectedPoint];

    recob::SpacePoint const& spacePointRef = *pointProxy;

    BOOST_TEST
      (std::addressof(spacePointRef) == std::addressof(expectedSpacePoint));
    BOOST_TEST
      (std::addressof(pointProxy.point()) == std::addressof(expectedSpacePoint));
    BOOST_TEST(pointProxy.position() == geo::vect::makePointFromCoords(expectedSpacePoint.XYZ()));
    BOOST_TEST(pointProxy.ID() == expectedSpacePoint.ID());
    BOOST_TEST(pointProxy.hasCharge() == expectedChargeInfo.hasCharge());
    BOOST_TEST(pointProxy.charge() == expectedChargeInfo.charge());

    decltype(auto) chargeInfo = pointProxy.get<recob::PointCharge>();
    static_assert(
      std::is_lvalue_reference<decltype(chargeInfo)>(),
      "Copy of parallel data element!"
      );
    BOOST_TEST
      (std::addressof(chargeInfo) == std::addressof(expectedChargeInfo));

    ++iExpectedPoint;
  } // for
  BOOST_TEST(iExpectedPoint == expectedSpacePoints.size());

} // ChargedSpacePointProxyTest::testChargedSpacePoints()


//------------------------------------------------------------------------------
void ChargedSpacePointProxyTest::analyze(art::Event const& event) {

  // usage example (supposed to be educational)
  proxyUsageExample(event);

  // actual test
  testChargedSpacePoints(event);

} // ChargedSpacePointProxyTest::analyze()


//------------------------------------------------------------------------------

DEFINE_ART_MODULE(ChargedSpacePointProxyTest)
