###############################################################################
###  ProxyBase tests
###
add_subdirectory(ProxyBase)

cet_build_plugin(ProxyBaseTest art::EDAnalyzer USE_BOOST_UNIT NO_INSTALL
  LIBRARIES PRIVATE
  lardata_RecoBaseProxy
  lardataobj::RecoBase
  lardataalg::UtilitiesHeaders
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  ROOT::GenVector
)

cet_test(ProxyBase_test HANDBUILT
  DATAFILES test_proxybase.fcl
  TEST_EXEC lar_ut
  TEST_ARGS -- --rethrow-all -c ./test_proxybase.fcl
)


###############################################################################
###  TrackProxy tests
###
cet_build_plugin(TrackProxyTest art::EDAnalyzer USE_BOOST_UNIT NO_INSTALL
  LIBRARIES PRIVATE
  lardata_RecoBaseProxy
  lardataobj::RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  ROOT::GenVector
)

cet_build_plugin(TrackProxyHitMaker art::EDProducer NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
)

cet_build_plugin(TrackProxyTrackMaker art::EDProducer NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  ROOT::Core
)

cet_test(TrackProxy_test
  HANDBUILT
  DATAFILES test_trackproxy.fcl
  TEST_EXEC lar_ut
  TEST_ARGS -- --rethrow-all -c ./test_trackproxy.fcl
)

###############################################################################
###  ChargedSpacePointProxy tests
###
cet_build_plugin(ChargedSpacePointProxyTest art::module USE_BOOST_UNIT NO_INSTALL
  LIBRARIES PRIVATE
  lardata_RecoBaseProxy
  lardataobj::RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  ROOT::GenVector
)

cet_build_plugin(ChargedSpacePointProxyInputMaker art::module USE_BOOST_UNIT NO_INSTALL
  LIBRARIES PRIVATE
  lardata_ArtDataHelper
  lardataobj::RecoBase
  lardata_ArtDataHelper
  messagefacility::MF_MessageLogger
)

cet_test(ChargedSpacePointProxy_test
  HANDBUILT
  DATAFILES test_chargedspacepointproxy.fcl
  TEST_EXEC lar_ut
  TEST_ARGS -- --rethrow-all -c ./test_chargedspacepointproxy.fcl
)

###############################################################################
