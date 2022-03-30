###############################################################################
###  ProxyBase tests
###
add_subdirectory(ProxyBase)

cet_build_plugin(ProxyBaseTest art::module
  LIBRARIES PRIVATE
  lardata_RecoBaseProxy
  lardataobj::RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  ROOT::GenVector
  USE_BOOST_UNIT
  )

cet_test(ProxyBase_test
  HANDBUILT
  DATAFILES test_proxybase.fcl
  TEST_EXEC lar_ut
  TEST_ARGS -- --rethrow-all -c ./test_proxybase.fcl
  USE_BOOST_UNIT
  )


###############################################################################
###  TrackProxy tests
###
cet_build_plugin(TrackProxyTest art::module
  LIBRARIES PRIVATE
  lardata_RecoBaseProxy
  lardataobj::RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  ROOT::GenVector
  USE_BOOST_UNIT
  )

cet_build_plugin(TrackProxyHitMaker art::module
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  )

cet_build_plugin(TrackProxyTrackMaker art::module
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
  USE_BOOST_UNIT
  )

###############################################################################
###  ChargedSpacePointProxy tests
###
cet_build_plugin(ChargedSpacePointProxyTest art::module
  LIBRARIES PRIVATE
  lardata_RecoBaseProxy
  lardataobj::RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  ROOT::GenVector
  USE_BOOST_UNIT
  )

cet_build_plugin(ChargedSpacePointProxyInputMaker art::module
  LIBRARIES PRIVATE
  lardata_ArtDataHelper
  lardataobj::RecoBase
  lardata_ArtDataHelper
  messagefacility::MF_MessageLogger
  USE_BOOST_UNIT
  )

cet_test(ChargedSpacePointProxy_test
  HANDBUILT
  DATAFILES test_chargedspacepointproxy.fcl
  TEST_EXEC lar_ut
  TEST_ARGS -- --rethrow-all -c ./test_chargedspacepointproxy.fcl
  USE_BOOST_UNIT
  )

###############################################################################
