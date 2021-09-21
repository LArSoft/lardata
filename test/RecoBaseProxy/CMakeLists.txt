###############################################################################
###  ProxyBase tests
###
add_subdirectory(ProxyBase)

simple_plugin(ProxyBaseTest "module"
  lardata_RecoBaseProxy
  lardataobj_RecoBase
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
simple_plugin(TrackProxyTest "module"
  lardata_RecoBaseProxy
  lardataobj_RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger

  ROOT::GenVector
  USE_BOOST_UNIT
  )

simple_plugin(TrackProxyHitMaker "module"
  lardataobj_RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger

  )

simple_plugin(TrackProxyTrackMaker "module"
  lardataobj_RecoBase
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
simple_plugin(ChargedSpacePointProxyTest "module"
  lardata_RecoBaseProxy
  lardataobj_RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger

  ROOT::GenVector
  USE_BOOST_UNIT
  )

simple_plugin(ChargedSpacePointProxyInputMaker "module"
  lardataobj_RecoBase
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
