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
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

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
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(TrackProxyHitMaker "module"
  lardataobj_RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger

  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(TrackProxyTrackMaker "module"
  lardataobj_RecoBase
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  ROOT::Core
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

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
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

simple_plugin(ChargedSpacePointProxyInputMaker "module"
  lardataobj_RecoBase
  lardata_ArtDataHelper
  messagefacility::MF_MessageLogger

  USE_BOOST_UNIT
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

cet_test(ChargedSpacePointProxy_test
  HANDBUILT
  DATAFILES test_chargedspacepointproxy.fcl
  TEST_EXEC lar_ut
  TEST_ARGS -- --rethrow-all -c ./test_chargedspacepointproxy.fcl
  USE_BOOST_UNIT
  )

###############################################################################
