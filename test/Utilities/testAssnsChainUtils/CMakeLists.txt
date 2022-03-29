simple_plugin(
  AssnsChainHitMaker "module"
    lardataobj_RecoBase
    art_Framework_Core
    art::Framework_Principal
    art::Persistency_Provenance
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
  NO_INSTALL
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
simple_plugin(
  AssnsChainClusterMaker "module"
    lardataobj_RecoBase
    art_Framework_Core
    art::Framework_Principal
    art::Persistency_Provenance
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
  NO_INSTALL
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
simple_plugin(
  AssnsChainPFParticleMaker "module"
    lardataobj_RecoBase
    art_Framework_Core
    art::Framework_Principal
    art::Persistency_Provenance
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
  NO_INSTALL
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
simple_plugin(
  AssnsChainShowerMaker "module"
    lardataobj_RecoBase
    art_Framework_Core
    art::Framework_Principal
    art::Persistency_Provenance
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
  NO_INSTALL
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
simple_plugin(
  AssnsChainTest "module"
    lardataobj_RecoBase
    art_Framework_Core
    art::Framework_Principal
    art::Persistency_Provenance
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
    ROOT::Core
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

cet_test(AssnsChainUtil_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config assnschainutils_test.fcl
  DATAFILES assnschainutils_test.fcl
  )
