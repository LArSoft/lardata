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
  )
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
  )
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
  )
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
  )
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
  )

cet_test(AssnsChainUtil_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config assnschainutils_test.fcl
  DATAFILES assnschainutils_test.fcl
  )
