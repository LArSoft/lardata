cet_build_plugin(AssnsChainHitMaker art::EDProducer NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Framework_Core
  art::Framework_Principal
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)
cet_build_plugin(AssnsChainClusterMaker art::EDProducer NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Framework_Core
  art::Framework_Principal
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)
cet_build_plugin(AssnsChainPFParticleMaker art::EDProducer NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Framework_Core
  art::Framework_Principal
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)
cet_build_plugin(AssnsChainShowerMaker art::EDProducer NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Framework_Core
  art::Framework_Principal
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)
cet_build_plugin(AssnsChainTest art::EDAnalyzer NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  larcorealg::CoreUtils
  art::Framework_Core
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
