cet_build_plugin(AssnsChainHitMaker art::module NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Framework_Core
  art::Framework_Principal
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)
cet_build_plugin(AssnsChainClusterMaker art::module NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Framework_Core
  art::Framework_Principal
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)
cet_build_plugin(AssnsChainPFParticleMaker art::module NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Framework_Core
  art::Framework_Principal
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)
cet_build_plugin(AssnsChainShowerMaker art::module NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  art::Framework_Core
  art::Framework_Principal
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)
cet_build_plugin(AssnsChainTest art::module NO_INSTALL
  LIBRARIES PRIVATE
  lardataobj::RecoBase
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
