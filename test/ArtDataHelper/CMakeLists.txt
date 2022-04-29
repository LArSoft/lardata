
cet_build_plugin(HitCollectionCreatorTest art::EDProducer NO_INSTALL
  LIBRARIES PRIVATE
  lardata_ArtDataHelper
  lardataobj::RecoBase
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

cet_build_plugin(HitDataProductChecker art::EDAnalyzer NO_INSTALL
  LIBRARIES PRIVATE
  lardata_ArtDataHelper
  lardataobj::RecoBase
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

cet_test(HitCollectorTest HANDBUILT
  DATAFILES hitcollectioncreator_test.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./hitcollectioncreator_test.fcl
)

install_fhicl()
install_source()

