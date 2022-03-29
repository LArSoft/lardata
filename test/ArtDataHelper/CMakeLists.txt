art_make(
  MODULE_LIBRARIES
    lardata_ArtDataHelper
    lardataobj_RecoBase
    messagefacility::MF_MessageLogger
    
    fhiclcpp::fhiclcpp
    cetlib_except::cetlib_except
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

cet_test(HitCollectorTest HANDBUILT
  DATAFILES hitcollectioncreator_test.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./hitcollectioncreator_test.fcl
  )

install_headers()
install_fhicl()
install_source()

