art_make(
  MODULE_LIBRARIES
    larcorealg_Geometry
    lardataalg_DetectorInfo
    art::Framework_Services_Registry
    messagefacility::MF_MessageLogger
    ROOT::Core
)

cet_test(LarProp HANDBUILT
  DATAFILES lartest.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./lartest.fcl
  )

cet_test(ServicePackTest HANDBUILT
  DATAFILES servicepacktest.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./servicepacktest.fcl
  )

# ------------------------------------------------------------------------------
# ---  LArProperties service
# ---
cet_test( LArPropertiesServiceStandard_test
  HANDBUILT
  DATAFILES larservicetest_standard.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./larservicetest_standard.fcl
)

cet_test( LArPropertiesServiceLArTPCdetector_test
  HANDBUILT
  DATAFILES larservicetest_lartpcdetector.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./larservicetest_lartpcdetector.fcl
)

cet_test( LArPropertiesServiceBo_test
  HANDBUILT
  DATAFILES larservicetest_bo.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./larservicetest_bo.fcl
)

# ------------------------------------------------------------------------------
# ---  DetectorProperties service
# ---
cet_test( DetectorPropertiesServiceStandard_test
  HANDBUILT
  DATAFILES detectorpropertiesservicetest_standard.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./detectorpropertiesservicetest_standard.fcl
)

cet_test( DetectorPropertiesServiceLArTPCdetector_test
  HANDBUILT
  DATAFILES detectorpropertiesservicetest_lartpcdetector.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./detectorpropertiesservicetest_lartpcdetector.fcl
)

cet_test( DetectorPropertiesServiceBo_test
  HANDBUILT
  DATAFILES detectorpropertiesservicetest_bo.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./detectorpropertiesservicetest_bo.fcl
)

# ------------------------------------------------------------------------------


install_headers()
install_fhicl()
install_source()
