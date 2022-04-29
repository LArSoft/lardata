
# this test requires larcore/Geometry/geometry_bo.fcl
# DetectorPropertiesStandard_test is built by lardataalg
cet_test( DetectorPropertiesBo_test
  HANDBUILT
  DATAFILES dettest_bo.fcl
  TEST_EXEC lardataalg::DetectorPropertiesStandard_test
  TEST_ARGS ./dettest_bo.fcl
)

install_source()
install_fhicl()
