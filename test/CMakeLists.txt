# ======================================================================
#
# Testing
#
# ======================================================================

include(CetTest)
cet_enable_asserts()

cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR})

cet_transitive_paths(GDML_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_GDML_DIR})

cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})

add_subdirectory( RecoObjects )
add_subdirectory( Utilities )
add_subdirectory( DetectorInfo )
add_subdirectory( DetectorInfoServices )
add_subdirectory( ArtDataHelper )
add_subdirectory( RecoBaseProxy )

# various integration tests

# these tests run a FCL file and fail only if lar exits with a bad exit code;
# they are fast and always run.
# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
# OPTIONAL_GROUPS is set to execute the test in the "BROKEN" test suite;
# that is because this test is currently failing for MicroBooNE
cet_test(donothing HANDBUILT
  DATAFILES test_donothing.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./test_donothing.fcl
)
