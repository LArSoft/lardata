simple_plugin(IntVectorProducer          "module"  NO_INSTALL cetlib_except art_Persistency_Provenance) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
simple_plugin(PtrMakerProducer           "module"  NO_INSTALL cetlib_except art_Persistency_Provenance) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
simple_plugin(PtrMakerProducer2          "module"  NO_INSTALL cetlib_except art_Persistency_Provenance) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types
simple_plugin(PtrMakerAnalyzer           "module"  NO_INSTALL cetlib_except art_Persistency_Provenance) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.09.01) - deprecated: use cet_build_plugin() with explicit source lists and plugin base types

cet_test(test_ptrmaker HANDBUILT
  TEST_EXEC lar 
  TEST_ARGS --rethrow-all --config test_ptrmaker.fcl
  DATAFILES test_ptrmaker.fcl
)

cet_test(test_ptrmaker2 HANDBUILT
  TEST_EXEC lar 
  TEST_ARGS --rethrow-all --config test_ptrmaker2.fcl
  DATAFILES test_ptrmaker2.fcl
)

cet_test(test_ptrmaker3 HANDBUILT
  TEST_EXEC lar 
  TEST_ARGS --rethrow-all --config test_ptrmaker3.fcl
  DATAFILES test_ptrmaker3.fcl
  TEST_PROPERTIES DEPENDS test_ptrmaker
)

cet_test(test_ptrmaker4 HANDBUILT
  TEST_EXEC lar 
  TEST_ARGS --rethrow-all --config test_ptrmaker4.fcl
  DATAFILES test_ptrmaker4.fcl
  TEST_PROPERTIES DEPENDS test_ptrmaker2
)
