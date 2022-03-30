cet_build_plugin(IntVectorProducer          art::module  NO_INSTALL
   LIBRARIES PRIVATE cetlib_except::cetlib_except art::Persistency_Provenance)
cet_build_plugin(PtrMakerProducer           art::module  NO_INSTALL 
   LIBRARIES PRIVATE cetlib_except::cetlib_except art::Persistency_Provenance)
cet_build_plugin(PtrMakerProducer2          art::module  NO_INSTALL 
   LIBRARIES PRIVATE cetlib_except::cetlib_except art::Persistency_Provenance)
cet_build_plugin(PtrMakerAnalyzer           art::module  NO_INSTALL
   LIBRARIES PRIVATE cetlib_except::cetlib_except art::Persistency_Provenance)

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
