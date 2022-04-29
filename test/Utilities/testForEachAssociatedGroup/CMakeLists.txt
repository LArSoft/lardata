cet_build_plugin(AssnProducer1 art::module NO_INSTALL
  LIBRARIES PRIVATE
  cetlib_except::cetlib_except
  art::Persistency_Provenance
)
cet_build_plugin(AssnProducer2 art::module NO_INSTALL
  LIBRARIES PRIVATE
  cetlib_except::cetlib_except
  art::Persistency_Provenance
)
cet_build_plugin(AssnAnalyzer art::module NO_INSTALL
  LIBRARIES PRIVATE
  cetlib_except::cetlib_except
  art::Persistency_Provenance
  range-v3::range-v3
)

cet_test(associated_groups_test USE_BOOST_UNIT
  LIBRARIES
  canvas::canvas
  cetlib::cetlib
  cetlib_except::cetlib_except
  range-v3::range-v3
  ROOT::Core
)

cet_test(associated_groups_with_left_test USE_BOOST_UNIT
  LIBRARIES
  larcorealg::CoreUtils
  canvas::canvas
  cetlib::cetlib
  cetlib_except::cetlib_except
  range-v3::range-v3
  ROOT::Core
)

cet_test(test_feag HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_feag.fcl
  DATAFILES test_feag.fcl
)

cet_test(test_feag2 HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_feag2.fcl
  DATAFILES test_feag2.fcl
  TEST_PROPERTIES DEPENDS test_feag
)

cet_test(test_associatedgroups_full HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_associatedgroups_full.fcl
  DATAFILES test_associatedgroups_full.fcl
)

cet_test(test_associatedgroups HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_associatedgroups.fcl
  DATAFILES test_associatedgroups.fcl
  TEST_PROPERTIES DEPENDS test_feag2
)

cet_test(test_associatedgroupswithleft_full HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_associatedgroupswithleft_full.fcl
  DATAFILES test_associatedgroupswithleft_full.fcl
)

cet_test(test_associatedgroupswithleft HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_associatedgroupswithleft.fcl
  DATAFILES test_associatedgroupswithleft.fcl
  TEST_PROPERTIES DEPENDS test_feag2
)

art_dictionary(NO_INSTALL
  DICTIONARY_LIBRARIES
  cetlib::cetlib
  cetlib_except::cetlib_except
)
