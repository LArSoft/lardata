add_subdirectory( testPtrMaker )
add_subdirectory( testForEachAssociatedGroup )
add_subdirectory( testAssnsChainUtils )

# BulkAllocator_test, NestedIterator_test, CountersMap_test 
# and test pure header libraries (they are templates)

# test removed per issue #19494
# cet_test(BulkAllocator_test USE_BOOST_UNIT)

cet_test(NestedIterator_test USE_BOOST_UNIT)
cet_test(CountersMap_test USE_BOOST_UNIT)
cet_test(FastMatrixMath_test LIBRARIES lardataalg::UtilitiesHeaders USE_BOOST_UNIT)
cet_test(SimpleFits_test LIBRARIES lardataalg::UtilitiesHeaders USE_BOOST_UNIT)
cet_test(ChiSquareAccumulator_test USE_BOOST_UNIT)
cet_test(Dereference_test USE_BOOST_UNIT)
cet_test(TensorIndices_test USE_BOOST_UNIT)
cet_test(TensorIndicesStress_test)
cet_test(GridContainers_test USE_BOOST_UNIT)
cet_test(RangeForWrapper_test USE_BOOST_UNIT)
cet_test(filterRangeFor_test USE_BOOST_UNIT)
cet_test(CollectionView_test USE_BOOST_UNIT)
cet_test(TupleLookupByTag_test LIBRARIES larcorealg::CoreUtils)

# run a FHiCL file with only ComputePi inside
cet_test(timingreference_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config timingreference_test.fcl
  DATAFILES timingreference_test.fcl
  OPTIONAL_GROUPS ONLY_SLF
)


install_fhicl()
install_source()
